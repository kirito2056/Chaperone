use crate::forcefield::native::NativeContacts;
use crate::system::{Real, System};

pub const CONTACT_TOLERANCE: Real = 1.2;

pub fn fraction_of_native_contacts(
    sys: &System,
    contacts: &NativeContacts,
    tol: Real,
) -> Option<Real> {
    if contacts.is_empty() {
        return None;
    }

    let mut formed = 0usize;
    for c in 0..contacts.len() {
        let i = contacts.i[c] as usize;
        let j = contacts.j[c] as usize;

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r2 = dx * dx + dy * dy + dz * dz;

        let threshold = tol * contacts.sigma[c];
        if r2 < threshold * threshold {
            formed += 1;
        }
    }

    Some(formed as Real / contacts.len() as Real)
}

pub struct PeriodTracker {
    prev_offset: Real,
    prev_time: Real,
    first: Option<Real>,
    last: Option<Real>,
    count: usize,
}

impl PeriodTracker {
    pub fn new(initial_offset: Real, initial_time: Real) -> Self {
        PeriodTracker {
            prev_offset: initial_offset,
            prev_time: initial_time,
            first: None,
            last: None,
            count: 0,
        }
    }

    pub fn push(&mut self, offset: Real, time: Real) {
        if self.prev_offset < 0.0 && offset >= 0.0 {
            let frac = -self.prev_offset / (offset - self.prev_offset);
            let crossing = self.prev_time + frac * (time - self.prev_time);
            if self.first.is_none() {
                self.first = Some(crossing);
            }
            self.last = Some(crossing);
            self.count += 1;
        }
        self.prev_offset = offset;
        self.prev_time = time;
    }

    pub fn crossings(&self) -> usize {
        self.count
    }

    pub fn period(&self) -> Option<Real> {
        match (self.first, self.last) {
            (Some(first), Some(last)) if self.count >= 2 => {
                Some((last - first) / (self.count - 1) as Real)
            }
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct EnergySummary {
    pub e_initial: Real,
    pub drift_scale: Real,
    pub max_angular_momentum_drift: Real,
    pub max_abs_drift: Real,
    pub head_mean: Real,
    pub tail_mean: Real,
    pub max_force: Real,
    pub max_speed: Real,
    pub first_nonfinite: Option<usize>,
}

impl EnergySummary {
    pub fn secular_drift(&self) -> Real {
        (self.tail_mean - self.head_mean).abs()
    }

    pub fn is_finite(&self) -> bool {
        self.first_nonfinite.is_none()
    }
}

pub struct EnergyMonitor {
    e_initial: Real,
    drift_scale: Real,
    l_initial: (Real, Real, Real),
    max_l_drift: Real,
    steps: usize,
    window: usize,
    head_sum: Real,
    head_count: usize,
    tail_sum: Real,
    tail_count: usize,
    max_abs_drift: Real,
    max_force: Real,
    max_speed: Real,
    first_nonfinite: Option<usize>,
}

impl EnergyMonitor {
    pub fn new(sys: &System, e_initial: Real, steps: usize) -> Self {
        EnergyMonitor {
            e_initial,
            drift_scale: e_initial.abs().max(1.0),
            l_initial: sys.angular_momentum(),
            max_l_drift: 0.0,
            steps,
            window: (steps / 10).max(1),
            head_sum: 0.0,
            head_count: 0,
            tail_sum: 0.0,
            tail_count: 0,
            max_abs_drift: 0.0,
            max_force: 0.0,
            max_speed: 0.0,
            first_nonfinite: None,
        }
    }

    pub fn update(&mut self, step: usize, sys: &System, e_total: Real) -> Real {
        let force = sys.max_force();
        let speed = sys.max_speed();
        let drift = (e_total - self.e_initial) / self.drift_scale;

        let healthy = e_total.is_finite() && force.is_finite() && speed.is_finite();
        if !healthy {
            if self.first_nonfinite.is_none() {
                self.first_nonfinite = Some(step);
            }
            return drift;
        }

        let (lx, ly, lz) = sys.angular_momentum();
        let dlx = lx - self.l_initial.0;
        let dly = ly - self.l_initial.1;
        let dlz = lz - self.l_initial.2;
        let l_drift = (dlx * dlx + dly * dly + dlz * dlz).sqrt();
        if l_drift.is_finite() && l_drift > self.max_l_drift {
            self.max_l_drift = l_drift;
        }

        if force > self.max_force {
            self.max_force = force;
        }
        if speed > self.max_speed {
            self.max_speed = speed;
        }
        if drift.abs() > self.max_abs_drift {
            self.max_abs_drift = drift.abs();
        }

        if step < self.window {
            self.head_sum += drift;
            self.head_count += 1;
        } else if step + self.window >= self.steps {
            self.tail_sum += drift;
            self.tail_count += 1;
        }

        drift
    }

    pub fn summary(&self) -> EnergySummary {
        EnergySummary {
            e_initial: self.e_initial,
            drift_scale: self.drift_scale,
            max_angular_momentum_drift: self.max_l_drift,
            max_abs_drift: self.max_abs_drift,
            head_mean: mean(self.head_sum, self.head_count),
            tail_mean: mean(self.tail_sum, self.tail_count),
            max_force: self.max_force,
            max_speed: self.max_speed,
            first_nonfinite: self.first_nonfinite,
        }
    }
}

fn mean(sum: Real, count: usize) -> Real {
    if count == 0 {
        Real::NAN
    } else {
        sum / count as Real
    }
}
