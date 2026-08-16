use crate::analysis::{
    fraction_of_local_contacts, fraction_of_native_contacts, fraction_of_tertiary_contacts,
    CONTACT_TOLERANCE,
};
use crate::forcefield::native::NativeContacts;
use crate::forcefield::{Energies, ForceField};
use crate::system::{Real, System};
use crate::thermostat::{instantaneous_temperature, Langevin};

#[derive(Debug, Clone, Copy)]
pub struct FoldingSummary {
    pub q_initial: Real,
    pub q_min: Real,
    pub q_max: Real,
    pub q_final: Real,
    pub q_local_final: Real,
    pub q_tertiary_final: Real,
    pub rg_initial: Real,
    pub rg_min: Real,
    pub rg_max: Real,
    pub rg_final: Real,
    pub temperature_mean: Real,
    pub max_force: Real,
    pub max_speed: Real,
    pub first_nonfinite: Option<usize>,
}

impl FoldingSummary {
    pub fn is_finite(&self) -> bool {
        self.first_nonfinite.is_none()
    }
}

pub struct FoldingMonitor {
    tail_start: usize,
    q_initial: Real,
    q_min: Real,
    q_max: Real,
    q_last: Real,
    q_local_last: Real,
    q_tertiary_last: Real,
    rg_initial: Real,
    rg_min: Real,
    rg_max: Real,
    rg_last: Real,
    temperature_sum: Real,
    temperature_count: usize,
    max_force: Real,
    max_speed: Real,
    first_nonfinite: Option<usize>,
}

fn observables(sys: &System, contacts: &NativeContacts) -> (Real, Real, Real, Real) {
    (
        fraction_of_native_contacts(sys, contacts, CONTACT_TOLERANCE).unwrap_or(Real::NAN),
        fraction_of_local_contacts(sys, contacts, CONTACT_TOLERANCE).unwrap_or(Real::NAN),
        fraction_of_tertiary_contacts(sys, contacts, CONTACT_TOLERANCE).unwrap_or(Real::NAN),
        sys.radius_of_gyration(),
    )
}

impl FoldingMonitor {
    pub fn new(sys: &System, contacts: &NativeContacts, steps: usize) -> Self {
        let (q, q_local, q_tertiary, rg) = observables(sys, contacts);
        FoldingMonitor {
            tail_start: steps / 2,
            q_initial: q,
            q_min: q,
            q_max: q,
            q_last: q,
            q_local_last: q_local,
            q_tertiary_last: q_tertiary,
            rg_initial: rg,
            rg_min: rg,
            rg_max: rg,
            rg_last: rg,
            temperature_sum: 0.0,
            temperature_count: 0,
            max_force: 0.0,
            max_speed: 0.0,
            first_nonfinite: None,
        }
    }

    pub fn update(&mut self, step: usize, sys: &System, contacts: &NativeContacts) {
        let force = sys.max_force();
        let speed = sys.max_speed();
        let (q, q_local, q_tertiary, rg) = observables(sys, contacts);

        if !(force.is_finite() && speed.is_finite() && q.is_finite() && rg.is_finite()) {
            if self.first_nonfinite.is_none() {
                self.first_nonfinite = Some(step);
            }
            return;
        }

        if force > self.max_force {
            self.max_force = force;
        }
        if speed > self.max_speed {
            self.max_speed = speed;
        }

        self.q_last = q;
        self.q_local_last = q_local;
        self.q_tertiary_last = q_tertiary;
        if q < self.q_min {
            self.q_min = q;
        }
        if q > self.q_max {
            self.q_max = q;
        }

        self.rg_last = rg;
        if rg < self.rg_min {
            self.rg_min = rg;
        }
        if rg > self.rg_max {
            self.rg_max = rg;
        }

        if step >= self.tail_start {
            self.temperature_sum += instantaneous_temperature(sys);
            self.temperature_count += 1;
        }
    }

    pub fn summary(&self) -> FoldingSummary {
        FoldingSummary {
            q_initial: self.q_initial,
            q_min: self.q_min,
            q_max: self.q_max,
            q_final: self.q_last,
            q_local_final: self.q_local_last,
            q_tertiary_final: self.q_tertiary_last,
            rg_initial: self.rg_initial,
            rg_min: self.rg_min,
            rg_max: self.rg_max,
            rg_final: self.rg_last,
            temperature_mean: if self.temperature_count == 0 {
                Real::NAN
            } else {
                self.temperature_sum / self.temperature_count as Real
            },
            max_force: self.max_force,
            max_speed: self.max_speed,
            first_nonfinite: self.first_nonfinite,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Equipartition {
    bond_sum: Real,
    bond_count: usize,
    angle_sum: Real,
    angle_count: usize,
}

impl Equipartition {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn update(&mut self, sys: &System, ff: &ForceField) {
        for b in 0..ff.bonds.len() {
            let i = ff.bonds.i[b] as usize;
            let j = ff.bonds.j[b] as usize;
            let dr = sys.distance(i, j) - ff.bonds.r0[b];
            self.bond_sum += dr * dr;
            self.bond_count += 1;
        }
        for t in 0..ff.angles.len() {
            let i = ff.angles.i[t] as usize;
            let j = ff.angles.j[t] as usize;
            let k = ff.angles.k[t] as usize;
            let dtheta = sys.angle(i, j, k) - ff.angles.theta0[t];
            self.angle_sum += dtheta * dtheta;
            self.angle_count += 1;
        }
    }

    pub fn bond_mean_square(&self) -> Real {
        if self.bond_count == 0 {
            Real::NAN
        } else {
            self.bond_sum / self.bond_count as Real
        }
    }

    pub fn angle_mean_square(&self) -> Real {
        if self.angle_count == 0 {
            Real::NAN
        } else {
            self.angle_sum / self.angle_count as Real
        }
    }

    pub fn expected_bond_mean_square(temperature: Real, bond_k: Real) -> Real {
        temperature / (2.0 * bond_k)
    }

    pub fn expected_angle_mean_square(temperature: Real, angle_k: Real) -> Real {
        temperature / (2.0 * angle_k)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Stage {
    pub name: &'static str,
    pub temperature: Real,
    pub gamma: Real,
    pub steps: usize,
    pub seed: u64,
}

pub struct Sample<'a> {
    pub stage: &'a Stage,
    pub step: usize,
    pub sys: &'a System,
    pub energies: Energies,
    pub q: Real,
    pub q_local: Real,
    pub q_tertiary: Real,
    pub rg: Real,
    pub temperature: Real,
}

pub fn run_stage<F>(
    sys: &mut System,
    ff: &ForceField,
    contacts: &NativeContacts,
    stage: &Stage,
    dt: Real,
    observe_every: usize,
    mut on_sample: F,
) -> FoldingSummary
where
    F: FnMut(&Sample),
{
    let observe_every = observe_every.max(1);
    let mut bath = Langevin::new(stage.gamma, stage.temperature, dt, stage.seed);
    let mut monitor = FoldingMonitor::new(sys, contacts, stage.steps);

    crate::integrator::initialize(sys, ff);

    for step in 0..stage.steps {
        let energies = bath.step(sys, ff);

        if (step + 1) % observe_every == 0 || step + 1 == stage.steps {
            monitor.update(step, sys, contacts);
            let (q, q_local, q_tertiary, rg) = observables(sys, contacts);
            on_sample(&Sample {
                stage,
                step: step + 1,
                sys,
                energies,
                q,
                q_local,
                q_tertiary,
                rg,
                temperature: instantaneous_temperature(sys),
            });
        }
    }

    monitor.summary()
}
