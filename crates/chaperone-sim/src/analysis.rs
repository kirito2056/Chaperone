use crate::system::Real;

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
