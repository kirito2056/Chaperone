use crate::system::{Real, System};

pub const MAX_EXTENSION: Real = 5.0;
const MIN_EXTENSION: Real = 1e-12;

pub struct Pull {
    pub index: Option<u32>,
    pub target: [Real; 3],
    pub k: Real,
}

impl Pull {
    pub fn new(k: Real) -> Self {
        assert!(k > 0.0, "pull stiffness must be positive, got {k}");
        Pull {
            index: None,
            target: [0.0; 3],
            k,
        }
    }

    pub fn is_active(&self) -> bool {
        self.index.is_some()
    }

    pub fn release(&mut self) {
        self.index = None;
    }

    pub fn validate(&self, n: usize) {
        if let Some(index) = self.index {
            assert!(
                (index as usize) < n,
                "pulled bead {index} out of range for system of {n} atoms"
            );
        }
    }

    pub fn extension(&self, sys: &System) -> Real {
        match self.offset(sys) {
            Some((_, _, _, ext)) => ext,
            None => 0.0,
        }
    }

    pub fn force_magnitude(&self, sys: &System) -> Real {
        match self.offset(sys) {
            Some((_, _, _, ext)) => self.k * ext.min(MAX_EXTENSION),
            None => 0.0,
        }
    }

    fn offset(&self, sys: &System) -> Option<(Real, Real, Real, Real)> {
        let i = self.index? as usize;
        if i >= sys.n {
            return None;
        }
        let dx = self.target[0] - sys.pos_x[i];
        let dy = self.target[1] - sys.pos_y[i];
        let dz = self.target[2] - sys.pos_z[i];
        let ext = (dx * dx + dy * dy + dz * dz).sqrt();
        if ext < MIN_EXTENSION {
            return None;
        }
        Some((dx, dy, dz, ext))
    }
}

impl Default for Pull {
    fn default() -> Self {
        Pull {
            index: None,
            target: [0.0; 3],
            k: 1.0,
        }
    }
}

fn potential(k: Real, ext: Real) -> Real {
    if ext <= MAX_EXTENSION {
        0.5 * k * ext * ext
    } else {
        0.5 * k * MAX_EXTENSION * MAX_EXTENSION + k * MAX_EXTENSION * (ext - MAX_EXTENSION)
    }
}

pub fn accumulate(sys: &mut System, pull: &Pull) -> Real {
    let Some((dx, dy, dz, ext)) = pull.offset(sys) else {
        return 0.0;
    };
    let i = pull.index.expect("offset implies an index") as usize;

    let scale = pull.k * ext.min(MAX_EXTENSION) / ext;
    sys.frc_x[i] += dx * scale;
    sys.frc_y[i] += dy * scale;
    sys.frc_z[i] += dz * scale;

    potential(pull.k, ext)
}

pub fn energy(sys: &System, pull: &Pull) -> Real {
    match pull.offset(sys) {
        Some((_, _, _, ext)) => potential(pull.k, ext),
        None => 0.0,
    }
}
