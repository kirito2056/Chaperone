use crate::system::{Real, System};

const MIN_R: Real = 1e-12;

pub struct Bonds {
    pub i: Vec<u32>,
    pub j: Vec<u32>,
    pub r0: Vec<Real>,
}

impl Bonds {
    pub fn new() -> Self {
        Bonds {
            i: Vec::new(),
            j: Vec::new(),
            r0: Vec::new(),
        }
    }

    pub fn push(&mut self, i: u32, j: u32, r0: Real) {
        assert!(i != j, "self-bond ({i}, {j}) is not a valid interaction");
        assert!(r0 > 0.0, "bond ({i}, {j}) has non-positive r0 = {r0}");
        self.i.push(i);
        self.j.push(j);
        self.r0.push(r0);
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    pub fn validate(&self, n: usize) {
        for b in 0..self.len() {
            let (i, j) = (self.i[b] as usize, self.j[b] as usize);
            assert!(
                i < n && j < n,
                "bond ({i}, {j}) out of range for system of {n} atoms"
            );
        }
    }
}

impl Default for Bonds {
    fn default() -> Self {
        Self::new()
    }
}

pub fn accumulate(sys: &mut System, bonds: &Bonds, k: Real) -> Real {
    let mut e_pot = 0.0;

    for b in 0..bonds.len() {
        let i = bonds.i[b] as usize;
        let j = bonds.j[b] as usize;
        let r0 = bonds.r0[b];

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        if r < MIN_R {
            e_pot += k * r0 * r0;
            continue;
        }

        let dr = r - r0;
        e_pot += k * dr * dr;

        let scale = 2.0 * k * dr / r;
        sys.frc_x[i] += dx * scale;
        sys.frc_y[i] += dy * scale;
        sys.frc_z[i] += dz * scale;
        sys.frc_x[j] -= dx * scale;
        sys.frc_y[j] -= dy * scale;
        sys.frc_z[j] -= dz * scale;
    }

    e_pot
}

pub fn energy(sys: &System, bonds: &Bonds, k: Real) -> Real {
    let mut e_pot = 0.0;

    for b in 0..bonds.len() {
        let i = bonds.i[b] as usize;
        let j = bonds.j[b] as usize;
        let r0 = bonds.r0[b];

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        if r < MIN_R {
            e_pot += k * r0 * r0;
            continue;
        }

        let dr = r - r0;
        e_pot += k * dr * dr;
    }

    e_pot
}
