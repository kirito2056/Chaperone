use crate::forcefield::pairlist::PairList;
use crate::system::{Real, System};

const MIN_R2: Real = 1e-24;

pub struct NativeContacts {
    pub i: Vec<u32>,
    pub j: Vec<u32>,
    pub sigma: Vec<Real>,
}

impl NativeContacts {
    pub fn new() -> Self {
        NativeContacts {
            i: Vec::new(),
            j: Vec::new(),
            sigma: Vec::new(),
        }
    }

    pub fn push(&mut self, i: u32, j: u32, sigma: Real) {
        assert!(i != j, "self-contact ({i}, {j}) is not a valid interaction");
        assert!(
            sigma > 0.0,
            "contact ({i}, {j}) has non-positive sigma = {sigma}"
        );
        self.i.push(i);
        self.j.push(j);
        self.sigma.push(sigma);
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    pub fn validate(&self, n: usize) {
        for c in 0..self.len() {
            let (i, j) = (self.i[c] as usize, self.j[c] as usize);
            assert!(
                i < n && j < n,
                "contact ({i}, {j}) out of range for system of {n} atoms"
            );
        }
    }

    pub fn from_structure(native: &System, cutoff: Real, min_sep: usize) -> Self {
        let mut contacts = NativeContacts::new();
        let candidates = PairList::all_pairs(native.n, min_sep);

        for p in 0..candidates.len() {
            let (i, j) = (candidates.i[p], candidates.j[p]);
            let r = native.distance(i as usize, j as usize);
            if r < cutoff {
                contacts.push(i, j, r);
            }
        }

        contacts
    }
}

impl Default for NativeContacts {
    fn default() -> Self {
        Self::new()
    }
}

pub fn accumulate(sys: &mut System, contacts: &NativeContacts, eps: Real) -> Real {
    let mut e_pot = 0.0;

    for c in 0..contacts.len() {
        let i = contacts.i[c] as usize;
        let j = contacts.j[c] as usize;
        let sigma = contacts.sigma[c];

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r2 = dx * dx + dy * dy + dz * dz;

        if r2 < MIN_R2 {
            continue;
        }

        let s2 = sigma * sigma / r2;
        let s4 = s2 * s2;
        let s8 = s4 * s4;
        let s10 = s8 * s2;
        let s12 = s10 * s2;

        e_pot += eps * (5.0 * s12 - 6.0 * s10);

        let scale = 60.0 * eps * (s10 - s12) / r2;
        sys.frc_x[i] += dx * scale;
        sys.frc_y[i] += dy * scale;
        sys.frc_z[i] += dz * scale;
        sys.frc_x[j] -= dx * scale;
        sys.frc_y[j] -= dy * scale;
        sys.frc_z[j] -= dz * scale;
    }

    e_pot
}

pub fn energy(sys: &System, contacts: &NativeContacts, eps: Real) -> Real {
    let mut e_pot = 0.0;

    for c in 0..contacts.len() {
        let i = contacts.i[c] as usize;
        let j = contacts.j[c] as usize;
        let sigma = contacts.sigma[c];

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r2 = dx * dx + dy * dy + dz * dz;

        if r2 < MIN_R2 {
            continue;
        }

        let s2 = sigma * sigma / r2;
        let s4 = s2 * s2;
        let s8 = s4 * s4;
        let s10 = s8 * s2;
        let s12 = s10 * s2;

        e_pot += eps * (5.0 * s12 - 6.0 * s10);
    }

    e_pot
}
