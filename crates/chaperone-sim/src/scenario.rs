use crate::forcefield::pairlist::PairList;
use crate::forcefield::ForceField;
use crate::system::{Real, System};

pub const BOND_K: Real = 100.0;
pub const R0: Real = 3.8;
pub const EPS: Real = 1.0;
pub const SIGMA: Real = 4.0;

pub fn spring(initial_separation: Real) -> (System, ForceField) {
    let mut sys = System::new(2);
    sys.pos_x[1] = initial_separation;

    let mut ff = ForceField::new(BOND_K, EPS, SIGMA);
    ff.bonds.push(0, 1, R0);

    (sys, ff)
}

pub fn chain4() -> (System, ForceField) {
    let mut sys = System::new(4);
    let corners = [(0.0, 0.0), (R0, 0.0), (R0, R0), (0.0, R0)];
    for (i, (x, y)) in corners.iter().enumerate() {
        sys.pos_x[i] = *x;
        sys.pos_y[i] = *y;
    }

    let mut ff = ForceField::new(BOND_K, EPS, SIGMA);
    for i in 0..3 {
        ff.bonds.push(i as u32, (i + 1) as u32, R0);
    }
    ff.repulsion_pairs = PairList::all_pairs(4, 3);

    (sys, ff)
}
