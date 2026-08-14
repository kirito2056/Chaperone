use crate::forcefield::angle::Angles;
use crate::forcefield::native::NativeContacts;
use crate::forcefield::pairlist::PairList;
use crate::forcefield::ForceField;
use crate::system::{Real, System};

pub const BOND_K: Real = 100.0;
pub const ANGLE_K: Real = 20.0;
pub const R0: Real = 3.8;
pub const EPS: Real = 1.0;
pub const SIGMA: Real = 4.0;
pub const CONTACT_CUTOFF: Real = 4.0;
pub const MIN_SEQUENCE_SEPARATION: usize = 3;

pub fn spring(initial_separation: Real) -> (System, ForceField) {
    let mut sys = System::new(2);
    sys.pos_x[1] = initial_separation;

    let mut ff = ForceField::new(BOND_K, ANGLE_K, EPS, SIGMA);
    ff.bonds.push(0, 1, R0);

    (sys, ff)
}

pub fn chain4() -> (System, ForceField) {
    chain4_with_gap(R0)
}

pub fn chain4_with_gap(gap: Real) -> (System, ForceField) {
    assert!(
        (gap - R0).abs() < 2.0 * R0,
        "gap {gap} is unreachable with bond length {R0}"
    );
    let mut sys = System::new(4);
    let b = (gap - R0) / 2.0;
    let a = (R0 * R0 - b * b).sqrt();
    let corners = [(0.0, 0.0), (a, b), (a, gap - b), (0.0, gap)];
    for (i, (x, y)) in corners.iter().enumerate() {
        sys.pos_x[i] = *x;
        sys.pos_y[i] = *y;
    }

    let mut ff = ForceField::new(BOND_K, ANGLE_K, EPS, SIGMA);
    for i in 0..3 {
        ff.bonds.push(i as u32, (i + 1) as u32, R0);
    }
    ff.repulsion_pairs = PairList::all_pairs(4, MIN_SEQUENCE_SEPARATION);

    (sys, ff)
}

pub fn native_pair(sigma: Real, initial_separation: Real) -> (System, ForceField) {
    let mut sys = System::new(2);
    sys.pos_x[1] = initial_separation;

    let mut ff = ForceField::new(BOND_K, ANGLE_K, EPS, SIGMA);
    ff.native.push(0, 1, sigma);

    (sys, ff)
}

pub fn chain5() -> (System, ForceField) {
    let mut sys = System::new(5);

    let d04 = 4.2;
    let y4 = d04 * d04 / (2.0 * R0);
    let z4 = (d04 * d04 - y4 * y4).sqrt();

    let coords = [
        (0.0, 0.0, 0.0),
        (R0, 0.0, 0.0),
        (R0, R0, 0.0),
        (0.0, R0, 0.0),
        (0.0, y4, z4),
    ];
    for (i, (x, y, z)) in coords.iter().enumerate() {
        sys.pos_x[i] = *x;
        sys.pos_y[i] = *y;
        sys.pos_z[i] = *z;
    }

    let mut ff = ForceField::new(BOND_K, ANGLE_K, EPS, SIGMA);
    for i in 0..4 {
        ff.bonds.push(i as u32, (i + 1) as u32, R0);
    }

    ff.angles = Angles::from_chain(&sys);
    ff.native = NativeContacts::from_structure(&sys, CONTACT_CUTOFF, MIN_SEQUENCE_SEPARATION);
    ff.repulsion_pairs =
        PairList::all_pairs(5, MIN_SEQUENCE_SEPARATION).exclude(&ff.native.i, &ff.native.j);

    (sys, ff)
}
