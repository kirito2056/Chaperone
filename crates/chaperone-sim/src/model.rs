use chaperone_pdb::Structure;

use crate::forcefield::angle::Angles;
use crate::forcefield::dihedral::Dihedrals;
use crate::forcefield::native::NativeContacts;
use crate::forcefield::pairlist::PairList;
use crate::forcefield::ForceField;
use crate::system::{Real, System};

pub const CONTACT_CUTOFF: Real = 4.5;
pub const MIN_SEQUENCE_SEPARATION: usize = 3;
pub const PLAXCO_CUTOFF: Real = 6.0;

pub fn system_from_structure(structure: &Structure) -> System {
    let mut sys = System::new(structure.len());
    for (i, residue) in structure.residues.iter().enumerate() {
        sys.pos_x[i] = residue.ca[0];
        sys.pos_y[i] = residue.ca[1];
        sys.pos_z[i] = residue.ca[2];
    }
    sys
}

pub fn go_model(
    structure: &Structure,
    bond_k: Real,
    angle_k: Real,
    k_phi1: Real,
    k_phi3: Real,
    eps: Real,
    sigma_nn: Real,
) -> (System, ForceField) {
    let sys = system_from_structure(structure);
    let n = sys.n;

    let mut ff = ForceField::new(bond_k, angle_k, k_phi1, k_phi3, eps, sigma_nn);

    for i in 1..n {
        ff.bonds
            .push((i - 1) as u32, i as u32, sys.distance(i - 1, i));
    }
    ff.angles = Angles::from_chain(&sys);
    ff.dihedrals = Dihedrals::from_chain(&sys);

    let mut native = NativeContacts::new();
    for (i, j) in structure.heavy_atom_contacts(CONTACT_CUTOFF, MIN_SEQUENCE_SEPARATION) {
        native.push(i as u32, j as u32, sys.distance(i, j));
    }
    ff.repulsion_pairs =
        PairList::all_pairs(n, MIN_SEQUENCE_SEPARATION).exclude(&native.i, &native.j);
    ff.native = native;

    (sys, ff)
}
