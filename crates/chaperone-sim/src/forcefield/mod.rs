pub mod angle;
pub mod bond;
pub mod dihedral;
pub mod native;
pub mod pairlist;
pub mod repulsion;

use crate::system::{Real, System};
use angle::Angles;
use bond::Bonds;
use dihedral::Dihedrals;
use native::NativeContacts;
use pairlist::PairList;

#[derive(Debug, Clone, Copy, Default)]
pub struct Energies {
    pub bond: Real,
    pub angle: Real,
    pub dihedral: Real,
    pub native: Real,
    pub repulsion: Real,
    pub kinetic: Real,
}

impl Energies {
    pub fn potential(&self) -> Real {
        self.bond + self.angle + self.dihedral + self.native + self.repulsion
    }

    pub fn total(&self) -> Real {
        self.potential() + self.kinetic
    }
}

pub struct ForceField {
    pub bonds: Bonds,
    pub bond_k: Real,
    pub angles: Angles,
    pub angle_k: Real,
    pub dihedrals: Dihedrals,
    pub k_phi1: Real,
    pub k_phi3: Real,
    pub native: NativeContacts,
    pub repulsion_pairs: PairList,
    pub eps: Real,
    pub sigma: Real,
}

impl ForceField {
    pub fn new(
        bond_k: Real,
        angle_k: Real,
        k_phi1: Real,
        k_phi3: Real,
        eps: Real,
        sigma: Real,
    ) -> Self {
        ForceField {
            bonds: Bonds::new(),
            bond_k,
            angles: Angles::new(),
            angle_k,
            dihedrals: Dihedrals::new(),
            k_phi1,
            k_phi3,
            native: NativeContacts::new(),
            repulsion_pairs: PairList::new(),
            eps,
            sigma,
        }
    }

    pub fn validate(&self, n: usize) {
        self.bonds.validate(n);
        self.angles.validate(n);
        self.dihedrals.validate(n);
        self.native.validate(n);
        self.repulsion_pairs.validate(n);
        assert!(
            self.sigma > 0.0,
            "sigma must be positive, got {}",
            self.sigma
        );
    }

    pub fn accumulate(&self, sys: &mut System) -> Energies {
        Energies {
            bond: bond::accumulate(sys, &self.bonds, self.bond_k),
            angle: angle::accumulate(sys, &self.angles, self.angle_k),
            dihedral: dihedral::accumulate(sys, &self.dihedrals, self.k_phi1, self.k_phi3),
            native: native::accumulate(sys, &self.native, self.eps),
            repulsion: repulsion::accumulate(sys, &self.repulsion_pairs, self.eps, self.sigma),
            kinetic: 0.0,
        }
    }

    pub fn potential_energy(&self, sys: &System) -> Real {
        bond::energy(sys, &self.bonds, self.bond_k)
            + angle::energy(sys, &self.angles, self.angle_k)
            + dihedral::energy(sys, &self.dihedrals, self.k_phi1, self.k_phi3)
            + native::energy(sys, &self.native, self.eps)
            + repulsion::energy(sys, &self.repulsion_pairs, self.eps, self.sigma)
    }
}
