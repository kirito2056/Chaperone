pub mod bond;
pub mod native;
pub mod pairlist;
pub mod repulsion;

use crate::system::{Real, System};
use bond::Bonds;
use native::NativeContacts;
use pairlist::PairList;

#[derive(Debug, Clone, Copy, Default)]
pub struct Energies {
    pub bond: Real,
    pub native: Real,
    pub repulsion: Real,
    pub kinetic: Real,
}

impl Energies {
    pub fn potential(&self) -> Real {
        self.bond + self.native + self.repulsion
    }

    pub fn total(&self) -> Real {
        self.potential() + self.kinetic
    }
}

pub struct ForceField {
    pub bonds: Bonds,
    pub bond_k: Real,
    pub native: NativeContacts,
    pub repulsion_pairs: PairList,
    pub eps: Real,
    pub sigma: Real,
}

impl ForceField {
    pub fn new(bond_k: Real, eps: Real, sigma: Real) -> Self {
        ForceField {
            bonds: Bonds::new(),
            bond_k,
            native: NativeContacts::new(),
            repulsion_pairs: PairList::new(),
            eps,
            sigma,
        }
    }

    pub fn validate(&self, n: usize) {
        self.bonds.validate(n);
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
            native: native::accumulate(sys, &self.native, self.eps),
            repulsion: repulsion::accumulate(sys, &self.repulsion_pairs, self.eps, self.sigma),
            kinetic: 0.0,
        }
    }

    pub fn potential_energy(&self, sys: &System) -> Real {
        bond::energy(sys, &self.bonds, self.bond_k)
            + native::energy(sys, &self.native, self.eps)
            + repulsion::energy(sys, &self.repulsion_pairs, self.eps, self.sigma)
    }
}
