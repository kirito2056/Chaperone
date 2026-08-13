pub mod bond;
pub mod pairlist;
pub mod repulsion;

use crate::system::{Real, System};
use bond::Bonds;
use pairlist::PairList;

#[derive(Debug, Clone, Copy, Default)]
pub struct Energies {
    pub bond: Real,
    pub repulsion: Real,
    pub kinetic: Real,
}

impl Energies {
    pub fn potential(&self) -> Real {
        self.bond + self.repulsion
    }

    pub fn total(&self) -> Real {
        self.potential() + self.kinetic
    }
}

pub struct ForceField {
    pub bonds: Bonds,
    pub bond_k: Real,
    pub repulsion_pairs: PairList,
    pub eps: Real,
    pub sigma: Real,
}

impl ForceField {
    pub fn new(bond_k: Real, eps: Real, sigma: Real) -> Self {
        ForceField {
            bonds: Bonds::new(),
            bond_k,
            repulsion_pairs: PairList::new(),
            eps,
            sigma,
        }
    }

    pub fn validate(&self, n: usize) {
        self.bonds.validate(n);
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
            repulsion: repulsion::accumulate(sys, &self.repulsion_pairs, self.eps, self.sigma),
            kinetic: 0.0,
        }
    }

    pub fn potential_energy(&self, sys: &System) -> Real {
        bond::energy(sys, &self.bonds, self.bond_k)
            + repulsion::energy(sys, &self.repulsion_pairs, self.eps, self.sigma)
    }
}
