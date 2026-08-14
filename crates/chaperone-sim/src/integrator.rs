use crate::forcefield::{Energies, ForceField};
use crate::system::{Real, System};

pub fn kick(sys: &mut System, h: Real) {
    for i in 0..sys.n {
        let inv_m = 1.0 / sys.mass[i];
        sys.vel_x[i] += h * sys.frc_x[i] * inv_m;
        sys.vel_y[i] += h * sys.frc_y[i] * inv_m;
        sys.vel_z[i] += h * sys.frc_z[i] * inv_m;
    }
}

pub fn drift(sys: &mut System, h: Real) {
    for i in 0..sys.n {
        sys.pos_x[i] += h * sys.vel_x[i];
        sys.pos_y[i] += h * sys.vel_y[i];
        sys.pos_z[i] += h * sys.vel_z[i];
    }
}

pub fn initialize(sys: &mut System, ff: &ForceField) -> Energies {
    ff.validate(sys.n);
    sys.clear_forces();
    let mut energies = ff.accumulate(sys);
    energies.kinetic = sys.kinetic_energy();
    energies
}

pub fn step(sys: &mut System, ff: &ForceField, dt: Real) -> Energies {
    let half = 0.5 * dt;
    kick(sys, half);
    drift(sys, dt);
    sys.clear_forces();
    let mut energies = ff.accumulate(sys);
    kick(sys, half);
    energies.kinetic = sys.kinetic_energy();
    energies
}
