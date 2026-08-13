use crate::forcefield::{Energies, ForceField};
use crate::system::{Real, System};

pub fn kick_drift(sys: &mut System, dt: Real) {
    let half = 0.5 * dt;
    for i in 0..sys.n {
        let inv_m = 1.0 / sys.mass[i];
        sys.vel_x[i] += half * sys.frc_x[i] * inv_m;
        sys.vel_y[i] += half * sys.frc_y[i] * inv_m;
        sys.vel_z[i] += half * sys.frc_z[i] * inv_m;
        sys.pos_x[i] += dt * sys.vel_x[i];
        sys.pos_y[i] += dt * sys.vel_y[i];
        sys.pos_z[i] += dt * sys.vel_z[i];
    }
}

pub fn kick(sys: &mut System, dt: Real) {
    let half = 0.5 * dt;
    for i in 0..sys.n {
        let inv_m = 1.0 / sys.mass[i];
        sys.vel_x[i] += half * sys.frc_x[i] * inv_m;
        sys.vel_y[i] += half * sys.frc_y[i] * inv_m;
        sys.vel_z[i] += half * sys.frc_z[i] * inv_m;
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
    kick_drift(sys, dt);
    sys.clear_forces();
    let mut energies = ff.accumulate(sys);
    kick(sys, dt);
    energies.kinetic = sys.kinetic_energy();
    energies
}
