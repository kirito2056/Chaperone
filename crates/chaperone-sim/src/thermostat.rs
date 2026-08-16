use crate::forcefield::{Energies, ForceField};
use crate::integrator::{drift, kick};
use crate::rng::Noise;
use crate::system::{Real, System};

pub struct Langevin {
    pub gamma: Real,
    pub temperature: Real,
    pub dt: Real,
    noise: Noise,
    step_count: u64,
    c1: Real,
    c2: Real,
}

impl Langevin {
    pub fn new(gamma: Real, temperature: Real, dt: Real, seed: u64) -> Self {
        assert!(gamma >= 0.0, "friction must not be negative");
        assert!(temperature > 0.0, "temperature must be positive");
        assert!(dt > 0.0, "timestep must be positive");

        let c1 = (-gamma * dt).exp();
        Langevin {
            gamma,
            temperature,
            dt,
            noise: Noise::new(seed),
            step_count: 0,
            c1,
            c2: (1.0 - c1 * c1).sqrt(),
        }
    }

    pub fn steps_taken(&self) -> u64 {
        self.step_count
    }

    fn o_step(&self, sys: &mut System) {
        let step = self.step_count;
        for i in 0..sys.n {
            let scale = (self.temperature / sys.mass[i]).sqrt();
            sys.vel_x[i] =
                self.c1 * sys.vel_x[i] + self.c2 * scale * self.noise.gaussian(step, i, 0);
            sys.vel_y[i] =
                self.c1 * sys.vel_y[i] + self.c2 * scale * self.noise.gaussian(step, i, 1);
            sys.vel_z[i] =
                self.c1 * sys.vel_z[i] + self.c2 * scale * self.noise.gaussian(step, i, 2);
        }
    }

    pub fn step(&mut self, sys: &mut System, ff: &ForceField) -> Energies {
        self.step_count += 1;
        let half = 0.5 * self.dt;
        kick(sys, half);
        drift(sys, half);
        self.o_step(sys);
        drift(sys, half);
        sys.clear_forces();
        let mut energies = ff.accumulate(sys);
        kick(sys, half);
        energies.kinetic = sys.kinetic_energy();
        energies
    }
}

pub fn instantaneous_temperature(sys: &System) -> Real {
    2.0 * sys.kinetic_energy() / (3.0 * sys.n as Real)
}

pub fn center_of_mass_velocity(sys: &System) -> (Real, Real, Real) {
    let mut m_total = 0.0;
    let (mut px, mut py, mut pz) = (0.0, 0.0, 0.0);
    for i in 0..sys.n {
        let m = sys.mass[i];
        m_total += m;
        px += m * sys.vel_x[i];
        py += m * sys.vel_y[i];
        pz += m * sys.vel_z[i];
    }
    (px / m_total, py / m_total, pz / m_total)
}
