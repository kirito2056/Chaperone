pub type Real = f64;

pub struct System {
    pub n: usize,
    pub pos_x: Vec<Real>,
    pub pos_y: Vec<Real>,
    pub pos_z: Vec<Real>,
    pub vel_x: Vec<Real>,
    pub vel_y: Vec<Real>,
    pub vel_z: Vec<Real>,
    pub frc_x: Vec<Real>,
    pub frc_y: Vec<Real>,
    pub frc_z: Vec<Real>,
    pub mass: Vec<Real>,
}

impl System {
    pub fn new(n: usize) -> Self {
        System {
            n,
            pos_x: vec![0.0; n],
            pos_y: vec![0.0; n],
            pos_z: vec![0.0; n],
            vel_x: vec![0.0; n],
            vel_y: vec![0.0; n],
            vel_z: vec![0.0; n],
            frc_x: vec![0.0; n],
            frc_y: vec![0.0; n],
            frc_z: vec![0.0; n],
            mass: vec![1.0; n],
        }
    }

    pub fn clear_forces(&mut self) {
        self.frc_x.fill(0.0);
        self.frc_y.fill(0.0);
        self.frc_z.fill(0.0);
    }

    pub fn kinetic_energy(&self) -> Real {
        let mut e = 0.0;
        for i in 0..self.n {
            let v2 = self.vel_x[i] * self.vel_x[i]
                + self.vel_y[i] * self.vel_y[i]
                + self.vel_z[i] * self.vel_z[i];
            e += 0.5 * self.mass[i] * v2;
        }
        e
    }

    pub fn temperature(&self) -> Real {
        if self.n == 0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (3.0 * self.n as Real)
    }

    pub fn distance(&self, i: usize, j: usize) -> Real {
        let dx = self.pos_x[j] - self.pos_x[i];
        let dy = self.pos_y[j] - self.pos_y[i];
        let dz = self.pos_z[j] - self.pos_z[i];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    pub fn total_force(&self) -> (Real, Real, Real) {
        (
            self.frc_x.iter().sum(),
            self.frc_y.iter().sum(),
            self.frc_z.iter().sum(),
        )
    }
}
