pub type Real = f64;

pub const PI: Real = std::f64::consts::PI;

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

    pub fn distance(&self, i: usize, j: usize) -> Real {
        let dx = self.pos_x[j] - self.pos_x[i];
        let dy = self.pos_y[j] - self.pos_y[i];
        let dz = self.pos_z[j] - self.pos_z[i];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    pub fn angle(&self, i: usize, j: usize, k: usize) -> Real {
        let ax = self.pos_x[i] - self.pos_x[j];
        let ay = self.pos_y[i] - self.pos_y[j];
        let az = self.pos_z[i] - self.pos_z[j];
        let bx = self.pos_x[k] - self.pos_x[j];
        let by = self.pos_y[k] - self.pos_y[j];
        let bz = self.pos_z[k] - self.pos_z[j];

        let ra = (ax * ax + ay * ay + az * az).sqrt();
        let rb = (bx * bx + by * by + bz * bz).sqrt();
        let cos = ((ax * bx + ay * by + az * bz) / (ra * rb)).clamp(-1.0, 1.0);
        cos.acos()
    }

    pub fn dihedral(&self, i: usize, j: usize, k: usize, l: usize) -> Real {
        let b1 = [
            self.pos_x[j] - self.pos_x[i],
            self.pos_y[j] - self.pos_y[i],
            self.pos_z[j] - self.pos_z[i],
        ];
        let b2 = [
            self.pos_x[k] - self.pos_x[j],
            self.pos_y[k] - self.pos_y[j],
            self.pos_z[k] - self.pos_z[j],
        ];
        let b3 = [
            self.pos_x[l] - self.pos_x[k],
            self.pos_y[l] - self.pos_y[k],
            self.pos_z[l] - self.pos_z[k],
        ];

        let cross = |u: [Real; 3], v: [Real; 3]| {
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        let dot = |u: [Real; 3], v: [Real; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];

        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let b2_len = dot(b2, b2).sqrt();

        (b2_len * dot(b1, n2)).atan2(dot(n1, n2))
    }

    pub fn angular_momentum(&self) -> (Real, Real, Real) {
        let mut lx = 0.0;
        let mut ly = 0.0;
        let mut lz = 0.0;
        for i in 0..self.n {
            let m = self.mass[i];
            lx += m * (self.pos_y[i] * self.vel_z[i] - self.pos_z[i] * self.vel_y[i]);
            ly += m * (self.pos_z[i] * self.vel_x[i] - self.pos_x[i] * self.vel_z[i]);
            lz += m * (self.pos_x[i] * self.vel_y[i] - self.pos_y[i] * self.vel_x[i]);
        }
        (lx, ly, lz)
    }

    pub fn total_torque(&self) -> (Real, Real, Real) {
        let mut tx = 0.0;
        let mut ty = 0.0;
        let mut tz = 0.0;
        for i in 0..self.n {
            tx += self.pos_y[i] * self.frc_z[i] - self.pos_z[i] * self.frc_y[i];
            ty += self.pos_z[i] * self.frc_x[i] - self.pos_x[i] * self.frc_z[i];
            tz += self.pos_x[i] * self.frc_y[i] - self.pos_y[i] * self.frc_x[i];
        }
        (tx, ty, tz)
    }

    pub fn total_force(&self) -> (Real, Real, Real) {
        (
            self.frc_x.iter().sum(),
            self.frc_y.iter().sum(),
            self.frc_z.iter().sum(),
        )
    }

    pub fn max_force(&self) -> Real {
        let mut m: Real = 0.0;
        for i in 0..self.n {
            let f = (self.frc_x[i] * self.frc_x[i]
                + self.frc_y[i] * self.frc_y[i]
                + self.frc_z[i] * self.frc_z[i])
                .sqrt();
            if f.is_finite() && f > m {
                m = f;
            } else if !f.is_finite() {
                return Real::NAN;
            }
        }
        m
    }

    pub fn max_speed(&self) -> Real {
        let mut m: Real = 0.0;
        for i in 0..self.n {
            let v = (self.vel_x[i] * self.vel_x[i]
                + self.vel_y[i] * self.vel_y[i]
                + self.vel_z[i] * self.vel_z[i])
                .sqrt();
            if v.is_finite() && v > m {
                m = v;
            } else if !v.is_finite() {
                return Real::NAN;
            }
        }
        m
    }
}
