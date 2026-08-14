use crate::system::{Real, System};

const MIN_LEN2: Real = 1e-24;
const MIN_SIN: Real = 1e-8;

pub struct Angles {
    pub i: Vec<u32>,
    pub j: Vec<u32>,
    pub k: Vec<u32>,
    pub theta0: Vec<Real>,
}

impl Angles {
    pub fn new() -> Self {
        Angles {
            i: Vec::new(),
            j: Vec::new(),
            k: Vec::new(),
            theta0: Vec::new(),
        }
    }

    pub fn push(&mut self, i: u32, j: u32, k: u32, theta0: Real) {
        assert!(
            i != j && j != k && i != k,
            "degenerate triple ({i}, {j}, {k})"
        );
        assert!(
            theta0 > 0.0 && theta0 < std::f64::consts::PI,
            "theta0 = {theta0} is outside (0, PI)"
        );
        self.i.push(i);
        self.j.push(j);
        self.k.push(k);
        self.theta0.push(theta0);
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    pub fn validate(&self, n: usize) {
        for t in 0..self.len() {
            let (i, j, k) = (self.i[t] as usize, self.j[t] as usize, self.k[t] as usize);
            assert!(
                i < n && j < n && k < n,
                "triple ({i}, {j}, {k}) out of range for system of {n} atoms"
            );
        }
    }

    pub fn from_chain(native: &System) -> Self {
        let mut angles = Angles::new();
        for t in 0..native.n.saturating_sub(2) {
            let theta0 = native.angle(t, t + 1, t + 2);
            angles.push(t as u32, (t + 1) as u32, (t + 2) as u32, theta0);
        }
        angles
    }
}

impl Default for Angles {
    fn default() -> Self {
        Self::new()
    }
}

pub fn accumulate(sys: &mut System, angles: &Angles, k_theta: Real) -> Real {
    let mut e_pot = 0.0;

    for t in 0..angles.len() {
        let i = angles.i[t] as usize;
        let j = angles.j[t] as usize;
        let k = angles.k[t] as usize;

        let ax = sys.pos_x[i] - sys.pos_x[j];
        let ay = sys.pos_y[i] - sys.pos_y[j];
        let az = sys.pos_z[i] - sys.pos_z[j];
        let bx = sys.pos_x[k] - sys.pos_x[j];
        let by = sys.pos_y[k] - sys.pos_y[j];
        let bz = sys.pos_z[k] - sys.pos_z[j];

        let ra2 = ax * ax + ay * ay + az * az;
        let rb2 = bx * bx + by * by + bz * bz;
        if ra2 < MIN_LEN2 || rb2 < MIN_LEN2 {
            continue;
        }

        let ra = ra2.sqrt();
        let rb = rb2.sqrt();
        let cos = ((ax * bx + ay * by + az * bz) / (ra * rb)).clamp(-1.0, 1.0);
        let sin = (1.0 - cos * cos).sqrt();
        if sin < MIN_SIN {
            continue;
        }

        let dtheta = cos.acos() - angles.theta0[t];
        e_pot += k_theta * dtheta * dtheta;

        let pref = 2.0 * k_theta * dtheta / sin;
        let inv_rarb = 1.0 / (ra * rb);
        let ca = cos / ra2;
        let cb = cos / rb2;

        let fix = pref * (bx * inv_rarb - ca * ax);
        let fiy = pref * (by * inv_rarb - ca * ay);
        let fiz = pref * (bz * inv_rarb - ca * az);
        let fkx = pref * (ax * inv_rarb - cb * bx);
        let fky = pref * (ay * inv_rarb - cb * by);
        let fkz = pref * (az * inv_rarb - cb * bz);

        sys.frc_x[i] += fix;
        sys.frc_y[i] += fiy;
        sys.frc_z[i] += fiz;
        sys.frc_x[k] += fkx;
        sys.frc_y[k] += fky;
        sys.frc_z[k] += fkz;
        sys.frc_x[j] -= fix + fkx;
        sys.frc_y[j] -= fiy + fky;
        sys.frc_z[j] -= fiz + fkz;
    }

    e_pot
}

pub fn energy(sys: &System, angles: &Angles, k_theta: Real) -> Real {
    let mut e_pot = 0.0;

    for t in 0..angles.len() {
        let i = angles.i[t] as usize;
        let j = angles.j[t] as usize;
        let k = angles.k[t] as usize;

        let ax = sys.pos_x[i] - sys.pos_x[j];
        let ay = sys.pos_y[i] - sys.pos_y[j];
        let az = sys.pos_z[i] - sys.pos_z[j];
        let bx = sys.pos_x[k] - sys.pos_x[j];
        let by = sys.pos_y[k] - sys.pos_y[j];
        let bz = sys.pos_z[k] - sys.pos_z[j];

        let ra2 = ax * ax + ay * ay + az * az;
        let rb2 = bx * bx + by * by + bz * bz;
        if ra2 < MIN_LEN2 || rb2 < MIN_LEN2 {
            continue;
        }

        let ra = ra2.sqrt();
        let rb = rb2.sqrt();
        let cos = ((ax * bx + ay * by + az * bz) / (ra * rb)).clamp(-1.0, 1.0);
        let sin = (1.0 - cos * cos).sqrt();
        if sin < MIN_SIN {
            continue;
        }

        let dtheta = cos.acos() - angles.theta0[t];
        e_pot += k_theta * dtheta * dtheta;
    }

    e_pot
}
