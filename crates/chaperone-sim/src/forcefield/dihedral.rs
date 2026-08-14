use crate::system::{Real, System};

const MIN_SIN: Real = 1e-8;
const MIN_BUILD_SIN: Real = 1e-3;
const MIN_LEN2: Real = 1e-24;

pub struct Dihedrals {
    pub i: Vec<u32>,
    pub j: Vec<u32>,
    pub k: Vec<u32>,
    pub l: Vec<u32>,
    pub phi0: Vec<Real>,
}

impl Dihedrals {
    pub fn new() -> Self {
        Dihedrals {
            i: Vec::new(),
            j: Vec::new(),
            k: Vec::new(),
            l: Vec::new(),
            phi0: Vec::new(),
        }
    }

    pub fn push(&mut self, i: u32, j: u32, k: u32, l: u32, phi0: Real) {
        let idx = [i, j, k, l];
        for a in 0..4 {
            for b in (a + 1)..4 {
                assert!(
                    idx[a] != idx[b],
                    "degenerate quadruple ({i}, {j}, {k}, {l})"
                );
            }
        }
        assert!(phi0.is_finite(), "phi0 = {phi0} is not finite");
        self.i.push(i);
        self.j.push(j);
        self.k.push(k);
        self.l.push(l);
        self.phi0.push(phi0);
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    pub fn validate(&self, n: usize) {
        for t in 0..self.len() {
            let idx = [
                self.i[t] as usize,
                self.j[t] as usize,
                self.k[t] as usize,
                self.l[t] as usize,
            ];
            for &a in &idx {
                assert!(
                    a < n,
                    "quadruple {idx:?} out of range for system of {n} atoms"
                );
            }
        }
    }

    pub fn from_chain(native: &System) -> Self {
        let mut dihedrals = Dihedrals::new();
        for t in 0..native.n.saturating_sub(3) {
            for (a, b, c) in [(t, t + 1, t + 2), (t + 1, t + 2, t + 3)] {
                let sin = native.angle(a, b, c).sin();
                assert!(
                    sin >= MIN_BUILD_SIN,
                    "triple ({a}, {b}, {c}) is nearly collinear (sin = {sin:.3e}); phi0 for \
                     the quadruple at {t} would be ill-conditioned and the term would be \
                     silently skipped at runtime"
                );
            }
            let phi0 = native.dihedral(t, t + 1, t + 2, t + 3);
            dihedrals.push(
                t as u32,
                (t + 1) as u32,
                (t + 2) as u32,
                (t + 3) as u32,
                phi0,
            );
        }
        dihedrals
    }
}

impl Default for Dihedrals {
    fn default() -> Self {
        Self::new()
    }
}

struct Geometry {
    b1: [Real; 3],
    b2: [Real; 3],
    b3: [Real; 3],
    n1: [Real; 3],
    n2: [Real; 3],
    n1_sq: Real,
    n2_sq: Real,
    b2_len: Real,
    b2_sq: Real,
    phi: Real,
}

fn cross(u: [Real; 3], v: [Real; 3]) -> [Real; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn dot(u: [Real; 3], v: [Real; 3]) -> Real {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

fn geometry(sys: &System, i: usize, j: usize, k: usize, l: usize) -> Option<Geometry> {
    let b1 = [
        sys.pos_x[j] - sys.pos_x[i],
        sys.pos_y[j] - sys.pos_y[i],
        sys.pos_z[j] - sys.pos_z[i],
    ];
    let b2 = [
        sys.pos_x[k] - sys.pos_x[j],
        sys.pos_y[k] - sys.pos_y[j],
        sys.pos_z[k] - sys.pos_z[j],
    ];
    let b3 = [
        sys.pos_x[l] - sys.pos_x[k],
        sys.pos_y[l] - sys.pos_y[k],
        sys.pos_z[l] - sys.pos_z[k],
    ];

    let b1_sq = dot(b1, b1);
    let b2_sq = dot(b2, b2);
    let b3_sq = dot(b3, b3);
    if b1_sq < MIN_LEN2 || b2_sq < MIN_LEN2 || b3_sq < MIN_LEN2 {
        return None;
    }

    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let n1_sq = dot(n1, n1);
    let n2_sq = dot(n2, n2);

    let floor1 = MIN_SIN * MIN_SIN * b1_sq * b2_sq;
    let floor2 = MIN_SIN * MIN_SIN * b2_sq * b3_sq;
    if n1_sq < floor1 || n2_sq < floor2 {
        return None;
    }

    let b2_len = b2_sq.sqrt();
    let phi = (b2_len * dot(b1, n2)).atan2(dot(n1, n2));

    Some(Geometry {
        b1,
        b2,
        b3,
        n1,
        n2,
        n1_sq,
        n2_sq,
        b2_len,
        b2_sq,
        phi,
    })
}

fn derivative(dphi: Real, k1: Real, k3: Real) -> Real {
    k1 * dphi.sin() + 3.0 * k3 * (3.0 * dphi).sin()
}

fn potential(dphi: Real, k1: Real, k3: Real) -> Real {
    k1 * (1.0 - dphi.cos()) + k3 * (1.0 - (3.0 * dphi).cos())
}

pub fn accumulate(sys: &mut System, dihedrals: &Dihedrals, k1: Real, k3: Real) -> Real {
    let mut e_pot = 0.0;

    for t in 0..dihedrals.len() {
        let i = dihedrals.i[t] as usize;
        let j = dihedrals.j[t] as usize;
        let k = dihedrals.k[t] as usize;
        let l = dihedrals.l[t] as usize;

        let Some(g) = geometry(sys, i, j, k, l) else {
            continue;
        };

        let dphi = g.phi - dihedrals.phi0[t];
        e_pot += potential(dphi, k1, k3);

        let dv = derivative(dphi, k1, k3);
        let si = dv * g.b2_len / g.n1_sq;
        let sl = -dv * g.b2_len / g.n2_sq;

        let fi = [g.n1[0] * si, g.n1[1] * si, g.n1[2] * si];
        let fl = [g.n2[0] * sl, g.n2[1] * sl, g.n2[2] * sl];

        let m = dot(g.b1, g.b2) / g.b2_sq;
        let n = dot(g.b3, g.b2) / g.b2_sq;

        let fj = [
            -fi[0] - m * fi[0] + n * fl[0],
            -fi[1] - m * fi[1] + n * fl[1],
            -fi[2] - m * fi[2] + n * fl[2],
        ];
        let fk = [
            -fl[0] + m * fi[0] - n * fl[0],
            -fl[1] + m * fi[1] - n * fl[1],
            -fl[2] + m * fi[2] - n * fl[2],
        ];

        sys.frc_x[i] += fi[0];
        sys.frc_y[i] += fi[1];
        sys.frc_z[i] += fi[2];
        sys.frc_x[j] += fj[0];
        sys.frc_y[j] += fj[1];
        sys.frc_z[j] += fj[2];
        sys.frc_x[k] += fk[0];
        sys.frc_y[k] += fk[1];
        sys.frc_z[k] += fk[2];
        sys.frc_x[l] += fl[0];
        sys.frc_y[l] += fl[1];
        sys.frc_z[l] += fl[2];
    }

    e_pot
}

pub fn energy(sys: &System, dihedrals: &Dihedrals, k1: Real, k3: Real) -> Real {
    let mut e_pot = 0.0;

    for t in 0..dihedrals.len() {
        let i = dihedrals.i[t] as usize;
        let j = dihedrals.j[t] as usize;
        let k = dihedrals.k[t] as usize;
        let l = dihedrals.l[t] as usize;

        if let Some(g) = geometry(sys, i, j, k, l) {
            e_pot += potential(g.phi - dihedrals.phi0[t], k1, k3);
        }
    }

    e_pot
}
