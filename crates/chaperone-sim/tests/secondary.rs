use chaperone_sim::secondary::{window, Assignment, Ss};
use chaperone_sim::system::{Real, System, PI};

mod common;
use common::UBIQUITIN_CA;

const DEG: Real = PI / 180.0;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> Real {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as Real) / ((1u64 << 53) as Real)
    }
    fn range(&mut self, lo: Real, hi: Real) -> Real {
        lo + (hi - lo) * self.next()
    }
}

fn build(thetas: &[Real], phis: &[Real]) -> System {
    let b = 3.8;
    let n = thetas.len() + 2;
    let mut p: Vec<[Real; 3]> = Vec::with_capacity(n);
    p.push([0.0, 0.0, 0.0]);
    p.push([b, 0.0, 0.0]);
    let t0 = PI - thetas[0];
    p.push([b + b * t0.cos(), b * t0.sin(), 0.0]);

    let sub = |u: [Real; 3], v: [Real; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    let cross = |u: [Real; 3], v: [Real; 3]| {
        [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ]
    };
    let unit = |u: [Real; 3]| {
        let m = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        [u[0] / m, u[1] / m, u[2] / m]
    };

    for k in 3..n {
        let th = thetas[k - 2];
        let ph = phis[k - 3];
        let (p0, p1, p2) = (p[k - 3], p[k - 2], p[k - 1]);
        let bc = unit(sub(p2, p1));
        let nrm = unit(cross(sub(p1, p0), bc));
        let m = cross(nrm, bc);
        let d = [
            -b * th.cos(),
            b * th.sin() * ph.cos(),
            b * th.sin() * ph.sin(),
        ];
        p.push([
            p2[0] + d[0] * bc[0] + d[1] * m[0] + d[2] * nrm[0],
            p2[1] + d[0] * bc[1] + d[1] * m[1] + d[2] * nrm[1],
            p2[2] + d[0] * bc[2] + d[1] * m[2] + d[2] * nrm[2],
        ]);
    }

    let mut sys = System::new(n);
    for (i, q) in p.iter().enumerate() {
        sys.pos_x[i] = q[0];
        sys.pos_y[i] = q[1];
        sys.pos_z[i] = q[2];
    }
    sys
}

fn chain(points: &[[Real; 3]]) -> System {
    let mut sys = System::new(points.len());
    for (i, p) in points.iter().enumerate() {
        sys.pos_x[i] = p[0];
        sys.pos_y[i] = p[1];
        sys.pos_z[i] = p[2];
    }
    sys
}

fn ubiquitin() -> System {
    chain(&UBIQUITIN_CA)
}

fn ideal_helix(n: usize) -> System {
    let pts: Vec<[Real; 3]> = (0..n)
        .map(|k| {
            let a = k as Real * 100.0 * DEG;
            [2.3 * a.cos(), 2.3 * a.sin(), 1.5 * k as Real]
        })
        .collect();
    chain(&pts)
}

fn ideal_strand(n: usize) -> System {
    let pts: Vec<[Real; 3]> = (0..n)
        .map(|k| {
            let s = if k % 2 == 0 { 1.0 } else { -1.0 };
            [0.0, 0.942 * s, 3.3 * k as Real]
        })
        .collect();
    chain(&pts)
}

fn mirrored(sys: &System) -> System {
    let mut out = System::new(sys.n);
    for i in 0..sys.n {
        out.pos_x[i] = sys.pos_x[i];
        out.pos_y[i] = sys.pos_y[i];
        out.pos_z[i] = -sys.pos_z[i];
    }
    out
}

fn assign(sys: &System) -> String {
    let mut a = Assignment::new(sys.n);
    a.update(sys);
    a.code()
}

fn segments(code: &str, sym: char) -> Vec<(usize, usize)> {
    let b: Vec<char> = code.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] != sym {
            i += 1;
            continue;
        }
        let start = i;
        while i < b.len() && b[i] == sym {
            i += 1;
        }
        out.push((start + 1, i));
    }
    out
}

#[test]
fn ubiquitin_matches_the_golden_labels() {
    let code = assign(&ubiquitin());
    assert_eq!(
        segments(&code, 'H'),
        vec![(23, 33)],
        "helix segments, got {code}"
    );
    assert_eq!(
        segments(&code, 'E'),
        vec![(2, 6), (11, 17), (41, 44), (64, 73)],
        "strand segments, got {code}"
    );
}

#[test]
fn ubiquitin_agrees_with_the_deposited_records() {
    const HELIX: [(usize, usize); 1] = [(23, 34)];
    const SHEET: [(usize, usize); 5] = [(1, 7), (10, 17), (40, 45), (48, 50), (64, 72)];

    let mut reference = ['C'; 76];
    for (a, b) in HELIX {
        for r in a..=b {
            reference[r - 1] = 'H';
        }
    }
    for (a, b) in SHEET {
        for r in a..=b {
            reference[r - 1] = 'E';
        }
    }

    let code: Vec<char> = assign(&ubiquitin()).chars().collect();
    let agree = code.iter().zip(&reference).filter(|(a, b)| a == b).count();
    let q3 = agree as Real / 76.0;
    assert!(
        q3 >= 0.85,
        "Q3 = {agree}/76 = {q3:.3}; the 3_10 helix at 56-59 is scored as coil on purpose"
    );
}

#[test]
fn an_ideal_helix_sits_at_the_centre_of_every_range() {
    let sys = ideal_helix(12);
    let w = window(&sys, 3);
    assert!((w.d2 - 5.43).abs() < 0.01, "d2 {}", w.d2);
    assert!((w.d3 - 5.05).abs() < 0.01, "d3 {}", w.d3);
    assert!((w.d4 - 6.20).abs() < 0.01, "d4 {}", w.d4);
    assert!(
        (w.theta / DEG - 90.4).abs() < 0.1,
        "theta {}",
        w.theta / DEG
    );
}

#[test]
fn a_right_handed_helix_has_a_positive_dihedral() {
    let sys = ideal_helix(12);
    let tau = window(&sys, 3).tau / DEG;
    assert!(
        (tau - 50.0).abs() < 1.0,
        "right-handed alpha helix must sit at the paper's +50 deg, got {tau:.2}"
    );

    let flipped = window(&mirrored(&sys), 3).tau / DEG;
    assert!(
        (flipped + 50.0).abs() < 1.0,
        "mirroring must flip the sign, got {flipped:.2}"
    );
}

#[test]
fn an_ideal_helix_labels_its_interior() {
    assert_eq!(assign(&ideal_helix(12)), "CHHHHHHHHHHC");
}

#[test]
fn an_ideal_strand_labels_its_interior() {
    assert_eq!(assign(&ideal_strand(10)), "CEEEEEEECC");
}

#[test]
fn labels_are_invariant_under_a_rigid_motion() {
    let sys = ubiquitin();
    let before = assign(&sys);

    let (s, c) = (0.7_f64.sin(), 0.7_f64.cos());
    let mut moved = System::new(sys.n);
    for i in 0..sys.n {
        let (x, y, z) = (sys.pos_x[i], sys.pos_y[i], sys.pos_z[i]);
        moved.pos_x[i] = c * x - s * y + 137.0;
        moved.pos_y[i] = s * x + c * y - 42.0;
        moved.pos_z[i] = z + 9.5;
    }

    assert_eq!(before, assign(&moved));
}

#[test]
fn coincident_beads_stay_coil_and_do_not_panic() {
    let mut sys = ideal_helix(12);
    sys.pos_x[5] = sys.pos_x[4];
    sys.pos_y[5] = sys.pos_y[4];
    sys.pos_z[5] = sys.pos_z[4];

    let code = assign(&sys);
    assert_eq!(code.len(), 12);
    assert!(code.chars().all(|c| c == 'C' || c == 'H' || c == 'E'));
}

#[test]
fn a_chain_shorter_than_a_window_is_all_coil() {
    for n in 0..5 {
        let sys = ideal_helix(n);
        assert_eq!(assign(&sys), "C".repeat(n), "n = {n}");
    }
}

#[test]
fn the_labels_agree_with_the_fractions() {
    let mut a = Assignment::new(76);
    a.update(&ubiquitin());
    let h = a.labels().iter().filter(|&&s| s == Ss::Helix).count();
    let e = a.labels().iter().filter(|&&s| s == Ss::Strand).count();
    assert!((a.helix_fraction() - h as Real / 76.0).abs() < 1e-12);
    assert!((a.strand_fraction() - e as Real / 76.0).abs() < 1e-12);
}

#[test]
fn a_random_coil_is_not_called_structured() {
    let mut rng = Lcg(20260819);
    let mut a = Assignment::new(76);
    let (mut h, mut e) = (0.0, 0.0);
    let chains = 300;
    for _ in 0..chains {
        let thetas: Vec<Real> = (0..74)
            .map(|_| rng.range(80.0 * DEG, 150.0 * DEG))
            .collect();
        let phis: Vec<Real> = (0..73).map(|_| rng.range(-PI, PI)).collect();
        a.update(&build(&thetas, &phis));
        h += a.helix_fraction();
        e += a.strand_fraction();
    }
    let (h, e) = (100.0 * h / chains as Real, 100.0 * e / chains as Real);
    assert!(
        h < 1.0,
        "helix on random coil = {h:.2}%, measured baseline 0.00%"
    );
    assert!(
        e < 4.5,
        "strand on random coil = {e:.2}%; baseline 2.48%, the 4/3 relaxation gives 5.96%"
    );
}
