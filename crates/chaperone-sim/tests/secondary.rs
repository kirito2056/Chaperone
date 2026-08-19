use chaperone_sim::secondary::{window, Assignment, Ss};
use chaperone_sim::system::{Real, System, PI};

const DEG: Real = PI / 180.0;

const UBIQUITIN_CA: [[Real; 3]; 76] = [
    [26.266, 25.413, 2.842],
    [26.85, 29.021, 3.898],
    [26.235, 30.058, 7.497],
    [26.772, 33.436, 9.197],
    [28.605, 33.965, 12.503],
    [27.691, 37.315, 14.143],
    [30.225, 38.643, 16.662],
    [29.607, 41.18, 19.467],
    [31.422, 43.94, 17.553],
    [28.978, 43.96, 14.678],
    [31.191, 42.012, 12.331],
    [29.542, 39.02, 10.653],
    [31.72, 36.289, 9.176],
    [30.505, 33.884, 6.512],
    [31.677, 30.275, 6.639],
    [31.22, 27.341, 4.275],
    [30.288, 24.245, 6.193],
    [28.468, 20.94, 5.98],
    [25.829, 19.825, 8.494],
    [28.054, 16.835, 9.21],
    [30.796, 19.083, 10.566],
    [31.398, 19.064, 14.286],
    [31.288, 22.201, 16.417],
    [35.031, 21.722, 17.069],
    [35.59, 21.945, 13.302],
    [33.533, 25.097, 12.978],
    [35.596, 26.715, 15.736],
    [38.794, 25.761, 13.88],
    [37.471, 27.391, 10.668],
    [36.731, 30.57, 12.645],
    [40.269, 30.508, 14.115],
    [41.718, 30.022, 10.643],
    [39.808, 32.994, 9.233],
    [39.676, 35.547, 12.072],
    [42.345, 34.269, 14.431],
    [40.226, 33.716, 17.509],
    [41.461, 30.751, 19.594],
    [38.817, 28.02, 19.889],
    [39.063, 28.063, 23.695],
    [37.738, 31.637, 23.712],
    [34.738, 30.875, 21.473],
    [31.2, 30.329, 22.78],
    [28.762, 29.573, 19.906],
    [25.034, 30.17, 20.401],
    [22.126, 29.062, 18.183],
    [18.443, 29.143, 19.083],
    [19.399, 29.894, 22.655],
    [21.55, 26.796, 23.133],
    [25.349, 26.872, 23.643],
    [26.826, 24.521, 21.012],
    [29.015, 21.657, 22.288],
    [32.262, 20.67, 20.514],
    [31.568, 16.962, 19.825],
    [28.108, 17.439, 18.276],
    [27.574, 18.192, 14.563],
    [25.594, 21.109, 13.072],
    [22.924, 18.583, 12.025],
    [22.418, 17.638, 15.693],
    [21.079, 21.149, 16.251],
    [19.065, 21.352, 12.999],
    [21.184, 24.263, 11.69],
    [20.081, 24.773, 8.033],
    [21.656, 26.847, 5.24],
    [21.907, 30.563, 5.881],
    [21.419, 30.253, 9.62],
    [23.212, 32.762, 11.891],
    [25.149, 31.609, 14.98],
    [26.179, 34.127, 17.65],
    [29.801, 34.145, 18.829],
    [30.479, 35.369, 22.374],
    [34.145, 35.472, 23.481],
    [35.161, 34.174, 26.896],
    [38.668, 35.502, 27.68],
    [40.873, 33.802, 30.253],
    [41.845, 36.55, 32.686],
    [40.373, 39.813, 33.944],
];

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
