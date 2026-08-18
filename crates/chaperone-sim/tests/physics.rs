use std::sync::OnceLock;

use chaperone_sim::analysis::{
    fraction_of_native_contacts, EnergyMonitor, EnergySummary, PeriodTracker, CONTACT_TOLERANCE,
};
use chaperone_sim::forcefield::angle::{self, Angles};
use chaperone_sim::forcefield::bond::{self, Bonds};
use chaperone_sim::forcefield::dihedral::{self, Dihedrals};
use chaperone_sim::forcefield::native::{self, NativeContacts};
use chaperone_sim::forcefield::pairlist::PairList;
use chaperone_sim::forcefield::pull::{self, Pull, MAX_EXTENSION};
use chaperone_sim::forcefield::repulsion;
use chaperone_sim::forcefield::ForceField;
use chaperone_sim::integrator;
use chaperone_sim::scenario::{
    self, ANCHOR_K, ANGLE_K, BOND_K, EPS, K_PHI1, K_PHI3, MIN_SEQUENCE_SEPARATION, PULL_K, R0,
    SIGMA,
};
use chaperone_sim::system::{Real, System, PI};

const MIN_SEPARATION: Real = 0.8 * SIGMA;
const MIN_TRIPLE_SIN: Real = 0.2;

struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed)
    }

    fn next_unit(&mut self) -> Real {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as Real) / ((1u64 << 53) as Real)
    }
}

fn coord(sys: &mut System, i: usize, dim: usize) -> &mut Real {
    match dim {
        0 => &mut sys.pos_x[i],
        1 => &mut sys.pos_y[i],
        _ => &mut sys.pos_z[i],
    }
}

fn shake(sys: &mut System, rng: &mut Lcg, delta: Real) {
    for i in 0..sys.n {
        sys.pos_x[i] += (rng.next_unit() - 0.5) * 2.0 * delta;
        sys.pos_y[i] += (rng.next_unit() - 0.5) * 2.0 * delta;
        sys.pos_z[i] += (rng.next_unit() - 0.5) * 2.0 * delta;
    }
}

fn random_chain(seed: u64, n: usize, box_size: Real, min_separation: Real) -> (System, Bonds, Lcg) {
    let mut rng = Lcg::new(seed);
    let mut sys = System::new(n);
    let min_sq = min_separation * min_separation;

    for i in 0..n {
        let mut placed = false;
        for _ in 0..100_000 {
            let x = rng.next_unit() * box_size;
            let y = rng.next_unit() * box_size;
            let z = rng.next_unit() * box_size;

            let clear = (0..i).all(|k| {
                let dx = x - sys.pos_x[k];
                let dy = y - sys.pos_y[k];
                let dz = z - sys.pos_z[k];
                dx * dx + dy * dy + dz * dz >= min_sq
            });

            if clear {
                sys.pos_x[i] = x;
                sys.pos_y[i] = y;
                sys.pos_z[i] = z;
                placed = true;
                break;
            }
        }
        assert!(
            placed,
            "seed {seed}: could not place bead {i} of {n} with min separation \
             {min_separation} in box {box_size}"
        );
    }

    let mut bonds = Bonds::new();
    for i in 1..n {
        bonds.push((i - 1) as u32, i as u32, R0);
    }

    (sys, bonds, rng)
}

fn assert_numerical_gradient<A, E>(sys: &mut System, accumulate: A, energy: E, label: &str)
where
    A: Fn(&mut System) -> Real,
    E: Fn(&System) -> Real,
{
    const H: Real = 1e-5;
    const ABS_TOL: Real = 1e-6;
    const REL_TOL: Real = 1e-4;

    sys.clear_forces();
    accumulate(sys);
    let analytic: Vec<[Real; 3]> = (0..sys.n)
        .map(|i| [sys.frc_x[i], sys.frc_y[i], sys.frc_z[i]])
        .collect();

    #[allow(clippy::needless_range_loop)]
    for i in 0..sys.n {
        for dim in 0..3 {
            let orig = *coord(sys, i, dim);

            *coord(sys, i, dim) = orig + H;
            let e_plus = energy(sys);

            *coord(sys, i, dim) = orig - H;
            let e_minus = energy(sys);

            *coord(sys, i, dim) = orig;

            let numerical = -(e_plus - e_minus) / (2.0 * H);
            let a = analytic[i][dim];
            let tol = ABS_TOL + REL_TOL * a.abs().max(numerical.abs());

            assert!(
                (a - numerical).abs() <= tol,
                "{label} atom {i} dim {dim}: analytic {a:.12e} vs numerical {numerical:.12e} \
                 (diff {:.3e}, tol {tol:.3e})",
                (a - numerical).abs()
            );
        }
    }
}

fn run_spring(dt: Real, steps: usize) -> EnergySummary {
    let (mut sys, ff) = scenario::spring(5.0);
    let energies = integrator::initialize(&mut sys, &ff);
    let mut monitor = EnergyMonitor::new(&sys, energies.total(), steps);

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, dt);
        monitor.update(step, &sys, energies.total());
    }

    monitor.summary()
}

struct Chain4Result {
    summary: EnergySummary,
    gap_min: Real,
    gap_max: Real,
}

fn run_chain4(dt: Real, steps: usize) -> Chain4Result {
    let (mut sys, ff) = scenario::chain4();
    let energies = integrator::initialize(&mut sys, &ff);
    let mut monitor = EnergyMonitor::new(&sys, energies.total(), steps);

    let mut gap_min = sys.distance(0, 3);
    let mut gap_max = gap_min;

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, dt);
        monitor.update(step, &sys, energies.total());

        let gap = sys.distance(0, 3);
        gap_min = gap_min.min(gap);
        gap_max = gap_max.max(gap);
    }

    Chain4Result {
        summary: monitor.summary(),
        gap_min,
        gap_max,
    }
}

fn spring_1m() -> &'static EnergySummary {
    static CACHE: OnceLock<EnergySummary> = OnceLock::new();
    CACHE.get_or_init(|| run_spring(1e-3, 1_000_000))
}

fn chain4_1m() -> &'static Chain4Result {
    static CACHE: OnceLock<Chain4Result> = OnceLock::new();
    CACHE.get_or_init(|| run_chain4(1e-3, 1_000_000))
}

struct Chain5Result {
    summary: EnergySummary,
    max_bond: Real,
    max_angle: Real,
    max_dihedral: Real,
    max_native_abs: Real,
    max_repulsion: Real,
}

fn run_chain5(dt: Real, steps: usize) -> Chain5Result {
    let (mut sys, ff) = scenario::chain5();
    let initial = integrator::initialize(&mut sys, &ff);
    let mut monitor = EnergyMonitor::new(&sys, initial.total(), steps);

    let mut max_bond = initial.bond;
    let mut max_angle = initial.angle;
    let mut max_dihedral = initial.dihedral;
    let mut max_native_abs = initial.native.abs();
    let mut max_repulsion = initial.repulsion;

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, dt);
        monitor.update(step, &sys, energies.total());

        max_bond = max_bond.max(energies.bond);
        max_angle = max_angle.max(energies.angle);
        max_dihedral = max_dihedral.max(energies.dihedral);
        max_native_abs = max_native_abs.max(energies.native.abs());
        max_repulsion = max_repulsion.max(energies.repulsion);
    }

    Chain5Result {
        summary: monitor.summary(),
        max_bond,
        max_angle,
        max_dihedral,
        max_native_abs,
        max_repulsion,
    }
}

fn chain5_1m() -> &'static Chain5Result {
    static CACHE: OnceLock<Chain5Result> = OnceLock::new();
    CACHE.get_or_init(|| run_chain5(1e-3, 1_000_000))
}

#[test]
fn bond_force_matches_numerical_gradient() {
    for seed in 0..20u64 {
        let (mut sys, bonds, _) = random_chain(seed, 5, 12.0, MIN_SEPARATION);
        assert_numerical_gradient(
            &mut sys,
            |s| bond::accumulate(s, &bonds, BOND_K),
            |s| bond::energy(s, &bonds, BOND_K),
            &format!("bond seed {seed}"),
        );
    }
}

#[test]
fn repulsion_force_matches_numerical_gradient() {
    for seed in 0..20u64 {
        let (mut sys, _, _) = random_chain(seed, 5, 12.0, MIN_SEPARATION);
        let pairs = PairList::all_pairs(sys.n, 3);
        assert_numerical_gradient(
            &mut sys,
            |s| repulsion::accumulate(s, &pairs, EPS, SIGMA),
            |s| repulsion::energy(s, &pairs, EPS, SIGMA),
            &format!("repulsion seed {seed}"),
        );
    }
}

#[test]
fn repulsion_is_always_repulsive() {
    let mut pairs = PairList::new();
    pairs.push(0, 1);

    for step in 1..=20 {
        let r = 0.5 * SIGMA + 0.15 * SIGMA * step as Real;
        let mut sys = System::new(2);
        sys.pos_x[1] = r;

        sys.clear_forces();
        repulsion::accumulate(&mut sys, &pairs, EPS, SIGMA);

        let dot = sys.frc_x[0] * r;
        assert!(
            dot < 0.0,
            "r = {r:.3}: F_0 . d = {dot:.6e}, expected negative (repulsive)"
        );
    }
}

#[test]
fn newton_third_law() {
    for seed in 0..20u64 {
        let (mut sys, bonds, _) = random_chain(seed, 8, 20.0, MIN_SEPARATION);
        let pairs = PairList::all_pairs(sys.n, 3);

        sys.clear_forces();
        bond::accumulate(&mut sys, &bonds, BOND_K);
        repulsion::accumulate(&mut sys, &pairs, EPS, SIGMA);

        let (fx, fy, fz) = sys.total_force();
        assert!(
            fx.abs() < 1e-9 && fy.abs() < 1e-9 && fz.abs() < 1e-9,
            "seed {seed}: net force ({fx:.3e}, {fy:.3e}, {fz:.3e}) is not zero"
        );
    }
}

#[test]
fn pairlist_respects_sequence_separation() {
    let n = 10;
    let pairs = PairList::all_pairs(n, 3);

    assert_eq!(pairs.len(), (n - 3) * (n - 2) / 2);

    let mut seen = std::collections::HashSet::new();
    for p in 0..pairs.len() {
        let (i, j) = (pairs.i[p], pairs.j[p]);
        assert!(i < j, "pair ({i}, {j}) is not ordered");
        assert!(j - i >= 3, "pair ({i}, {j}) violates separation >= 3");
        assert!(seen.insert((i, j)), "duplicate pair ({i}, {j})");
    }
}

#[test]
fn pairlist_never_emits_self_pairs() {
    for sep in 0..3usize {
        let pairs = PairList::all_pairs(6, sep);
        for p in 0..pairs.len() {
            assert_ne!(pairs.i[p], pairs.j[p], "self-pair from separation {sep}");
        }
    }
}

#[test]
fn spring_energy_stays_bounded() {
    let s = spring_1m();
    assert!(s.is_finite(), "blew up at step {:?}", s.first_nonfinite);
    assert!(
        s.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        s.max_abs_drift
    );
}

#[test]
fn spring_has_no_secular_drift() {
    let s = spring_1m();
    assert!(
        s.secular_drift() < 1e-6,
        "secular drift = {:.3e} (head {:.3e}, tail {:.3e}), expected < 1e-6",
        s.secular_drift(),
        s.head_mean,
        s.tail_mean
    );
}

#[test]
fn chain4_energy_is_conserved() {
    let r = chain4_1m();
    assert!(
        r.summary.is_finite(),
        "blew up at step {:?}",
        r.summary.first_nonfinite
    );
    assert!(
        r.summary.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        r.summary.max_abs_drift
    );
    assert!(
        r.summary.secular_drift() < 1e-6,
        "secular drift = {:.3e}, expected < 1e-6",
        r.summary.secular_drift()
    );
}

#[test]
fn chain4_repulsion_actually_does_work() {
    let r = chain4_1m();

    assert!(
        r.gap_max > R0 + 1e-3,
        "d(0,3) max = {:.6}, never exceeded r0 = {R0}: only repulsion can push 0 and 3 \
         apart, so a static gap means the term is inert",
        r.gap_max
    );

    assert!(
        r.gap_min < SIGMA,
        "d(0,3) min = {:.6}, never entered the repulsive range sigma = {SIGMA}",
        r.gap_min
    );
}

#[test]
fn spring_matches_analytic_solution() {
    const DT: Real = 1e-3;
    const STEPS: usize = 1_000_000;
    const INITIAL_SEPARATION: Real = 5.0;

    let (mut sys, ff) = scenario::spring(INITIAL_SEPARATION);
    integrator::initialize(&mut sys, &ff);

    let mut tracker = PeriodTracker::new(sys.distance(0, 1) - R0, 0.0);
    let mut r_min = INITIAL_SEPARATION;
    let mut r_max = INITIAL_SEPARATION;

    for step in 0..STEPS {
        integrator::step(&mut sys, &ff, DT);

        let r = sys.distance(0, 1);
        r_min = r_min.min(r);
        r_max = r_max.max(r);
        tracker.push(r - R0, (step + 1) as Real * DT);
    }

    let reduced_mass = 0.5;
    let omega = (2.0 * BOND_K / reduced_mass).sqrt();
    let theoretical = 2.0 * PI / omega;
    let measured = tracker.period().expect("no full period observed");
    let rel_error = (measured - theoretical).abs() / theoretical;

    assert!(
        rel_error < 5e-5,
        "period {measured:.7} vs theory {theoretical:.7}, relative error {rel_error:.3e}"
    );

    let amplitude = INITIAL_SEPARATION - R0;
    assert!(
        (r_max - (R0 + amplitude)).abs() < 1e-4 && (r_min - (R0 - amplitude)).abs() < 1e-4,
        "r range [{r_min:.6}, {r_max:.6}], expected symmetric about r0 = {R0}"
    );
}

#[test]
fn native_force_matches_numerical_gradient() {
    let mut beyond_sigma = 0usize;
    let mut within_sigma = 0usize;

    for seed in 0..20u64 {
        let (mut sys, _, mut rng) = random_chain(seed, 5, 12.0, MIN_SEPARATION);
        let contacts = NativeContacts::from_structure(&sys, 12.0, MIN_SEQUENCE_SEPARATION);
        assert!(
            !contacts.is_empty(),
            "seed {seed}: no native contacts found"
        );

        shake(&mut sys, &mut rng, 0.3);

        for c in 0..contacts.len() {
            let r = sys.distance(contacts.i[c] as usize, contacts.j[c] as usize);
            if r > contacts.sigma[c] {
                beyond_sigma += 1;
            } else {
                within_sigma += 1;
            }
        }

        assert_numerical_gradient(
            &mut sys,
            |s| native::accumulate(s, &contacts, EPS),
            |s| native::energy(s, &contacts, EPS),
            &format!("native seed {seed}"),
        );
    }

    assert!(
        beyond_sigma > 0 && within_sigma > 0,
        "shaken configurations covered only one branch \
         (r > sigma: {beyond_sigma}, r < sigma: {within_sigma}); \
         both attraction and repulsion must be exercised"
    );
}

#[test]
fn native_minimum_is_at_sigma() {
    let (mut sys, ff) = scenario::native_pair(SIGMA, SIGMA);
    sys.clear_forces();
    let e = native::accumulate(&mut sys, &ff.native, EPS);

    assert!(
        (e + EPS).abs() < 1e-12,
        "V(sigma) = {e:.15}, expected exactly -eps = {:.15}",
        -EPS
    );

    let (fx, fy, fz) = (sys.frc_x[0], sys.frc_y[0], sys.frc_z[0]);
    let f = (fx * fx + fy * fy + fz * fz).sqrt();
    assert!(f < 1e-9, "|F(sigma)| = {f:.3e}, expected 0");
}

#[test]
fn native_switches_sign_at_sigma() {
    let mut contacts = NativeContacts::new();
    contacts.push(0, 1, SIGMA);

    for step in 0..=26 {
        let r = 0.7 * SIGMA + 0.05 * SIGMA * step as Real;
        if (r - SIGMA).abs() < 1e-6 {
            continue;
        }

        let mut sys = System::new(2);
        sys.pos_x[1] = r;
        sys.clear_forces();
        native::accumulate(&mut sys, &contacts, EPS);

        let dot = sys.frc_x[0] * r;
        if r < SIGMA {
            assert!(
                dot < 0.0,
                "r = {r:.3} < sigma: F.d = {dot:.6e}, expected repulsive"
            );
        } else {
            assert!(
                dot > 0.0,
                "r = {r:.3} > sigma: F.d = {dot:.6e}, expected attractive"
            );
        }
    }
}

#[test]
fn native_oscillation_matches_harmonic_approximation() {
    const DT: Real = 1e-3;
    const STEPS: usize = 1_000_000;
    const AMPLITUDE: Real = 0.02;

    let (mut sys, ff) = scenario::native_pair(SIGMA, SIGMA + AMPLITUDE);
    integrator::initialize(&mut sys, &ff);

    let mut tracker = PeriodTracker::new(sys.distance(0, 1) - SIGMA, 0.0);
    for step in 0..STEPS {
        integrator::step(&mut sys, &ff, DT);
        tracker.push(sys.distance(0, 1) - SIGMA, (step + 1) as Real * DT);
    }

    let k_eff = 120.0 * EPS / (SIGMA * SIGMA);
    let reduced_mass = 0.5;
    let omega = (k_eff / reduced_mass).sqrt();
    let theoretical = 2.0 * PI / omega;
    let measured = tracker.period().expect("no full period observed");
    let rel_error = (measured - theoretical).abs() / theoretical;

    assert!(
        rel_error < 1e-2,
        "period {measured:.6} vs harmonic {theoretical:.6}, relative error {rel_error:.3e}; \
         anharmonic shift is expected to be ~2.1 * A^2 = {:.3e}",
        2.1 * AMPLITUDE * AMPLITUDE
    );
}

#[test]
fn pair_lists_partition_all_pairs() {
    let (sys, ff) = scenario::chain5();
    let all = PairList::all_pairs(sys.n, MIN_SEQUENCE_SEPARATION);

    assert_eq!(
        ff.native.len() + ff.repulsion_pairs.len(),
        all.len(),
        "native ({}) + non-native ({}) != all pairs ({})",
        ff.native.len(),
        ff.repulsion_pairs.len(),
        all.len()
    );

    let key = |i: u32, j: u32| (i.min(j), i.max(j));
    let native: std::collections::HashSet<_> = (0..ff.native.len())
        .map(|c| key(ff.native.i[c], ff.native.j[c]))
        .collect();

    for p in 0..ff.repulsion_pairs.len() {
        let k = key(ff.repulsion_pairs.i[p], ff.repulsion_pairs.j[p]);
        assert!(!native.contains(&k), "pair {k:?} is in both lists");
    }
}

#[test]
fn chain5_energy_is_conserved() {
    let r = chain5_1m();

    assert!(
        r.summary.is_finite(),
        "blew up at step {:?}",
        r.summary.first_nonfinite
    );
    assert!(
        r.summary.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        r.summary.max_abs_drift
    );
    assert!(
        r.summary.secular_drift() < 3e-5,
        "secular drift = {:.3e}, expected < 3e-5; a coupled three-term system is chaotic, \
         so the shadow energy diffuses (~sqrt(T)) instead of oscillating and the 1e-6 bound \
         used for the integrable spring and U-chain does not apply",
        r.summary.secular_drift()
    );
}

#[test]
fn chain5_exercises_every_term() {
    let r = chain5_1m();

    assert!(
        r.max_bond > 1e-3,
        "bond energy peaked at {:.3e}, term looks inert",
        r.max_bond
    );
    assert!(
        r.max_angle > 0.1,
        "angle energy peaked at {:.3e}, term looks inert",
        r.max_angle
    );
    assert!(
        r.max_dihedral > 0.1,
        "dihedral energy peaked at {:.3e}, term looks inert",
        r.max_dihedral
    );
    assert!(
        r.max_native_abs > 1e-3,
        "native energy peaked at {:.3e}, term looks inert",
        r.max_native_abs
    );
    assert!(
        r.max_repulsion > 1e-3,
        "repulsion energy peaked at {:.3e}, term looks inert",
        r.max_repulsion
    );
}

#[test]
fn q_is_one_for_the_native_structure() {
    let (sys, ff) = scenario::chain5();
    let q = fraction_of_native_contacts(&sys, &ff.native, CONTACT_TOLERANCE)
        .expect("chain5 has native contacts");
    assert_eq!(q, 1.0, "native structure should have Q = 1");
}

#[test]
fn q_drops_when_the_structure_is_scrambled() {
    let (mut sys, ff) = scenario::chain5();
    for i in 0..sys.n {
        sys.pos_x[i] *= 5.0;
        sys.pos_y[i] *= 5.0;
        sys.pos_z[i] *= 5.0;
    }

    let q = fraction_of_native_contacts(&sys, &ff.native, CONTACT_TOLERANCE)
        .expect("chain5 has native contacts");
    assert_eq!(q, 0.0, "expanded structure should have no contacts formed");
}

#[test]
fn q_counts_the_threshold_exactly() {
    let mut sys = System::new(4);
    sys.pos_x[1] = 1.19 * SIGMA;
    sys.pos_y[2] = 20.0;
    sys.pos_y[3] = 20.0;
    sys.pos_x[3] = 1.21 * SIGMA;

    let mut contacts = NativeContacts::new();
    contacts.push(0, 1, SIGMA);
    contacts.push(2, 3, SIGMA);

    let q = fraction_of_native_contacts(&sys, &contacts, CONTACT_TOLERANCE)
        .expect("two contacts defined");
    assert_eq!(q, 0.5, "one contact inside 1.2 sigma, one outside");
}

#[test]
fn q_is_none_without_contacts() {
    let sys = System::new(2);
    let contacts = NativeContacts::new();
    assert!(fraction_of_native_contacts(&sys, &contacts, CONTACT_TOLERANCE).is_none());
}

fn triples_are_well_conditioned(sys: &System, min_sin: Real) -> bool {
    (0..sys.n.saturating_sub(2)).all(|t| sys.angle(t, t + 1, t + 2).sin() >= min_sin)
}

fn random_conditioned_chain(seed: u64, n: usize, delta: Real) -> (System, System) {
    for attempt in 0..64u64 {
        let (mut sys, _, mut rng) = random_chain(seed * 64 + attempt, n, 12.0, MIN_SEPARATION);
        if !triples_are_well_conditioned(&sys, MIN_TRIPLE_SIN) {
            continue;
        }

        let mut native = System::new(n);
        native.pos_x.copy_from_slice(&sys.pos_x);
        native.pos_y.copy_from_slice(&sys.pos_y);
        native.pos_z.copy_from_slice(&sys.pos_z);

        for _ in 0..64 {
            sys.pos_x.copy_from_slice(&native.pos_x);
            sys.pos_y.copy_from_slice(&native.pos_y);
            sys.pos_z.copy_from_slice(&native.pos_z);
            shake(&mut sys, &mut rng, delta);
            if triples_are_well_conditioned(&sys, MIN_TRIPLE_SIN) {
                return (sys, native);
            }
        }
    }
    panic!("seed {seed}: could not build a chain with all triple sines >= {MIN_TRIPLE_SIN}");
}

#[test]
fn angle_force_matches_numerical_gradient() {
    let mut wider = 0usize;
    let mut narrower = 0usize;

    for seed in 0..20u64 {
        let (mut sys, native) = random_conditioned_chain(seed, 5, 0.3);
        let angles = Angles::from_chain(&native);

        for t in 0..angles.len() {
            let theta = sys.angle(
                angles.i[t] as usize,
                angles.j[t] as usize,
                angles.k[t] as usize,
            );
            if theta > angles.theta0[t] {
                wider += 1;
            } else {
                narrower += 1;
            }
        }

        assert_numerical_gradient(
            &mut sys,
            |s| angle::accumulate(s, &angles, ANGLE_K),
            |s| angle::energy(s, &angles, ANGLE_K),
            &format!("angle seed {seed}"),
        );
    }

    assert!(
        wider > 0 && narrower > 0,
        "shaken chains covered only one branch (theta > theta0: {wider}, theta < theta0: \
         {narrower}); both restoring directions must be exercised"
    );
}

#[test]
fn angle_produces_no_net_force_or_torque() {
    for seed in 0..20u64 {
        let (mut sys, native) = random_conditioned_chain(seed, 5, 0.3);
        let angles = Angles::from_chain(&native);
        sys.clear_forces();
        angle::accumulate(&mut sys, &angles, ANGLE_K);

        let (fx, fy, fz) = sys.total_force();
        assert!(
            fx.abs() < 1e-9 && fy.abs() < 1e-9 && fz.abs() < 1e-9,
            "seed {seed}: net force ({fx:.3e}, {fy:.3e}, {fz:.3e})"
        );

        let (tx, ty, tz) = sys.total_torque();
        assert!(
            tx.abs() < 1e-9 && ty.abs() < 1e-9 && tz.abs() < 1e-9,
            "seed {seed}: net torque ({tx:.3e}, {ty:.3e}, {tz:.3e})"
        );

        for i in 0..sys.n {
            sys.pos_x[i] += 100.0;
            sys.pos_y[i] -= 37.0;
            sys.pos_z[i] += 5.0;
        }
        let (sx, sy, sz) = sys.total_torque();
        assert!(
            (sx - tx).abs() < 1e-8 && (sy - ty).abs() < 1e-8 && (sz - tz).abs() < 1e-8,
            "seed {seed}: torque depends on origin choice"
        );
    }
}

#[test]
fn angle_force_is_perpendicular_to_its_own_arm() {
    let mut rng = Lcg::new(9);

    for _ in 0..50 {
        let mut sys = System::new(3);
        for i in [0usize, 2] {
            sys.pos_x[i] = (rng.next_unit() - 0.5) * 8.0;
            sys.pos_y[i] = (rng.next_unit() - 0.5) * 8.0;
            sys.pos_z[i] = (rng.next_unit() - 0.5) * 8.0;
        }

        let theta = sys.angle(0, 1, 2);
        if theta.sin() < MIN_TRIPLE_SIN {
            continue;
        }

        let theta0 = theta - 0.15;
        if theta0 <= 0.0 {
            continue;
        }
        let mut angles = Angles::new();
        angles.push(0, 1, 2, theta0);

        sys.clear_forces();
        angle::accumulate(&mut sys, &angles, ANGLE_K);

        let arms = [(0usize, 1usize), (2usize, 1usize)];
        for (end, vertex) in arms {
            let ax = sys.pos_x[end] - sys.pos_x[vertex];
            let ay = sys.pos_y[end] - sys.pos_y[vertex];
            let az = sys.pos_z[end] - sys.pos_z[vertex];
            let ra = (ax * ax + ay * ay + az * az).sqrt();

            let dot = sys.frc_x[end] * ax + sys.frc_y[end] * ay + sys.frc_z[end] * az;
            assert!(
                dot.abs() < 1e-9,
                "F.arm = {dot:.3e} for bead {end}: the projection term is missing"
            );

            let f = (sys.frc_x[end] * sys.frc_x[end]
                + sys.frc_y[end] * sys.frc_y[end]
                + sys.frc_z[end] * sys.frc_z[end])
                .sqrt();
            let expected = (2.0 * ANGLE_K * (theta - theta0)).abs() / ra;
            assert!(
                (f - expected).abs() < 1e-9 * expected.max(1.0),
                "|F| = {f:.9} but |dV/dtheta|/r = {expected:.9}: the 1/sin factor is missing"
            );
        }
    }
}

#[test]
fn angle_force_matches_closed_form_at_60_degrees() {
    let root3 = (3.0 as Real).sqrt();
    let mut sys = System::new(3);
    sys.pos_x[0] = 1.0;
    sys.pos_x[2] = 0.5;
    sys.pos_y[2] = root3 / 2.0;

    let mut angles = Angles::new();
    angles.push(0, 1, 2, PI / 3.0 - 0.1);

    sys.clear_forces();
    let e = angle::accumulate(&mut sys, &angles, 20.0);

    assert!(
        (e - 0.2).abs() < 1e-12,
        "V = {e:.15}, expected K * 0.1^2 = 0.2"
    );

    let expected: [[Real; 3]; 3] = [
        [0.0, 4.0, 0.0],
        [-2.0 * root3, -2.0, 0.0],
        [2.0 * root3, -2.0, 0.0],
    ];
    for (i, want) in expected.iter().enumerate() {
        let got = [sys.frc_x[i], sys.frc_y[i], sys.frc_z[i]];
        for dim in 0..3 {
            assert!(
                (got[dim] - want[dim]).abs() < 1e-12,
                "bead {i} dim {dim}: got {:.15}, expected {:.15}",
                got[dim],
                want[dim]
            );
        }
    }
}

#[test]
fn angle_force_matches_closed_form_at_90_degrees() {
    let ra = 2.0;
    let rb = 3.0;
    let mut sys = System::new(3);
    sys.pos_x[0] = ra;
    sys.pos_y[2] = rb;

    let mut angles = Angles::new();
    angles.push(0, 1, 2, PI / 2.0 - 0.1);

    sys.clear_forces();
    angle::accumulate(&mut sys, &angles, 20.0);

    let dv = 2.0 * 20.0 * 0.1;
    assert!((sys.frc_y[0] - dv / ra).abs() < 1e-12 && sys.frc_x[0].abs() < 1e-12);
    assert!((sys.frc_x[2] - dv / rb).abs() < 1e-12 && sys.frc_y[2].abs() < 1e-12);
}

#[test]
fn angle_is_skipped_at_a_linear_geometry() {
    let mut sys = System::new(3);
    sys.pos_x[0] = -3.0;
    sys.pos_x[2] = 4.0;

    let mut angles = Angles::new();
    angles.push(0, 1, 2, PI / 2.0);

    sys.clear_forces();
    let e = angle::accumulate(&mut sys, &angles, ANGLE_K);

    assert_eq!(e, 0.0, "a collinear triple must be skipped, not evaluated");
    for i in 0..3 {
        assert!(
            sys.frc_x[i].is_finite() && sys.frc_y[i].is_finite() && sys.frc_z[i].is_finite(),
            "bead {i} force is not finite at a linear geometry"
        );
        assert_eq!(sys.frc_x[i], 0.0);
        assert_eq!(sys.frc_y[i], 0.0);
        assert_eq!(sys.frc_z[i], 0.0);
    }
}

#[test]
fn chain5_conserves_angular_momentum() {
    let r = chain5_1m();
    assert!(
        r.summary.max_angular_momentum_drift < 1e-9,
        "max |L - L0| = {:.3e}, expected < 1e-9; velocity Verlet conserves L exactly when \
         the total torque vanishes, so any growth means a force term produces net torque",
        r.summary.max_angular_momentum_drift
    );
}

fn torsion_frame(i: [Real; 3], j: [Real; 3], k: [Real; 3], l: [Real; 3]) -> System {
    let mut sys = System::new(4);
    for (idx, p) in [i, j, k, l].iter().enumerate() {
        sys.pos_x[idx] = p[0];
        sys.pos_y[idx] = p[1];
        sys.pos_z[idx] = p[2];
    }
    sys
}

fn assert_forces_match(sys: &System, expected: &[[Real; 3]], tol: Real, label: &str) {
    for (i, want) in expected.iter().enumerate() {
        let got = [sys.frc_x[i], sys.frc_y[i], sys.frc_z[i]];
        for dim in 0..3 {
            assert!(
                (got[dim] - want[dim]).abs() < tol,
                "{label} bead {i} dim {dim}: got {:.15}, expected {:.15}",
                got[dim],
                want[dim]
            );
        }
    }
}

#[test]
fn dihedral_force_matches_numerical_gradient() {
    let mut wider = 0usize;
    let mut narrower = 0usize;

    for seed in 0..20u64 {
        let (mut sys, native) = random_conditioned_chain(seed, 6, 0.3);
        let dihedrals = Dihedrals::from_chain(&native);
        assert!(!dihedrals.is_empty(), "seed {seed}: no dihedrals built");

        for t in 0..dihedrals.len() {
            let phi = sys.dihedral(
                dihedrals.i[t] as usize,
                dihedrals.j[t] as usize,
                dihedrals.k[t] as usize,
                dihedrals.l[t] as usize,
            );
            if phi > dihedrals.phi0[t] {
                wider += 1;
            } else {
                narrower += 1;
            }
        }

        assert_numerical_gradient(
            &mut sys,
            |s| dihedral::accumulate(s, &dihedrals, K_PHI1, K_PHI3),
            |s| dihedral::energy(s, &dihedrals, K_PHI1, K_PHI3),
            &format!("dihedral seed {seed}"),
        );
    }

    assert!(
        wider > 0 && narrower > 0,
        "shaken chains covered only one branch (phi > phi0: {wider}, phi < phi0: {narrower})"
    );
}

#[test]
fn dihedral_matches_closed_form_perpendicular() {
    let psi: Real = 0.9;
    let mut sys = torsion_frame(
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, psi.cos(), psi.sin()],
    );

    let phi = sys.dihedral(0, 1, 2, 3);
    assert!(
        (phi - psi).abs() < 1e-14,
        "phi = {phi:.15}, expected {psi:.15}"
    );

    let mut dihedrals = Dihedrals::new();
    dihedrals.push(0, 1, 2, 3, psi - 0.4);

    sys.clear_forces();
    let e = dihedral::accumulate(&mut sys, &dihedrals, 1.0, 0.5);
    assert!(
        (e - 0.397760128758778).abs() < 1e-12,
        "V = {e:.15}, expected 0.397760128758778"
    );

    let fi = [0.0, 0.0, 1.78747697125949];
    let fk = [0.0, -1.40017881192699, 1.111113503389155];
    assert_forces_match(
        &sys,
        &[fi, [-fi[0], -fi[1], -fi[2]], fk, [-fk[0], -fk[1], -fk[2]]],
        1e-12,
        "perpendicular",
    );
}

#[test]
fn dihedral_matches_closed_form_skewed() {
    let mut sys = torsion_frame(
        [-0.5, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.6, 0.8, 0.9],
    );

    let phi = sys.dihedral(0, 1, 2, 3);
    assert!(
        (phi - 0.844153986113171).abs() < 1e-12,
        "phi = {phi:.15}, expected 0.844153986113171"
    );

    let mut dihedrals = Dihedrals::new();
    dihedrals.push(0, 1, 2, 3, 0.444153986113171);

    sys.clear_forces();
    let e = dihedral::accumulate(&mut sys, &dihedrals, 1.0, 0.5);
    assert!(
        (e - 0.397760128758778).abs() < 1e-12,
        "V = {e:.15}, expected 0.397760128758778"
    );

    assert_forces_match(
        &sys,
        &[
            [0.0, 0.0, 1.78747697125949],
            [0.0, 0.332840539475905, -2.530204471386278],
            [0.0, -1.442309004395589, 1.728921691166507],
            [0.0, 1.109468464919684, -0.986194191039719],
        ],
        1e-11,
        "skewed",
    );
}

#[test]
fn dihedral_force_magnitude_identity_holds() {
    for seed in 0..20u64 {
        let (mut sys, native) = random_conditioned_chain(seed, 4, 0.3);
        let dihedrals = Dihedrals::from_chain(&native);

        let phi = sys.dihedral(0, 1, 2, 3);
        let dphi = phi - dihedrals.phi0[0];
        let dv = (K_PHI1 * dphi.sin() + 3.0 * K_PHI3 * (3.0 * dphi).sin()).abs();

        sys.clear_forces();
        dihedral::accumulate(&mut sys, &dihedrals, K_PHI1, K_PHI3);

        let b1 = [
            sys.pos_x[1] - sys.pos_x[0],
            sys.pos_y[1] - sys.pos_y[0],
            sys.pos_z[1] - sys.pos_z[0],
        ];
        let b2 = [
            sys.pos_x[2] - sys.pos_x[1],
            sys.pos_y[2] - sys.pos_y[1],
            sys.pos_z[2] - sys.pos_z[1],
        ];
        let n1 = [
            b1[1] * b2[2] - b1[2] * b2[1],
            b1[2] * b2[0] - b1[0] * b2[2],
            b1[0] * b2[1] - b1[1] * b2[0],
        ];
        let n1_len = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
        let b2_len = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();

        let f_i = (sys.frc_x[0] * sys.frc_x[0]
            + sys.frc_y[0] * sys.frc_y[0]
            + sys.frc_z[0] * sys.frc_z[0])
            .sqrt();

        let lhs = f_i * n1_len;
        let rhs = dv * b2_len;
        assert!(
            (lhs - rhs).abs() < 1e-9 * rhs.max(1.0),
            "seed {seed}: |F_i|.|n1| = {lhs:.12} but |dV/dphi|.|b2| = {rhs:.12}; \
             the |b2| factor is missing"
        );

        let dot_b1 = sys.frc_x[0] * b1[0] + sys.frc_y[0] * b1[1] + sys.frc_z[0] * b1[2];
        let dot_b2 = sys.frc_x[0] * b2[0] + sys.frc_y[0] * b2[1] + sys.frc_z[0] * b2[2];
        assert!(
            dot_b1.abs() < 1e-9 && dot_b2.abs() < 1e-9,
            "seed {seed}: F_i is not perpendicular to the i-j-k plane normal"
        );
    }
}

#[test]
fn dihedral_is_periodic_in_phi0() {
    let psi: Real = 0.9;
    let build = |phi0: Real| {
        let mut sys = torsion_frame(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, psi.cos(), psi.sin()],
        );
        let mut dihedrals = Dihedrals::new();
        dihedrals.push(0, 1, 2, 3, phi0);
        sys.clear_forces();
        let e = dihedral::accumulate(&mut sys, &dihedrals, K_PHI1, K_PHI3);
        (sys, e)
    };

    let (a, ea) = build(0.5);
    let (b, eb) = build(0.5 + 2.0 * PI);

    assert!(
        (ea - eb).abs() < 1e-13,
        "V differs by 2*pi shift: {ea:.15} vs {eb:.15}"
    );
    for i in 0..4 {
        assert!(
            (a.frc_x[i] - b.frc_x[i]).abs() < 1e-12
                && (a.frc_y[i] - b.frc_y[i]).abs() < 1e-12
                && (a.frc_z[i] - b.frc_z[i]).abs() < 1e-12,
            "bead {i} force differs by 2*pi shift"
        );
    }
}

#[test]
fn dihedral_profile_matches_one_and_three_fold_terms() {
    for step in 0..40 {
        let psi = -PI + (step as Real + 0.5) * (2.0 * PI / 40.0);
        let mut sys = torsion_frame(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, psi.cos(), psi.sin()],
        );

        let mut dihedrals = Dihedrals::new();
        dihedrals.push(0, 1, 2, 3, 0.0);

        sys.clear_forces();
        let e = dihedral::accumulate(&mut sys, &dihedrals, K_PHI1, K_PHI3);

        let expected_e = K_PHI1 * (1.0 - psi.cos()) + K_PHI3 * (1.0 - (3.0 * psi).cos());
        assert!(
            (e - expected_e).abs() < 1e-12,
            "psi = {psi:.6}: V = {e:.12}, expected {expected_e:.12}"
        );

        let expected_dv = (K_PHI1 * psi.sin() + 3.0 * K_PHI3 * (3.0 * psi).sin()).abs();
        let f_i = (sys.frc_x[0] * sys.frc_x[0]
            + sys.frc_y[0] * sys.frc_y[0]
            + sys.frc_z[0] * sys.frc_z[0])
            .sqrt();
        assert!(
            (f_i - expected_dv).abs() < 1e-12,
            "psi = {psi:.6}: |F_i| = {f_i:.12}, expected |dV/dphi| = {expected_dv:.12}; \
             the 3-fold term factor of 3 is wrong"
        );
    }
}

#[test]
fn dihedral_produces_no_net_force_or_torque() {
    for seed in 0..20u64 {
        let (mut sys, native) = random_conditioned_chain(seed, 6, 0.3);
        let dihedrals = Dihedrals::from_chain(&native);

        sys.clear_forces();
        dihedral::accumulate(&mut sys, &dihedrals, K_PHI1, K_PHI3);

        let (fx, fy, fz) = sys.total_force();
        assert!(
            fx.abs() < 1e-9 && fy.abs() < 1e-9 && fz.abs() < 1e-9,
            "seed {seed}: net force ({fx:.3e}, {fy:.3e}, {fz:.3e})"
        );

        let (tx, ty, tz) = sys.total_torque();
        assert!(
            tx.abs() < 1e-9 && ty.abs() < 1e-9 && tz.abs() < 1e-9,
            "seed {seed}: net torque ({tx:.3e}, {ty:.3e}, {tz:.3e}); the F_j / F_k \
             redistribution coefficients are wrong"
        );

        for i in 0..sys.n {
            sys.pos_x[i] += 100.0;
            sys.pos_y[i] -= 37.0;
            sys.pos_z[i] += 5.0;
        }
        let (sx, sy, sz) = sys.total_torque();
        assert!(
            (sx - tx).abs() < 1e-7 && (sy - ty).abs() < 1e-7 && (sz - tz).abs() < 1e-7,
            "seed {seed}: torque depends on origin choice"
        );
    }
}

#[test]
#[should_panic(expected = "nearly collinear")]
fn dihedral_from_chain_rejects_a_collinear_native_structure() {
    let sys = torsion_frame(
        [0.0, 0.0, 0.0],
        [3.8, 0.0, 0.0],
        [7.6, 0.0, 0.0],
        [7.6, 3.8, 0.0],
    );
    Dihedrals::from_chain(&sys);
}

#[test]
fn dihedral_from_chain_accepts_a_bent_native_structure() {
    let (_, ff) = scenario::chain5();
    assert_eq!(ff.dihedrals.len(), 2);
}

#[test]
fn dihedral_is_skipped_at_a_collinear_geometry() {
    let mut sys = torsion_frame(
        [-3.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.0, 1.0, 1.0],
    );

    let mut dihedrals = Dihedrals::new();
    dihedrals.push(0, 1, 2, 3, 0.7);

    sys.clear_forces();
    let e = dihedral::accumulate(&mut sys, &dihedrals, K_PHI1, K_PHI3);

    assert_eq!(e, 0.0, "a collinear i-j-k must be skipped, not evaluated");
    for i in 0..4 {
        assert_eq!(sys.frc_x[i], 0.0);
        assert_eq!(sys.frc_y[i], 0.0);
        assert_eq!(sys.frc_z[i], 0.0);
    }
}

#[test]
fn native_structure_is_a_bonded_force_free_state() {
    let (mut sys, ff) = scenario::chain5();

    type TermFn = fn(&mut System, &chaperone_sim::ForceField) -> Real;
    let terms: [(&str, TermFn); 4] = [
        ("bond", |s, f| bond::accumulate(s, &f.bonds, f.bond_k)),
        ("angle", |s, f| angle::accumulate(s, &f.angles, f.angle_k)),
        ("dihedral", |s, f| {
            dihedral::accumulate(s, &f.dihedrals, f.k_phi1, f.k_phi3)
        }),
        ("native", |s, f| native::accumulate(s, &f.native, f.eps)),
    ];

    for (name, term) in terms {
        sys.clear_forces();
        let e = term(&mut sys, &ff);
        let f = sys.max_force();
        assert!(
            f < 1e-9,
            "{name}: max |F| = {f:.3e} at the native structure (E = {e:.6}); the parameters \
             recovered by from_structure do not round-trip through the force routine"
        );
    }

    sys.clear_forces();
    repulsion::accumulate(&mut sys, &ff.repulsion_pairs, ff.eps, ff.sigma);
    let f = sys.max_force();
    assert!(
        f.is_finite() && f > 0.0,
        "non-native repulsion should be the only term with force at the native structure, \
         got max |F| = {f:.3e}"
    );
}

#[test]
fn integrator_is_second_order() {
    let coarse = run_spring(1e-3, 200_000).max_abs_drift;
    let fine = run_spring(1e-4, 2_000_000).max_abs_drift;
    let ratio = coarse / fine;

    assert!(
        (50.0..=200.0).contains(&ratio),
        "drift ratio {ratio:.1} (coarse {coarse:.3e}, fine {fine:.3e}), expected ~100"
    );
}

// ---------------------------------------------------------------------------
// H. steered MD — the only external force in the model
// ---------------------------------------------------------------------------

fn pulled(index: u32, target: [Real; 3]) -> Pull {
    Pull {
        index: Some(index),
        target,
        k: PULL_K,
    }
}

fn single_bead(at: [Real; 3]) -> System {
    let mut sys = System::new(1);
    sys.pos_x[0] = at[0];
    sys.pos_y[0] = at[1];
    sys.pos_z[0] = at[2];
    sys
}

#[test]
fn pull_matches_a_closed_form() {
    // Inside the clamp: |d| = 3, V = k/2 * 9 = 90, F = k * d
    let mut sys = single_bead([0.0, 0.0, 0.0]);
    let p = pulled(0, [1.0, 2.0, 2.0]);
    sys.clear_forces();
    let e = pull::accumulate(&mut sys, &p);
    assert!((e - 90.0).abs() < 1e-12, "V = {e:.12}, expected 90");
    assert!((sys.frc_x[0] - 20.0).abs() < 1e-12);
    assert!((sys.frc_y[0] - 40.0).abs() < 1e-12);
    assert!((sys.frc_z[0] - 40.0).abs() < 1e-12);

    // Exactly at the clamp: |d| = 5, V = 250, F = (60, 80, 0)
    let mut sys = single_bead([1.0, 2.0, 2.0]);
    let p = pulled(0, [4.0, 6.0, 2.0]);
    sys.clear_forces();
    let e = pull::accumulate(&mut sys, &p);
    assert!((e - 250.0).abs() < 1e-12, "V = {e:.12}, expected 250");
    assert!((sys.frc_x[0] - 60.0).abs() < 1e-12);
    assert!((sys.frc_y[0] - 80.0).abs() < 1e-12);
    assert!(sys.frc_z[0].abs() < 1e-12);

    // Beyond the clamp: |d| = 10, V = 250 + k*5*5 = 750, |F| = k*5 = 100
    let mut sys = single_bead([0.0, 0.0, 0.0]);
    let p = pulled(0, [0.0, 0.0, 10.0]);
    sys.clear_forces();
    let e = pull::accumulate(&mut sys, &p);
    assert!((e - 750.0).abs() < 1e-12, "V = {e:.12}, expected 750");
    assert!((sys.frc_z[0] - 100.0).abs() < 1e-12);
}

#[test]
fn pull_force_is_capped_beyond_the_clamp() {
    let p = pulled(0, [0.0, 0.0, 0.0]);
    for reach in [1.0, 4.9, MAX_EXTENSION, 20.0, 200.0] {
        let sys = single_bead([0.0, 0.0, -reach]);
        let f = p.force_magnitude(&sys);
        let expected = PULL_K * reach.min(MAX_EXTENSION);
        assert!(
            (f - expected).abs() < 1e-9,
            "reach {reach}: |F| = {f:.6}, expected {expected:.6}; without the clamp a \
             20 A mouse flick injects k/2 * 400 = 4000 eps, about twenty times the \
             whole folding well"
        );
    }
}

#[test]
fn pull_force_matches_numerical_gradient() {
    for (seed, offset) in [
        (0u64, [1.7, -0.9, 0.4]),
        (1, [-2.2, 1.1, 3.0]),
        (2, [0.3, 0.2, -0.1]),
        (3, [9.0, -7.0, 4.0]),
    ] {
        let (mut sys, _, _) = random_chain(seed, 5, 12.0, MIN_SEPARATION);
        let target = [
            sys.pos_x[2] + offset[0],
            sys.pos_y[2] + offset[1],
            sys.pos_z[2] + offset[2],
        ];
        let p = pulled(2, target);
        assert_numerical_gradient(
            &mut sys,
            |s| pull::accumulate(s, &p),
            |s| pull::energy(s, &p),
            &format!("pull seed {seed}"),
        );
    }
}

#[test]
fn pull_restores_toward_the_target() {
    let mut rng = Lcg::new(808);
    for _ in 0..40 {
        let at = [
            (rng.next_unit() - 0.5) * 20.0,
            (rng.next_unit() - 0.5) * 20.0,
            (rng.next_unit() - 0.5) * 20.0,
        ];
        let target = [
            (rng.next_unit() - 0.5) * 20.0,
            (rng.next_unit() - 0.5) * 20.0,
            (rng.next_unit() - 0.5) * 20.0,
        ];
        let mut sys = single_bead(at);
        let p = pulled(0, target);
        sys.clear_forces();
        pull::accumulate(&mut sys, &p);

        let dot = sys.frc_x[0] * (target[0] - at[0])
            + sys.frc_y[0] * (target[1] - at[1])
            + sys.frc_z[0] * (target[2] - at[2]);
        assert!(dot > 0.0, "F . (target - r) = {dot:.6e}, expected positive");
    }
}

#[test]
fn pull_equipartition_matches_three_t_over_k() {
    const T: Real = 0.5;
    const GAMMA: Real = 10.0;
    const DT: Real = 0.005;
    const STEPS: usize = 2_000_000;

    let mut sys = single_bead([0.0, 0.0, 0.0]);
    let mut ff = ForceField::new(BOND_K, ANGLE_K, K_PHI1, K_PHI3, EPS, SIGMA);
    ff.pull = pulled(0, [0.0, 0.0, 0.0]);

    let mut bath = chaperone_sim::thermostat::Langevin::new(GAMMA, T, DT, 555);
    integrator::initialize(&mut sys, &ff);

    let mut sum = 0.0;
    for _ in 0..STEPS {
        bath.step(&mut sys, &ff);
        sum +=
            sys.pos_x[0] * sys.pos_x[0] + sys.pos_y[0] * sys.pos_y[0] + sys.pos_z[0] * sys.pos_z[0];
    }

    let measured = sum / STEPS as Real;
    let expected = 3.0 * T / PULL_K;
    assert!(
        (measured / expected - 1.0).abs() < 0.08,
        "<|r - target|^2> = {measured:.6} against 3T/k = {expected:.6}; dropping the \
         one half from V = k/2 d^2 reads exactly half here"
    );
}

#[test]
fn pulling_accounts_for_the_entire_net_force() {
    let (mut sys, ff_base) = scenario::chain5();
    let mut ff = ff_base;
    let target = [sys.pos_x[1] + 2.0, sys.pos_y[1] - 1.5, sys.pos_z[1] + 0.8];
    ff.pull = pulled(1, target);

    sys.clear_forces();
    ff.accumulate(&mut sys);

    let (fx, fy, fz) = sys.total_force();
    let expected = [
        PULL_K * (target[0] - sys.pos_x[1]),
        PULL_K * (target[1] - sys.pos_y[1]),
        PULL_K * (target[2] - sys.pos_z[1]),
    ];
    for (got, want) in [(fx, expected[0]), (fy, expected[1]), (fz, expected[2])] {
        assert!(
            (got - want).abs() < 1e-9,
            "net force {got:.9} vs the pull alone {want:.9}; every other term must \
             still cancel, and the pull must be the whole remainder"
        );
    }
}

#[test]
fn pull_with_a_static_target_conserves_energy() {
    const STEPS: usize = 100_000;
    let (mut sys, ff_base) = scenario::chain5();
    let mut ff = ff_base;
    ff.pull = pulled(
        0,
        [sys.pos_x[0] + 1.5, sys.pos_y[0] + 1.0, sys.pos_z[0] - 0.7],
    );

    let energies = integrator::initialize(&mut sys, &ff);
    assert!(
        energies.pull > 0.0,
        "the static target must actually stretch the spring"
    );

    let mut monitor = EnergyMonitor::new(&sys, energies.total(), STEPS);
    for step in 0..STEPS {
        let energies = integrator::step(&mut sys, &ff, 1e-3);
        monitor.update(step, &sys, energies.total());
    }

    let s = monitor.summary();
    assert!(s.is_finite(), "blew up at {:?}", s.first_nonfinite);
    assert!(
        s.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}; a static target makes V_pull time independent, so the \
         total energy is conserved even though the net force and the angular momentum \
         are not",
        s.max_abs_drift
    );
    assert!(
        s.secular_drift() < 3e-5,
        "secular drift = {:.3e}",
        s.secular_drift()
    );
}

#[test]
fn an_empty_pull_matches_the_golden_trajectory() {
    const GOLDEN_TOLERANCE: Real = 1e-7;
    const GOLDEN: [(u64, u64, u64); 5] = [
        (0x402AB6837970BE4E, 0x40098C9873558125, 0xC02213779CA658FF),
        (0x402AA46CE9868146, 0x3FE6DA08A8E88730, 0xC018D091BB861E9C),
        (0x4030DDF2148C8828, 0x3FF84DC48BDD5630, 0xC01429E5D064A36E),
        (0x4030F9571E31956A, 0x401072DE4D0F1DE7, 0xC01F0ECD57530092),
        (0x402C1469954558C6, 0x4016FA6B56256D7E, 0xC0179B0867EEC1CE),
    ];

    let (mut sys, ff) = scenario::chain5();
    assert!(!ff.pull.is_active(), "chain5 must not come pre-grabbed");

    chaperone_sim::thermostat::sample_initial_velocities(&mut sys, 0.5, 90210);
    integrator::initialize(&mut sys, &ff);
    let mut bath = chaperone_sim::thermostat::Langevin::new(0.2, 0.5, 0.005, 31415);
    for _ in 0..10_000 {
        bath.step(&mut sys, &ff);
    }

    for (i, (x, y, z)) in GOLDEN.iter().enumerate() {
        let dx = sys.pos_x[i] - Real::from_bits(*x);
        let dy = sys.pos_y[i] - Real::from_bits(*y);
        let dz = sys.pos_z[i] - Real::from_bits(*z);
        let deviation = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(
            deviation < GOLDEN_TOLERANCE,
            "bead {i} drifted {deviation:.3e} from the trajectory recorded before the \
             pull term existed"
        );
    }
}

#[test]
fn an_anchor_and_a_pull_sum_in_the_net_force() {
    let (mut sys, mut ff) = scenario::chain5();

    let pull_target = [sys.pos_x[4] + 1.8, sys.pos_y[4] - 1.1, sys.pos_z[4] + 0.6];
    let anchor_target = [sys.pos_x[0] - 0.4, sys.pos_y[0] + 0.9, sys.pos_z[0] - 0.3];
    ff.pull = Pull {
        index: Some(4),
        target: pull_target,
        k: PULL_K,
    };
    ff.anchor = Pull {
        index: Some(0),
        target: anchor_target,
        k: ANCHOR_K,
    };

    sys.clear_forces();
    ff.accumulate(&mut sys);

    let (fx, fy, fz) = sys.total_force();
    let expected = [
        PULL_K * (pull_target[0] - sys.pos_x[4]) + ANCHOR_K * (anchor_target[0] - sys.pos_x[0]),
        PULL_K * (pull_target[1] - sys.pos_y[4]) + ANCHOR_K * (anchor_target[1] - sys.pos_y[0]),
        PULL_K * (pull_target[2] - sys.pos_z[4]) + ANCHOR_K * (anchor_target[2] - sys.pos_z[0]),
    ];
    for (got, want) in [(fx, expected[0]), (fy, expected[1]), (fz, expected[2])] {
        assert!(
            (got - want).abs() < 1e-9,
            "net force {got:.9} vs pull plus anchor {want:.9}; the two restraints are the \
             only external forces and every internal term must still cancel"
        );
    }
}

#[test]
fn an_anchor_holds_its_bead_against_a_pull() {
    const STEPS: usize = 40_000;
    let (mut sys, mut ff) = scenario::chain5();

    let held = [sys.pos_x[0], sys.pos_y[0], sys.pos_z[0]];
    ff.anchor = Pull {
        index: Some(0),
        target: held,
        k: ANCHOR_K,
    };
    ff.pull = Pull {
        index: Some(4),
        target: [sys.pos_x[4] + 30.0, sys.pos_y[4], sys.pos_z[4]],
        k: PULL_K,
    };

    integrator::initialize(&mut sys, &ff);
    let mut bath = chaperone_sim::thermostat::Langevin::new(0.2, 0.3, 1e-3, 4242);
    for _ in 0..STEPS {
        bath.step(&mut sys, &ff);
    }

    let held_drift = ((sys.pos_x[0] - held[0]).powi(2)
        + (sys.pos_y[0] - held[1]).powi(2)
        + (sys.pos_z[0] - held[2]).powi(2))
    .sqrt();

    assert!(
        held_drift < 5.0,
        "the anchored bead wandered {held_drift:.3} A; the clamp caps the restoring force \
         at k * 5 A, so it must not exceed that reach"
    );
    assert!(
        sys.distance(0, 4) > 6.0,
        "d(0,4) = {:.3}; the pull should stretch the chain against the anchor",
        sys.distance(0, 4)
    );
}
