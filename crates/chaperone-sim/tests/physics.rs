use std::sync::OnceLock;

use chaperone_sim::analysis::{
    fraction_of_native_contacts, EnergyMonitor, EnergySummary, PeriodTracker, CONTACT_TOLERANCE,
};
use chaperone_sim::forcefield::bond::{self, Bonds};
use chaperone_sim::forcefield::native::{self, NativeContacts};
use chaperone_sim::forcefield::pairlist::PairList;
use chaperone_sim::forcefield::repulsion;
use chaperone_sim::integrator;
use chaperone_sim::scenario::{self, BOND_K, EPS, MIN_SEQUENCE_SEPARATION, R0, SIGMA};
use chaperone_sim::system::{Real, System, PI};

const MIN_SEPARATION: Real = 0.8 * SIGMA;

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
    let mut monitor = EnergyMonitor::new(energies.total(), steps);

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
    let mut monitor = EnergyMonitor::new(energies.total(), steps);

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
    max_native_abs: Real,
    max_repulsion: Real,
}

fn run_chain5(dt: Real, steps: usize) -> Chain5Result {
    let (mut sys, ff) = scenario::chain5();
    let initial = integrator::initialize(&mut sys, &ff);
    let mut monitor = EnergyMonitor::new(initial.total(), steps);

    let mut max_bond = initial.bond;
    let mut max_native_abs = initial.native.abs();
    let mut max_repulsion = initial.repulsion;

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, dt);
        monitor.update(step, &sys, energies.total());

        max_bond = max_bond.max(energies.bond);
        max_native_abs = max_native_abs.max(energies.native.abs());
        max_repulsion = max_repulsion.max(energies.repulsion);
    }

    Chain5Result {
        summary: monitor.summary(),
        max_bond,
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
fn chain5_exercises_all_three_terms() {
    let r = chain5_1m();

    assert!(
        r.max_bond > 1e-3,
        "bond energy peaked at {:.3e}, term looks inert",
        r.max_bond
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
