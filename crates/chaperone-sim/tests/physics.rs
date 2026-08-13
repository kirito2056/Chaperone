use chaperone_sim::analysis::PeriodTracker;
use chaperone_sim::forcefield::bond::{self, Bonds};
use chaperone_sim::forcefield::pairlist::PairList;
use chaperone_sim::forcefield::repulsion;
use chaperone_sim::integrator;
use chaperone_sim::system::{Real, System};

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

const K: Real = 100.0;
const R0: Real = 3.8;
const EPS: Real = 1.0;
const SIGMA: Real = 4.0;
const MIN_SEPARATION: Real = 0.8 * SIGMA;

fn coord(sys: &mut System, i: usize, dim: usize) -> &mut Real {
    match dim {
        0 => &mut sys.pos_x[i],
        1 => &mut sys.pos_y[i],
        _ => &mut sys.pos_z[i],
    }
}

fn random_chain(seed: u64, n: usize, box_size: Real, min_separation: Real) -> (System, Bonds) {
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
    for i in 0..n - 1 {
        bonds.push(i as u32, (i + 1) as u32, R0);
    }

    (sys, bonds)
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

fn two_bead_system(initial_separation: Real) -> (System, Bonds) {
    let mut sys = System::new(2);
    sys.pos_x[1] = initial_separation;

    let mut bonds = Bonds::new();
    bonds.push(0, 1, R0);

    (sys, bonds)
}

fn chain4_system() -> (System, Bonds, PairList) {
    let mut sys = System::new(4);
    let corners = [(0.0, 0.0), (R0, 0.0), (R0, R0), (0.0, R0)];
    for (i, (x, y)) in corners.iter().enumerate() {
        sys.pos_x[i] = *x;
        sys.pos_y[i] = *y;
    }

    let mut bonds = Bonds::new();
    for i in 0..3 {
        bonds.push(i as u32, (i + 1) as u32, R0);
    }

    (sys, bonds, PairList::all_pairs(4, 3))
}

struct EnergyStats {
    max_abs_drift: Real,
    head_mean: Real,
    tail_mean: Real,
}

fn run_two_bead(dt: Real, steps: usize) -> EnergyStats {
    let (mut sys, bonds) = two_bead_system(5.0);

    sys.clear_forces();
    let mut e_pot = bond::accumulate(&mut sys, &bonds, K);
    let e_initial = e_pot + sys.kinetic_energy();

    let window = steps / 10;
    let mut head_sum = 0.0;
    let mut tail_sum = 0.0;
    let mut max_abs_drift: Real = 0.0;

    for step in 0..steps {
        integrator::kick_drift(&mut sys, dt);
        sys.clear_forces();
        e_pot = bond::accumulate(&mut sys, &bonds, K);
        integrator::kick(&mut sys, dt);

        let e_total = e_pot + sys.kinetic_energy();
        assert!(e_total.is_finite(), "step {step}: E_total is not finite");

        let drift = (e_total - e_initial) / e_initial;
        max_abs_drift = max_abs_drift.max(drift.abs());

        if step < window {
            head_sum += drift;
        } else if step >= steps - window {
            tail_sum += drift;
        }
    }

    EnergyStats {
        max_abs_drift,
        head_mean: head_sum / window as Real,
        tail_mean: tail_sum / window as Real,
    }
}

fn run_chain4(dt: Real, steps: usize) -> (EnergyStats, Real) {
    let (mut sys, bonds, pairs) = chain4_system();

    sys.clear_forces();
    let mut e_pot = bond::accumulate(&mut sys, &bonds, K)
        + repulsion::accumulate(&mut sys, &pairs, EPS, SIGMA);
    let e_initial = e_pot + sys.kinetic_energy();

    let window = steps / 10;
    let mut head_sum = 0.0;
    let mut tail_sum = 0.0;
    let mut max_abs_drift: Real = 0.0;
    let mut min_gap = sys.distance(0, 3);

    for step in 0..steps {
        integrator::kick_drift(&mut sys, dt);
        sys.clear_forces();
        e_pot = bond::accumulate(&mut sys, &bonds, K)
            + repulsion::accumulate(&mut sys, &pairs, EPS, SIGMA);
        integrator::kick(&mut sys, dt);

        let e_total = e_pot + sys.kinetic_energy();
        assert!(e_total.is_finite(), "step {step}: E_total is not finite");

        let drift = (e_total - e_initial) / e_initial;
        max_abs_drift = max_abs_drift.max(drift.abs());
        min_gap = min_gap.min(sys.distance(0, 3));

        if step < window {
            head_sum += drift;
        } else if step >= steps - window {
            tail_sum += drift;
        }
    }

    (
        EnergyStats {
            max_abs_drift,
            head_mean: head_sum / window as Real,
            tail_mean: tail_sum / window as Real,
        },
        min_gap,
    )
}

#[test]
fn bond_force_matches_numerical_gradient() {
    for seed in 0..20u64 {
        let (mut sys, bonds) = random_chain(seed, 5, 12.0, MIN_SEPARATION);
        assert_numerical_gradient(
            &mut sys,
            |s| bond::accumulate(s, &bonds, K),
            |s| bond::energy(s, &bonds, K),
            &format!("bond seed {seed}"),
        );
    }
}

#[test]
fn repulsion_force_matches_numerical_gradient() {
    for seed in 0..20u64 {
        let (mut sys, _) = random_chain(seed, 5, 12.0, MIN_SEPARATION);
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
        let (mut sys, bonds) = random_chain(seed, 8, 20.0, MIN_SEPARATION);
        let pairs = PairList::all_pairs(sys.n, 3);

        sys.clear_forces();
        bond::accumulate(&mut sys, &bonds, K);
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
fn energy_oscillation_stays_bounded() {
    let stats = run_two_bead(1e-3, 1_000_000);
    assert!(
        stats.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        stats.max_abs_drift
    );
}

#[test]
fn energy_has_no_secular_drift() {
    let stats = run_two_bead(1e-3, 1_000_000);
    let secular = (stats.tail_mean - stats.head_mean).abs();
    assert!(
        secular < 1e-6,
        "secular drift = {:.3e} (head {:.3e}, tail {:.3e}), expected < 1e-6",
        secular,
        stats.head_mean,
        stats.tail_mean
    );
}

#[test]
fn energy_conserved_with_bond_and_repulsion() {
    let (stats, min_gap) = run_chain4(1e-3, 1_000_000);

    assert!(
        stats.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        stats.max_abs_drift
    );

    let secular = (stats.tail_mean - stats.head_mean).abs();
    assert!(
        secular < 1e-6,
        "secular drift = {:.3e} (head {:.3e}, tail {:.3e}), expected < 1e-6",
        secular,
        stats.head_mean,
        stats.tail_mean
    );

    assert!(
        min_gap <= R0 + 1e-9,
        "min d(0,3) = {min_gap:.6}, repulsion never engaged"
    );
}

#[test]
fn oscillation_matches_analytic_solution() {
    const DT: Real = 1e-3;
    const STEPS: usize = 1_000_000;
    const INITIAL_SEPARATION: Real = 5.0;

    let (mut sys, bonds) = two_bead_system(INITIAL_SEPARATION);

    sys.clear_forces();
    bond::accumulate(&mut sys, &bonds, K);

    let mut tracker = PeriodTracker::new(sys.distance(0, 1) - R0, 0.0);
    let mut r_min = INITIAL_SEPARATION;
    let mut r_max = INITIAL_SEPARATION;

    for step in 0..STEPS {
        integrator::kick_drift(&mut sys, DT);
        sys.clear_forces();
        bond::accumulate(&mut sys, &bonds, K);
        integrator::kick(&mut sys, DT);

        let r = sys.distance(0, 1);
        r_min = r_min.min(r);
        r_max = r_max.max(r);
        tracker.push(r - R0, (step + 1) as Real * DT);
    }

    let reduced_mass = 0.5;
    let omega = (2.0 * K / reduced_mass).sqrt();
    let theoretical = 2.0 * std::f64::consts::PI / omega;
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
fn integrator_is_second_order() {
    let coarse = run_two_bead(1e-3, 200_000).max_abs_drift;
    let fine = run_two_bead(1e-4, 2_000_000).max_abs_drift;
    let ratio = coarse / fine;

    assert!(
        (50.0..=200.0).contains(&ratio),
        "drift ratio {ratio:.1} (coarse {coarse:.3e}, fine {fine:.3e}), expected ~100"
    );
}
