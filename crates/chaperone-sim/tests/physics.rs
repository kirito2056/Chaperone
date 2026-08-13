use chaperone_sim::analysis::PeriodTracker;
use chaperone_sim::forcefield::bond::{self, Bonds};
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

fn random_chain(seed: u64, n: usize, box_size: Real) -> (System, Bonds) {
    let mut rng = Lcg::new(seed);
    let mut sys = System::new(n);
    for i in 0..n {
        sys.pos_x[i] = rng.next_unit() * box_size;
        sys.pos_y[i] = rng.next_unit() * box_size;
        sys.pos_z[i] = rng.next_unit() * box_size;
    }

    let mut bonds = Bonds::new();
    for i in 0..n - 1 {
        bonds.push(i as u32, (i + 1) as u32, R0);
    }

    (sys, bonds)
}

fn two_bead_system(initial_separation: Real) -> (System, Bonds) {
    let mut sys = System::new(2);
    sys.pos_x[1] = initial_separation;

    let mut bonds = Bonds::new();
    bonds.push(0, 1, R0);

    (sys, bonds)
}

struct EnergyStats {
    max_abs_drift: Real,
    head_mean: Real,
    tail_mean: Real,
}

fn run_energy_stats(dt: Real, steps: usize) -> EnergyStats {
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

        let drift = (e_pot + sys.kinetic_energy() - e_initial) / e_initial;
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

#[test]
fn bond_force_matches_numerical_gradient() {
    const H: Real = 1e-5;
    const ABS_TOL: Real = 1e-6;
    const REL_TOL: Real = 1e-4;

    for seed in 0..20u64 {
        let (mut sys, bonds) = random_chain(seed, 5, 10.0);

        for b in 0..bonds.len() {
            let r = sys.distance(bonds.i[b] as usize, bonds.j[b] as usize);
            assert!(r > 0.01, "seed {seed}: degenerate pair, r = {r}");
        }

        sys.clear_forces();
        bond::accumulate(&mut sys, &bonds, K);
        let analytic: Vec<[Real; 3]> = (0..sys.n)
            .map(|i| [sys.frc_x[i], sys.frc_y[i], sys.frc_z[i]])
            .collect();

        #[allow(clippy::needless_range_loop)]
        for i in 0..sys.n {
            for dim in 0..3 {
                let slot = match dim {
                    0 => &mut sys.pos_x[i],
                    1 => &mut sys.pos_y[i],
                    _ => &mut sys.pos_z[i],
                };
                let orig = *slot;

                *slot = orig + H;
                let e_plus = bond::energy(&sys, &bonds, K);

                let slot = match dim {
                    0 => &mut sys.pos_x[i],
                    1 => &mut sys.pos_y[i],
                    _ => &mut sys.pos_z[i],
                };
                *slot = orig - H;
                let e_minus = bond::energy(&sys, &bonds, K);

                let slot = match dim {
                    0 => &mut sys.pos_x[i],
                    1 => &mut sys.pos_y[i],
                    _ => &mut sys.pos_z[i],
                };
                *slot = orig;

                let numerical = -(e_plus - e_minus) / (2.0 * H);
                let a = analytic[i][dim];
                let tol = ABS_TOL + REL_TOL * a.abs().max(numerical.abs());

                assert!(
                    (a - numerical).abs() <= tol,
                    "seed {seed} atom {i} dim {dim}: analytic {a:.12e} vs numerical \
                     {numerical:.12e} (diff {:.3e}, tol {tol:.3e})",
                    (a - numerical).abs()
                );
            }
        }
    }
}

#[test]
fn newton_third_law() {
    for seed in 0..20u64 {
        let (mut sys, bonds) = random_chain(seed, 8, 10.0);
        sys.clear_forces();
        bond::accumulate(&mut sys, &bonds, K);

        let (fx, fy, fz) = sys.total_force();
        assert!(
            fx.abs() < 1e-9 && fy.abs() < 1e-9 && fz.abs() < 1e-9,
            "seed {seed}: net force ({fx:.3e}, {fy:.3e}, {fz:.3e}) is not zero"
        );
    }
}

#[test]
fn energy_oscillation_stays_bounded() {
    let stats = run_energy_stats(1e-3, 1_000_000);
    assert!(
        stats.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        stats.max_abs_drift
    );
}

#[test]
fn energy_has_no_secular_drift() {
    let stats = run_energy_stats(1e-3, 1_000_000);
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
    let coarse = run_energy_stats(1e-3, 200_000).max_abs_drift;
    let fine = run_energy_stats(1e-4, 2_000_000).max_abs_drift;
    let ratio = coarse / fine;

    assert!(
        (50.0..=200.0).contains(&ratio),
        "drift ratio {ratio:.1} (coarse {coarse:.3e}, fine {fine:.3e}), expected ~100"
    );
}
