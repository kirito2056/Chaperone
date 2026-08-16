use chaperone_sim::analysis::EnergyMonitor;
use chaperone_sim::forcefield::ForceField;
use chaperone_sim::integrator;
use chaperone_sim::scenario::{self, R0};
use chaperone_sim::system::{Real, System};
use chaperone_sim::thermostat::{center_of_mass_velocity, instantaneous_temperature, Langevin};

const T_TEST: Real = 0.8;
const M_TEST: Real = 1.3;
const GAMMA_TEST: Real = 2.1;
const DT_TEST: Real = 0.4;
const SEED: u64 = 20260814;

fn empty_field() -> ForceField {
    ForceField::new(0.0, 0.0, 0.0, 0.0, 1.0, 4.0)
}

fn free_particles(n: usize) -> (System, ForceField) {
    let mut sys = System::new(n);
    sys.mass.fill(M_TEST);
    (sys, empty_field())
}

fn zigzag_chain5() -> (System, ForceField) {
    scenario::chain5()
}

// ---------------------------------------------------------------------------
// A. γ=0 환원
// ---------------------------------------------------------------------------

#[test]
fn langevin_at_zero_friction_tracks_velocity_verlet_over_a_short_horizon() {
    const DT: Real = 5e-3;
    const STEPS: usize = 1_000;

    let (mut reference, ff) = zigzag_chain5();
    let (mut thermostatted, _) = zigzag_chain5();

    integrator::initialize(&mut reference, &ff);
    integrator::initialize(&mut thermostatted, &ff);
    let mut langevin = Langevin::new(0.0, T_TEST, DT, SEED);

    for _ in 0..STEPS {
        integrator::step(&mut reference, &ff, DT);
        langevin.step(&mut thermostatted, &ff);
    }

    for i in 0..reference.n {
        for (a, b) in [
            (reference.pos_x[i], thermostatted.pos_x[i]),
            (reference.pos_y[i], thermostatted.pos_y[i]),
            (reference.pos_z[i], thermostatted.pos_z[i]),
            (reference.vel_x[i], thermostatted.vel_x[i]),
            (reference.vel_y[i], thermostatted.vel_y[i]),
            (reference.vel_z[i], thermostatted.vel_z[i]),
        ] {
            let scale = a.abs().max(1.0);
            assert!(
                (a - b).abs() / scale < 1e-12,
                "bead {i}: {a:.17e} vs {b:.17e}"
            );
        }
    }
}

#[test]
fn langevin_at_zero_friction_conserves_energy_and_angular_momentum() {
    const DT: Real = 5e-3;
    const STEPS: usize = 200_000;

    let (mut sys, ff) = zigzag_chain5();
    let initial = integrator::initialize(&mut sys, &ff);
    let mut monitor = EnergyMonitor::new(&sys, initial.total(), STEPS);
    let mut langevin = Langevin::new(0.0, T_TEST, DT, SEED);

    for step in 0..STEPS {
        let energies = langevin.step(&mut sys, &ff);
        monitor.update(step, &sys, energies.total());
    }

    let s = monitor.summary();
    assert!(s.is_finite(), "trajectory left the finite range");
    assert!(
        s.max_abs_drift < 5e-4,
        "energy drift {:.3e}",
        s.max_abs_drift
    );
    assert!(
        s.secular_drift() < 3e-5,
        "secular drift {:.3e}",
        s.secular_drift()
    );
    assert!(
        s.max_angular_momentum_drift < 1e-9,
        "angular momentum drift {:.3e}",
        s.max_angular_momentum_drift
    );
}

// ---------------------------------------------------------------------------
// B. 자유입자 Ornstein-Uhlenbeck
// ---------------------------------------------------------------------------

#[test]
fn free_particle_velocity_relaxes_at_the_ornstein_uhlenbeck_rate() {
    const N: usize = 40_000;
    const STEPS: usize = 8;
    let v0 = 5.0;

    let (mut sys, ff) = free_particles(N);
    sys.vel_x.fill(v0);
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, DT_TEST, SEED);

    for step in 1..=STEPS {
        langevin.step(&mut sys, &ff);
        let mean: Real = sys.vel_x.iter().sum::<Real>() / N as Real;
        let t = step as Real * DT_TEST;
        let expected = v0 * (-GAMMA_TEST * t).exp();
        let sigma = (T_TEST / M_TEST / N as Real).sqrt();
        assert!(
            (mean - expected).abs() < 4.0 * sigma,
            "step {step}: mean {mean:.6} vs {expected:.6}, tolerance {:.6}",
            4.0 * sigma
        );
    }
}

#[test]
fn free_particle_velocity_equilibrates_to_the_bath_temperature() {
    const N: usize = 40_000;
    const STEPS: usize = 40;

    let (mut sys, ff) = free_particles(N);
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, DT_TEST, SEED);
    for _ in 0..STEPS {
        langevin.step(&mut sys, &ff);
    }

    let expected = T_TEST / M_TEST;
    let measured: Real = sys
        .vel_x
        .iter()
        .zip(&sys.vel_y)
        .zip(&sys.vel_z)
        .map(|((x, y), z)| (x * x + y * y + z * z) / 3.0)
        .sum::<Real>()
        / N as Real;
    let relative_error = 4.0 * (2.0 / (3.0 * N as Real)).sqrt();
    assert!(
        (measured / expected - 1.0).abs() < relative_error,
        "<v^2> {measured:.6} vs {expected:.6}"
    );
}

// ---------------------------------------------------------------------------
// C. 속도 모멘트
// ---------------------------------------------------------------------------

#[test]
fn velocity_second_and_fourth_moments_match_maxwell_boltzmann() {
    const N: usize = 20_000;
    const BURN_IN: usize = 40;
    const SAMPLES: usize = 20;

    let (mut sys, ff) = free_particles(N);
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, DT_TEST, SEED);
    for _ in 0..BURN_IN {
        langevin.step(&mut sys, &ff);
    }

    let s2 = T_TEST / M_TEST;
    let (mut m2, mut m4) = (0.0, 0.0);
    let mut count = 0.0;
    for _ in 0..SAMPLES {
        for _ in 0..2 {
            langevin.step(&mut sys, &ff);
        }
        for i in 0..N {
            for v in [sys.vel_x[i], sys.vel_y[i], sys.vel_z[i]] {
                m2 += v * v;
                m4 += v * v * v * v;
                count += 1.0;
            }
        }
    }
    m2 /= count;
    m4 /= count;

    let independent = (N * SAMPLES) as Real;
    assert!(
        (m2 / s2 - 1.0).abs() < 4.0 * (2.0 / (3.0 * independent)).sqrt(),
        "<v^2> {m2:.6} vs {s2:.6}"
    );
    assert!(
        (m4 / (3.0 * s2 * s2) - 1.0).abs() < 4.0 * (24.0 / (3.0 * independent)).sqrt(),
        "<v^4> {m4:.6} vs {:.6}",
        3.0 * s2 * s2
    );
}

#[test]
fn every_particle_is_coupled_to_the_bath() {
    const N: usize = 76;
    const STEPS: usize = 400;

    let (mut sys, ff) = free_particles(N);
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, DT_TEST, SEED);

    let mut per_bead = vec![0.0; N];
    for _ in 0..STEPS {
        langevin.step(&mut sys, &ff);
        for (i, acc) in per_bead.iter_mut().enumerate() {
            *acc += (sys.vel_x[i] * sys.vel_x[i]
                + sys.vel_y[i] * sys.vel_y[i]
                + sys.vel_z[i] * sys.vel_z[i])
                / 3.0;
        }
    }

    let expected = T_TEST / M_TEST;
    for (i, sum) in per_bead.iter().enumerate() {
        let mean = sum / STEPS as Real;
        assert!(
            (mean / expected - 1.0).abs() < 0.5,
            "bead {i} has <v^2> {mean:.6}, expected near {expected:.6}"
        );
    }
}

// ---------------------------------------------------------------------------
// D. 조화진동자 닫힌 형태 (stiff bond — 곡률 보정을 통계오차 아래로)
// ---------------------------------------------------------------------------

const PAIRS: usize = 400;
const K_EFF: Real = 2.0 * scenario::BOND_K;
const MU: Real = M_TEST / 2.0;

fn bond_dt() -> Real {
    0.675 / (K_EFF / MU).sqrt()
}

fn bond_pairs() -> (System, ForceField) {
    let mut sys = System::new(2 * PAIRS);
    sys.mass.fill(M_TEST);
    let mut ff = ForceField::new(scenario::BOND_K, 0.0, 0.0, 0.0, 1.0, 4.0);
    for p in 0..PAIRS {
        let (a, b) = (2 * p, 2 * p + 1);
        sys.pos_x[a] = 40.0 * p as Real;
        sys.pos_x[b] = 40.0 * p as Real + R0;
        ff.bonds.push(a as u32, b as u32, R0);
    }
    (sys, ff)
}

struct BondSample {
    stretch_second: Real,
    stretch_fourth: Real,
    radial_velocity_second: Real,
    stretch_velocity_cross: Real,
    samples: Real,
}

fn sample_bonds(interval: usize, rounds: usize) -> BondSample {
    let (mut sys, ff) = bond_pairs();
    let dt = bond_dt();
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, dt, SEED);
    integrator::initialize(&mut sys, &ff);
    for _ in 0..200 {
        langevin.step(&mut sys, &ff);
    }

    let mut out = BondSample {
        stretch_second: 0.0,
        stretch_fourth: 0.0,
        radial_velocity_second: 0.0,
        stretch_velocity_cross: 0.0,
        samples: 0.0,
    };
    for _ in 0..rounds {
        for _ in 0..interval {
            langevin.step(&mut sys, &ff);
        }
        for p in 0..PAIRS {
            let (a, b) = (2 * p, 2 * p + 1);
            let d = [
                sys.pos_x[b] - sys.pos_x[a],
                sys.pos_y[b] - sys.pos_y[a],
                sys.pos_z[b] - sys.pos_z[a],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let u = [d[0] / r, d[1] / r, d[2] / r];
            let vr = u[0] * (sys.vel_x[b] - sys.vel_x[a])
                + u[1] * (sys.vel_y[b] - sys.vel_y[a])
                + u[2] * (sys.vel_z[b] - sys.vel_z[a]);
            let s = r - R0;
            out.stretch_second += s * s;
            out.stretch_fourth += s * s * s * s;
            out.radial_velocity_second += vr * vr;
            out.stretch_velocity_cross += s * vr;
            out.samples += 1.0;
        }
    }
    out.stretch_second /= out.samples;
    out.stretch_fourth /= out.samples;
    out.radial_velocity_second /= out.samples;
    out.stretch_velocity_cross /= out.samples;
    out
}

fn relaxation_interval() -> usize {
    (2.0 / (GAMMA_TEST * bond_dt())).ceil() as usize
}

#[test]
fn bond_length_variance_is_exactly_kt_over_the_effective_stiffness() {
    let s = sample_bonds(relaxation_interval(), 30);
    let expected = T_TEST / K_EFF;
    let statistical = 4.0 * (2.0 / s.samples).sqrt();
    let tolerance = statistical.max(1e-3);
    assert!(
        (s.stretch_second / expected - 1.0).abs() < tolerance,
        "<(r-r0)^2> {:.8} vs {expected:.8}, tolerance {tolerance:.4}",
        s.stretch_second
    );
}

#[test]
fn bond_velocity_variance_carries_the_dt_squared_signature() {
    let dt = bond_dt();
    let epsilon = dt * dt * K_EFF / (4.0 * MU);
    let s = sample_bonds(relaxation_interval(), 100);
    let expected = (T_TEST / MU) * (1.0 - epsilon);
    let naive = T_TEST / MU;
    let tolerance = (4.0 * (2.0 / s.samples).sqrt()).max(1e-3);

    assert!(
        (s.radial_velocity_second / expected - 1.0).abs() < tolerance,
        "<v_r^2> {:.8} vs {expected:.8}",
        s.radial_velocity_second
    );
    assert!(
        epsilon > 3.0 * tolerance,
        "dt^2 signature {epsilon:.5} is not resolvable above tolerance {tolerance:.5}"
    );
    assert!(
        (s.radial_velocity_second / naive - 1.0).abs() > 2.0 * tolerance,
        "the epsilon-free value would also pass; the signature is not being tested"
    );
}

#[test]
fn bond_length_is_gaussian_in_the_fourth_moment() {
    let s = sample_bonds(relaxation_interval(), 30);
    let variance = T_TEST / K_EFF;
    let expected = 3.0 * variance * variance;
    let tolerance = (4.0 * (24.0 / s.samples).sqrt()).max(2e-3);
    assert!(
        (s.stretch_fourth / expected - 1.0).abs() < tolerance,
        "<(r-r0)^4> {:.10} vs {expected:.10}",
        s.stretch_fourth
    );
}

#[test]
fn position_and_velocity_are_uncorrelated_in_the_stationary_state() {
    let s = sample_bonds(relaxation_interval(), 30);
    let scale = (T_TEST / K_EFF).sqrt() * (T_TEST / MU).sqrt();
    let tolerance = 4.0 * scale / s.samples.sqrt();
    assert!(
        s.stretch_velocity_cross.abs() < tolerance,
        "<(r-r0) v_r> {:.3e}, tolerance {tolerance:.3e}",
        s.stretch_velocity_cross
    );
}

#[test]
fn reported_kinetic_energy_is_measured_after_the_final_kick() {
    let (mut sys, ff) = zigzag_chain5();
    integrator::initialize(&mut sys, &ff);
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, 5e-3, SEED);
    for _ in 0..50 {
        let energies = langevin.step(&mut sys, &ff);
        assert_eq!(
            energies.kinetic.to_bits(),
            sys.kinetic_energy().to_bits(),
            "reported kinetic energy does not match the post-step state"
        );
    }
}

// ---------------------------------------------------------------------------
// E. 자유도 규약
// ---------------------------------------------------------------------------

#[test]
fn instantaneous_temperature_uses_three_n_degrees_of_freedom() {
    const N: usize = 76;
    let mut sys = System::new(N);
    sys.mass.fill(M_TEST);
    let speed = (T_TEST / M_TEST).sqrt();
    sys.vel_x.fill(speed);
    sys.vel_y.fill(speed);
    sys.vel_z.fill(speed);

    let measured = instantaneous_temperature(&sys);
    assert!(
        (measured - T_TEST).abs() < 1e-12,
        "temperature {measured:.12} vs {T_TEST:.12}"
    );

    let with_removed_modes = 2.0 * sys.kinetic_energy() / (3.0 * N as Real - 6.0);
    assert!(
        (with_removed_modes / T_TEST - 1.0) > 0.02,
        "the 3N-6 convention must be visibly different, got {with_removed_modes:.12}"
    );
}

// ---------------------------------------------------------------------------
// F. 재현성
// ---------------------------------------------------------------------------

#[test]
fn the_same_seed_reproduces_the_trajectory_bit_for_bit() {
    let run = |seed: u64| {
        let (mut sys, ff) = zigzag_chain5();
        integrator::initialize(&mut sys, &ff);
        let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, 5e-3, seed);
        for _ in 0..500 {
            langevin.step(&mut sys, &ff);
        }
        sys
    };
    let a = run(SEED);
    let b = run(SEED);
    for i in 0..a.n {
        assert_eq!(a.pos_x[i].to_bits(), b.pos_x[i].to_bits());
        assert_eq!(a.vel_z[i].to_bits(), b.vel_z[i].to_bits());
    }
}

#[test]
fn different_seeds_diverge_but_agree_on_the_temperature() {
    const N: usize = 20_000;
    const STEPS: usize = 60;

    let run = |seed: u64| {
        let (mut sys, ff) = free_particles(N);
        let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, DT_TEST, seed);
        for _ in 0..STEPS {
            langevin.step(&mut sys, &ff);
        }
        sys
    };
    let a = run(SEED);
    let b = run(SEED + 1);

    let identical = (0..N).filter(|&i| a.vel_x[i] == b.vel_x[i]).count();
    assert_eq!(
        identical, 0,
        "{identical} velocities are shared across seeds"
    );

    let temperature = |s: &System| instantaneous_temperature(s);
    let (ta, tb) = (temperature(&a), temperature(&b));
    let tolerance = 4.0 * T_TEST * (2.0 / (3.0 * N as Real)).sqrt();
    assert!(
        (ta - tb).abs() < 2.0 * tolerance,
        "temperatures {ta:.6} and {tb:.6} disagree"
    );
    assert!((ta - T_TEST).abs() < tolerance, "temperature {ta:.6}");
}

// ---------------------------------------------------------------------------
// G. 노이즈 독립성
// ---------------------------------------------------------------------------

#[test]
fn noise_is_independent_across_particles() {
    const N: usize = 76;
    const BURN_IN: usize = 40;
    const SAMPLES: usize = 4_000;

    let (mut sys, ff) = free_particles(N);
    let mut langevin = Langevin::new(GAMMA_TEST, T_TEST, DT_TEST, SEED);
    for _ in 0..BURN_IN {
        langevin.step(&mut sys, &ff);
    }

    let mut sum = 0.0;
    for _ in 0..SAMPLES {
        for _ in 0..2 {
            langevin.step(&mut sys, &ff);
        }
        let (vx, vy, vz) = center_of_mass_velocity(&sys);
        sum += (vx * vx + vy * vy + vz * vz) / 3.0;
    }
    let measured = sum / SAMPLES as Real;
    let expected = T_TEST / (M_TEST * N as Real);

    let tolerance = 4.0 * (2.0 / (3.0 * SAMPLES as Real)).sqrt();
    assert!(
        (measured / expected - 1.0).abs() < tolerance,
        "<V_com^2> {measured:.8} vs {expected:.8}; shared noise would give {:.8}",
        T_TEST / M_TEST
    );
}
