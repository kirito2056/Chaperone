use chaperone_sim::analysis::{
    fraction_of_local_contacts, fraction_of_native_contacts, fraction_of_tertiary_contacts,
    EnergyMonitor, CONTACT_TOLERANCE,
};
use chaperone_sim::folding::{run_stage, Equipartition, FoldingMonitor, Stage};
use chaperone_sim::forcefield::angle::Angles;
use chaperone_sim::forcefield::dihedral::Dihedrals;
use chaperone_sim::forcefield::native::NativeContacts;
use chaperone_sim::forcefield::ForceField;
use chaperone_sim::model;
use chaperone_sim::scenario::{self, ANGLE_K, BOND_K, EPS, K_PHI1, K_PHI3, R0, SIGMA};
use chaperone_sim::system::{Real, System, PI};
use chaperone_sim::{integrator, thermostat};

const DT: Real = 0.005;
const GAMMA: Real = 0.2;
const PDB_PATH: &str = "../../data/pdb/1UBQ.pdb";

// ---------------------------------------------------------------------------
// A. build guards — the three layouts of version_11 §2
// ---------------------------------------------------------------------------

#[test]
#[should_panic(expected = "is outside (0, PI)")]
fn extended_chain_rejects_a_collinear_layout() {
    Angles::from_chain(&scenario::extended_chain(6, PI));
}

#[test]
#[should_panic(expected = "nearly collinear")]
fn extended_chain_rejects_a_nearly_collinear_layout() {
    let sys = scenario::extended_chain(6, PI - 5e-4);
    Angles::from_chain(&sys);
    Dihedrals::from_chain(&sys);
}

#[test]
fn extended_chain_accepts_a_zigzag_layout() {
    let sys = scenario::extended_chain(6, 120.0_f64.to_radians());
    assert_eq!(Angles::from_chain(&sys).len(), 4);
    assert_eq!(Dihedrals::from_chain(&sys).len(), 3);
    for i in 1..sys.n {
        assert!((sys.distance(i - 1, i) - R0).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// B. radius of gyration — duplicate-consistency against the pdb crate
// ---------------------------------------------------------------------------

fn synthetic_pdb(n: usize) -> String {
    let mut lines = Vec::new();
    for i in 0..n {
        let x = 3.7 * i as Real;
        let y = if i % 2 == 0 { 0.0 } else { 0.866 };
        let seq = i as i32 + 1;
        for (serial, name, element, dz) in
            [(2 * seq - 1, " N  ", "N", 1.0), (2 * seq, " CA ", "C", 0.0)]
        {
            lines.push(format!(
                "ATOM  {serial:>5} {name:<4} ALA A{seq:>4}    \
                 {x:8.3}{y:8.3}{dz:8.3}  1.00  0.00          {element:>2}"
            ));
        }
    }
    lines.join("\n")
}

#[test]
fn system_and_structure_agree_on_the_radius_of_gyration() {
    let text = synthetic_pdb(12);
    let structure = chaperone_pdb::parse(&text, None).expect("fixture parses");
    let sys = model::system_from_structure(&structure);

    let from_structure = structure.radius_of_gyration();
    let from_system = sys.radius_of_gyration();
    assert!(
        (from_structure - from_system).abs() < 1e-12,
        "Rg {from_structure:.12} vs {from_system:.12}; the System copy dropped the \
         centre-of-mass subtraction or a residue"
    );
    assert!(from_system > 0.0);
}

// ---------------------------------------------------------------------------
// C. initial velocities — the reserved step = 0 stream
// ---------------------------------------------------------------------------

const T_INIT: Real = 0.7;

fn two_mass_system(n: usize) -> System {
    let mut sys = System::new(n);
    for i in 0..n {
        sys.mass[i] = if i % 2 == 0 { 1.0 } else { 3.0 };
    }
    sys
}

#[test]
fn initial_velocities_match_the_bath_temperature() {
    let mut sys = two_mass_system(20_000);
    thermostat::sample_initial_velocities(&mut sys, T_INIT, 4242);

    for (class, mass) in [(0usize, 1.0), (1, 3.0)] {
        let mut sum = 0.0;
        let mut count = 0usize;
        for i in (class..sys.n).step_by(2) {
            sum += sys.vel_x[i] * sys.vel_x[i]
                + sys.vel_y[i] * sys.vel_y[i]
                + sys.vel_z[i] * sys.vel_z[i];
            count += 3;
        }
        let measured = sum / count as Real;
        let expected = T_INIT / mass;
        assert!(
            (measured / expected - 1.0).abs() < 0.03,
            "mass {mass}: <v^2> per component {measured:.6} vs T/m {expected:.6}; \
             a missing sqrt or a shared m[0] would land here"
        );
    }
}

#[test]
fn initial_velocities_come_from_the_reserved_zeroth_stream() {
    let seed = 90210;
    let mut sys = two_mass_system(64);
    thermostat::sample_initial_velocities(&mut sys, T_INIT, seed);

    let noise = chaperone_sim::rng::Noise::new(seed);
    for i in 0..sys.n {
        let scale = (T_INIT / sys.mass[i]).sqrt();
        assert_eq!(sys.vel_x[i], scale * noise.gaussian(0, i, 0));
        assert_eq!(sys.vel_y[i], scale * noise.gaussian(0, i, 1));
        assert_eq!(sys.vel_z[i], scale * noise.gaussian(0, i, 2));
    }

    let mut again = two_mass_system(64);
    thermostat::sample_initial_velocities(&mut again, T_INIT, seed);
    assert_eq!(again.vel_x, sys.vel_x, "same seed must replay bit for bit");
}

// ---------------------------------------------------------------------------
// D. FoldingMonitor — deterministic synthetic trajectory
// ---------------------------------------------------------------------------

#[test]
fn folding_monitor_summarises_a_synthetic_trajectory() {
    let n = 16;
    let mut contacts = NativeContacts::new();
    contacts.push(0, 3, 4.0);
    contacts.push(0, 15, 4.0);

    let place = |sys: &mut System, d_local: Real, d_tertiary: Real| {
        for i in 0..n {
            sys.pos_x[i] = 0.0;
            sys.pos_y[i] = 0.0;
            sys.pos_z[i] = 0.0;
        }
        sys.pos_x[3] = d_local;
        sys.pos_x[15] = d_tertiary;
    };

    let mut sys = System::new(n);
    place(&mut sys, 1.0, 1.0);
    let mut monitor = FoldingMonitor::new(&sys, &contacts, 4);

    place(&mut sys, 1.0, 99.0);
    monitor.update(0, &sys, &contacts);
    place(&mut sys, 99.0, 99.0);
    monitor.update(1, &sys, &contacts);
    place(&mut sys, 1.0, 99.0);
    monitor.update(2, &sys, &contacts);
    place(&mut sys, 1.0, 1.0);
    monitor.update(3, &sys, &contacts);

    let s = monitor.summary();
    assert_eq!(s.q_initial, 1.0);
    assert_eq!(s.q_min, 0.0);
    assert_eq!(s.q_max, 1.0);
    assert_eq!(s.q_final, 1.0);
    assert_eq!(s.q_local_final, 1.0);
    assert_eq!(s.q_tertiary_final, 1.0);
    assert!(s.is_finite());
    assert!(
        s.rg_max > s.rg_min,
        "Rg must move when a bead is thrown out to 99"
    );
    assert!(
        s.temperature_mean.is_finite(),
        "the tail window (steps >= 2) must contribute"
    );
}

// ---------------------------------------------------------------------------
// E. equipartition — the formula and the bonded wiring, without a PDB
// ---------------------------------------------------------------------------

fn bonded_only_chain(n: usize) -> (System, ForceField) {
    let sys = scenario::extended_chain(n, 120.0_f64.to_radians());
    let mut ff = ForceField::new(BOND_K, ANGLE_K, K_PHI1, K_PHI3, EPS, SIGMA);
    for i in 1..n {
        ff.bonds
            .push((i - 1) as u32, i as u32, sys.distance(i - 1, i));
    }
    ff.angles = Angles::from_chain(&sys);
    (sys, ff)
}

#[test]
fn equipartition_recovers_the_bond_and_angle_temperatures() {
    const T: Real = 0.3;
    let (mut sys, ff) = bonded_only_chain(40);
    thermostat::sample_initial_velocities(&mut sys, T, 5150);

    let empty = NativeContacts::new();
    let stage = Stage {
        name: "equipartition",
        temperature: T,
        gamma: GAMMA,
        steps: 60_000,
        seed: 77,
    };
    let mut equi = Equipartition::new();
    run_stage(&mut sys, &ff, &empty, &stage, DT, 20, |s| {
        equi.update(s.sys, &ff)
    });

    let bond = equi.bond_mean_square() / Equipartition::expected_bond_mean_square(T, BOND_K);
    let angle = equi.angle_mean_square() / Equipartition::expected_angle_mean_square(T, ANGLE_K);

    assert!(
        (0.9..=1.1).contains(&bond),
        "<(r-r0)^2> ratio {bond:.4}, expected ~1 for an uncoupled harmonic bond"
    );
    assert!(
        (0.9..=1.1).contains(&angle),
        "<(th-th0)^2> ratio {angle:.4}, expected ~1 when only bond and angle act"
    );
}

// ---------------------------------------------------------------------------
// F. quench — the invariant is a non-increasing potential at reset instants
// ---------------------------------------------------------------------------

#[test]
fn quench_never_increases_the_potential_at_reset_instants() {
    let (mut sys, ff) = scenario::chain5();
    thermostat::sample_initial_velocities(&mut sys, 0.5, 31337);
    integrator::initialize(&mut sys, &ff);

    let mut previous: Option<Real> = None;
    for cycle in 0..40 {
        for v in sys
            .vel_x
            .iter_mut()
            .chain(sys.vel_y.iter_mut())
            .chain(sys.vel_z.iter_mut())
        {
            *v = 0.0;
        }

        let mut energies = integrator::initialize(&mut sys, &ff);
        let e_pot = energies.potential();
        if let Some(before) = previous {
            assert!(
                e_pot <= before + 1e-6,
                "cycle {cycle}: E_pot rose from {before:.9} to {e_pot:.9} across a quench"
            );
        }
        previous = Some(e_pot);

        for _ in 0..200 {
            energies = integrator::step(&mut sys, &ff, 1e-3);
        }
        assert!(energies.total().is_finite());
    }
}

// ---------------------------------------------------------------------------
// G. 1UBQ — the whole pipeline
// ---------------------------------------------------------------------------

fn ubiquitin(bond_k: Real, angle_k: Real, eps: Real, sigma: Real) -> (System, ForceField) {
    let text = std::fs::read_to_string(PDB_PATH).unwrap_or_else(|e| {
        panic!("{PDB_PATH}: {e}\nfetch it with: curl -o data/pdb/1UBQ.pdb https://files.rcsb.org/download/1UBQ.pdb")
    });
    let structure = chaperone_pdb::parse(&text, None).expect("1UBQ parses");
    model::go_model(&structure, bond_k, angle_k, K_PHI1, K_PHI3, eps, sigma)
}

fn ubiquitin_reference() -> (System, ForceField) {
    ubiquitin(BOND_K, ANGLE_K, EPS, SIGMA)
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn ubiquitin_conserves_energy_and_angular_momentum_in_nve() {
    let steps = 100_000;
    let (mut sys, ff) = ubiquitin_reference();
    thermostat::sample_initial_velocities(&mut sys, 0.05, 20260816);

    let energies = integrator::initialize(&mut sys, &ff);
    let mut monitor = EnergyMonitor::new(&sys, energies.total(), steps);
    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, DT);
        monitor.update(step, &sys, energies.total());
    }

    let s = monitor.summary();
    assert!(s.is_finite(), "blew up at step {:?}", s.first_nonfinite);
    assert!(
        s.max_abs_drift < 5e-4,
        "max |drift| = {:.3e}, expected < 5e-4",
        s.max_abs_drift
    );
    assert!(
        s.secular_drift() < 1e-6,
        "secular drift = {:.3e}, expected < 1e-6",
        s.secular_drift()
    );
    assert!(
        s.max_angular_momentum_drift < 1e-9,
        "max |L - L0| = {:.3e}, expected < 1e-9",
        s.max_angular_momentum_drift
    );
}

fn hold(temperature: Real, index: u64) -> (chaperone_sim::folding::FoldingSummary, Real, Real) {
    let (mut sys, ff) = ubiquitin_reference();
    thermostat::sample_initial_velocities(&mut sys, temperature, 1000 + index);

    let stage = Stage {
        name: "hold",
        temperature,
        gamma: GAMMA,
        steps: 20_000,
        seed: 100 + index,
    };
    let mut equi = Equipartition::new();
    let summary = run_stage(&mut sys, &ff, &ff.native, &stage, DT, 50, |s| {
        equi.update(s.sys, &ff)
    });

    (
        summary,
        equi.bond_mean_square() / Equipartition::expected_bond_mean_square(temperature, BOND_K),
        equi.angle_mean_square() / Equipartition::expected_angle_mean_square(temperature, ANGLE_K),
    )
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn the_native_state_stays_folded_across_a_temperature_ladder() {
    for (index, temperature) in [0.1, 0.3, 0.5].into_iter().enumerate() {
        let (s, bond, angle) = hold(temperature, index as u64);

        assert!(s.is_finite(), "T = {temperature}: blew up");
        assert!(
            s.q_final > 0.95,
            "T = {temperature}: Q fell to {:.4}; the native state is a static force minimum \
             (ubiquitin.rs), so losing it at a low rung points at the dynamics, not at Tf",
            s.q_final
        );
        assert!(
            (s.temperature_mean / temperature - 1.0).abs() < 0.15,
            "T = {temperature}: measured <T> = {:.4}; a gamma <-> temperature swap in \
             Langevin::new lands here",
            s.temperature_mean
        );

        assert!(
            (0.9..=1.1).contains(&bond),
            "T = {temperature}: <(r-r0)^2> ratio {bond:.4} against the intended BOND_K; \
             a bond_k <-> angle_k swap reads about 4.6 here"
        );
        assert!(
            (0.55..=0.80).contains(&angle),
            "T = {temperature}: <(th-th0)^2> ratio {angle:.4}; the dihedral and contact \
             terms stiffen the angles to about 0.63-0.73, and an eps <-> sigma_nn swap \
             drops it to about 0.47"
        );
    }
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn equipartition_detects_swapped_force_field_arguments() {
    const T: Real = 0.3;

    let run = |bond_k, angle_k, eps, sigma| {
        let (mut sys, ff) = ubiquitin(bond_k, angle_k, eps, sigma);
        thermostat::sample_initial_velocities(&mut sys, T, 7);
        let stage = Stage {
            name: "swap",
            temperature: T,
            gamma: GAMMA,
            steps: 40_000,
            seed: 11,
        };
        let mut equi = Equipartition::new();
        run_stage(&mut sys, &ff, &ff.native, &stage, DT, 50, |s| {
            equi.update(s.sys, &ff)
        });
        (
            equi.bond_mean_square() / Equipartition::expected_bond_mean_square(T, BOND_K),
            equi.angle_mean_square() / Equipartition::expected_angle_mean_square(T, ANGLE_K),
        )
    };

    let (bond, angle) = run(ANGLE_K, BOND_K, EPS, SIGMA);
    assert!(
        !(0.9..=1.1).contains(&bond) || !(0.55..=0.80).contains(&angle),
        "a bond_k <-> angle_k swap slipped through: bond {bond:.4}, angle {angle:.4}"
    );

    let (bond, angle) = run(BOND_K, ANGLE_K, SIGMA, EPS);
    assert!(
        !(0.9..=1.1).contains(&bond) || !(0.55..=0.80).contains(&angle),
        "an eps <-> sigma_nn swap slipped through: bond {bond:.4}, angle {angle:.4}"
    );
}

fn melt() -> (System, ForceField, chaperone_sim::folding::FoldingSummary) {
    let (mut sys, ff) = ubiquitin_reference();
    thermostat::sample_initial_velocities(&mut sys, 1.8, 2001);
    let stage = Stage {
        name: "melt",
        temperature: 1.8,
        gamma: GAMMA,
        steps: 200_000,
        seed: 201,
    };
    let summary = run_stage(&mut sys, &ff, &ff.native, &stage, DT, 100, |_| {});
    (sys, ff, summary)
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn heating_expands_the_native_state_to_a_coil() {
    let (_, _, s) = melt();

    assert!(s.is_finite(), "melt blew up at {:?}", s.first_nonfinite);
    assert!(
        s.rg_final > 18.0,
        "Rg only reached {:.3} from the native 11.49; the chain did not open up",
        s.rg_final
    );
    assert!(
        s.q_tertiary_final < 0.10,
        "Q_tertiary is still {:.4}; the tertiary fold survived the melt",
        s.q_tertiary_final
    );
    assert!(
        (s.temperature_mean / 1.8 - 1.0).abs() < 0.15,
        "measured <T> = {:.4} against the 1.8 set point",
        s.temperature_mean
    );
}

#[test]
#[ignore = "requires data/pdb/1UBQ.pdb"]
fn the_melted_chain_recovers_its_tertiary_contacts() {
    let (mut sys, ff, melted) = melt();
    assert!(
        melted.q_tertiary_final < 0.10,
        "the handoff must start from an unfolded chain, got Q_tertiary {:.4}",
        melted.q_tertiary_final
    );

    let stage = Stage {
        name: "refold",
        temperature: 0.6,
        gamma: GAMMA,
        steps: 400_000,
        seed: 301,
    };
    assert_ne!(stage.seed, 201, "each stage needs its own noise stream");

    let q_before = fraction_of_tertiary_contacts(&sys, &ff.native, CONTACT_TOLERANCE).unwrap();
    let s = run_stage(&mut sys, &ff, &ff.native, &stage, DT, 100, |_| {});

    assert!(s.is_finite(), "refold blew up at {:?}", s.first_nonfinite);
    assert!(
        s.q_initial < 0.30,
        "Q started at {:.4}; the contact map must come from the PDB, never from a snapshot",
        s.q_initial
    );
    assert!(
        s.q_tertiary_final > 0.50,
        "Q_tertiary went {q_before:.4} -> {:.4}; the chain did not recover its fold",
        s.q_tertiary_final
    );
    assert!(
        s.rg_final < 14.0,
        "Rg stayed at {:.3}; the chain did not collapse",
        s.rg_final
    );

    let q = fraction_of_native_contacts(&sys, &ff.native, CONTACT_TOLERANCE).unwrap();
    let q_local = fraction_of_local_contacts(&sys, &ff.native, CONTACT_TOLERANCE).unwrap();
    assert!(q > 0.6 && q_local > 0.8, "Q {q:.4}, Q_local {q_local:.4}");
}
