use chaperone_sim::analysis::{
    fraction_of_native_contacts, EnergyMonitor, EnergySummary, PeriodTracker, CONTACT_TOLERANCE,
};
use chaperone_sim::integrator;
use chaperone_sim::model;
use chaperone_sim::scenario::{self, BOND_K, R0, SIGMA};
use chaperone_sim::system::{Real, PI};

const DT: Real = 1e-3;
const DEFAULT_STEPS: usize = 1_000_000;
const SAMPLE_EVERY: usize = 1000;

fn steps() -> usize {
    std::env::args()
        .nth(2)
        .map(|s| {
            s.parse()
                .unwrap_or_else(|_| panic!("steps must be an integer, got {s:?}"))
        })
        .unwrap_or(DEFAULT_STEPS)
}

fn report(summary: &EnergySummary, steps: usize, elapsed: std::time::Duration) {
    eprintln!("steps                {steps}");
    eprintln!("dt                   {DT:.1e}");
    eprintln!("wall time            {:.3} s", elapsed.as_secs_f64());
    eprintln!(
        "throughput           {:.2e} steps/s",
        steps as Real / elapsed.as_secs_f64()
    );
    eprintln!();
    eprintln!("E_total(0)           {:.9}", summary.e_initial);
    eprintln!("drift scale          {:.9}", summary.drift_scale);
    eprintln!(
        "max |drift|          {:.3e}   (limit 5e-4)",
        summary.max_abs_drift
    );
    eprintln!("head 10% mean drift  {:.3e}", summary.head_mean);
    eprintln!("tail 10% mean drift  {:.3e}", summary.tail_mean);
    eprintln!(
        "secular drift        {:.3e}   (limit 1e-6)",
        summary.secular_drift()
    );
    eprintln!();
    eprintln!("max |F|              {:.3e}", summary.max_force);
    eprintln!(
        "max |L - L0|         {:.3e}   (limit 1e-9)",
        summary.max_angular_momentum_drift
    );
    eprintln!("max |v|              {:.3e}", summary.max_speed);
    match summary.first_nonfinite {
        None => eprintln!("finite               ok"),
        Some(step) => eprintln!("finite               BLEW UP at step {}", step + 1),
    }
}

fn run_spring() {
    let steps = steps();
    let initial_separation = 5.0;
    let (mut sys, ff) = scenario::spring(initial_separation);
    let energies = integrator::initialize(&mut sys, &ff);
    let e_initial = energies.total();

    let mut monitor = EnergyMonitor::new(&sys, e_initial, steps);
    let mut tracker = PeriodTracker::new(sys.distance(0, 1) - R0, 0.0);
    let mut r_min = initial_separation;
    let mut r_max = initial_separation;

    println!("step,time,r,e_pot,e_kin,e_total,drift");
    println!(
        "0,0.000000,{:.6},{:.9},{:.9},{e_initial:.9},{:.3e}",
        sys.distance(0, 1),
        energies.potential(),
        energies.kinetic,
        0.0
    );

    let start = std::time::Instant::now();

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, DT);
        let e_total = energies.total();
        let drift = monitor.update(step, &sys, e_total);

        let time = (step + 1) as Real * DT;
        let r = sys.distance(0, 1);
        r_min = r_min.min(r);
        r_max = r_max.max(r);
        tracker.push(r - R0, time);

        if (step + 1) % SAMPLE_EVERY == 0 {
            println!(
                "{},{time:.6},{r:.6},{:.9},{:.9},{e_total:.9},{drift:.3e}",
                step + 1,
                energies.potential(),
                energies.kinetic
            );
        }
    }

    let elapsed = start.elapsed();
    report(&monitor.summary(), steps, elapsed);

    let reduced_mass = 0.5;
    let omega = (2.0 * BOND_K / reduced_mass).sqrt();
    let theoretical_period = 2.0 * PI / omega;
    let measured_period = tracker.period().unwrap_or(Real::NAN);

    eprintln!();
    eprintln!(
        "r range              [{r_min:.6}, {r_max:.6}]   center {:.6}   (r0 = {R0})",
        (r_min + r_max) / 2.0
    );
    eprintln!("crossings            {}", tracker.crossings());
    eprintln!("period measured      {measured_period:.7}");
    eprintln!("period theoretical   {theoretical_period:.7}");
    eprintln!(
        "relative error       {:.3e}   (Verlet phase error ~{:.3e})",
        (measured_period - theoretical_period).abs() / theoretical_period,
        (omega * DT).powi(2) / 24.0
    );
}

fn run_chain4(gap: Real) {
    let steps = steps();
    let (mut sys, ff) = scenario::chain4_with_gap(gap);
    let energies = integrator::initialize(&mut sys, &ff);
    let e_initial = energies.total();

    let mut monitor = EnergyMonitor::new(&sys, e_initial, steps);
    let mut gap_min = sys.distance(0, 3);
    let mut gap_max = gap_min;

    println!("step,time,gap,e_bond,e_rep,e_kin,e_total,drift");
    println!(
        "0,0.000000,{:.6},{:.9},{:.9},{:.9},{e_initial:.9},{:.3e}",
        sys.distance(0, 3),
        energies.bond,
        energies.repulsion,
        energies.kinetic,
        0.0
    );

    let start = std::time::Instant::now();

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, DT);
        let e_total = energies.total();
        let drift = monitor.update(step, &sys, e_total);

        let time = (step + 1) as Real * DT;
        let gap = sys.distance(0, 3);
        gap_min = gap_min.min(gap);
        gap_max = gap_max.max(gap);

        if (step + 1) % SAMPLE_EVERY == 0 {
            println!(
                "{},{time:.6},{gap:.6},{:.9},{:.9},{:.9},{e_total:.9},{drift:.3e}",
                step + 1,
                energies.bond,
                energies.repulsion,
                energies.kinetic
            );
        }
    }

    let elapsed = start.elapsed();
    report(&monitor.summary(), steps, elapsed);

    eprintln!();
    eprintln!("pairs                {}", ff.repulsion_pairs.len());
    eprintln!("d(0,3) range         [{gap_min:.6}, {gap_max:.6}]   (sigma = {SIGMA})");
    eprintln!("closest approach     {:.3} sigma", gap_min / SIGMA);
}

fn run_chain5() {
    let steps = steps();
    let (mut sys, ff) = scenario::chain5();
    let initial = integrator::initialize(&mut sys, &ff);
    let e_initial = initial.total();

    let mut monitor = EnergyMonitor::new(&sys, e_initial, steps);
    let mut max_bond = initial.bond;
    let mut max_angle = initial.angle;
    let mut max_dihedral = initial.dihedral;
    let mut max_native_abs = initial.native.abs();
    let mut max_repulsion = initial.repulsion;
    let mut q_min = 1.0;

    println!("step,time,q,e_bond,e_angle,e_dih,e_native,e_rep,e_kin,e_total,drift");

    let start = std::time::Instant::now();

    for step in 0..steps {
        let energies = integrator::step(&mut sys, &ff, DT);
        let e_total = energies.total();
        let drift = monitor.update(step, &sys, e_total);

        max_bond = max_bond.max(energies.bond);
        max_angle = max_angle.max(energies.angle);
        max_dihedral = max_dihedral.max(energies.dihedral);
        max_native_abs = max_native_abs.max(energies.native.abs());
        max_repulsion = max_repulsion.max(energies.repulsion);

        let q = fraction_of_native_contacts(&sys, &ff.native, CONTACT_TOLERANCE).unwrap_or(0.0);
        q_min = q.min(q_min);

        if (step + 1) % SAMPLE_EVERY == 0 {
            println!(
                "{},{:.6},{q:.4},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{e_total:.9},{drift:.3e}",
                step + 1,
                (step + 1) as Real * DT,
                energies.bond,
                energies.angle,
                energies.dihedral,
                energies.native,
                energies.repulsion,
                energies.kinetic
            );
        }
    }

    let elapsed = start.elapsed();
    report(&monitor.summary(), steps, elapsed);

    eprintln!();
    eprintln!("dihedrals            {}", ff.dihedrals.len());
    eprintln!("native contacts      {}", ff.native.len());
    eprintln!("non-native pairs     {}", ff.repulsion_pairs.len());
    eprintln!("peak |E_bond|        {max_bond:.6}");
    eprintln!("peak |E_angle|       {max_angle:.6}");
    eprintln!("peak |E_dihedral|    {max_dihedral:.6}");
    eprintln!("peak |E_native|      {max_native_abs:.6}");
    eprintln!("peak |E_rep|         {max_repulsion:.6}");
    eprintln!("Q min                {q_min:.4}");
}

fn run_pdb(path: &str) {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            std::process::exit(1);
        }
    };

    let structure = match chaperone_pdb::parse(&text, None) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{path}: {e}");
            std::process::exit(1);
        }
    };

    let n = structure.len();
    let bonds = structure.peptide_bond_lengths();
    let bond_min = bonds.iter().cloned().fold(Real::INFINITY, Real::min);
    let bond_max = bonds.iter().cloned().fold(Real::NEG_INFINITY, Real::max);

    let heavy =
        structure.heavy_atom_contacts(model::CONTACT_CUTOFF, model::MIN_SEQUENCE_SEPARATION);
    let ca = structure.ca_contacts(7.5, model::MIN_SEQUENCE_SEPARATION);
    let plaxco = structure.heavy_atom_contacts(model::PLAXCO_CUTOFF, 1);

    let heavy_set: std::collections::HashSet<_> = heavy.iter().copied().collect();
    let ca_set: std::collections::HashSet<_> = ca.iter().copied().collect();
    let shared = heavy_set.intersection(&ca_set).count();

    let (sys, ff) = model::go_model(
        &structure,
        scenario::BOND_K,
        scenario::ANGLE_K,
        scenario::K_PHI1,
        scenario::K_PHI3,
        scenario::EPS,
        scenario::SIGMA,
    );

    let sigma_min = ff
        .native
        .sigma
        .iter()
        .cloned()
        .fold(Real::INFINITY, Real::min);
    let sigma_max = ff
        .native
        .sigma
        .iter()
        .cloned()
        .fold(Real::NEG_INFINITY, Real::max);
    let sigma_mean = ff.native.sigma.iter().sum::<Real>() / ff.native.len() as Real;

    println!("chain                {}", structure.chain as char);
    println!("residues             {n}");
    println!(
        "first / last         {} / {}",
        structure.residues[0].seq,
        structure.residues[n - 1].seq
    );
    println!("peptide CA-CA        [{bond_min:.3}, {bond_max:.3}]");
    println!("radius of gyration   {:.3}", structure.radius_of_gyration());
    println!();
    println!("contacts heavy 4.5   {}   (sep >= 3)", heavy.len());
    println!("contacts CA 7.5      {}   (sep >= 3)", ca.len());
    println!(
        "  shared / heavy-only / CA-only   {} / {} / {}",
        shared,
        heavy.len() - shared,
        ca.len() - shared
    );
    println!();
    println!(
        "model CO             {:.4}   (heavy 4.5, sep >= 3)",
        chaperone_pdb::relative_contact_order(&heavy, n).unwrap_or(Real::NAN)
    );
    println!(
        "Plaxco CO            {:.4}   (heavy 6.0, sep >= 1)",
        chaperone_pdb::relative_contact_order(&plaxco, n).unwrap_or(Real::NAN)
    );
    println!();
    println!("sigma range          [{sigma_min:.3}, {sigma_max:.3}]   mean {sigma_mean:.3}");
    println!(
        "bonds / angles       {} / {}",
        ff.bonds.len(),
        ff.angles.len()
    );
    println!("dihedrals            {}", ff.dihedrals.len());
    println!(
        "native / non-native  {} / {}",
        ff.native.len(),
        ff.repulsion_pairs.len()
    );

    let mut probe = sys;
    let mut worst: Real = 0.0;
    for (name, force) in [
        ("bond", 0usize),
        ("angle", 1),
        ("dihedral", 2),
        ("native", 3),
    ] {
        probe.clear_forces();
        match force {
            0 => {
                chaperone_sim::forcefield::bond::accumulate(&mut probe, &ff.bonds, ff.bond_k);
            }
            1 => {
                chaperone_sim::forcefield::angle::accumulate(&mut probe, &ff.angles, ff.angle_k);
            }
            2 => {
                chaperone_sim::forcefield::dihedral::accumulate(
                    &mut probe,
                    &ff.dihedrals,
                    ff.k_phi1,
                    ff.k_phi3,
                );
            }
            _ => {
                chaperone_sim::forcefield::native::accumulate(&mut probe, &ff.native, ff.eps);
            }
        }
        let f = probe.max_force();
        worst = worst.max(f);
        println!("native-state |F| {name:<9} {f:.3e}");
    }
    println!("worst bonded |F|     {worst:.3e}   (limit 1e-6)");
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        None | Some("spring") => run_spring(),
        Some("chain4") => run_chain4(R0),
        Some("chain4-open") => run_chain4(4.5),
        Some("chain5") => run_chain5(),
        Some("pdb") => match std::env::args().nth(2) {
            Some(path) => run_pdb(&path),
            None => {
                eprintln!("usage: chaperone pdb <file.pdb>");
                std::process::exit(1);
            }
        },
        Some(other) => {
            eprintln!("unknown scenario: {other}");
            eprintln!("usage: chaperone [spring|chain4|chain4-open|chain5] [steps]");
            eprintln!("       chaperone pdb <file.pdb>");
            std::process::exit(1);
        }
    }
}
