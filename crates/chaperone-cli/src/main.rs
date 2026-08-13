use chaperone_sim::analysis::{EnergyMonitor, EnergySummary, PeriodTracker};
use chaperone_sim::integrator;
use chaperone_sim::scenario::{self, BOND_K, R0, SIGMA};
use chaperone_sim::system::{Real, PI};

const DT: Real = 1e-3;
const STEPS: usize = 1_000_000;
const SAMPLE_EVERY: usize = 1000;

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
    eprintln!("max |v|              {:.3e}", summary.max_speed);
    match summary.first_nonfinite {
        None => eprintln!("finite               ok"),
        Some(step) => eprintln!("finite               BLEW UP at step {}", step + 1),
    }
}

fn run_spring() {
    let initial_separation = 5.0;
    let (mut sys, ff) = scenario::spring(initial_separation);
    let energies = integrator::initialize(&mut sys, &ff);
    let e_initial = energies.total();

    let mut monitor = EnergyMonitor::new(e_initial, STEPS);
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

    for step in 0..STEPS {
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
    report(&monitor.summary(), STEPS, elapsed);

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

fn run_chain4() {
    let (mut sys, ff) = scenario::chain4();
    let energies = integrator::initialize(&mut sys, &ff);
    let e_initial = energies.total();

    let mut monitor = EnergyMonitor::new(e_initial, STEPS);
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

    for step in 0..STEPS {
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
    report(&monitor.summary(), STEPS, elapsed);

    eprintln!();
    eprintln!("pairs                {}", ff.repulsion_pairs.len());
    eprintln!("d(0,3) range         [{gap_min:.6}, {gap_max:.6}]   (sigma = {SIGMA})");
    eprintln!("closest approach     {:.3} sigma", gap_min / SIGMA);
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        None | Some("spring") => run_spring(),
        Some("chain4") => run_chain4(),
        Some(other) => {
            eprintln!("unknown scenario: {other}");
            eprintln!("usage: chaperone [spring|chain4]");
            std::process::exit(1);
        }
    }
}
