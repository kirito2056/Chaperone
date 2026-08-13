use chaperone_sim::analysis::PeriodTracker;
use chaperone_sim::forcefield::bond::{self, Bonds};
use chaperone_sim::integrator;
use chaperone_sim::system::{Real, System};

const K: Real = 100.0;
const R0: Real = 3.8;
const DT: Real = 1e-3;
const STEPS: usize = 1_000_000;
const SAMPLE_EVERY: usize = 1000;
const INITIAL_SEPARATION: Real = 5.0;

fn main() {
    let mut sys = System::new(2);
    sys.pos_x[1] = INITIAL_SEPARATION;

    let mut bonds = Bonds::new();
    bonds.push(0, 1, R0);

    sys.clear_forces();
    let mut e_pot = bond::accumulate(&mut sys, &bonds, K);
    let e_initial = e_pot + sys.kinetic_energy();

    println!("step,time,r,e_pot,e_kin,e_total,drift");
    println!(
        "0,0.000000,{:.6},{:.9},{:.9},{:.9},{:.3e}",
        sys.distance(0, 1),
        e_pot,
        sys.kinetic_energy(),
        e_initial,
        0.0
    );

    let window = STEPS / 10;
    let mut head_sum = 0.0;
    let mut tail_sum = 0.0;
    let mut max_abs_drift: Real = 0.0;

    let mut r_min = INITIAL_SEPARATION;
    let mut r_max = INITIAL_SEPARATION;
    let mut tracker = PeriodTracker::new(sys.distance(0, 1) - R0, 0.0);

    let start = std::time::Instant::now();

    for step in 0..STEPS {
        integrator::kick_drift(&mut sys, DT);
        sys.clear_forces();
        e_pot = bond::accumulate(&mut sys, &bonds, K);
        integrator::kick(&mut sys, DT);

        let time = (step + 1) as Real * DT;
        let r = sys.distance(0, 1);
        let e_kin = sys.kinetic_energy();
        let e_total = e_pot + e_kin;
        let drift = (e_total - e_initial) / e_initial;

        max_abs_drift = max_abs_drift.max(drift.abs());
        r_min = r_min.min(r);
        r_max = r_max.max(r);
        tracker.push(r - R0, time);

        if step < window {
            head_sum += drift;
        } else if step >= STEPS - window {
            tail_sum += drift;
        }

        if (step + 1) % SAMPLE_EVERY == 0 {
            println!(
                "{},{:.6},{:.6},{:.9},{:.9},{:.9},{:.3e}",
                step + 1,
                time,
                r,
                e_pot,
                e_kin,
                e_total,
                drift
            );
        }
    }

    let elapsed = start.elapsed();

    let head_mean = head_sum / window as Real;
    let tail_mean = tail_sum / window as Real;

    let reduced_mass = 0.5;
    let omega = (2.0 * K / reduced_mass).sqrt();
    let theoretical_period = 2.0 * std::f64::consts::PI / omega;
    let measured_period = tracker.period().unwrap_or(Real::NAN);

    eprintln!("steps                {STEPS}");
    eprintln!("dt                   {DT:.1e}");
    eprintln!("wall time            {:.3} s", elapsed.as_secs_f64());
    eprintln!(
        "throughput           {:.2e} steps/s",
        STEPS as Real / elapsed.as_secs_f64()
    );
    eprintln!();
    eprintln!("E_total(0)           {e_initial:.9}");
    eprintln!("max |drift|          {max_abs_drift:.3e}   (limit 5e-4)");
    eprintln!("head 10% mean drift  {head_mean:.3e}");
    eprintln!("tail 10% mean drift  {tail_mean:.3e}");
    eprintln!(
        "secular drift        {:.3e}   (limit 1e-6)",
        (tail_mean - head_mean).abs()
    );
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
