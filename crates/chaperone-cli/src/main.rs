use chaperone_sim::analysis::PeriodTracker;
use chaperone_sim::forcefield::bond::{self, Bonds};
use chaperone_sim::forcefield::pairlist::PairList;
use chaperone_sim::forcefield::repulsion;
use chaperone_sim::integrator;
use chaperone_sim::system::{Real, System};

const K: Real = 100.0;
const R0: Real = 3.8;
const EPS: Real = 1.0;
const SIGMA: Real = 4.0;
const DT: Real = 1e-3;
const STEPS: usize = 1_000_000;
const SAMPLE_EVERY: usize = 1000;

struct Monitor {
    e_initial: Real,
    steps: usize,
    window: usize,
    head_sum: Real,
    tail_sum: Real,
    max_abs_drift: Real,
    max_force: Real,
    max_speed: Real,
    first_nonfinite: Option<usize>,
}

impl Monitor {
    fn new(e_initial: Real, steps: usize) -> Self {
        Monitor {
            e_initial,
            steps,
            window: steps / 10,
            head_sum: 0.0,
            tail_sum: 0.0,
            max_abs_drift: 0.0,
            max_force: 0.0,
            max_speed: 0.0,
            first_nonfinite: None,
        }
    }

    fn update(&mut self, step: usize, sys: &System, e_total: Real) -> Real {
        let mut nonfinite = !e_total.is_finite();

        for i in 0..sys.n {
            let f = (sys.frc_x[i] * sys.frc_x[i]
                + sys.frc_y[i] * sys.frc_y[i]
                + sys.frc_z[i] * sys.frc_z[i])
                .sqrt();
            let v = (sys.vel_x[i] * sys.vel_x[i]
                + sys.vel_y[i] * sys.vel_y[i]
                + sys.vel_z[i] * sys.vel_z[i])
                .sqrt();

            if f.is_finite() && v.is_finite() {
                if f > self.max_force {
                    self.max_force = f;
                }
                if v > self.max_speed {
                    self.max_speed = v;
                }
            } else {
                nonfinite = true;
            }
        }

        let drift = (e_total - self.e_initial) / self.e_initial;
        if drift.is_finite() {
            if drift.abs() > self.max_abs_drift {
                self.max_abs_drift = drift.abs();
            }
        } else {
            nonfinite = true;
        }

        if nonfinite && self.first_nonfinite.is_none() {
            self.first_nonfinite = Some(step);
        }

        if step < self.window {
            self.head_sum += drift;
        } else if step >= self.steps - self.window {
            self.tail_sum += drift;
        }

        drift
    }

    fn report(&self, elapsed: std::time::Duration) {
        let head_mean = self.head_sum / self.window as Real;
        let tail_mean = self.tail_sum / self.window as Real;

        eprintln!("steps                {}", self.steps);
        eprintln!("dt                   {DT:.1e}");
        eprintln!("wall time            {:.3} s", elapsed.as_secs_f64());
        eprintln!(
            "throughput           {:.2e} steps/s",
            self.steps as Real / elapsed.as_secs_f64()
        );
        eprintln!();
        eprintln!("E_total(0)           {:.9}", self.e_initial);
        eprintln!("max |drift|          {:.3e}   (limit 5e-4)", self.max_abs_drift);
        eprintln!("head 10% mean drift  {head_mean:.3e}");
        eprintln!("tail 10% mean drift  {tail_mean:.3e}");
        eprintln!(
            "secular drift        {:.3e}   (limit 1e-6)",
            (tail_mean - head_mean).abs()
        );
        eprintln!();
        eprintln!("max |F|              {:.3e}", self.max_force);
        eprintln!("max |v|              {:.3e}", self.max_speed);
        match self.first_nonfinite {
            None => eprintln!("finite               ok"),
            Some(step) => eprintln!("finite               BLEW UP at step {step}"),
        }
    }
}

fn run_spring() {
    let mut sys = System::new(2);
    let initial_separation = 5.0;
    sys.pos_x[1] = initial_separation;

    let mut bonds = Bonds::new();
    bonds.push(0, 1, R0);

    sys.clear_forces();
    let mut e_pot = bond::accumulate(&mut sys, &bonds, K);
    let e_initial = e_pot + sys.kinetic_energy();

    let mut monitor = Monitor::new(e_initial, STEPS);
    let mut tracker = PeriodTracker::new(sys.distance(0, 1) - R0, 0.0);
    let mut r_min = initial_separation;
    let mut r_max = initial_separation;

    println!("step,time,r,e_pot,e_kin,e_total,drift");
    println!(
        "0,0.000000,{:.6},{:.9},{:.9},{:.9},{:.3e}",
        sys.distance(0, 1),
        e_pot,
        sys.kinetic_energy(),
        e_initial,
        0.0
    );

    let start = std::time::Instant::now();

    for step in 0..STEPS {
        integrator::kick_drift(&mut sys, DT);
        sys.clear_forces();
        e_pot = bond::accumulate(&mut sys, &bonds, K);
        integrator::kick(&mut sys, DT);

        let time = (step + 1) as Real * DT;
        let e_kin = sys.kinetic_energy();
        let e_total = e_pot + e_kin;
        let drift = monitor.update(step, &sys, e_total);

        let r = sys.distance(0, 1);
        r_min = r_min.min(r);
        r_max = r_max.max(r);
        tracker.push(r - R0, time);

        if (step + 1) % SAMPLE_EVERY == 0 {
            println!("{},{time:.6},{r:.6},{e_pot:.9},{e_kin:.9},{e_total:.9},{drift:.3e}", step + 1);
        }
    }

    let elapsed = start.elapsed();
    monitor.report(elapsed);

    let reduced_mass = 0.5;
    let omega = (2.0 * K / reduced_mass).sqrt();
    let theoretical_period = 2.0 * std::f64::consts::PI / omega;
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
    let pairs = PairList::all_pairs(4, 3);

    sys.clear_forces();
    let mut e_bond = bond::accumulate(&mut sys, &bonds, K);
    let mut e_rep = repulsion::accumulate(&mut sys, &pairs, EPS, SIGMA);
    let e_initial = e_bond + e_rep + sys.kinetic_energy();

    let mut monitor = Monitor::new(e_initial, STEPS);
    let mut gap_min = sys.distance(0, 3);
    let mut gap_max = gap_min;

    println!("step,time,gap,e_bond,e_rep,e_kin,e_total,drift");
    println!(
        "0,0.000000,{:.6},{e_bond:.9},{e_rep:.9},{:.9},{e_initial:.9},{:.3e}",
        sys.distance(0, 3),
        sys.kinetic_energy(),
        0.0
    );

    let start = std::time::Instant::now();

    for step in 0..STEPS {
        integrator::kick_drift(&mut sys, DT);
        sys.clear_forces();
        e_bond = bond::accumulate(&mut sys, &bonds, K);
        e_rep = repulsion::accumulate(&mut sys, &pairs, EPS, SIGMA);
        integrator::kick(&mut sys, DT);

        let time = (step + 1) as Real * DT;
        let e_kin = sys.kinetic_energy();
        let e_total = e_bond + e_rep + e_kin;
        let drift = monitor.update(step, &sys, e_total);

        let gap = sys.distance(0, 3);
        gap_min = gap_min.min(gap);
        gap_max = gap_max.max(gap);

        if (step + 1) % SAMPLE_EVERY == 0 {
            println!(
                "{},{time:.6},{gap:.6},{e_bond:.9},{e_rep:.9},{e_kin:.9},{e_total:.9},{drift:.3e}",
                step + 1
            );
        }
    }

    let elapsed = start.elapsed();
    monitor.report(elapsed);

    eprintln!();
    eprintln!("pairs                {}", pairs.len());
    eprintln!("d(0,3) range         [{gap_min:.6}, {gap_max:.6}]   (sigma = {SIGMA})");
    eprintln!(
        "closest approach     {:.3} sigma",
        gap_min / SIGMA
    );
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
