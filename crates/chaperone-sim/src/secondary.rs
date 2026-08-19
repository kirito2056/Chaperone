use crate::system::{Real, System, PI};

const DEG: Real = PI / 180.0;

const H_D3: (Real, Real) = (4.8, 5.8);
const H_D4: (Real, Real) = (5.8, 7.0);
const H_TH: (Real, Real) = (77.0 * DEG, 101.0 * DEG);
const H_TA: (Real, Real) = (30.0 * DEG, 70.0 * DEG);

const S_D2: (Real, Real) = (6.1, 7.3);
const S_D3: (Real, Real) = (9.0, 10.8);
const S_D4: (Real, Real) = (11.3, 13.5);
const S_TH: (Real, Real) = (110.0 * DEG, 138.0 * DEG);
const S_TA_LO: Real = -125.0 * DEG;
const S_TA_HI: Real = 145.0 * DEG;

pub const HELIX_RUN: usize = 5;
pub const STRAND_RUN: usize = 4;
pub const SHORT_STRAND_RUN: usize = 3;
pub const SHORT_STRAND_CONTACTS: usize = 5;
const CONTACT_MIN: Real = 4.2;
const CONTACT_MAX: Real = 5.2;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ss {
    Coil,
    Helix,
    Strand,
}

impl Ss {
    pub fn code(self) -> char {
        match self {
            Ss::Coil => 'C',
            Ss::Helix => 'H',
            Ss::Strand => 'E',
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Window {
    pub d2: Real,
    pub d3: Real,
    pub d4: Real,
    pub theta: Real,
    pub tau: Real,
}

pub fn window(sys: &System, k: usize) -> Window {
    let n = sys.n;
    let has2 = k >= 1 && k + 1 < n;
    let has3 = k >= 1 && k + 2 < n;
    let has4 = k >= 1 && k + 3 < n;
    Window {
        d2: if has2 {
            sys.distance(k - 1, k + 1)
        } else {
            Real::NAN
        },
        d3: if has3 {
            sys.distance(k - 1, k + 2)
        } else {
            Real::NAN
        },
        d4: if has4 {
            sys.distance(k - 1, k + 3)
        } else {
            Real::NAN
        },
        theta: if has2 {
            sys.angle(k - 1, k, k + 1)
        } else {
            Real::NAN
        },
        tau: if has3 {
            sys.dihedral(k - 1, k, k + 1, k + 2)
        } else {
            Real::NAN
        },
    }
}

fn inside(x: Real, r: (Real, Real)) -> bool {
    r.0 <= x && x <= r.1
}

fn strand_tau(x: Real) -> bool {
    x <= S_TA_LO || x >= S_TA_HI
}

fn mask_consecutive(src: &[bool], run: usize, out: &mut [bool]) {
    out.fill(false);
    let n = src.len();
    if run == 0 || n < run {
        return;
    }
    for q in 0..=(n - run) {
        if src[q..q + run].iter().all(|&b| b) {
            out[q..q + run].fill(true);
        }
    }
}

fn extend_region(base: &[bool], ext: &[bool], out: &mut [bool]) {
    out.copy_from_slice(base);
    let n = base.len();
    for i in 0..n {
        if base[i] {
            continue;
        }
        let touches = (i + 1 < n && base[i + 1]) || (i > 0 && base[i - 1]);
        if touches && ext[i] {
            out[i] = true;
        }
    }
}

fn regions_with_contacts(sys: &System, cand: &[bool], contacts: &mut [usize], out: &mut [bool]) {
    let n = cand.len();
    out.fill(false);
    for a in 0..n {
        contacts[a] = 0;
        if !cand[a] {
            continue;
        }
        let mut c = 0;
        for b in 0..n {
            let d = sys.distance(a, b);
            if d > CONTACT_MIN && d <= CONTACT_MAX {
                c += 1;
            }
        }
        contacts[a] = c;
    }

    let mut i = 0;
    while i < n {
        if !cand[i] {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < n && cand[j] {
            j += 1;
        }
        if contacts[i..j].iter().sum::<usize>() >= SHORT_STRAND_CONTACTS {
            out[i..j].fill(true);
        }
        i = j;
    }
}

pub struct Assignment {
    labels: Vec<Ss>,
    strict_h: Vec<bool>,
    relaxed_h: Vec<bool>,
    strict_e: Vec<bool>,
    relaxed_e: Vec<bool>,
    run_h: Vec<bool>,
    run_e: Vec<bool>,
    short_e: Vec<bool>,
    helix: Vec<bool>,
    strand: Vec<bool>,
    contacts: Vec<usize>,
}

impl Assignment {
    pub fn new(n: usize) -> Self {
        Assignment {
            labels: vec![Ss::Coil; n],
            strict_h: vec![false; n],
            relaxed_h: vec![false; n],
            strict_e: vec![false; n],
            relaxed_e: vec![false; n],
            run_h: vec![false; n],
            run_e: vec![false; n],
            short_e: vec![false; n],
            helix: vec![false; n],
            strand: vec![false; n],
            contacts: vec![0; n],
        }
    }

    fn resize(&mut self, n: usize) {
        self.labels.resize(n, Ss::Coil);
        self.strict_h.resize(n, false);
        self.relaxed_h.resize(n, false);
        self.strict_e.resize(n, false);
        self.relaxed_e.resize(n, false);
        self.run_h.resize(n, false);
        self.run_e.resize(n, false);
        self.short_e.resize(n, false);
        self.helix.resize(n, false);
        self.strand.resize(n, false);
        self.contacts.resize(n, 0);
    }

    pub fn update(&mut self, sys: &System) {
        let n = sys.n;
        if self.labels.len() != n {
            self.resize(n);
        }

        let Assignment {
            labels,
            strict_h,
            relaxed_h,
            strict_e,
            relaxed_e,
            run_h,
            run_e,
            short_e,
            helix,
            strand,
            contacts,
        } = self;

        for k in 0..n {
            let w = window(sys, k);
            strict_h[k] = (inside(w.d3, H_D3) && inside(w.d4, H_D4))
                || (inside(w.theta, H_TH) && inside(w.tau, H_TA));
            relaxed_h[k] = inside(w.d3, H_D3) || inside(w.theta, H_TH);
            strict_e[k] = (inside(w.d2, S_D2) && inside(w.d3, S_D3) && inside(w.d4, S_D4))
                || (inside(w.theta, S_TH) && strand_tau(w.tau));
            relaxed_e[k] = inside(w.d3, S_D3);
        }

        mask_consecutive(strict_h, HELIX_RUN, run_h);
        extend_region(run_h, relaxed_h, helix);

        mask_consecutive(strict_e, STRAND_RUN, run_e);
        mask_consecutive(strict_e, SHORT_STRAND_RUN, short_e);
        regions_with_contacts(sys, short_e, contacts, strand);
        for k in 0..n {
            run_e[k] |= strand[k];
        }
        extend_region(run_e, relaxed_e, strand);

        for k in 0..n {
            labels[k] = if strand[k] {
                Ss::Strand
            } else if helix[k] {
                Ss::Helix
            } else {
                Ss::Coil
            };
        }
    }

    pub fn labels(&self) -> &[Ss] {
        &self.labels
    }

    pub fn code(&self) -> String {
        self.labels.iter().map(|s| s.code()).collect()
    }

    fn fraction(&self, want: Ss) -> Real {
        if self.labels.is_empty() {
            return 0.0;
        }
        let hits = self.labels.iter().filter(|&&s| s == want).count();
        hits as Real / self.labels.len() as Real
    }

    pub fn helix_fraction(&self) -> Real {
        self.fraction(Ss::Helix)
    }

    pub fn strand_fraction(&self) -> Real {
        self.fraction(Ss::Strand)
    }
}
