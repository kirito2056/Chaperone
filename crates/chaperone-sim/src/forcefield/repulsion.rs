use crate::forcefield::pairlist::PairList;
use crate::system::{Real, System};

const MIN_R2: Real = 1e-24;

pub fn accumulate(sys: &mut System, pairs: &PairList, eps: Real, sigma: Real) -> Real {
    let sigma_sq = sigma * sigma;
    let mut e_pot = 0.0;

    for p in 0..pairs.len() {
        let i = pairs.i[p] as usize;
        let j = pairs.j[p] as usize;

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r2 = dx * dx + dy * dy + dz * dz;

        if r2 < MIN_R2 {
            continue;
        }

        let s2 = sigma_sq / r2;
        let s6 = s2 * s2 * s2;
        let s12 = s6 * s6;

        e_pot += eps * s12;

        let scale = -12.0 * eps * s12 / r2;
        sys.frc_x[i] += dx * scale;
        sys.frc_y[i] += dy * scale;
        sys.frc_z[i] += dz * scale;
        sys.frc_x[j] -= dx * scale;
        sys.frc_y[j] -= dy * scale;
        sys.frc_z[j] -= dz * scale;
    }

    e_pot
}

pub fn energy(sys: &System, pairs: &PairList, eps: Real, sigma: Real) -> Real {
    let sigma_sq = sigma * sigma;
    let mut e_pot = 0.0;

    for p in 0..pairs.len() {
        let i = pairs.i[p] as usize;
        let j = pairs.j[p] as usize;

        let dx = sys.pos_x[j] - sys.pos_x[i];
        let dy = sys.pos_y[j] - sys.pos_y[i];
        let dz = sys.pos_z[j] - sys.pos_z[i];
        let r2 = dx * dx + dy * dy + dz * dz;

        if r2 < MIN_R2 {
            continue;
        }

        let s2 = sigma_sq / r2;
        let s6 = s2 * s2 * s2;
        e_pot += eps * s6 * s6;
    }

    e_pot
}
