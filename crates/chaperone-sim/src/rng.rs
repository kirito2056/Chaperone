use crate::system::Real;

const DIM_BITS: u32 = 2;
const I_BITS: u32 = 20;
const STEP_BITS: u32 = 41;

const MAX_I: u64 = 1 << I_BITS;
const MAX_STEP: u64 = 1 << STEP_BITS;

const TWO_POW_53: Real = 9007199254740992.0;

fn pack(flag: u64, step: u64, i: usize, dim: u32) -> u64 {
    debug_assert!(flag < 2);
    debug_assert!(step < MAX_STEP);
    debug_assert!((i as u64) < MAX_I);
    debug_assert!(dim < 3);
    (flag << (STEP_BITS + I_BITS + DIM_BITS))
        | (step << (I_BITS + DIM_BITS))
        | ((i as u64) << DIM_BITS)
        | dim as u64
}

fn splitmix64(counter: u64) -> u64 {
    let mut z = counter.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn uniform_open_unit(bits: u64) -> Real {
    ((bits >> 11) + 1) as Real / TWO_POW_53
}

fn box_muller(u1: Real, u2: Real) -> (Real, Real) {
    let r = (-2.0 * u1.ln()).sqrt();
    let theta = 2.0 * crate::system::PI * u2;
    (r * theta.cos(), r * theta.sin())
}

pub struct Noise {
    mixed: u64,
}

impl Noise {
    pub fn new(seed: u64) -> Self {
        Noise {
            mixed: splitmix64(seed),
        }
    }

    fn uniform(&self, flag: u64, step: u64, i: usize, dim: u32) -> Real {
        uniform_open_unit(splitmix64(self.mixed ^ pack(flag, step, i, dim)))
    }

    pub fn gaussian(&self, step: u64, i: usize, dim: u32) -> Real {
        let u1 = self.uniform(0, step, i, dim);
        let u2 = self.uniform(1, step, i, dim);
        box_muller(u1, u2).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn uniform_never_returns_zero_and_never_exceeds_one() {
        assert_eq!(uniform_open_unit(0), 1.0 / TWO_POW_53);
        assert_eq!(uniform_open_unit(u64::MAX), 1.0);
        for bits in [0, 1, 2047, 2048, 1 << 40, u64::MAX - 1, u64::MAX] {
            let u = uniform_open_unit(bits);
            assert!(u > 0.0, "u must be strictly positive, got {u} from {bits}");
            assert!(u <= 1.0, "u must not exceed one, got {u} from {bits}");
            assert!((-2.0 * u.ln()).sqrt().is_finite());
        }
    }

    #[test]
    fn pack_is_injective() {
        let mut seen = HashSet::new();
        for flag in 0..2 {
            for step in [0u64, 1, 2, 1000, MAX_STEP - 1] {
                for i in [0usize, 1, 2, 75, (MAX_I - 1) as usize] {
                    for dim in 0..3 {
                        assert!(seen.insert(pack(flag, step, i, dim)));
                    }
                }
            }
        }
    }

    #[test]
    fn the_same_key_always_gives_the_same_value() {
        let noise = Noise::new(7);
        for step in [0u64, 1, 12345] {
            for i in [0usize, 3, 75] {
                for dim in 0..3 {
                    assert_eq!(
                        noise.gaussian(step, i, dim).to_bits(),
                        noise.gaussian(step, i, dim).to_bits()
                    );
                }
            }
        }
    }

    #[test]
    fn distinct_keys_give_distinct_values() {
        let noise = Noise::new(7);
        let mut seen = HashSet::new();
        for step in 0..64u64 {
            for i in 0..64usize {
                for dim in 0..3 {
                    assert!(seen.insert(noise.gaussian(step, i, dim).to_bits()));
                }
            }
        }
    }

    // 인접 시드가 스트림을 공유하면 Tf 앙상블의 오차막대만 조용히 죽는다.
    // langevin.rs 의 어떤 테스트도 이걸 못 잡는다.
    #[test]
    fn adjacent_seeds_share_no_values() {
        for (a, b) in [(0u64, 1u64), (1, 2), (7, 8), (100, 101)] {
            let (na, nb) = (Noise::new(a), Noise::new(b));
            let mut from_a = HashSet::new();
            for step in 0..40u64 {
                for i in 0..40usize {
                    for dim in 0..3 {
                        from_a.insert(na.gaussian(step, i, dim).to_bits());
                    }
                }
            }
            for step in 0..40u64 {
                for i in 0..40usize {
                    for dim in 0..3 {
                        assert!(
                            !from_a.contains(&nb.gaussian(step, i, dim).to_bits()),
                            "seeds {a} and {b} share a stream"
                        );
                    }
                }
            }
        }
    }

    fn moments(samples: &[Real]) -> (Real, Real, Real) {
        let n = samples.len() as Real;
        let mean = samples.iter().sum::<Real>() / n;
        let m2 = samples.iter().map(|x| (x - mean).powi(2)).sum::<Real>() / n;
        let m4 = samples.iter().map(|x| (x - mean).powi(4)).sum::<Real>() / n;
        (mean, m2, m4 / (m2 * m2))
    }

    #[test]
    fn samples_are_standard_normal_through_the_fourth_moment() {
        let noise = Noise::new(20260814);
        let mut samples = Vec::with_capacity(1_000_000);
        for step in 0..10_000u64 {
            for i in 0..34usize {
                for dim in 0..3 {
                    samples.push(noise.gaussian(step, i, dim));
                }
            }
        }
        let (mean, var, kurtosis) = moments(&samples);
        let n = samples.len() as Real;
        assert!(mean.abs() < 4.0 / n.sqrt(), "mean {mean}");
        assert!((var - 1.0).abs() < 4.0 * (2.0 / n).sqrt(), "variance {var}");
        assert!(
            (kurtosis - 3.0).abs() < 4.0 * (24.0 / n).sqrt(),
            "kurtosis {kurtosis}"
        );
    }

    // (1) 정수/uniform 레벨 — 정수 연산 + 정확한 2^53 나눗셈이라 플랫폼 불변.
    //     깨지면 무조건 규약 회귀다. 값을 갱신하는 것이 아니라
    //     기존 runs/ 가 무효가 됐다는 신호로 읽는다.
    #[test]
    fn integer_layer_golden_values_are_stable() {
        assert_eq!(pack(0, 0, 0, 0), 0x0000_0000_0000_0000);
        assert_eq!(pack(1, 0, 0, 0), 0x8000_0000_0000_0000);
        assert_eq!(pack(0, 1, 0, 0), 0x0000_0000_0040_0000);
        assert_eq!(pack(0, 0, 1, 0), 0x0000_0000_0000_0004);
        assert_eq!(pack(0, 0, 0, 1), 0x0000_0000_0000_0001);
        assert_eq!(pack(0, 12345, 76, 2), 0x0000_000C_0E40_0132);

        assert_eq!(splitmix64(0), 0xE220_A839_7B1D_CDAF);
        assert_eq!(splitmix64(1), 0x910A_2DEC_8902_5CC1);
        assert_eq!(splitmix64(0xDEAD_BEEF), 0x4ADF_B90F_68C9_EB9B);

        assert_eq!(Noise::new(0).mixed, 0xE220_A839_7B1D_CDAF);
        assert_eq!(Noise::new(20260814).mixed, 0xA1CD_F6F6_A186_D52D);

        assert_eq!(uniform_open_unit(0).to_bits(), 0x3CA0_0000_0000_0000);
        assert_eq!(uniform_open_unit(u64::MAX).to_bits(), 0x3FF0_0000_0000_0000);
        assert_eq!(
            Noise::new(20260814).uniform(0, 1, 0, 0).to_bits(),
            0x3FEE_B57B_0E08_80C5
        );
    }

    // (2) gaussian 레벨 — ln·cos 는 correctly-rounded 보장이 없으므로 libm 까지 핀한다.
    //     x86_64 linux glibc 에서 채취. 깨지면 (1) 을 먼저 본다:
    //     (1) 이 살아있으면 규약이 아니라 libm 차이다.
    #[test]
    fn gaussian_layer_golden_values_are_stable() {
        let noise = Noise::new(20260814);
        assert_eq!(noise.gaussian(1, 0, 0).to_bits(), 0xBFD1_DFA1_471F_26AE);
        assert_eq!(noise.gaussian(1, 0, 1).to_bits(), 0xBFF4_AD93_7C6A_BDB1);
        assert_eq!(noise.gaussian(1, 75, 2).to_bits(), 0x3FD6_BC8C_BC9A_E5E3);
        assert_eq!(noise.gaussian(999, 12, 1).to_bits(), 0xC000_06A6_DFB4_2340);
    }
}

