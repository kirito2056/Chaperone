use chaperone_sim::ribbon::{Frame, Frames};
use chaperone_sim::system::{Real, PI};

mod common;
use common::UBIQUITIN_CA;

const DEG: Real = PI / 180.0;

fn dot(u: [Real; 3], v: [Real; 3]) -> Real {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}
fn sub(u: [Real; 3], v: [Real; 3]) -> [Real; 3] {
    [u[0] - v[0], u[1] - v[1], u[2] - v[2]]
}
fn cross(u: [Real; 3], v: [Real; 3]) -> [Real; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}
fn unit(u: [Real; 3]) -> [Real; 3] {
    let m = dot(u, u).sqrt();
    [u[0] / m, u[1] / m, u[2] / m]
}

fn ubiquitin() -> Vec<[Real; 3]> {
    UBIQUITIN_CA.to_vec()
}

fn ideal_helix(n: usize) -> Vec<[Real; 3]> {
    (0..n)
        .map(|k| {
            let t = k as Real * 100.0 * DEG;
            [2.3 * t.cos(), 2.3 * t.sin(), 1.5 * k as Real]
        })
        .collect()
}

fn ideal_strand(n: usize) -> Vec<[Real; 3]> {
    (0..n)
        .map(|k| {
            let s = if k % 2 == 0 { 1.0 } else { -1.0 };
            [0.0, 0.942 * s, 3.3 * k as Real]
        })
        .collect()
}

fn framed(points: &[[Real; 3]]) -> Frames {
    let mut f = Frames::new(points.len());
    f.update(points);
    f
}

fn assert_orthonormal(frames: &[Frame], what: &str) {
    for (k, f) in frames.iter().enumerate() {
        for (name, v) in [("tangent", f.tangent), ("width", f.width), ("face", f.face)] {
            assert!(
                v.iter().all(|c| c.is_finite()),
                "{what} bead {k}: {name} is not finite"
            );
            assert!(
                (dot(v, v).sqrt() - 1.0).abs() < 1e-12,
                "{what} bead {k}: |{name}| = {}",
                dot(v, v).sqrt()
            );
        }
        assert!(
            dot(f.tangent, f.width).abs() < 1e-12,
            "{what} bead {k}: tangent . width = {}",
            dot(f.tangent, f.width)
        );
        let expect = cross(f.tangent, f.width);
        for (e, g) in expect.iter().zip(f.face.iter()) {
            assert!(
                (e - g).abs() < 1e-12,
                "{what} bead {k}: face is not right handed"
            );
        }
    }
}

#[test]
fn every_frame_is_right_handed_and_orthonormal() {
    assert_orthonormal(framed(&ubiquitin()).frames(), "1UBQ");
    assert_orthonormal(framed(&ideal_helix(12)).frames(), "ideal helix");
    assert_orthonormal(framed(&ideal_strand(10)).frames(), "ideal strand");
}

#[test]
fn the_raw_normal_really_does_flip_along_a_strand() {
    let p = ideal_strand(10);
    let mut flips = 0;
    for k in 1..p.len() - 2 {
        let a0 = cross(sub(p[k], p[k - 1]), sub(p[k + 1], p[k]));
        let a1 = cross(sub(p[k + 1], p[k]), sub(p[k + 2], p[k + 1]));
        assert!(dot(a0, a1) < 0.0, "window {k} did not flip");
        flips += 1;
    }
    assert!(flips >= 6, "expected the whole strand to flip, got {flips}");

    let ca = ubiquitin();
    let mut negative = 0;
    let mut total = 0;
    for r in [
        2usize, 3, 4, 5, 12, 13, 14, 15, 16, 41, 42, 43, 65, 66, 67, 68,
    ] {
        let k = r - 1;
        let a0 = cross(sub(ca[k], ca[k - 1]), sub(ca[k + 1], ca[k]));
        let a1 = cross(sub(ca[k + 1], ca[k]), sub(ca[k + 2], ca[k + 1]));
        total += 1;
        if dot(a0, a1) < 0.0 {
            negative += 1;
        }
    }
    assert_eq!(
        negative, total,
        "every strand step of 1UBQ must flip before correction"
    );
}

#[test]
fn the_corrected_frame_never_twists_by_half_a_turn() {
    for (what, p) in [
        ("1UBQ", ubiquitin()),
        ("ideal helix", ideal_helix(12)),
        ("ideal strand", ideal_strand(10)),
    ] {
        let f = framed(&p);
        let mut worst: Real = 0.0;
        for k in 0..p.len() - 1 {
            let c = dot(f.frames()[k].width, f.frames()[k + 1].width).clamp(-1.0, 1.0);
            worst = worst.max(c.acos() / DEG);
        }
        assert!(
            worst < 100.0,
            "{what}: worst neighbour twist {worst:.1} deg; without the sign fix a strand \
             reaches 155 deg on average"
        );
    }
}

#[test]
fn a_strand_points_its_width_at_the_paired_strand() {
    const PAIRS: [(usize, usize, usize, usize); 4] = [
        (1, 7, 10, 17),
        (1, 7, 64, 72),
        (64, 72, 40, 45),
        (40, 45, 48, 50),
    ];
    let ca = ubiquitin();
    let f = framed(&ca);

    let mut total = 0.0;
    let mut count = 0usize;
    for (lo1, hi1, lo2, hi2) in PAIRS {
        for (x0, x1, y0, y1) in [(lo1, hi1, lo2, hi2), (lo2, hi2, lo1, hi1)] {
            for r in x0..=x1 {
                let i = r - 1;
                let mut best = (Real::INFINITY, 0usize);
                for q in y0..=y1 {
                    let d = dot(sub(ca[q - 1], ca[i]), sub(ca[q - 1], ca[i]));
                    if d < best.0 {
                        best = (d, q - 1);
                    }
                }
                let p = sub(ca[best.1], ca[i]);
                let t = f.frames()[i].tangent;
                let s = dot(p, t);
                let perp = sub(p, [t[0] * s, t[1] * s, t[2] * s]);
                total += dot(f.frames()[i].width, unit(perp)).abs();
                count += 1;
            }
        }
    }

    let mean = total / count as Real;
    assert!(
        mean >= 0.8,
        "mean |width . partner| = {mean:.3} over {count} residues; the binormal gives 0.872 \
         and the normal gives 0.336"
    );
}

#[test]
fn an_ideal_helix_wraps_its_face_around_the_axis() {
    let p = ideal_helix(12);
    let f = framed(&p);
    let axis = [0.0, 0.0, 1.0];

    for (k, q) in p.iter().enumerate().take(p.len() - 2).skip(2) {
        let tilt = dot(f.frames()[k].width, axis).abs().acos() / DEG;
        assert!(
            (tilt - 33.5).abs() < 0.5,
            "bead {k}: width sits {tilt:.2} deg off the axis, expected 33.5"
        );
        let radial = unit([q[0], q[1], 0.0]);
        assert!(
            dot(f.frames()[k].face, radial).abs() > 0.999,
            "bead {k}: the face must look along the radius"
        );
    }
}

#[test]
fn frames_rotate_with_the_structure() {
    let ca = ubiquitin();
    let before = framed(&ca);

    let (s, c) = (0.7_f64.sin(), 0.7_f64.cos());
    let rot = |v: [Real; 3]| [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]];
    let moved: Vec<[Real; 3]> = ca
        .iter()
        .map(|p| {
            let q = rot(*p);
            [q[0] + 137.0, q[1] - 42.0, q[2] + 9.5]
        })
        .collect();
    let after = framed(&moved);

    for k in 0..ca.len() {
        let want = rot(before.frames()[k].width);
        let got = after.frames()[k].width;
        for (a, b) in want.iter().zip(got.iter()) {
            assert!(
                (a - b).abs() < 1e-9,
                "bead {k}: frames must be equivariant, not invariant"
            );
        }
    }
}

#[test]
fn a_mirror_flips_the_frame_by_one_consistent_sign() {
    let ca = ubiquitin();
    let before = framed(&ca);
    let mirrored: Vec<[Real; 3]> = ca.iter().map(|p| [p[0], p[1], -p[2]]).collect();
    let after = framed(&mirrored);

    let m = |v: [Real; 3]| [v[0], v[1], -v[2]];
    let mut sign: Option<Real> = None;
    for k in 1..ca.len() - 1 {
        let d = dot(after.frames()[k].width, m(before.frames()[k].width));
        assert!(
            (d.abs() - 1.0).abs() < 1e-9,
            "bead {k}: mirrored width must stay parallel, |dot| = {}",
            d.abs()
        );
        let s = d.signum();
        match sign {
            None => sign = Some(s),
            Some(first) => assert!(
                (first - s).abs() < 1e-9,
                "bead {k}: the mirror sign must be global, not per bead"
            ),
        }
    }
}

#[test]
fn degenerate_geometry_still_yields_a_frame() {
    let mut collinear = ideal_helix(12);
    collinear[5] = [
        0.5 * (collinear[4][0] + collinear[6][0]),
        0.5 * (collinear[4][1] + collinear[6][1]),
        0.5 * (collinear[4][2] + collinear[6][2]),
    ];
    assert_orthonormal(framed(&collinear).frames(), "collinear triple");

    let mut coincident = ideal_helix(12);
    coincident[7] = coincident[6];
    assert_orthonormal(framed(&coincident).frames(), "coincident beads");

    let straight: Vec<[Real; 3]> = (0..8).map(|k| [0.0, 0.0, 3.8 * k as Real]).collect();
    assert_orthonormal(framed(&straight).frames(), "straight line");

    for n in 0..3 {
        assert_orthonormal(framed(&ideal_helix(n)).frames(), "short chain");
    }
}

#[test]
fn interpolated_frames_stay_orthonormal_and_smooth() {
    let residues = ideal_helix(12);
    let base = framed(&residues);

    const SUB: usize = 3;
    let mut points = Vec::new();
    for i in 0..residues.len() - 1 {
        for s in 0..SUB {
            let u = s as Real / SUB as Real;
            points.push([
                residues[i][0] * (1.0 - u) + residues[i + 1][0] * u,
                residues[i][1] * (1.0 - u) + residues[i + 1][1] * u,
                residues[i][2] * (1.0 - u) + residues[i + 1][2] * u,
            ]);
        }
    }
    points.push(residues[residues.len() - 1]);

    let mut fine = Frames::new(points.len());
    fine.interpolate(base.frames(), &points, SUB);
    assert_orthonormal(fine.frames(), "interpolated");

    let mut worst: Real = 0.0;
    for k in 0..points.len() - 1 {
        let c = dot(fine.frames()[k].width, fine.frames()[k + 1].width).clamp(-1.0, 1.0);
        worst = worst.max(c.acos() / DEG);
    }
    assert!(
        worst < 40.0,
        "interpolated twist {worst:.1} deg per point; subdividing by {SUB} must reduce the \
         per-residue 50 deg"
    );
}

#[test]
fn frames_do_not_depend_on_the_scale_of_the_coordinates() {
    let ca = ubiquitin();
    let base = framed(&ca);

    for factor in [1e-4, 1e4] {
        let scaled: Vec<[Real; 3]> = ca
            .iter()
            .map(|p| [p[0] * factor, p[1] * factor, p[2] * factor])
            .collect();
        let other = framed(&scaled);
        for k in 0..ca.len() {
            for (a, b) in base.frames()[k]
                .width
                .iter()
                .zip(other.frames()[k].width.iter())
            {
                assert!(
                    (a - b).abs() < 1e-9,
                    "bead {k} at scale {factor} (angstrom to micrometre): a degeneracy threshold \
                     must be relative, not absolute"
                );
            }
        }
    }
}

#[test]
fn interpolation_survives_opposed_source_frames() {
    const SUB: usize = 3;
    let source = [
        Frame {
            tangent: [0.0, 0.0, 1.0],
            width: [1.0, 0.0, 0.0],
            face: [0.0, 1.0, 0.0],
        },
        Frame {
            tangent: [0.0, 0.0, 1.0],
            width: [-0.9397, 0.3420, 0.0],
            face: [-0.3420, -0.9397, 0.0],
        },
    ];
    assert!(
        dot(source[0].width, source[1].width) < 0.0,
        "this test is pointless unless the two frames really do oppose"
    );

    let points: Vec<[Real; 3]> = (0..=SUB).map(|k| [0.0, 0.0, 3.8 * k as Real]).collect();
    let mut fine = Frames::new(points.len());
    fine.interpolate(&source, &points, SUB);
    assert_orthonormal(fine.frames(), "opposed sources");

    let mut worst: Real = 0.0;
    for k in 0..points.len() - 1 {
        let c = dot(fine.frames()[k].width, fine.frames()[k + 1].width).clamp(-1.0, 1.0);
        worst = worst.max(c.acos() / DEG);
    }
    assert!(
        worst < 30.0,
        "opposed source frames must be matched before blending; twist {worst:.1} deg"
    );
}
