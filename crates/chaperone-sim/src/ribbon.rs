use crate::system::Real;

const DEGENERATE: Real = 1e-6;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    pub tangent: [Real; 3],
    pub width: [Real; 3],
    pub face: [Real; 3],
}

impl Frame {
    fn identity() -> Self {
        Frame {
            tangent: [0.0, 0.0, 1.0],
            width: [1.0, 0.0, 0.0],
            face: [0.0, 1.0, 0.0],
        }
    }
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

fn dot(u: [Real; 3], v: [Real; 3]) -> Real {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

fn scale(u: [Real; 3], s: Real) -> [Real; 3] {
    [u[0] * s, u[1] * s, u[2] * s]
}

fn norm2(u: [Real; 3]) -> Real {
    dot(u, u)
}

fn normalize(u: [Real; 3]) -> Option<[Real; 3]> {
    let n2 = norm2(u);
    if n2 > 0.0 && n2.is_finite() {
        Some(scale(u, 1.0 / n2.sqrt()))
    } else {
        None
    }
}

fn any_perpendicular(t: [Real; 3]) -> [Real; 3] {
    let axis = if t[0].abs() <= t[1].abs() && t[0].abs() <= t[2].abs() {
        [1.0, 0.0, 0.0]
    } else if t[1].abs() <= t[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalize(cross(t, axis)).unwrap_or([1.0, 0.0, 0.0])
}

fn orthonormalize(v: [Real; 3], t: [Real; 3]) -> Option<[Real; 3]> {
    normalize(sub(v, scale(t, dot(v, t))))
}

pub struct Frames {
    frames: Vec<Frame>,
    raw: Vec<Option<[Real; 3]>>,
    aligned: Vec<[Real; 3]>,
}

impl Frames {
    pub fn new(n: usize) -> Self {
        Frames {
            frames: vec![Frame::identity(); n],
            raw: vec![None; n],
            aligned: Vec::new(),
        }
    }

    fn resize(&mut self, n: usize) {
        self.frames.resize(n, Frame::identity());
        self.raw.resize(n, None);
    }

    pub fn frames(&self) -> &[Frame] {
        &self.frames
    }

    pub fn raw(&self) -> &[Option<[Real; 3]>] {
        &self.raw
    }

    fn tangent(points: &[[Real; 3]], k: usize) -> [Real; 3] {
        let n = points.len();
        let raw = if n < 2 {
            [0.0, 0.0, 1.0]
        } else if k == 0 {
            sub(points[1], points[0])
        } else if k + 1 == n {
            sub(points[n - 1], points[n - 2])
        } else {
            sub(points[k + 1], points[k - 1])
        };
        normalize(raw).unwrap_or([0.0, 0.0, 1.0])
    }

    pub fn update(&mut self, points: &[[Real; 3]]) {
        let n = points.len();
        if self.frames.len() != n {
            self.resize(n);
        }
        if n == 0 {
            return;
        }

        for k in 0..n {
            self.raw[k] = None;
            if k >= 1 && k + 1 < n {
                let a = sub(points[k], points[k - 1]);
                let b = sub(points[k + 1], points[k]);
                let c = cross(a, b);
                if norm2(c) > DEGENERATE * DEGENERATE * norm2(a) * norm2(b) {
                    self.raw[k] = Some(c);
                }
            }
        }

        let mut last: Option<[Real; 3]> = None;
        for k in 0..n {
            if let Some(c) = self.raw[k] {
                let c = match last {
                    Some(prev) if dot(c, prev) < 0.0 => scale(c, -1.0),
                    _ => c,
                };
                self.raw[k] = Some(c);
                last = Some(c);
            }
        }

        let mut carried: Option<[Real; 3]> = None;
        for k in 0..n {
            let t = Self::tangent(points, k);
            let width = self.raw[k]
                .or(carried)
                .and_then(|v| orthonormalize(v, t))
                .unwrap_or_else(|| any_perpendicular(t));

            carried = Some(width);
            self.frames[k] = Frame {
                tangent: t,
                width,
                face: cross(t, width),
            };
        }

        if n >= 3 {
            for (k, source) in [(0usize, 1usize), (n - 1, n - 2)] {
                let t = self.frames[k].tangent;
                let width = orthonormalize(self.frames[source].width, t)
                    .unwrap_or_else(|| any_perpendicular(t));
                self.frames[k] = Frame {
                    tangent: t,
                    width,
                    face: cross(t, width),
                };
            }
        }
    }

    pub fn interpolate(&mut self, source: &[Frame], points: &[[Real; 3]], subdivisions: usize) {
        let n = points.len();
        if self.frames.len() != n {
            self.resize(n);
        }
        if n == 0 || source.is_empty() || subdivisions == 0 {
            return;
        }

        self.aligned.clear();
        let mut prev: Option<[Real; 3]> = None;
        for f in source {
            let w = match prev {
                Some(p) if dot(f.width, p) < 0.0 => scale(f.width, -1.0),
                _ => f.width,
            };
            self.aligned.push(w);
            prev = Some(w);
        }

        for k in 0..n {
            let t = Self::tangent(points, k);
            let i = (k / subdivisions).min(source.len() - 1);
            let j = (i + 1).min(source.len() - 1);
            let u = (k % subdivisions) as Real / subdivisions as Real;

            let wi = self.aligned[i];
            let wj = self.aligned[j];
            let blend = [
                wi[0] * (1.0 - u) + wj[0] * u,
                wi[1] * (1.0 - u) + wj[1] * u,
                wi[2] * (1.0 - u) + wj[2] * u,
            ];

            let width = orthonormalize(blend, t)
                .or_else(|| orthonormalize(wi, t))
                .unwrap_or_else(|| any_perpendicular(t));

            self.raw[k] = None;
            self.frames[k] = Frame {
                tangent: t,
                width,
                face: cross(t, width),
            };
        }
    }
}
