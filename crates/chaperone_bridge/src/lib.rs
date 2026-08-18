use core::pin::Pin;
use cxx_qt::CxxQtType;
use cxx_qt_lib::{QQuaternion, QString, QVector3D};

use chaperone_sim::analysis::{
    fraction_of_native_contacts, fraction_of_tertiary_contacts, CONTACT_TOLERANCE,
};
use chaperone_sim::forcefield::ForceField;
use chaperone_sim::scenario::{ANGLE_K, BOND_K, EPS, K_PHI1, K_PHI3, SIGMA};
use chaperone_sim::system::{Real, System};
use chaperone_sim::thermostat::{sample_initial_velocities, Langevin};
use chaperone_sim::{integrator, model};

const DT: Real = 0.005;
const GAMMA: Real = 0.2;
const SEED: u64 = 20260818;
const DEFAULT_TEMPERATURE: f32 = 0.5;

#[cxx_qt::bridge]
pub mod qobject {
    unsafe extern "C++" {
        include!("cxx-qt-lib/qstring.h");
        type QString = cxx_qt_lib::QString;
        include!("cxx-qt-lib/qvector3d.h");
        type QVector3D = cxx_qt_lib::QVector3D;
        include!("cxx-qt-lib/qquaternion.h");
        type QQuaternion = cxx_qt_lib::QQuaternion;
    }

    extern "RustQt" {
        #[qobject]
        #[qml_element]
        #[qproperty(i32, atom_count, cxx_name = "atomCount")]
        #[qproperty(i32, bond_count, cxx_name = "bondCount")]
        #[qproperty(f32, bounding_radius, cxx_name = "boundingRadius")]
        #[qproperty(f32, view_radius, cxx_name = "viewRadius")]
        #[qproperty(QString, status)]
        #[qproperty(bool, running)]
        #[qproperty(f32, temperature)]
        #[qproperty(i32, frame)]
        #[qproperty(f32, q)]
        #[qproperty(f32, q_tertiary, cxx_name = "qTertiary")]
        #[qproperty(f32, rg)]
        #[qproperty(f32, steps_per_second, cxx_name = "stepsPerSecond")]
        type Simulation = super::SimulationRust;

        #[qinvokable]
        #[cxx_name = "loadPdb"]
        fn load_pdb(self: Pin<&mut Simulation>, path: &QString) -> bool;

        #[qinvokable]
        fn advance(self: Pin<&mut Simulation>, steps: i32);

        #[qinvokable]
        fn reset(self: Pin<&mut Simulation>);

        #[qinvokable]
        #[cxx_name = "setBathTemperature"]
        fn set_bath_temperature(self: Pin<&mut Simulation>, temperature: f32);

        #[qinvokable]
        #[cxx_name = "positionAt"]
        fn position_at(self: &Simulation, index: i32) -> QVector3D;

        #[qinvokable]
        #[cxx_name = "hueAt"]
        fn hue_at(self: &Simulation, index: i32) -> f32;

        #[qinvokable]
        #[cxx_name = "bondMidpoint"]
        fn bond_midpoint(self: &Simulation, index: i32) -> QVector3D;

        #[qinvokable]
        #[cxx_name = "bondRotation"]
        fn bond_rotation(self: &Simulation, index: i32) -> QQuaternion;

        #[qinvokable]
        #[cxx_name = "bondLength"]
        fn bond_length(self: &Simulation, index: i32) -> f32;
    }
}

pub struct SimulationRust {
    atom_count: i32,
    bond_count: i32,
    bounding_radius: f32,
    view_radius: f32,
    status: QString,
    running: bool,
    temperature: f32,
    frame: i32,
    q: f32,
    q_tertiary: f32,
    rg: f32,
    steps_per_second: f32,

    sys: Option<System>,
    ff: Option<ForceField>,
    bath: Option<Langevin>,
    native: Vec<[Real; 3]>,
    positions: Vec<[f32; 3]>,
}

impl Default for SimulationRust {
    fn default() -> Self {
        SimulationRust {
            atom_count: 0,
            bond_count: 0,
            bounding_radius: 1.0,
            view_radius: 1.0,
            status: QString::default(),
            running: false,
            temperature: DEFAULT_TEMPERATURE,
            frame: 0,
            q: 0.0,
            q_tertiary: 0.0,
            rg: 0.0,
            steps_per_second: 0.0,
            sys: None,
            ff: None,
            bath: None,
            native: Vec::new(),
            positions: Vec::new(),
        }
    }
}

struct Observables {
    q: f32,
    q_tertiary: f32,
    rg: f32,
    radius: f32,
    view_radius: f32,
}

impl SimulationRust {
    fn refresh(&mut self) -> Observables {
        let (Some(sys), Some(ff)) = (self.sys.as_ref(), self.ff.as_ref()) else {
            return Observables {
                q: 0.0,
                q_tertiary: 0.0,
                rg: 0.0,
                radius: 1.0,
                view_radius: 1.0,
            };
        };

        let n = sys.n as Real;
        let mut centre = [0.0 as Real; 3];
        for i in 0..sys.n {
            centre[0] += sys.pos_x[i];
            centre[1] += sys.pos_y[i];
            centre[2] += sys.pos_z[i];
        }
        for c in centre.iter_mut() {
            *c /= n;
        }

        self.positions.clear();
        let mut radius: f32 = 0.0;
        for i in 0..sys.n {
            let p = [
                (sys.pos_x[i] - centre[0]) as f32,
                (sys.pos_y[i] - centre[1]) as f32,
                (sys.pos_z[i] - centre[2]) as f32,
            ];
            radius = radius.max((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt());
            self.positions.push(p);
        }

        // 필드를 직접 쓰면 뒤따르는 set_view_radius 가 변화 없음으로 판단해
        // NOTIFY 를 생략한다. 읽기만 하고 쓰기는 setter 에 맡긴다.
        let smoothed = if self.view_radius <= 1.0 {
            radius
        } else {
            self.view_radius * 0.97 + radius * 0.03
        };

        Observables {
            q: fraction_of_native_contacts(sys, &ff.native, CONTACT_TOLERANCE).unwrap_or(0.0)
                as f32,
            q_tertiary: fraction_of_tertiary_contacts(sys, &ff.native, CONTACT_TOLERANCE)
                .unwrap_or(0.0) as f32,
            rg: sys.radius_of_gyration() as f32,
            radius,
            view_radius: smoothed,
        }
    }
}

impl qobject::Simulation {
    fn publish(mut self: Pin<&mut Self>, o: Observables) {
        self.as_mut().set_q(o.q);
        self.as_mut().set_q_tertiary(o.q_tertiary);
        self.as_mut().set_rg(o.rg);
        self.as_mut().set_bounding_radius(o.radius);
        self.as_mut().set_view_radius(o.view_radius);
    }

    pub fn load_pdb(mut self: Pin<&mut Self>, path: &QString) -> bool {
        let path = path.to_string();

        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => {
                self.as_mut()
                    .set_status(QString::from(&format!("cannot read {path}: {error}")));
                return false;
            }
        };

        let structure = match chaperone_pdb::parse(&text, None) {
            Ok(structure) => structure,
            Err(error) => {
                self.as_mut()
                    .set_status(QString::from(&format!("{path}: {error}")));
                return false;
            }
        };

        let temperature = *self.as_ref().temperature() as Real;
        let (mut sys, ff) =
            model::go_model(&structure, BOND_K, ANGLE_K, K_PHI1, K_PHI3, EPS, SIGMA);
        sample_initial_velocities(&mut sys, temperature, SEED);
        integrator::initialize(&mut sys, &ff);

        let native: Vec<[Real; 3]> = (0..sys.n)
            .map(|i| [sys.pos_x[i], sys.pos_y[i], sys.pos_z[i]])
            .collect();
        let count = sys.n as i32;
        let contacts = ff.native.len();
        let chain = structure.chain as char;

        {
            let rust = self.as_mut().rust_mut().get_mut();
            rust.sys = Some(sys);
            rust.ff = Some(ff);
            rust.bath = Some(Langevin::new(GAMMA, temperature, DT, SEED));
            rust.native = native;
        }

        let observables = self.as_mut().rust_mut().get_mut().refresh();
        self.as_mut().publish(observables);
        self.as_mut().set_atom_count(count);
        self.as_mut().set_bond_count((count - 1).max(0));
        self.as_mut().set_frame(0);
        self.as_mut().set_status(QString::from(&format!(
            "{count} residues, chain {chain}, {contacts} native contacts"
        )));
        true
    }

    pub fn advance(mut self: Pin<&mut Self>, steps: i32) {
        if steps <= 0 {
            return;
        }

        let start = std::time::Instant::now();
        {
            let rust = self.as_mut().rust_mut().get_mut();
            let (Some(sys), Some(ff), Some(bath)) =
                (rust.sys.as_mut(), rust.ff.as_ref(), rust.bath.as_mut())
            else {
                return;
            };
            for _ in 0..steps {
                bath.step(sys, ff);
            }
        }
        let elapsed = start.elapsed().as_secs_f64();

        let observables = self.as_mut().rust_mut().get_mut().refresh();
        self.as_mut().publish(observables);

        let frame = *self.as_ref().frame() + 1;
        self.as_mut().set_frame(frame);
        if elapsed > 0.0 {
            self.as_mut()
                .set_steps_per_second((steps as f64 / elapsed) as f32);
        }
    }

    pub fn reset(mut self: Pin<&mut Self>) {
        {
            let rust = self.as_mut().rust_mut().get_mut();
            let Some(sys) = rust.sys.as_mut() else {
                return;
            };
            for (i, p) in rust.native.iter().enumerate() {
                sys.pos_x[i] = p[0];
                sys.pos_y[i] = p[1];
                sys.pos_z[i] = p[2];
                sys.vel_x[i] = 0.0;
                sys.vel_y[i] = 0.0;
                sys.vel_z[i] = 0.0;
            }
        }

        let observables = self.as_mut().rust_mut().get_mut().refresh();
        self.as_mut().publish(observables);
        self.as_mut().set_frame(0);
    }

    pub fn set_bath_temperature(mut self: Pin<&mut Self>, temperature: f32) {
        if temperature <= 0.0 {
            return;
        }
        self.as_mut().set_temperature(temperature);
        if let Some(bath) = self.as_mut().rust_mut().get_mut().bath.as_mut() {
            bath.set_temperature(temperature as Real);
        }
    }

    pub fn position_at(&self, index: i32) -> QVector3D {
        match usize::try_from(index)
            .ok()
            .and_then(|i| self.positions.get(i))
        {
            Some(p) => QVector3D::new(p[0], p[1], p[2]),
            None => QVector3D::new(0.0, 0.0, 0.0),
        }
    }

    fn bond_ends(&self, index: i32) -> Option<([f32; 3], [f32; 3])> {
        let i = usize::try_from(index).ok()?;
        Some((*self.positions.get(i)?, *self.positions.get(i + 1)?))
    }

    pub fn bond_midpoint(&self, index: i32) -> QVector3D {
        match self.bond_ends(index) {
            Some((a, b)) => QVector3D::new(
                0.5 * (a[0] + b[0]),
                0.5 * (a[1] + b[1]),
                0.5 * (a[2] + b[2]),
            ),
            None => QVector3D::new(0.0, 0.0, 0.0),
        }
    }

    pub fn bond_length(&self, index: i32) -> f32 {
        match self.bond_ends(index) {
            Some((a, b)) => {
                let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            }
            None => 0.0,
        }
    }

    // #Cylinder 는 +Y 를 따라 서 있다. +Y 에서 결합 방향으로 가는 최단호 회전을 준다.
    pub fn bond_rotation(&self, index: i32) -> QQuaternion {
        let identity = QQuaternion::new(1.0, &QVector3D::new(0.0, 0.0, 0.0));
        let Some((a, b)) = self.bond_ends(index) else {
            return identity;
        };

        let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if len < 1e-6 {
            return identity;
        }
        let d = [d[0] / len, d[1] / len, d[2] / len];

        // w = 1 + Y·d,  v = Y × d = (dz, 0, -dx)
        let w = 1.0 + d[1];
        if w < 1e-6 {
            // 정확히 -Y 방향: 축이 사라지므로 X 축 180도로 고정한다
            return QQuaternion::new(0.0, &QVector3D::new(1.0, 0.0, 0.0));
        }
        let v = [d[2], 0.0, -d[0]];
        let norm = (w * w + v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        QQuaternion::new(
            w / norm,
            &QVector3D::new(v[0] / norm, v[1] / norm, v[2] / norm),
        )
    }

    pub fn hue_at(&self, index: i32) -> f32 {
        let n = self.positions.len().max(1) as f32;
        (index.max(0) as f32 / n) * 0.75
    }
}
