use core::pin::Pin;
use cxx_qt::CxxQtType;
use cxx_qt_lib::{QQuaternion, QString, QVector3D};

use chaperone_sim::analysis::{
    fraction_of_native_contacts, fraction_of_tertiary_contacts, CONTACT_TOLERANCE,
};
use chaperone_sim::forcefield::ForceField;
use chaperone_sim::scenario::{ANCHOR_K, PULL_K};
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
        #[qproperty(i32, grabbed_index, cxx_name = "grabbedIndex")]
        #[qproperty(f32, pull_force, cxx_name = "pullForce")]
        #[qproperty(f32, pull_extension, cxx_name = "pullExtension")]
        #[qproperty(i32, anchored_index, cxx_name = "anchoredIndex")]
        #[qproperty(f32, pull_coordinate, cxx_name = "pullCoordinate")]
        #[qproperty(i32, trace_length, cxx_name = "traceLength")]
        #[qproperty(f32, trace_min_coordinate, cxx_name = "traceMinCoordinate")]
        #[qproperty(f32, trace_max_coordinate, cxx_name = "traceMaxCoordinate")]
        #[qproperty(f32, trace_max_force, cxx_name = "traceMaxForce")]
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

        #[qinvokable]
        fn grab(self: Pin<&mut Simulation>, index: i32) -> bool;

        #[qinvokable]
        #[cxx_name = "dragTo"]
        fn drag_to(self: Pin<&mut Simulation>, x: f32, y: f32, z: f32);

        #[qinvokable]
        fn release(self: Pin<&mut Simulation>);

        #[qinvokable]
        #[cxx_name = "toggleAnchor"]
        fn toggle_anchor(self: Pin<&mut Simulation>, index: i32) -> bool;

        #[qinvokable]
        #[cxx_name = "traceCoordinateAt"]
        fn trace_coordinate_at(self: &Simulation, index: i32) -> f32;

        #[qinvokable]
        #[cxx_name = "traceForceAt"]
        fn trace_force_at(self: &Simulation, index: i32) -> f32;

        #[qinvokable]
        #[cxx_name = "clearTrace"]
        fn clear_trace(self: Pin<&mut Simulation>);

        #[qinvokable]
        #[cxx_name = "saveTrace"]
        fn save_trace(self: Pin<&mut Simulation>, path: &QString) -> bool;

        #[qinvokable]
        #[cxx_name = "pullMidpoint"]
        fn pull_midpoint(self: &Simulation) -> QVector3D;

        #[qinvokable]
        #[cxx_name = "pullRotation"]
        fn pull_rotation(self: &Simulation) -> QQuaternion;
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
    grabbed_index: i32,
    pull_force: f32,
    pull_extension: f32,
    anchored_index: i32,
    pull_coordinate: f32,
    trace_length: i32,
    trace_min_coordinate: f32,
    trace_max_coordinate: f32,
    trace_max_force: f32,

    sys: Option<System>,
    ff: Option<ForceField>,
    bath: Option<Langevin>,
    native: Vec<[Real; 3]>,
    positions: Vec<[f32; 3]>,
    centre: [Real; 3],
    anchor_local: [f32; 3],
    hold_local: [f32; 3],
    grab_origin: [f32; 3],
    trace: Vec<[f32; 4]>,
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
            grabbed_index: -1,
            pull_force: 0.0,
            pull_extension: 0.0,
            anchored_index: -1,
            pull_coordinate: 0.0,
            trace_length: 0,
            trace_min_coordinate: 0.0,
            trace_max_coordinate: 0.0,
            trace_max_force: 0.0,
            sys: None,
            ff: None,
            bath: None,
            native: Vec::new(),
            positions: Vec::new(),
            centre: [0.0; 3],
            anchor_local: [0.0; 3],
            hold_local: [0.0; 3],
            grab_origin: [0.0; 3],
            trace: Vec::new(),
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
        if self.sys.is_none() || self.ff.is_none() {
            return Observables {
                q: 0.0,
                q_tertiary: 0.0,
                rg: 0.0,
                radius: 1.0,
                view_radius: 1.0,
            };
        }

        // 1단계: 무게중심과 표시 좌표. sys 만 읽는다.
        let (centre, radius) = {
            let sys = self.sys.as_ref().expect("checked above");
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
            (centre, radius)
        };
        self.centre = centre;

        // 2단계: 앵커는 표시 좌표에 산다. 이 프레임의 무게중심으로 시뮬 좌표를 다시 만든다.
        let anchor = self.anchor_local;
        let hold = self.hold_local;
        if let Some(ff) = self.ff.as_mut() {
            if ff.pull.is_active() {
                ff.pull.target = [
                    anchor[0] as Real + centre[0],
                    anchor[1] as Real + centre[1],
                    anchor[2] as Real + centre[2],
                ];
            }
            if ff.anchor.is_active() {
                ff.anchor.target = [
                    hold[0] as Real + centre[0],
                    hold[1] as Real + centre[1],
                    hold[2] as Real + centre[2],
                ];
            }
        }

        // 필드를 직접 쓰면 뒤따르는 set_view_radius 가 변화 없음으로 판단해
        // NOTIFY 를 생략한다. 읽기만 하고 쓰기는 setter 에 맡긴다.
        let smoothed = if self.view_radius <= 1.0 {
            radius
        } else {
            self.view_radius * 0.97 + radius * 0.03
        };

        // 3단계: 관측량.
        let sys = self.sys.as_ref().expect("checked above");
        let ff = self.ff.as_ref().expect("checked above");
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
        self.as_mut().set_grabbed_index(-1);
        self.as_mut().set_anchored_index(-1);
        self.as_mut().set_pull_force(0.0);
        self.as_mut().set_pull_extension(0.0);
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
        self.as_mut().publish_pull();
        self.as_mut().record_trace();

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

    pub fn bond_rotation(&self, index: i32) -> QQuaternion {
        match self.bond_ends(index) {
            Some((a, b)) => arc_rotation(a, b),
            None => QQuaternion::new(1.0, &QVector3D::new(0.0, 0.0, 0.0)),
        }
    }

    pub fn pull_midpoint(&self) -> QVector3D {
        match self.grabbed_position() {
            Some(a) => QVector3D::new(
                0.5 * (a[0] + self.anchor_local[0]),
                0.5 * (a[1] + self.anchor_local[1]),
                0.5 * (a[2] + self.anchor_local[2]),
            ),
            None => QVector3D::new(0.0, 0.0, 0.0),
        }
    }

    pub fn pull_rotation(&self) -> QQuaternion {
        match self.grabbed_position() {
            Some(a) => arc_rotation(a, self.anchor_local),
            None => QQuaternion::new(1.0, &QVector3D::new(0.0, 0.0, 0.0)),
        }
    }

    fn grabbed_position(&self) -> Option<[f32; 3]> {
        let i = usize::try_from(self.grabbed_index).ok()?;
        self.positions.get(i).copied()
    }

    pub fn grab(mut self: Pin<&mut Self>, index: i32) -> bool {
        let Some(position) = usize::try_from(index)
            .ok()
            .and_then(|i| self.positions.get(i).copied())
        else {
            return false;
        };

        {
            let rust = self.as_mut().rust_mut().get_mut();
            rust.anchor_local = position;
            rust.grab_origin = position;
            let centre = rust.centre;
            if let Some(ff) = rust.ff.as_mut() {
                ff.pull.index = Some(index as u32);
                ff.pull.k = PULL_K;
                ff.pull.target = [
                    position[0] as Real + centre[0],
                    position[1] as Real + centre[1],
                    position[2] as Real + centre[2],
                ];
            } else {
                return false;
            }
        }

        self.as_mut().set_grabbed_index(index);
        self.publish_pull();
        true
    }

    // 일시정지 중에도 값이 바뀌어야 QML 바인딩이 살아난다. frame 은 안 늘어난다.
    pub fn drag_to(mut self: Pin<&mut Self>, x: f32, y: f32, z: f32) {
        if *self.as_ref().grabbed_index() < 0 {
            return;
        }
        {
            let rust = self.as_mut().rust_mut().get_mut();
            rust.anchor_local = [x, y, z];
            let centre = rust.centre;
            if let Some(ff) = rust.ff.as_mut() {
                ff.pull.target = [
                    x as Real + centre[0],
                    y as Real + centre[1],
                    z as Real + centre[2],
                ];
            }
        }
        self.publish_pull();
    }

    // 신장 좌표. 고정점이 있으면 끝-끝 거리(AFM 이 재는 것), 없으면 잡은 자리에서의 변위.
    fn coordinate(&self) -> f32 {
        let Some(grabbed) = self.grabbed_position() else {
            return 0.0;
        };
        let other = usize::try_from(self.anchored_index)
            .ok()
            .and_then(|i| self.positions.get(i).copied())
            .unwrap_or(self.grab_origin);
        let d = [
            grabbed[0] - other[0],
            grabbed[1] - other[1],
            grabbed[2] - other[2],
        ];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }

    pub fn trace_coordinate_at(&self, index: i32) -> f32 {
        usize::try_from(index)
            .ok()
            .and_then(|i| self.trace.get(i))
            .map(|s| s[0])
            .unwrap_or(0.0)
    }

    pub fn trace_force_at(&self, index: i32) -> f32 {
        usize::try_from(index)
            .ok()
            .and_then(|i| self.trace.get(i))
            .map(|s| s[1])
            .unwrap_or(0.0)
    }

    pub fn clear_trace(mut self: Pin<&mut Self>) {
        self.as_mut().rust_mut().get_mut().trace.clear();
        self.as_mut().set_trace_length(0);
        self.as_mut().set_trace_min_coordinate(0.0);
        self.as_mut().set_trace_max_coordinate(0.0);
        self.as_mut().set_trace_max_force(0.0);
    }

    pub fn save_trace(mut self: Pin<&mut Self>, path: &QString) -> bool {
        use std::fmt::Write as _;
        let path = path.to_string();
        let mut out = String::from("coordinate,force,q,q_tertiary\n");
        for s in &self.trace {
            let _ = writeln!(out, "{:.4},{:.4},{:.4},{:.4}", s[0], s[1], s[2], s[3]);
        }
        match std::fs::write(&path, out) {
            Ok(()) => {
                let n = self.trace.len();
                self.as_mut()
                    .set_status(QString::from(&format!("{n} samples -> {path}")));
                true
            }
            Err(error) => {
                self.as_mut()
                    .set_status(QString::from(&format!("cannot write {path}: {error}")));
                false
            }
        }
    }

    pub fn toggle_anchor(mut self: Pin<&mut Self>, index: i32) -> bool {
        if *self.as_ref().anchored_index() == index {
            if let Some(ff) = self.as_mut().rust_mut().get_mut().ff.as_mut() {
                ff.anchor.release();
            }
            self.as_mut().set_anchored_index(-1);
            return true;
        }

        let Some(position) = usize::try_from(index)
            .ok()
            .and_then(|i| self.positions.get(i).copied())
        else {
            return false;
        };

        {
            let rust = self.as_mut().rust_mut().get_mut();
            rust.hold_local = position;
            let centre = rust.centre;
            let Some(ff) = rust.ff.as_mut() else {
                return false;
            };
            ff.anchor.index = Some(index as u32);
            ff.anchor.k = ANCHOR_K;
            ff.anchor.target = [
                position[0] as Real + centre[0],
                position[1] as Real + centre[1],
                position[2] as Real + centre[2],
            ];
        }

        self.as_mut().set_anchored_index(index);
        true
    }

    pub fn release(mut self: Pin<&mut Self>) {
        if let Some(ff) = self.as_mut().rust_mut().get_mut().ff.as_mut() {
            ff.pull.release();
        }
        self.as_mut().set_grabbed_index(-1);
        self.as_mut().set_pull_force(0.0);
        self.as_mut().set_pull_extension(0.0);
    }

    fn publish_pull(mut self: Pin<&mut Self>) {
        let (force, extension) = {
            let rust = self.as_ref().get_ref();
            match (rust.sys.as_ref(), rust.ff.as_ref()) {
                (Some(sys), Some(ff)) => (
                    ff.pull.force_magnitude(sys) as f32,
                    ff.pull.extension(sys) as f32,
                ),
                _ => (0.0, 0.0),
            }
        };
        self.as_mut().set_pull_force(force);
        self.as_mut().set_pull_extension(extension);
        let coordinate = self.as_ref().get_ref().coordinate();
        self.as_mut().set_pull_coordinate(coordinate);
    }

    fn record_trace(mut self: Pin<&mut Self>) {
        const CAP: usize = 200_000;
        if *self.as_ref().grabbed_index() < 0 {
            return;
        }
        let sample = [
            *self.as_ref().pull_coordinate(),
            *self.as_ref().pull_force(),
            *self.as_ref().q(),
            *self.as_ref().q_tertiary(),
        ];
        let rust = self.as_mut().rust_mut().get_mut();
        if rust.trace.len() >= CAP {
            return;
        }
        rust.trace.push(sample);
        let n = rust.trace.len() as i32;

        let min_coordinate = if n == 1 {
            sample[0]
        } else {
            self.as_ref().trace_min_coordinate().min(sample[0])
        };
        let max_coordinate = self.as_ref().trace_max_coordinate().max(sample[0]);
        let max_force = self.as_ref().trace_max_force().max(sample[1]);
        self.as_mut().set_trace_length(n);
        self.as_mut().set_trace_min_coordinate(min_coordinate);
        self.as_mut().set_trace_max_coordinate(max_coordinate);
        self.as_mut().set_trace_max_force(max_force);
    }

    pub fn hue_at(&self, index: i32) -> f32 {
        let n = self.positions.len().max(1) as f32;
        (index.max(0) as f32 / n) * 0.75
    }
}

// #Cylinder 는 +Y 를 따라 서 있다. +Y 에서 a→b 방향으로 가는 최단호 회전.
fn arc_rotation(a: [f32; 3], b: [f32; 3]) -> QQuaternion {
    let identity = QQuaternion::new(1.0, &QVector3D::new(0.0, 0.0, 0.0));
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len < 1e-6 {
        return identity;
    }
    let d = [d[0] / len, d[1] / len, d[2] / len];

    let w = 1.0 + d[1];
    if w < 1e-6 {
        return QQuaternion::new(0.0, &QVector3D::new(1.0, 0.0, 0.0));
    }
    let v = [d[2], 0.0, -d[0]];
    let norm = (w * w + v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    QQuaternion::new(
        w / norm,
        &QVector3D::new(v[0] / norm, v[1] / norm, v[2] / norm),
    )
}
