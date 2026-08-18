use core::pin::Pin;
use cxx_qt::CxxQtType;
use cxx_qt_lib::{QString, QVector3D};

#[cxx_qt::bridge]
pub mod qobject {
    unsafe extern "C++" {
        include!("cxx-qt-lib/qstring.h");
        type QString = cxx_qt_lib::QString;
        include!("cxx-qt-lib/qvector3d.h");
        type QVector3D = cxx_qt_lib::QVector3D;
    }

    extern "RustQt" {
        #[qobject]
        #[qml_element]
        #[qproperty(i32, atom_count, cxx_name = "atomCount")]
        #[qproperty(f32, bounding_radius, cxx_name = "boundingRadius")]
        #[qproperty(QString, status)]
        type Simulation = super::SimulationRust;

        #[qinvokable]
        #[cxx_name = "loadPdb"]
        fn load_pdb(self: Pin<&mut Simulation>, path: &QString) -> bool;

        #[qinvokable]
        #[cxx_name = "positionAt"]
        fn position_at(self: &Simulation, index: i32) -> QVector3D;

        #[qinvokable]
        #[cxx_name = "hueAt"]
        fn hue_at(self: &Simulation, index: i32) -> f32;
    }
}

#[derive(Default)]
pub struct SimulationRust {
    atom_count: i32,
    bounding_radius: f32,
    status: QString,
    positions: Vec<[f32; 3]>,
}

impl qobject::Simulation {
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

        let mut positions: Vec<[f32; 3]> = structure
            .residues
            .iter()
            .map(|r| [r.ca[0] as f32, r.ca[1] as f32, r.ca[2] as f32])
            .collect();

        let n = positions.len() as f32;
        let mut centre = [0.0f32; 3];
        for p in &positions {
            for d in 0..3 {
                centre[d] += p[d];
            }
        }
        for c in centre.iter_mut() {
            *c /= n;
        }
        for p in positions.iter_mut() {
            for d in 0..3 {
                p[d] -= centre[d];
            }
        }

        let radius = positions
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
            .fold(0.0f32, f32::max);
        let count = positions.len() as i32;
        let chain = structure.chain as char;

        self.as_mut().rust_mut().get_mut().positions = positions;
        self.as_mut().set_atom_count(count);
        self.as_mut().set_bounding_radius(radius);
        self.as_mut().set_status(QString::from(&format!(
            "{count} residues, chain {chain}, radius {radius:.1} A"
        )));
        true
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

    pub fn hue_at(&self, index: i32) -> f32 {
        let n = self.positions.len().max(1) as f32;
        (index.max(0) as f32 / n) * 0.75
    }
}
