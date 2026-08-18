use cxx_qt_build::{CxxQtBuilder, QmlModule};

fn main() {
    CxxQtBuilder::new_qml_module(QmlModule::new("Chaperone").version(1, 0))
        .qt_module("Quick")
        .file("src/lib.rs")
        .build();
}
