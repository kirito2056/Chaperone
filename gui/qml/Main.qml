import QtQuick
import QtQuick3D
import QtQuick3D.Helpers
import Chaperone 1.0

Window {
    id: root
    width: 1100
    height: 760
    visible: true
    title: "Chaperone"
    color: "#0e1013"

    readonly property real beadRadius: 1.6
    readonly property real cameraDistance: Math.max(sim.boundingRadius, 1.0) * 2.3

    Simulation {
        id: sim
    }

    Component.onCompleted: {
        if (sim.loadPdb(defaultPdbPath))
            console.warn("chaperone:", sim.status)
        else
            console.warn("chaperone:", sim.status)
    }

    View3D {
        id: view
        anchors.fill: parent

        environment: SceneEnvironment {
            clearColor: root.color
            backgroundMode: SceneEnvironment.Color
            antialiasingMode: SceneEnvironment.MSAA
            antialiasingQuality: SceneEnvironment.High
        }

        Node {
            id: pivot

            PerspectiveCamera {
                id: camera
                z: root.cameraDistance
                clipNear: 0.5
                clipFar: 10000
            }
        }

        OrbitCameraController {
            anchors.fill: parent
            origin: pivot
            camera: camera
            panEnabled: true
        }

        DirectionalLight {
            eulerRotation: Qt.vector3d(-35, -60, 0)
            brightness: 1.5
        }

        DirectionalLight {
            eulerRotation: Qt.vector3d(25, 120, 0)
            brightness: 0.5
        }

        Repeater3D {
            model: sim.atomCount

            Model {
                required property int index

                source: "#Sphere"
                position: sim.positionAt(index)
                scale: Qt.vector3d(root.beadRadius / 50, root.beadRadius / 50, root.beadRadius / 50)

                materials: PrincipledMaterial {
                    baseColor: Qt.hsva(sim.hueAt(index), 0.65, 0.95, 1.0)
                    roughness: 0.35
                    metalness: 0.0
                }
            }
        }
    }

    Rectangle {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.margins: 14
        width: label.implicitWidth + 24
        height: label.implicitHeight + 16
        radius: 6
        color: "#c0181c22"

        Text {
            id: label
            anchors.centerIn: parent
            color: "#d7dde6"
            font.pixelSize: 13
            text: sim.status === "" ? "loading..." : sim.status
        }
    }

    Text {
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.margins: 14
        color: "#5c6672"
        font.pixelSize: 12
        text: "drag to orbit  ·  wheel to zoom  ·  middle-drag to pan"
    }
}
