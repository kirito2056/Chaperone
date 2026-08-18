import QtQuick
import QtQuick.Controls.Basic
import QtQuick3D
import QtQuick3D.Helpers
import Chaperone 1.0

Window {
    id: root
    width: 1180
    height: 800
    visible: true
    title: "Chaperone"
    color: "#0e1013"

    readonly property real beadRadius: 1.6
    property int stepsPerFrame: 120
    property real nativeRadius: 1.0

    Simulation {
        id: sim
    }

    Component.onCompleted: {
        sim.loadPdb(defaultPdbPath)
        console.warn("chaperone:", sim.status)
        root.nativeRadius = Math.max(sim.boundingRadius, 1.0)
        camera.z = root.nativeRadius * 2.8
    }

    FrameAnimation {
        running: sim.running
        onTriggered: sim.advance(root.stepsPerFrame)
    }

    Shortcut {
        sequence: "Space"
        onActivated: sim.running = !sim.running
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
                clipNear: 0.5
                clipFar: 20000
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

        Node {
            id: molecule
            readonly property real fit: root.nativeRadius / Math.max(sim.viewRadius, 1.0)
            scale: Qt.vector3d(fit, fit, fit)

            Repeater3D {
                model: sim.atomCount

                    Model {
                    required property int index

                    source: "#Sphere"
                    position: (sim.frame, sim.positionAt(index))
                    scale: Qt.vector3d(root.beadRadius / 50, root.beadRadius / 50, root.beadRadius / 50)

                    materials: PrincipledMaterial {
                        baseColor: Qt.hsva(sim.hueAt(index), 0.65, 0.95, 1.0)
                        roughness: 0.35
                    }
                }
            }
        }
    }

    Rectangle {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.margins: 14
        width: statusLabel.implicitWidth + 24
        height: statusLabel.implicitHeight + 16
        radius: 6
        color: "#c0181c22"

        Text {
            id: statusLabel
            anchors.centerIn: parent
            color: "#d7dde6"
            font.pixelSize: 13
            text: sim.status
        }
    }

    Rectangle {
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 14
        width: readout.implicitWidth + 28
        height: readout.implicitHeight + 20
        radius: 6
        color: "#c0181c22"

        Grid {
            id: readout
            anchors.centerIn: parent
            columns: 2
            columnSpacing: 16
            rowSpacing: 4

            component Key: Text {
                color: "#78828f"
                font.pixelSize: 12
            }
            component Value: Text {
                color: "#d7dde6"
                font.pixelSize: 12
                font.family: "monospace"
            }

            Key { text: "T" }
            Value { text: sim.temperature.toFixed(2) }
            Key { text: "Q" }
            Value { text: sim.q.toFixed(3) }
            Key { text: "Q tertiary" }
            Value { text: sim.qTertiary.toFixed(3) }
            Key { text: "Rg" }
            Value { text: sim.rg.toFixed(2) + " A" }
            Key { text: "steps/s" }
            Value { text: sim.stepsPerSecond.toFixed(0) }
            Key { text: "frame" }
            Value { text: sim.frame }
        }
    }

    Rectangle {
        anchors.horizontalCenter: parent.horizontalCenter
        anchors.bottom: parent.bottom
        anchors.bottomMargin: 18
        width: controls.implicitWidth + 32
        height: controls.implicitHeight + 22
        radius: 8
        color: "#d0181c22"

        Row {
            id: controls
            anchors.centerIn: parent
            spacing: 18

            Button {
                width: 84
                text: sim.running ? "pause" : "play"
                onClicked: sim.running = !sim.running
            }

            Button {
                width: 72
                text: "reset"
                onClicked: {
                    sim.running = false
                    sim.reset()
                }
            }

            Column {
                spacing: 2
                Text {
                    color: "#78828f"
                    font.pixelSize: 11
                    text: "temperature  " + sim.temperature.toFixed(2)
                }
                Slider {
                    width: 220
                    from: 0.05
                    to: 2.0
                    value: sim.temperature
                    onMoved: sim.setBathTemperature(value)
                }
            }

            Column {
                spacing: 2
                Text {
                    color: "#78828f"
                    font.pixelSize: 11
                    text: "steps / frame  " + root.stepsPerFrame
                }
                Slider {
                    width: 160
                    from: 10
                    to: 400
                    stepSize: 10
                    value: root.stepsPerFrame
                    onMoved: root.stepsPerFrame = value
                }
            }
        }
    }

    Text {
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.margins: 14
        color: "#4c5663"
        font.pixelSize: 11
        text: "space to play/pause  ·  drag to orbit  ·  wheel to zoom"
    }
}
