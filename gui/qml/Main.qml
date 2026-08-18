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

    readonly property real beadRadius: root.showSticks ? 1.25 : 1.6
    readonly property real bondRadius: 0.55
    readonly property real tubeRadius: 0.95
    property bool showSticks: true
    property bool showTube: false
    property int stepsPerFrame: 120
    property real nativeRadius: 1.0

    Simulation {
        id: sim
    }

    Component.onCompleted: {
        sim.loadPdb(defaultPdbPath)
        console.warn("chaperone:", sim.status)
        root.nativeRadius = Math.max(sim.boundingRadius, 1.0)
        sceneCamera.z = root.nativeRadius * 2.8
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
        camera: sceneCamera

        environment: SceneEnvironment {
            clearColor: root.color
            backgroundMode: SceneEnvironment.Color
            antialiasingMode: SceneEnvironment.MSAA
            antialiasingQuality: SceneEnvironment.High
        }

        Node {
            id: pivot

            PerspectiveCamera {
                id: sceneCamera
                clipNear: 0.5
                clipFar: 20000
            }
        }

        OrbitCameraController {
            id: orbit
            anchors.fill: parent
            origin: pivot
            camera: sceneCamera
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
                model: root.showTube ? 0 : sim.atomCount

                Model {
                    required property int index
                    property int beadIndex: index
                    readonly property bool grabbed: index === sim.grabbedIndex
                    readonly property bool anchored: index === sim.anchoredIndex

                    pickable: true
                    source: "#Sphere"
                    position: (sim.frame, sim.positionAt(index))
                    scale: {
                        var r = root.beadRadius * ((grabbed || anchored) ? 1.7 : 1.0) / 50
                        return Qt.vector3d(r, r, r)
                    }

                    materials: PrincipledMaterial {
                        baseColor: grabbed
                            ? "#ffffff"
                            : anchored
                                ? "#ff9d4d"
                                : Qt.hsva(sim.hueAt(index), 0.62, 0.95, 1.0)
                        roughness: 0.35
                    }
                }
            }

            Repeater3D {
                model: (root.showSticks && !root.showTube) ? sim.bondCount : 0

                Model {
                    required property int index

                    source: "#Cylinder"
                    position: (sim.frame, sim.bondMidpoint(index))
                    rotation: (sim.frame, sim.bondRotation(index))
                    scale: (sim.frame, Qt.vector3d(root.bondRadius / 50,
                                                   sim.bondLength(index) / 100,
                                                   root.bondRadius / 50))

                    materials: PrincipledMaterial {
                        baseColor: Qt.hsva(sim.hueAt(index), 0.45, 0.70, 1.0)
                        roughness: 0.5
                    }
                }
            }
            Repeater3D {
                model: root.showTube ? sim.splineCount : 0

                Model {
                    required property int index

                    source: "#Cylinder"
                    position: (sim.frame, sim.splineMidpoint(index))
                    rotation: (sim.frame, sim.splineRotation(index))
                    scale: (sim.frame, Qt.vector3d(root.tubeRadius / 50,
                                                   sim.splineLength(index) / 100 * 1.6,
                                                   root.tubeRadius / 50))

                    materials: PrincipledMaterial {
                        baseColor: Qt.hsva(sim.splineHue(index), 0.62, 0.95, 1.0)
                        roughness: 0.35
                    }
                }
            }

            Model {
                visible: sim.grabbedIndex >= 0
                source: "#Cylinder"
                position: (sim.frame, sim.pullExtension, sim.pullMidpoint())
                rotation: (sim.frame, sim.pullExtension, sim.pullRotation())
                scale: (sim.frame, Qt.vector3d(0.007, sim.pullExtension / 100, 0.007))

                materials: PrincipledMaterial {
                    baseColor: "#ffd479"
                    lighting: PrincipledMaterial.NoLighting
                }
            }
        }
    }

    MouseArea {
        id: grabArea
        anchors.fill: parent
        acceptedButtons: Qt.LeftButton | Qt.RightButton
        property real depth: 0

        function anchorFrom(px, py) {
            var scene = view.mapTo3DScene(Qt.vector3d(px, py, depth))
            return molecule.mapPositionFromScene(scene)
        }

        onPressed: (mouse) => {
            var hit = view.pick(mouse.x, mouse.y)

            if (mouse.button === Qt.RightButton) {
                if (hit.objectHit && hit.objectHit.beadIndex !== undefined)
                    sim.toggleAnchor(hit.objectHit.beadIndex)
                return
            }

            if (hit.objectHit && hit.objectHit.beadIndex !== undefined) {
                // 깊이는 반드시 mapFrom3DScene 의 z 로. pickResult.distance 는
                // 광선 길이라 mapTo3DScene 의 규약(근평면 거리)과 다르다.
                depth = view.mapFrom3DScene(hit.scenePosition).z
                if (sim.grab(hit.objectHit.beadIndex)) {
                    orbit.enabled = false
                    return
                }
            }
            mouse.accepted = false
        }

        onPositionChanged: (mouse) => {
            if (sim.grabbedIndex < 0)
                return
            var local = anchorFrom(mouse.x, mouse.y)
            sim.dragTo(local.x, local.y, local.z)
        }

        onReleased: {
            sim.release()
            orbit.enabled = true
        }

        onCanceled: {
            sim.release()
            orbit.enabled = true
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
            Key { text: "pull |F|" }
            Value { text: sim.grabbedIndex < 0 ? "-" : sim.pullForce.toFixed(1) }
            Key { text: "extension" }
            Value { text: sim.grabbedIndex < 0 ? "-" : sim.pullExtension.toFixed(2) + " A" }
            Key { text: "anchored" }
            Value { text: sim.anchoredIndex < 0 ? "-" : String(sim.anchoredIndex) }
            Key { text: "extent" }
            Value { text: sim.grabbedIndex < 0 ? "-" : sim.pullCoordinate.toFixed(2) + " A" }
        }
    }

    Rectangle {
        id: plotPanel
        visible: sim.traceLength > 1
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.bottomMargin: 100
        anchors.rightMargin: 14
        width: 320
        height: 200
        radius: 6
        color: "#d0181c22"

        Text {
            anchors.left: parent.left
            anchors.top: parent.top
            anchors.margins: 8
            color: "#78828f"
            font.pixelSize: 11
            text: "force vs extension  ·  " + sim.traceLength + " samples"
        }

        Canvas {
            id: plot
            anchors.fill: parent
            anchors.topMargin: 26
            anchors.margins: 10

            onPaint: {
                var ctx = getContext("2d")
                ctx.reset()

                var n = sim.traceLength
                if (n < 2)
                    return

                var x0 = sim.traceMinCoordinate
                var span = Math.max(sim.traceMaxCoordinate - x0, 1e-3)
                var ymax = Math.max(sim.traceMaxForce, 1e-3)

                ctx.strokeStyle = "#39424f"
                ctx.lineWidth = 1
                ctx.beginPath()
                ctx.moveTo(0, height)
                ctx.lineTo(width, height)
                ctx.moveTo(0, 0)
                ctx.lineTo(0, height)
                ctx.stroke()

                ctx.strokeStyle = "#ffd479"
                ctx.lineWidth = 1.4
                ctx.beginPath()
                var stride = Math.max(1, Math.floor(n / 900))
                for (var i = 0; i < n; i += stride) {
                    var px = (sim.traceCoordinateAt(i) - x0) / span * width
                    var py = height - sim.traceForceAt(i) / ymax * height
                    if (i === 0)
                        ctx.moveTo(px, py)
                    else
                        ctx.lineTo(px, py)
                }
                ctx.stroke()
            }
        }

        Text {
            anchors.right: parent.right
            anchors.bottom: parent.bottom
            anchors.margins: 8
            color: "#5c6672"
            font.pixelSize: 10
            text: sim.traceMinCoordinate.toFixed(0) + "-" + sim.traceMaxCoordinate.toFixed(0)
                  + " A  /  0-" + sim.traceMaxForce.toFixed(0)
        }
    }

    Timer {
        interval: 200
        repeat: true
        running: plotPanel.visible
        onTriggered: plot.requestPaint()
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
                    sim.clearTrace()
                }
            }

            Button {
                width: 78
                enabled: sim.traceLength > 1
                text: "save"
                onClicked: sim.saveTrace(runsDir + "/pull.csv")
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

            CheckBox {
                id: tubeBox
                anchors.verticalCenter: parent.verticalCenter
                text: "tube"
                checked: root.showTube
                onToggled: root.showTube = checked

                contentItem: Text {
                    text: tubeBox.text
                    color: "#a8b2bf"
                    font.pixelSize: 12
                    verticalAlignment: Text.AlignVCenter
                    leftPadding: tubeBox.indicator.width + tubeBox.spacing
                }
            }

            CheckBox {
                id: sticksBox
                enabled: !root.showTube
                anchors.verticalCenter: parent.verticalCenter
                text: "sticks"
                checked: root.showSticks
                onToggled: root.showSticks = checked

                contentItem: Text {
                    text: sticksBox.text
                    color: "#a8b2bf"
                    font.pixelSize: 12
                    verticalAlignment: Text.AlignVCenter
                    leftPadding: sticksBox.indicator.width + sticksBox.spacing
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
        text: "drag a bead to pull  ·  right-click a bead to anchor  ·  space to play/pause"
    }
}
