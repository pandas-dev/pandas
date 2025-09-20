/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

import QtQuick 2.15
import QtQuick.Layouts 1.15
import HelperWidgets 2.0
import StudioTheme 1.0 as StudioTheme

Column {
    width: parent.width

    Section {
        caption: qsTr("Scene Environment")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Antialiasing Mode")
                tooltip: qsTr("Sets the antialiasing mode applied to the scene.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "SceneEnvironment"
                    model: ["NoAA", "SSAA", "MSAA", "ProgressiveAA"]
                    backendValue: backendValues.antialiasingMode
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Antialiasing Quality")
                tooltip: qsTr("Sets the level of antialiasing applied to the scene.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "SceneEnvironment"
                    model: ["Medium", "High", "VeryHigh"]
                    backendValue: backendValues.antialiasingQuality
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Temporal AA")
                tooltip: qsTr("Enables temporal antialiasing using camera jittering and frame blending.")
            }

            SecondColumnLayout {
                CheckBox {
                    text: backendValues.temporalAAEnabled.valueToString
                    backendValue: backendValues.temporalAAEnabled
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Temporal AA Strength")
                tooltip: qsTr("Sets the amount of temporal antialiasing applied.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 2.0
                    minimumValue: 0.01
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.temporalAAStrength
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Background Mode")
                tooltip: qsTr("Controls if and how the background of the scene should be cleared.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "SceneEnvironment"
                    model: ["Transparent", "Unspecified", "Color", "SkyBox"]
                    backendValue: backendValues.backgroundMode
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Enable Depth Test")
                tooltip: qsTr("Enables depth testing. Disable to optimize render speed for layers with mostly transparent objects.")
            }

            SecondColumnLayout {
                CheckBox {
                    text: backendValues.depthTestEnabled.valueToString
                    backendValue: backendValues.depthTestEnabled
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Enable Depth Prepass")
                tooltip: qsTr("Draw depth buffer as a separate pass. Disable to optimize render speed for layers with low depth complexity.")
            }

            SecondColumnLayout {
                CheckBox {
                    text: backendValues.depthPrePassEnabled.valueToString
                    backendValue: backendValues.depthPrePassEnabled
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Effect")
                tooltip: qsTr("A post-processing effect applied to this scene.")
                Layout.alignment: Qt.AlignTop
                Layout.topMargin: 5
            }

            SecondColumnLayout {
                EditableListView {
                    backendValue: backendValues.effects
                    model: backendValues.effects.expressionAsList
                    Layout.fillWidth: true
                    typeFilter: "QtQuick3D.Effect"

                    onAdd: function(value) { backendValues.effects.idListAdd(value) }
                    onRemove: function(idx) { backendValues.effects.idListRemove(idx) }
                    onReplace: function (idx, value) { backendValues.effects.idListReplace(idx, value) }
                }

                ExpandingSpacer {}
            }

            PropertyLabel { text: qsTr("Clear Color") }

            ColorEditor {
                backendValue: backendValues.clearColor
                supportGradient: false
            }
        }
    }

    Section {
        caption: qsTr("Ambient Occlusion")
        width: parent.width

        SectionLayout {

            PropertyLabel {
                text: qsTr("AO Strength")
                tooltip: qsTr("Sets the amount of ambient occlusion applied.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 100
                    minimumValue: 0
                    decimals: 0
                    backendValue: backendValues.aoStrength
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("AO Distance")
                tooltip: qsTr("Sets how far ambient occlusion shadows spread away from objects.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 99999
                    minimumValue: 0
                    decimals: 0
                    backendValue: backendValues.aoDistance
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("AO Softness")
                tooltip: qsTr("Sets how smooth the edges of the ambient occlusion shading are.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 50
                    minimumValue: 0
                    decimals: 0
                    backendValue: backendValues.aoSoftness
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("AO Dither")
                tooltip: qsTr("Enables scattering of the ambient occlusion shadow band edges to improve smoothness (at the risk of sometimes producing obvious patterned artifacts).")
            }

            SecondColumnLayout {
                CheckBox {
                    text: backendValues.aoDither.valueToString
                    backendValue: backendValues.aoDither
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("AO Sample Rate")
                tooltip: qsTr("Sets the ambient occlusion quality (more shades of gray) at the expense of performance.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 4
                    minimumValue: 2
                    decimals: 0
                    backendValue: backendValues.aoSampleRate
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("AO Bias")
                tooltip: qsTr("Sets the cutoff distance preventing objects from exhibiting ambient occlusion at close distances.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 999999
                    minimumValue: -999999
                    realDragRange: 5000
                    decimals: 2
                    backendValue: backendValues.aoBias
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Image Based Lighting")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Light Probe")
                tooltip: qsTr("Defines a texture for overriding or setting an image based lighting texture for use with the skybox of this scene.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.lightProbe
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Probe Brightness")
                tooltip: qsTr("Sets the amount of light emitted by the light probe.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 999999
                    minimumValue: -999999
                    realDragRange: 5000
                    decimals: 0
                    backendValue: backendValues.probeBrightness
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Fast IBL")
                tooltip: qsTr("Use a faster approximation to image-based lighting.")
            }

            SecondColumnLayout {
                CheckBox {
                    text: backendValues.aoDither.valueToString
                    backendValue: backendValues.fastIBL
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Probe Horizon")
                tooltip: qsTr("Upper limit for horizon darkening of the light probe.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: -0.001
                    minimumValue: -1
                    decimals: 3
                    stepSize: 0.1
                    backendValue: backendValues.probeHorizon
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Probe FOV")
                tooltip: qsTr("Image source FOV for the case of using a camera-source as the IBL probe.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 180
                    minimumValue: 1.0
                    decimals: 1
                    backendValue: backendValues.probeFieldOfView
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }
}
