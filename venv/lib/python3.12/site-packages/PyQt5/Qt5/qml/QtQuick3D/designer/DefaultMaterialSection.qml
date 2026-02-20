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
        caption: qsTr("Default Material")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Lighting")
                tooltip: qsTr("Defines which lighting method is used when generating this material.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "DefaultMaterial"
                    model: ["NoLighting", "FragmentLighting"]
                    backendValue: backendValues.lighting
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Blend Mode")
                tooltip: qsTr("Determines how the colors of the model rendered blend with those behind it.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "DefaultMaterial"
                    model: ["SourceOver", "Screen", "Multiply", "Overlay", "ColorBurn", "ColorDodge"]
                    backendValue: backendValues.blendMode
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Enable Vertex Colors")
                tooltip: qsTr("Enables the use of vertex colors from the mesh.")
            }

            SecondColumnLayout {
                CheckBox {
                    text: backendValues.vertexColorsEnabled.valueToString
                    backendValue: backendValues.vertexColorsEnabled
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Diffuse")
        width: parent.width

        SectionLayout {
            PropertyLabel { text: qsTr("Diffuse Color") }

            ColorEditor {
                backendValue: backendValues.diffuseColor
                supportGradient: false
            }

            PropertyLabel {
                text: qsTr("Diffuse Map")
                tooltip: qsTr("Defines a texture to apply to the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.diffuseMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Emissive")
        width: parent.width

        SectionLayout {
            PropertyLabel { text: qsTr("Emissive Color") }

            ColorEditor {
                backendValue: backendValues.emissiveColor
                supportGradient: false
            }

            PropertyLabel {
                text: qsTr("Emissive Factor")
                tooltip: qsTr("Determines the amount of self-illumination from the material (will not light other objects).")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.emissiveFactor
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Emissive Map")
                tooltip: qsTr("Sets a texture to be used to set the emissive factor for different parts of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.emissiveMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Specular")
        width: parent.width

        SectionLayout {
            PropertyLabel { text: qsTr("Specular Tint") }

            ColorEditor {
                backendValue: backendValues.specularTint
                supportGradient: false
            }

            PropertyLabel {
                text: qsTr("Specular Amount")
                tooltip: qsTr("Controls the strength of specularity (highlights and reflections).")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.specularAmount
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Specular Map")
                tooltip: qsTr("Defines a RGB texture to modulate the amount and the color of specularity across the surface of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.specularMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Specular Model")
                tooltip: qsTr("Determines which functions are used to calculate specular highlights for lights in the scene.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "DefaultMaterial"
                    model: ["Default", "KGGX", "KWard"]
                    backendValue: backendValues.specularModel
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Reflection Map")
                tooltip: qsTr("Sets a texture used for specular highlights on the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.specularReflectionMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Index of Refraction")
                tooltip: qsTr("Controls what angles of reflections are affected by the Fresnel power.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 3
                    minimumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.indexOfRefraction
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Fresnel Power")
                tooltip: qsTr("Decreases head-on reflections (looking directly at the surface) while maintaining reflections seen at grazing angles.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 9999999
                    minimumValue: -9999999
                    realDragRange: 5000
                    decimals: 2
                    backendValue: backendValues.fresnelPower
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Specular Roughness")
                tooltip: qsTr("Controls the size of the specular highlight generated from lights and the clarity of reflections in general.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0.001
                    decimals: 3
                    backendValue: backendValues.specularRoughness
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Roughness Map")
                tooltip: qsTr("Defines a texture to control the specular roughness of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.roughnessMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Opacity")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Opacity")
                tooltip: qsTr("Sets the visibility of the geometry for this material.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.opacity
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Opacity Map")
                tooltip: qsTr("Defines a texture used to control the opacity differently for different parts of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.opacityMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Bump/Normal")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Bump Amount")
                tooltip: qsTr("Controls the amount of simulated displacement for the bump map or normal map.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.bumpAmount
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Bump Map")
                tooltip: qsTr("Defines a grayscale texture to simulate fine geometry displacement across the surface of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    id: bumpMapComboBox
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.bumpMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth

                    Connections {
                        target: normalMapComboBox.backendValue
                        function onExpressionChanged() {
                            if (normalMapComboBox.backendValue.expression !== "")
                                bumpMapComboBox.backendValue.resetValue()
                        }
                    }
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Normal Map")
                tooltip: qsTr("Defines a RGB image used to simulate fine geometry displacement across the surface of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    id: normalMapComboBox
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.normalMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth

                    Connections {
                        target: bumpMapComboBox.backendValue
                        function onExpressionChanged() {
                            if (bumpMapComboBox.backendValue.expression !== "")
                                normalMapComboBox.backendValue.resetValue()
                        }
                    }
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Translucency")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Translucency Falloff")
                tooltip: qsTr("Defines the amount of falloff for the translucency based on the angle of the normals of the object to the light source.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 999999
                    minimumValue: -999999
                    realDragRange: 5000
                    decimals: 2
                    backendValue: backendValues.translucentFalloff
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Diffuse Light Wrap")
                tooltip: qsTr("Determines the amount of light wrap for the translucency map.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.diffuseLightWrap
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Translucency Map")
                tooltip: qsTr("Defines a grayscale texture controlling how much light can pass through the material from behind.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.translucencyMap
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }
}
