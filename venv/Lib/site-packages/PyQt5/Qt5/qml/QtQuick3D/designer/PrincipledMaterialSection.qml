/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
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
import HelperWidgets 2.0
import QtQuick.Layouts 1.12

Column {
    width: parent.width

    Section {
        caption: qsTr("Principled Material")
        width: parent.width

        SectionLayout {
            Label {
                text: qsTr("Alpha Mode")
                tooltip: qsTr("Sets the mode for how the alpha channel of material color is used.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "PrincipledMaterial"
                    model: ["Opaque", "Mask", "Blend"]
                    backendValue: backendValues.alphaMode
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Alpha Cutoff")
                tooltip: qsTr("Specifies the cutoff value when using the Mask alphaMode.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.alphaCutoff
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Blend Mode")
                tooltip: qsTr("Determines how the colors of the model rendered blend with those behind it.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "PrincipledMaterial"
                    model: ["SourceOver", "Screen", "Multiply", "Overlay", "ColorBurn", "ColorDodge"]
                    backendValue: backendValues.blendMode
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Index Of Refraction")
                tooltip: qsTr("Controls how fast light travels through the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 1
                    maximumValue: 3
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.indexOfRefraction
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Lighting")
                tooltip: qsTr("Defines which lighting method is used when generating this material.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "PrincipledMaterial"
                    model: ["NoLighting", "FragmentLighting"]
                    backendValue: backendValues.lighting
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Metalness")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Metalness")
                tooltip: qsTr("Defines the metalness of the the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.metalness
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Metalness Map")
                tooltip: qsTr("Sets a texture to be used to set the metalness amount for the different parts of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.metalnessMap
                }
            }

            Label {
                text: qsTr("Metalness Channel")
                tooltip: qsTr("Defines the texture channel used to read the metalness value from metalnessMap.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "Material"
                    model: ["R", "G", "B", "A"]
                    backendValue: backendValues.metalnessChannel
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Normal")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Normal Map")
                tooltip: qsTr("Defines an RGB image used to simulate fine geometry displacement across the surface of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.normalMap
                }
            }

            Label {
                text: qsTr("Normal Strength")
                tooltip: qsTr("Controls the amount of simulated displacement for the normalMap.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.normalStrength
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Occlusion")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Occlusion Amount")
                tooltip: qsTr("Contains the factor used to modify the values from the occlusionMap texture.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.occlusionAmount
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Occlusion Map")
                tooltip: qsTr("Defines a texture used to determine how much indirect light the different areas of the material should receive.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.occlusionMap
                }
            }

            Label {
                text: qsTr("Occlusion Channel")
                tooltip: qsTr("Defines the texture channel used to read the occlusion value from occlusionMap.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "Material"
                    model: ["R", "G", "B", "A"]
                    backendValue: backendValues.occlusionChannel
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Opacity")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Opacity")
                tooltip: qsTr("Drops the opacity of just this material, separate from the model.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.opacity
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Opacity Map")
                tooltip: qsTr("Defines a texture used to control the opacity differently for different parts of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.opacityMap
                }
            }

            Label {
                text: qsTr("Opacity Channel")
                tooltip: qsTr("Defines the texture channel used to read the opacity value from opacityMap.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "Material"
                    model: ["R", "G", "B", "A"]
                    backendValue: backendValues.opacityChannel
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Roughness")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Roughness")
                tooltip: qsTr("Controls the size of the specular highlight generated from lights, and the clarity of reflections in general.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.roughness
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Roughness Map")
                tooltip: qsTr("Defines a texture to control the specular roughness of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.roughnessMap
                }
            }

            Label {
                text: qsTr("Roughness Channel")
                tooltip: qsTr("Defines the texture channel used to read the roughness value from roughnessMap.")
            }
            SecondColumnLayout {
                ComboBox {
                    scope: "Material"
                    model: ["R", "G", "B", "A"]
                    backendValue: backendValues.roughnessChannel
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Specular")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Specular Amount")
                tooltip: qsTr("Controls the strength of specularity (highlights and reflections).")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.specularAmount
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Specular Map")
                tooltip: qsTr("Defines a RGB Texture to modulate the amount and the color of specularity across the surface of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.specularMap
                }
            }

            Label {
                text: qsTr("Specular Reflection Map")
                tooltip: qsTr("Sets a texture used for specular highlights on the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.specularReflectionMap
                }
            }

            Label {
                text: qsTr("Specular Tint")
                tooltip: qsTr("Defines how much of the base color contributes to the specular reflections.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.specularTint
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        caption: qsTr("Base Color")
        width: parent.width

        Column {
            width: parent.width

            ColorEditor {
                caption: qsTr("Base Color")
                backendValue: backendValues.baseColor
                supportGradient: false
                Layout.fillWidth: true
            }
            SectionLayout {
                Label {
                    text: qsTr("Base Color Map")
                    tooltip: qsTr("Defines a texture used to set the base color of the material.")
                }
                SecondColumnLayout {
                    IdComboBox {
                        typeFilter: "QtQuick3D.Texture"
                        Layout.fillWidth: true
                        backendValue: backendValues.baseColorMap
                    }
                }
            }
        }
    }

    Section {
        caption: qsTr("Emissive Color")
        width: parent.width

        Column {
            width: parent.width

            ColorEditor {
                caption: qsTr("Emissive Color")
                backendValue: backendValues.emissiveColor
                supportGradient: false
                Layout.fillWidth: true
            }
            SectionLayout {
                Label {
                    text: qsTr("Emissive Map")
                    tooltip: qsTr("Sets a texture to be used to set the emissive factor for different parts of the material.")
                }
                SecondColumnLayout {
                    IdComboBox {
                        typeFilter: "QtQuick3D.Texture"
                        Layout.fillWidth: true
                        backendValue: backendValues.emissiveMap
                    }
                }
            }
        }
    }
}
