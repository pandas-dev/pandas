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
import StudioTheme 1.0 as StudioTheme

Column {
    id: materialRoot
    width: parent.width

    property int labelWidth: 10
    property int labelSpinBoxSpacing: 0
    property int spinBoxMinimumWidth: 120

    Section {
        caption: qsTr("Environment Map")
        width: parent.width

        SectionLayout {
            Label {
                text: qsTr("Enabled")
                tooltip: qsTr("Specifies if the environment map is enabled.")
            }
            SecondColumnLayout {
                CheckBox {
                    text: backendValues.uEnvironmentMappingEnabled.valueToString
                    backendValue: backendValues.uEnvironmentMappingEnabled
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Texture")
                tooltip: qsTr("Defines a texture for environment map.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.uEnvironmentTexture_texture
                    defaultItem: qsTr("Default")
                }
            }
        }
    }

    Section {
        caption: qsTr("Shadow Map")
        width: parent.width

        SectionLayout {
            Label {
                text: qsTr("Enabled")
                tooltip: qsTr("Specifies if the shadow map is enabled.")
            }
            SecondColumnLayout {
                CheckBox {
                    text: backendValues.uShadowMappingEnabled.valueToString
                    backendValue: backendValues.uShadowMappingEnabled
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Texture")
                tooltip: qsTr("Defines a texture for shadow map.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.uBakedShadowTexture_texture
                    defaultItem: qsTr("Default")
                }
            }
        }
    }


    Section {
        caption: qsTr("Glass Color")
        width: parent.width
        ColorEditor {
            caption: qsTr("Glass Color")
            backendValue: backendValues.glass_color
            supportGradient: false
            isVector3D: true
            Layout.fillWidth: true
        }
    }

    Section {
        caption: qsTr("Band Light Color")
        width: parent.width
        ColorEditor {
            caption: qsTr("Band Light Color")
            backendValue: backendValues.intLightCol
            supportGradient: false
            isVector3D: true
            Layout.fillWidth: true
        }
    }

    Section {
        caption: qsTr("Bump")
        width: parent.width
        ColumnLayout {
            width: parent.width - 16
            SectionLayout {
                Label {
                    text: qsTr("Scale")
                    tooltip: qsTr("Set the scale of the bump bands.")
                }
                SecondColumnLayout {
                    SpinBox {
                        maximumValue: 5
                        minimumValue: 0
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.bumpScale
                        Layout.fillWidth: true
                    }
                }
                Label {
                    text: qsTr("Bands")
                    tooltip: qsTr("Set the number of the bump bands.")
                }
                SecondColumnLayout {
                    SpinBox {
                        maximumValue: 10
                        minimumValue: 0
                        decimals: 0
                        backendValue: backendValues.bumpBands
                        Layout.fillWidth: true
                    }
                }
                Label {
                    text: qsTr("Strength")
                    tooltip: qsTr("Set the glass bump map strength.")
                }
                SecondColumnLayout {
                    SpinBox {
                        maximumValue: 1
                        minimumValue: 0
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.glass_bfactor
                        Layout.fillWidth: true
                    }
                }
                Label {
                    text: qsTr("Internal")
                    tooltip: qsTr("Specifies if the bump map be used only for the internal lighting.")
                }
                SecondColumnLayout {
                    CheckBox {
                        text: backendValues.glass_binside.valueToString
                        backendValue: backendValues.glass_binside
                        Layout.fillWidth: true
                    }
                }
            Label {
                text: qsTr("Texture")
                tooltip: qsTr("Defines a texture for bump map.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.glass_bump_texture
                    defaultItem: qsTr("Default")
                }
            }
            }
            ColumnLayout {
                width: parent.width
                Label {
                    width: 100
                    text: qsTr("Coordinates")
                    tooltip: qsTr("Sets the bump coordinates of the refraction.")
                }

                RowLayout {
                    spacing: materialRoot.labelSpinBoxSpacing

                    Label {
                        text: qsTr("X")
                        width: materialRoot.labelWidth
                        color: StudioTheme.Values.theme3DAxisXColor
                    }
                    SpinBox {
                        maximumValue: 10000
                        minimumValue: 0
                        realDragRange: 1000
                        decimals: 2
                        backendValue: backendValues.bumpCoords_x
                        Layout.fillWidth: true
                        Layout.minimumWidth: materialRoot.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: materialRoot.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Y")
                        width: materialRoot.labelWidth
                        color: StudioTheme.Values.theme3DAxisYColor
                    }
                    SpinBox {
                        maximumValue: 10000
                        minimumValue: 0
                        realDragRange: 1000
                        decimals: 2
                        backendValue: backendValues.bumpCoords_y
                        Layout.fillWidth: true
                        Layout.minimumWidth: materialRoot.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: materialRoot.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Z")
                        width: materialRoot.labelWidth
                        color: StudioTheme.Values.theme3DAxisZColor
                    }
                    SpinBox {
                        maximumValue: 10000
                        minimumValue: 0
                        realDragRange: 1000
                        decimals: 2
                        backendValue: backendValues.bumpCoords_z
                        Layout.fillWidth: true
                        Layout.minimumWidth: materialRoot.spinBoxMinimumWidth
                    }
                }
            }
        }
    }

    Section {
        caption: qsTr("General")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("Roughness")
                tooltip: qsTr("Set the material roughness.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.roughness
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Blur Size")
                tooltip: qsTr("Set the amount of blurring behind the glass.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 50
                    minimumValue: 0
                    decimals: 2
                    backendValue: backendValues.blur_size
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Refract Depth")
                tooltip: qsTr("Set the refract depth of the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 5
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.refract_depth
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Fresnel Power")
                tooltip: qsTr("Set the fresnel power of the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 100
                    minimumValue: 0
                    decimals: 2
                    backendValue: backendValues.uFresnelPower
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Reflectivity")
                tooltip: qsTr("Set the reflectivity of the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.reflectivity_amount
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Index of Refraction")
                tooltip: qsTr("Set the index of refraction for the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 2.1
                    minimumValue: 1.4
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.glass_ior
                    Layout.fillWidth: true
                }
            }
        }
    }
    Section {
        caption: qsTr("Band Light")
        width: parent.width
        ColumnLayout {
            width: parent.width - 16
            SectionLayout {
                Label {
                    text: qsTr("Falloff")
                    tooltip: qsTr("Set the light intensity falloff rate.")
                }
                SecondColumnLayout {
                    SpinBox {
                        maximumValue: 10
                        minimumValue: 0
                        decimals: 2
                        backendValue: backendValues.intLightFall
                        Layout.fillWidth: true
                    }
                }
                Label {
                    text: qsTr("Angle")
                    tooltip: qsTr("Set the angle of lightsource. Band is perpendicular to this.")
                }
                SecondColumnLayout {
                    SpinBox {
                        maximumValue: 360
                        minimumValue: 0
                        decimals: 2
                        backendValue: backendValues.intLightRot
                        Layout.fillWidth: true
                    }
                }
                Label {
                    text: qsTr("Brightness")
                    tooltip: qsTr("Set the brightness of the band light.")
                }
                SecondColumnLayout {
                    SpinBox {
                        maximumValue: 10000
                        minimumValue: 0
                        realDragRange: 1000
                        decimals: 2
                        backendValue: backendValues.intLightBrt
                        Layout.fillWidth: true
                    }
                }
            }
            ColumnLayout {
                width: parent.width
                Label {
                    width: 100
                    text: qsTr("Position")
                    tooltip: qsTr("Sets the Position of the band light in the UV space.")
                }

                RowLayout {
                    spacing: materialRoot.labelSpinBoxSpacing

                    Label {
                        text: qsTr("X")
                        width: materialRoot.labelWidth
                    }
                    SpinBox {
                        maximumValue: 1
                        minimumValue: 0
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.intLightPos_x
                        Layout.fillWidth: true
                        Layout.minimumWidth: materialRoot.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: materialRoot.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Y")
                        width: materialRoot.labelWidth
                    }
                    SpinBox {
                        maximumValue: 1
                        minimumValue: 0
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.intLightPos_y
                        Layout.fillWidth: true
                        Layout.minimumWidth: materialRoot.spinBoxMinimumWidth
                    }
                }
            }
        }
    }
    Section {
        caption: qsTr("Random Gradient Maps")
        width: parent.width
        SectionLayout {
            Label {
                text: qsTr("1D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.randomGradient1D_texture
                    defaultItem: qsTr("Default")
                }
            }
            Label {
                text: qsTr("2D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.randomGradient2D_texture
                    defaultItem: qsTr("Default")
                }
            }
            Label {
                text: qsTr("3D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.randomGradient3D_texture
                    defaultItem: qsTr("Default")
                }
            }
            Label {
                text: qsTr("4D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.randomGradient4D_texture
                    defaultItem: qsTr("Default")
                }
            }
        }
    }
}
