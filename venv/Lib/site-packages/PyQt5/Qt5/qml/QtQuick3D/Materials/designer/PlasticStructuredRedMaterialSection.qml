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
        caption: qsTr("Diffuse Color")
        width: parent.width
        ColorEditor {
            caption: qsTr("Diffuse Color")
            backendValue: backendValues.diffuse_color
            supportGradient: false
            isVector3D: true
            Layout.fillWidth: true
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
                text: qsTr("Index of Refraction")
                tooltip: qsTr("Set the index of refraction for the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1.6
                    minimumValue: 1.4
                    decimals: 3
                    stepSize: 0.01
                    backendValue: backendValues.material_ior
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Texture Scaling")
                tooltip: qsTr("Set the texture scaling of the material.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.texture_scaling
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Bump Factor")
                tooltip: qsTr("Set the strength of the bumpiness.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.bump_factor
                    Layout.fillWidth: true
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
