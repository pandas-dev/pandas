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
import QtQuick.Layouts 1.15
import HelperWidgets 2.0
import StudioTheme 1.0 as StudioTheme

Column {
    width: parent.width

    Section {
        caption: qsTr("Bump")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Scale")
                tooltip: qsTr("Set the scale of the bump bands.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 5
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.bumpScale
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Bands")
                tooltip: qsTr("Set the number of the bump bands.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10
                    decimals: 0
                    backendValue: backendValues.bumpBands
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Coordinates")
                tooltip: qsTr("Sets the bump coordinates of the refraction.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10000
                    decimals: 2
                    backendValue: backendValues.bumpCoords_x
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

                ControlLabel {
                    text: "X"
                    color: StudioTheme.Values.theme3DAxisXColor
                }

                ExpandingSpacer {}
            }

            PropertyLabel {}

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10000
                    decimals: 2
                    backendValue: backendValues.bumpCoords_y
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

                ControlLabel {
                    text: "Y"
                    color: StudioTheme.Values.theme3DAxisYColor
                }

                ExpandingSpacer {}
            }


            PropertyLabel {}

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10000
                    decimals: 2
                    backendValue: backendValues.bumpCoords_z
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

                ControlLabel {
                    text: "Z"
                    color: StudioTheme.Values.theme3DAxisZColor
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("General")
        width: parent.width

        SectionLayout {
            PropertyLabel { text: qsTr("Glass Color") }

            ColorEditor {
                backendValue: backendValues.glass_color
                supportGradient: false
                isVector3D: true
            }

            PropertyLabel {
                text: qsTr("Roughness")
                tooltip: qsTr("Set the material roughness.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.roughness
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Blur Size")
                tooltip: qsTr("Set the amount of blurring behind the glass.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 50
                    decimals: 2
                    backendValue: backendValues.blur_size
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Fresnel Power")
                tooltip: qsTr("Set the fresnel power of the material.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 100
                    decimals: 2
                    backendValue: backendValues.uFresnelPower
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Reflectivity")
                tooltip: qsTr("Set the reflectivity of the material.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.reflectivity_amount
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Index of Refraction")
                tooltip: qsTr("Set the index of refraction for the material.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 1.4
                    maximumValue: 2.1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.glass_ior
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Noise")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Scale")
                tooltip: qsTr("Set the noise scale.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 40
                    decimals: 2
                    backendValue: backendValues.noiseScale
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Coordinates")
                tooltip: qsTr("Sets the noise coordinates.")
            }


            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10000
                    decimals: 2
                    backendValue: backendValues.noiseCoords_x
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

                ControlLabel {
                    text: "X"
                    color: StudioTheme.Values.theme3DAxisXColor
                }

                ExpandingSpacer {}
            }

            PropertyLabel {}

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10000
                    decimals: 2
                    backendValue: backendValues.noiseCoords_y
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

                ControlLabel {
                    text: "Y"
                    color: StudioTheme.Values.theme3DAxisYColor
                }

                ExpandingSpacer {}
            }


            PropertyLabel {}

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 10000
                    decimals: 2
                    backendValue: backendValues.noiseCoords_z
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

                ControlLabel {
                    text: "Z"
                    color: StudioTheme.Values.theme3DAxisZColor
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Random Gradient Maps")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("1D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.randomGradient1D_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("2D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.randomGradient2D_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("3D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.randomGradient3D_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("4D")
                tooltip: qsTr("Defines a texture map used to create the random bumpiness of the material.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.randomGradient4D_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }
}
