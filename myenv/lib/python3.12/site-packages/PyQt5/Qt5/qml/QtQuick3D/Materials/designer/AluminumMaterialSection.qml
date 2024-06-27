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
        caption: qsTr("Roughness")
        width: parent.width

        SectionLayout {
            PropertyLabel { text: qsTr("Metal Color") }

            ColorEditor {
                backendValue: backendValues.metal_color
                supportGradient: false
                isVector3D: true
            }

            PropertyLabel {
                text: qsTr("Map Offset")
                tooltip: qsTr("Set the material roughness map offset.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.roughness_map_offset
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Map scale")
                tooltip: qsTr("Set the material roughness map scale.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.roughness_map_scale
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Texture")
                tooltip: qsTr("Defines a texture for roughness map.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.roughness_texture_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }

    Section {
        caption: qsTr("Reflection")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Map Offset")
                tooltip: qsTr("Set the material reclection map offset.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.reflection_map_offset
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Map scale")
                tooltip: qsTr("Set the material reclection map scale.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 1
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.reflection_map_scale
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Texture")
                tooltip: qsTr("Defines a texture for reflection map.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.reflection_texture_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Tiling")
                tooltip: qsTr("Sets the tiling repeat of the reflection map.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 1
                    maximumValue: 100
                    decimals: 0
                    backendValue: backendValues.tiling_x
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
                    minimumValue: 1
                    maximumValue: 100
                    decimals: 0
                    backendValue: backendValues.tiling_y
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
                    minimumValue: 1
                    maximumValue: 100
                    decimals: 0
                    backendValue: backendValues.tiling_z
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
        caption: qsTr("Bump")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Amount")
                tooltip: qsTr("Set the bump map bumpiness.")
            }

            SecondColumnLayout {
                SpinBox {
                    minimumValue: 0
                    maximumValue: 2
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.bump_amount
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Texture")
                tooltip: qsTr("Defines a texture for bump map.")
            }

            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    backendValue: backendValues.bump_texture_texture
                    defaultItem: qsTr("Default")
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }
}
