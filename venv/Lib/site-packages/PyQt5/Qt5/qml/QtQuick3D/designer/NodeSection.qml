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
import HelperWidgets 2.0
import QtQuick.Layouts 1.12
import StudioTheme 1.0 as StudioTheme

Column {
    width: parent.width

    Section {
        width: parent.width
        caption: qsTr("Node")

        SectionLayout {

            Label {
                text: qsTr("Opacity")
                tooltip: qsTr("Controls the local opacity value of the node.")
            }

            // ### should be a slider
            SpinBox {
                minimumValue: 0
                maximumValue: 1
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.opacity
                Layout.fillWidth: true
                sliderIndicatorVisible: true
            }

            Label {
                text: qsTr("Visibility")
                tooltip: qsTr("Sets the local visibility of the node.")
            }
            SecondColumnLayout {
                // ### should be a slider
                CheckBox {
                    text: qsTr("Is Visible")
                    backendValue: backendValues.visible
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        id: transformSection
        width: parent.width
        caption: qsTr("Transform")

        property int labelWidth: 10
        property int labelSpinBoxSpacing: 0
        property int spinBoxMinimumWidth: 120

        GridLayout {
            columns: 2
            rows: 2
            columnSpacing: 24
            rowSpacing: 12
            width: parent.width - 16

            ColumnLayout {

                Label {
                    width: 100
                    text: qsTr("Translation")
                    tooltip: qsTr("Sets the translation of the node.")
                }

                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("X")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisXColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.x
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Y")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisYColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.y
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Z")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisZColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.z
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
            }

            ColumnLayout {

                Label {
                    width: 100
                    text: qsTr("Rotation")
                    tooltip: qsTr("Sets the rotation of the node in degrees.")
                }

                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("X")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisXColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.eulerRotation_x
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Y")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisYColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.eulerRotation_y
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Z")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisZColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.eulerRotation_z
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
            }

            ColumnLayout {

                Label {
                    text: qsTr("Scale")
                    tooltip: qsTr("Sets the scale of the node.")
                }

                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("X")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisXColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 50
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.scale_x
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Y")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisYColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 50
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.scale_y
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Z")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisZColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 50
                        decimals: 2
                        stepSize: 0.1
                        backendValue: backendValues.scale_z
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
            }

            ColumnLayout {

                Label {
                    text: qsTr("Pivot")
                    tooltip: qsTr("Sets the pivot of the node.")
                }

                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("X")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisXColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.pivot_x
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Y")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisYColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.pivot_y
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
                RowLayout {
                    spacing: transformSection.labelSpinBoxSpacing

                    Label {
                        text: qsTr("Z")
                        width: transformSection.labelWidth
                        color: StudioTheme.Values.theme3DAxisZColor
                    }
                    SpinBox {
                        maximumValue: 9999999
                        minimumValue: -9999999
                        realDragRange: 5000
                        decimals: 2
                        backendValue: backendValues.pivot_z
                        Layout.fillWidth: true
                        Layout.minimumWidth: transformSection.spinBoxMinimumWidth
                    }
                }
            }

        }
    }
}
