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
        caption: qsTr("Distortion")
        width: parent.width

        SectionLayout {
            Label {
                text: qsTr("Radius")
                tooltip: qsTr("Radius of the effect.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.radius
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Strength")
                tooltip: qsTr("Strength of the distortion.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 10
                    minimumValue: -10
                    decimals: 2
                    backendValue: backendValues.distortionStrength
                    Layout.fillWidth: true
                }
            }
        }
    }

    Section {
        id: centerSection
        width: parent.width
        caption: qsTr("Position")

        property int labelWidth: 10
        property int labelSpinBoxSpacing: 0
        property int spinBoxMinimumWidth: 120

        ColumnLayout {
            width: parent.width - 16

            Label {
                width: 100
                text: qsTr("Center")
                tooltip: qsTr("Center of the distortion.")
            }
            RowLayout {
                spacing: centerSection.labelSpinBoxSpacing

                Label {
                    text: qsTr("X")
                    width: centerSection.labelWidth
                }
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.center_x
                    Layout.fillWidth: true
                    Layout.minimumWidth: centerSection.spinBoxMinimumWidth
                }
            }
            RowLayout {
                spacing: centerSection.labelSpinBoxSpacing

                Label {
                    text: qsTr("Y")
                    width: centerSection.labelWidth
                }
                SpinBox {
                    maximumValue: 1
                    minimumValue: 0
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.center_y
                    Layout.fillWidth: true
                    Layout.minimumWidth: centerSection.spinBoxMinimumWidth
                }
            }
        }
    }
}
