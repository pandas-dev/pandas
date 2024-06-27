/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Controls 2 module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

import QtQuick 2.12
import HelperWidgets 2.0
import QtQuick.Layouts 1.12

Column {
    width: parent.width

    Section {
        width: parent.width
        caption: qsTr("ProgressBar")

        SectionLayout {
            Label {
                text: qsTr("Indeterminate")
                tooltip: qsTr("Whether the progress is indeterminate.")
            }
            SecondColumnLayout {
                CheckBox {
                    text: backendValues.indeterminate.valueToString
                    backendValue: backendValues.indeterminate
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("Value")
                tooltip: qsTr("The current value of the progress.")
            }
            SecondColumnLayout {
                SpinBox {
                    minimumValue: Math.min(backendValues.from.value, backendValues.to.value)
                    maximumValue: Math.max(backendValues.from.value, backendValues.to.value)
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.value
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("From")
                tooltip: qsTr("The starting value for the progress.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 9999999
                    minimumValue: -9999999
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.from
                    Layout.fillWidth: true
                }
            }

            Label {
                text: qsTr("To")
                tooltip: qsTr("The ending value for the progress.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 9999999
                    minimumValue: -9999999
                    decimals: 2
                    stepSize: 0.1
                    backendValue: backendValues.to
                    Layout.fillWidth: true
                }
            }
        }
    }

    ControlSection {
        width: parent.width
    }

    FontSection {
        width: parent.width
    }

    PaddingSection {
        width: parent.width
    }
}
