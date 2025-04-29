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

Section {
    caption: qsTr("Control")

    SectionLayout {
        Label {
            text: qsTr("Enabled")
            tooltip: qsTr("Whether the control is enabled.")
        }
        SecondColumnLayout {
            CheckBox {
                text: backendValues.enabled.valueToString
                backendValue: backendValues.enabled
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("Focus Policy")
            tooltip: qsTr("Focus policy of the control.")
        }
        SecondColumnLayout {
            ComboBox {
                backendValue: backendValues.focusPolicy
                model: [ "TabFocus", "ClickFocus", "StrongFocus", "WheelFocus", "NoFocus" ]
                scope: "Qt"
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("Hover")
            tooltip: qsTr("Whether control accepts hover events.")
        }
        SecondColumnLayout {
            CheckBox {
                text: backendValues.hoverEnabled.valueToString
                backendValue: backendValues.hoverEnabled
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("Spacing")
            tooltip: qsTr("Spacing between internal elements of the control.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 9999999
                minimumValue: -9999999
                decimals: 0
                backendValue: backendValues.spacing
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("Wheel")
            tooltip: qsTr("Whether control accepts wheel events.")
        }
        SecondColumnLayout {
            CheckBox {
                text: backendValues.wheelEnabled.valueToString
                backendValue: backendValues.wheelEnabled
                Layout.fillWidth: true
            }
        }
    }
}
