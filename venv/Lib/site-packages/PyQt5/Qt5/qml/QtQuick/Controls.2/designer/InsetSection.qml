/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
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
    caption: qsTr("Inset")

    SectionLayout {
        Label {
            text: qsTr("Vertical")
        }
        SecondColumnLayout {
            Label {
                text: qsTr("Top")
                tooltip: qsTr("Top inset for the background.")
                width: 42
            }
            SpinBox {
                maximumValue: 10000
                minimumValue: -10000
                realDragRange: 5000
                decimals: 0
                backendValue: backendValues.topInset
                Layout.fillWidth: true
            }
            Item {
                width: 4
                height: 4
            }

            Label {
                text: qsTr("Bottom")
                tooltip: qsTr("Bottom inset for the background.")
                width: 42
            }
            SpinBox {
                maximumValue: 10000
                minimumValue: -10000
                realDragRange: 5000
                decimals: 0
                backendValue: backendValues.bottomInset
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("Horizontal")
        }
        SecondColumnLayout {
            Label {
                text: qsTr("Left")
                tooltip: qsTr("Left inset for the background.")
                width: 42
            }
            SpinBox {
                maximumValue: 10000
                minimumValue: -10000
                realDragRange: 5000
                decimals: 0
                backendValue: backendValues.leftInset
                Layout.fillWidth: true
            }
            Item {
                width: 4
                height: 4
            }

            Label {
                text: qsTr("Right")
                tooltip: qsTr("Right inset for the background.")
                width: 42
            }
            SpinBox {
                maximumValue: 10000
                minimumValue: -10000
                realDragRange: 5000
                decimals: 0
                backendValue: backendValues.rightInset
                Layout.fillWidth: true
            }
        }
    }
}
