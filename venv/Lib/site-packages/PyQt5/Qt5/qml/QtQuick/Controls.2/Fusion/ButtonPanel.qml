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
import QtQuick.Controls 2.12
import QtQuick.Controls.impl 2.12
import QtQuick.Controls.Fusion 2.12
import QtQuick.Controls.Fusion.impl 2.12

Rectangle {
    id: panel

    property Item control
    property bool highlighted: control.highlighted

    visible: !control.flat || control.down || control.checked

    color: Fusion.buttonColor(control.palette, panel.highlighted, control.down || control.checked, control.hovered)
    gradient: control.down || control.checked ? null : buttonGradient

    Gradient {
        id: buttonGradient
        GradientStop {
            position: 0
            color: Fusion.gradientStart(Fusion.buttonColor(panel.control.palette, panel.highlighted, panel.control.down, panel.control.hovered))
        }
        GradientStop {
            position: 1
            color: Fusion.gradientStop(Fusion.buttonColor(panel.control.palette, panel.highlighted, panel.control.down, panel.control.hovered))
        }
    }

    radius: 2
    border.color: Fusion.buttonOutline(control.palette, panel.highlighted || control.visualFocus, control.enabled)

    Rectangle {
        x: 1; y: 1
        width: parent.width - 2
        height: parent.height - 2
        border.color: Fusion.innerContrastLine
        color: "transparent"
        radius: 2
    }
}
