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
import QtQuick.Templates 2.12 as T
import QtQuick.Controls.Material 2.12

T.StackView {
    id: control

    popEnter: Transition {
        // slide_in_left
        NumberAnimation { property: "x"; from: (control.mirrored ? -0.5 : 0.5) *  -control.width; to: 0; duration: 200; easing.type: Easing.OutCubic }
        NumberAnimation { property: "opacity"; from: 0.0; to: 1.0; duration: 200; easing.type: Easing.OutCubic }
    }

    popExit: Transition {
        // slide_out_right
        NumberAnimation { property: "x"; from: 0; to: (control.mirrored ? -0.5 : 0.5) * control.width; duration: 200; easing.type: Easing.OutCubic }
        NumberAnimation { property: "opacity"; from: 1.0; to: 0.0; duration: 200; easing.type: Easing.OutCubic }
    }

    pushEnter: Transition {
        // slide_in_right
        NumberAnimation { property: "x"; from: (control.mirrored ? -0.5 : 0.5) * control.width; to: 0; duration: 200; easing.type: Easing.OutCubic }
        NumberAnimation { property: "opacity"; from: 0.0; to: 1.0; duration: 200; easing.type: Easing.OutCubic }
    }

    pushExit: Transition {
        // slide_out_left
        NumberAnimation { property: "x"; from: 0; to: (control.mirrored ? -0.5 : 0.5) * -control.width; duration: 200; easing.type: Easing.OutCubic }
        NumberAnimation { property: "opacity"; from: 1.0; to: 0.0; duration: 200; easing.type: Easing.OutCubic }
    }

    replaceEnter: Transition {
        // slide_in_right
        NumberAnimation { property: "x"; from: (control.mirrored ? -0.5 : 0.5) * control.width; to: 0; duration: 200; easing.type: Easing.OutCubic }
        NumberAnimation { property: "opacity"; from: 0.0; to: 1.0; duration: 200; easing.type: Easing.OutCubic }
    }

    replaceExit: Transition {
        // slide_out_left
        NumberAnimation { property: "x"; from: 0; to: (control.mirrored ? -0.5 : 0.5) * -control.width; duration: 200; easing.type: Easing.OutCubic }
        NumberAnimation { property: "opacity"; from: 1.0; to: 0.0; duration: 200; easing.type: Easing.OutCubic }
    }
}
