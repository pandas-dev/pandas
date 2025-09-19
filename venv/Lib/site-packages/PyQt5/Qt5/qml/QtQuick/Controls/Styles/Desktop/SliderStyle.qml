/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Quick Controls module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

import QtQuick 2.2
import QtQuick.Controls 1.2
import QtQuick.Controls.Private 1.0

Style {
    readonly property Item control: __control
    property Component panel: StyleItem {
        elementType: "slider"
        sunken: control.pressed
        implicitWidth: 200
        contentHeight: horizontal ? 22 : 200
        contentWidth: horizontal ? 200 : 22

        maximum: control.maximumValue*100
        minimum: control.minimumValue*100
        step: control.stepSize*100
        value: control.__handlePos*100
        horizontal: control.orientation === Qt.Horizontal
        enabled: control.enabled
        hasFocus: control.activeFocus
        hover: control.hovered
        hints: control.styleHints
        activeControl: control.tickmarksEnabled ? "ticks" : ""
        property int handleWidth: 15
        property int handleHeight: 15
    }
    padding { top: 0 ; left: 0 ; right: 0 ; bottom: 0 }
}
