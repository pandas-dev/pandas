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
import QtQuick.Controls 2.12
import QtQuick.Controls.Imagine 2.12
import QtQuick.Controls.Imagine.impl 2.12

T.GroupBox {
    id: control

    implicitWidth: Math.max(implicitBackgroundWidth + leftInset + rightInset,
                            contentWidth + leftPadding + rightPadding,
                            implicitLabelWidth + leftPadding + rightPadding)
    implicitHeight: Math.max(implicitBackgroundHeight + topInset + bottomInset,
                             contentHeight + topPadding + bottomPadding)

    topPadding: (background ? background.topPadding : 0) + (implicitLabelWidth > 0 ? implicitLabelHeight + spacing : 0)
    leftPadding: background ? background.leftPadding : 0
    rightPadding: background ? background.rightPadding : 0
    bottomPadding: background ? background.bottomPadding : 0

    label: Label {
        width: control.width

        topPadding: background.topPadding
        leftPadding: background.leftPadding
        rightPadding: background.rightPadding
        bottomPadding: background.bottomPadding

        text: control.title
        font: control.font
        elide: Text.ElideRight
        verticalAlignment: Text.AlignVCenter

        color: control.palette.windowText

        background: NinePatchImage {
            width: parent.width
            height: parent.height

            source: Imagine.url + "groupbox-title"
            NinePatchImageSelector on source {
                states: [
                    {"disabled": !control.enabled},
                    {"mirrored": control.mirrored}
                ]
            }
        }
    }

    background: NinePatchImage {
        x: -leftInset
        y: control.topPadding - control.bottomPadding - topInset
        width: control.width + leftInset + rightInset
        height: control.height + topInset + bottomInset - control.topPadding + control.bottomPadding

        source: Imagine.url + "groupbox-background"
        NinePatchImageSelector on source {
            states: [
                {"disabled": !control.enabled},
                {"mirrored": control.mirrored}
            ]
        }
    }
}
