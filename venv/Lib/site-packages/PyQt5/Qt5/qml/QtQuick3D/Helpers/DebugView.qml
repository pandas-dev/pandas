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
import QtQuick3D 1.15

Rectangle {
    property var source: null
    width: layout.width + 10
    height: layout.height + 10
    color: "#80778BA5"
    radius: 5

    Column {
        id: layout
        anchors.centerIn: parent

        Text {
            text: source.renderStats.fps + " FPS (" + (source.renderStats.frameTime).toFixed(3) + "ms)"
            font.pointSize: 13
            color: "white"
        }
        Text {
            text: "Sync: " + (source.renderStats.syncTime).toFixed(3) + "ms"
            font.pointSize: 9
            color: "white"
        }
        Text {
            text: "Render: " + (source.renderStats.renderTime).toFixed(3) + "ms"
            font.pointSize: 9
            color: "white"
        }
        Text {
            text: "Max: " + (source.renderStats.maxFrameTime).toFixed(3) + "ms"
            font.pointSize: 9
            color: "white"
        }
    }
}
