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

Node {
    id: axisGrid_obj

    property alias gridColor: gridMaterial.diffuseColor
    property alias gridOpacity: gridMaterial.opacity
    property alias enableXZGrid: gridXZ.visible
    property alias enableXYGrid: gridXY.visible
    property alias enableYZGrid: gridYZ.visible
    property bool enableAxisLines: true

    // Axis Lines
    Model {
        id: xAxis
        source: "#Cube"
        position: Qt.vector3d(5000, 0, 0)
        scale: Qt.vector3d(100, .05, .05)
        visible: enableAxisLines

        materials: DefaultMaterial {
            lighting: DefaultMaterial.NoLighting
            diffuseColor: "red"
        }
    }

    Model {
        id: yAxis
        source: "#Cube"
        position: Qt.vector3d(0, 5000, 0)
        scale: Qt.vector3d(0.05, 100, 0.05)
        visible: enableAxisLines
        materials: DefaultMaterial {
            lighting: DefaultMaterial.NoLighting
            diffuseColor: "green"
        }
    }

    Model {
        id: zAxis
        source: "#Cube"
        position: Qt.vector3d(0, 0, 5000)
        scale: Qt.vector3d(0.05, 0.05, 100)
        visible: enableAxisLines
        materials: DefaultMaterial {
            lighting: DefaultMaterial.NoLighting
            diffuseColor: "blue"
        }
    }

    // Grid Lines
    DefaultMaterial {
        id: gridMaterial
        lighting: DefaultMaterial.NoLighting
        opacity: 0.5
        diffuseColor: Qt.rgba(0.8, 0.8, 0.8, 1)
    }

    Model {
        id: gridXZ
        source: "meshes/axisGrid.mesh"
        scale: Qt.vector3d(100, 100, 100)
        materials: [
            gridMaterial
        ]
    }

    Model {
        id: gridXY
        visible: false
        source: "meshes/axisGrid.mesh"
        scale: Qt.vector3d(100, 100, 100)
        eulerRotation: Qt.vector3d(90, 0, 0)
        materials: [
            gridMaterial
        ]
    }

    Model {
        id: gridYZ
        visible: false
        source: "meshes/axisGrid.mesh"
        scale: Qt.vector3d(100, 100, 100)
        eulerRotation: Qt.vector3d(0, 0, 90)
        materials: [
            gridMaterial
        ]
    }
}
