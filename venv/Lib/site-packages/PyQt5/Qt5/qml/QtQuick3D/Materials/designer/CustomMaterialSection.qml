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
        caption: qsTr("Custom Material")
        width: parent.width

        SectionLayout {
            Label {
                text: qsTr("Transparency")
                tooltip: qsTr("Specifies if the material has transparency.")
            }
            SecondColumnLayout {
                CheckBox {
                    text: backendValues.hasTransparency.valueToString
                    backendValue: backendValues.hasTransparency
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Refraction")
                tooltip: qsTr("Specifies if the material has refraction.")
            }
            SecondColumnLayout {
                CheckBox {
                    text: backendValues.hasRefraction.valueToString
                    backendValue: backendValues.hasRefraction
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Always Dirty")
                tooltip: qsTr("Specifies if the material needs to be refreshed every time it is used.")
            }
            SecondColumnLayout {
                CheckBox {
                    text: backendValues.alwaysDirty.valueToString
                    backendValue: backendValues.alwaysDirty
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Shader Info")
                tooltip: qsTr("Shader info for the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.ShaderInfo"
                    Layout.fillWidth: true
                    backendValue: backendValues.shaderInfo
                }
            }
            Label {
                text: qsTr("Passes")
                tooltip: qsTr("Render passes of the material.")
            }
            SecondColumnLayout {
                EditableListView {
                    backendValue: backendValues.passes
                    model: backendValues.passes.expressionAsList
                    Layout.fillWidth: true
                    typeFilter: "QtQuick3D.Pass"

                    onAdd: function(value) { backendValues.passes.idListAdd(value) }
                    onRemove: function(idx) { backendValues.passes.idListRemove(idx) }
                    onReplace: function (idx, value) { backendValues.passes.idListReplace(idx, value) }
                }
            }
        }
    }

    Section {
        // Copied from quick3d's MaterialSection.qml
        caption: qsTr("Material")
        width: parent.width

        SectionLayout {
            Label {
                text: qsTr("Light Probe")
                tooltip: qsTr("Defines a texture for overriding or setting an image based lighting texture for use with this material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.lightProbe
                }
            }
            Label {
                text: qsTr("Displacement Map")
                tooltip: qsTr("Defines a grayscale image used to offset the vertices of geometry across the surface of the material.")
            }
            SecondColumnLayout {
                IdComboBox {
                    typeFilter: "QtQuick3D.Texture"
                    Layout.fillWidth: true
                    backendValue: backendValues.displacementMap
                }
            }
            Label {
                text: qsTr("Displacement Amount")
                tooltip: qsTr("Controls the offset amount for the displacement map.")
            }
            SecondColumnLayout {
                SpinBox {
                    maximumValue: 9999999
                    minimumValue: -9999999
                    realDragRange: 5000
                    decimals: 0
                    backendValue: backendValues.displacementAmount
                    Layout.fillWidth: true
                }
            }
            Label {
                text: qsTr("Culling Mode")
                tooltip: qsTr("Defines whether culling is enabled and which mode is actually enabled.")
            }
            ComboBox {
                scope: "Material"
                model: ["BackFaceCulling", "FrontFaceCulling", "NoCulling"]
                backendValue: backendValues.cullMode
                Layout.fillWidth: true
            }
        }
    }
}
