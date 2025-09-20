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
import QtQuick.Layouts 1.15
import HelperWidgets 2.0
import StudioTheme 1.0 as StudioTheme

Section {
    caption: qsTr("Model")

    SectionLayout {
        id: tessellationSection

        PropertyLabel {
            text: qsTr("Source")
            tooltip: qsTr("Defines the location of the mesh file containing the geometry of this model.")
        }

        SecondColumnLayout {
            UrlChooser {
                backendValue: backendValues.source
                filter: "*.mesh"
            }

            ExpandingSpacer {}
        }

        function hasTessellationMode(mode) {
            if (tessellationModeComboBox.backendValue.valueToString !== "" &&
                tessellationModeComboBox.backendValue.valueToString !== mode)
                return false

            if (tessellationModeComboBox.backendValue.enumeration !== "" &&
                tessellationModeComboBox.backendValue.enumeration !== mode)
                return false

            return true
        }

        PropertyLabel {
            text: qsTr("Tessellation Mode")
            tooltip: qsTr("Defines what method to use to dynamically generate additional geometry for the model.")
        }

        SecondColumnLayout {
            ComboBox {
                id: tessellationModeComboBox
                scope: "Model"
                model: ["NoTessellation", "Linear", "Phong", "NPatch"]
                backendValue: backendValues.tessellationMode
                implicitWidth: StudioTheme.Values.singleControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Edge Tessellation")
            tooltip: qsTr("Defines the edge multiplier to the tessellation generator.")
        }

        SecondColumnLayout {
            SpinBox {
                maximumValue: 64.0
                minimumValue: 0.0
                decimals: 0
                backendValue: backendValues.edgeTessellation
                enabled: !tessellationSection.hasTessellationMode("NoTessellation")
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Inner Tessellation")
            tooltip: qsTr("Defines the inner multiplier to the tessellation generator.")
        }

        SecondColumnLayout {
            SpinBox {
                maximumValue: 64.0
                minimumValue: 0.0
                decimals: 0
                backendValue: backendValues.innerTessellation
                enabled: !tessellationSection.hasTessellationMode("NoTessellation")
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Enable Wireframe Mode")
            tooltip: qsTr("Enables the wireframe mode if tesselation is enabled.")
        }

        SecondColumnLayout {
            CheckBox {
                text: backendValues.isWireframeMode.valueToString
                backendValue: backendValues.isWireframeMode
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Casts Shadows")
            tooltip: qsTr("Enables the geometry of this model to be rendered to the shadow maps.")
        }

        SecondColumnLayout {
            CheckBox {
                text: backendValues.castsShadows.valueToString
                backendValue: backendValues.castsShadows
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Receives Shadows")
            tooltip: qsTr("Enables the geometry of this model to receive shadows.")
        }

        SecondColumnLayout {
            CheckBox {
                text: backendValues.receivesShadows.valueToString
                backendValue: backendValues.receivesShadows
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Pickable")
            tooltip: qsTr("Controls whether the model is pickable or not.")
        }

        SecondColumnLayout {
            CheckBox {
                text: backendValues.pickable.valueToString
                backendValue: backendValues.pickable
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Materials")
        }

        SecondColumnLayout {
            EditableListView {
                backendValue: backendValues.materials
                model: backendValues.materials.expressionAsList
                Layout.fillWidth: true

                onAdd: function(value) { backendValues.materials.idListAdd(value) }
                onRemove: function(idx) { backendValues.materials.idListRemove(idx) }
                onReplace: function (idx, value) { backendValues.materials.idListReplace(idx, value) }
            }

            ExpandingSpacer {}
        }
    }
}
