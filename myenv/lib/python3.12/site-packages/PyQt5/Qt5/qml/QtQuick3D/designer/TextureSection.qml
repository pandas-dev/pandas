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
    caption: qsTr("Image")
    width: parent.width

    SectionLayout {
        PropertyLabel {
            text: qsTr("Source")
            tooltip: qsTr("Holds the location of an image file containing the data used by the texture.")
        }

        SecondColumnLayout {
            UrlChooser {
                backendValue: backendValues.source
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Source Item")
            tooltip: qsTr("Defines an item to be used as the source of the texture.")
        }

        SecondColumnLayout {
            IdComboBox {
                typeFilter: "QtQuick.Item"
                backendValue: backendValues.sourceItem
                implicitWidth: StudioTheme.Values.singleControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Scale")
        }

        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: 0
                realDragRange: 2
                decimals: 2
                backendValue: backendValues.scaleU
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

            ControlLabel {
                text: "U"
                tooltip: qsTr("Defines how to scale the U texture coordinate when mapping to UV coordinates of a mesh.")
            }

            Spacer { implicitWidth: StudioTheme.Values.controlGap }

            SpinBox {
                maximumValue: 999999
                minimumValue: 0
                realDragRange: 2
                decimals: 2
                backendValue: backendValues.scaleV
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

            ControlLabel {
                text: "V"
                tooltip: qsTr("Defines how to scale the V texture coordinate when mapping to UV coordinates of a mesh.")
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Texture Mapping")
            tooltip: qsTr("Defines which method of mapping to use when sampling this texture.")
        }

        SecondColumnLayout {
            ComboBox {
                scope: "Texture"
                model: ["UV", "Environment", "LightProbe"]
                backendValue: backendValues.mappingMode
                implicitWidth: StudioTheme.Values.singleControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("U Tiling")
            tooltip: qsTr("Controls how the texture is mapped when the U scaling value is greater than 1.")
        }

        SecondColumnLayout {
            ComboBox {
                scope: "Texture"
                model: ["Unknown", "ClampToEdge", "MirroredRepeat", "Repeat"]
                backendValue: backendValues.tilingModeHorizontal
                implicitWidth: StudioTheme.Values.singleControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("V Tiling")
            tooltip: qsTr("Controls how the texture is mapped when the V scaling value is greater than 1.")
        }

        SecondColumnLayout {
            ComboBox {
                scope: "Texture"
                model: ["Unknown", "ClampToEdge", "MirroredRepeat", "Repeat"]
                backendValue: backendValues.tilingModeVertical
                implicitWidth: StudioTheme.Values.singleControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("UV Rotation")
            tooltip: qsTr("Rotates the texture around the pivot point.")
        }

        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 360
                decimals: 0
                backendValue: backendValues.rotationUV
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Position")
        }

        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.positionU
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

            ControlLabel {
                text: "U"
                tooltip: qsTr("Offsets the U coordinate mapping from left to right.")
            }

            Spacer { implicitWidth: StudioTheme.Values.controlGap }

            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.positionV
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

            ControlLabel {
                text: "V"
                tooltip: qsTr("Offsets the V coordinate mapping from bottom to top.")
            }

            ExpandingSpacer {}
        }

        PropertyLabel {
            text: qsTr("Pivot")
        }

        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.pivotU
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

            ControlLabel {
                text: "U"
                tooltip: qsTr("Sets the pivot U position.")
            }

            Spacer { implicitWidth: StudioTheme.Values.controlGap }

            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.pivotV
                implicitWidth: StudioTheme.Values.twoControlColumnWidth
                               + StudioTheme.Values.actionIndicatorWidth
            }

            Spacer { implicitWidth: StudioTheme.Values.controlLabelGap }

            ControlLabel {
                text: "V"
                tooltip: qsTr("Sets the pivot V position.")
            }

            ExpandingSpacer {}
        }
    }
}
