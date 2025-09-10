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
import HelperWidgets 2.0
import QtQuick.Layouts 1.12

Section {
    caption: qsTr("Image")
    width: parent.width
    SectionLayout {
        Label {
            text: qsTr("Source")
            tooltip: qsTr("Holds the location of an image file containing the data used by the texture.")
        }
        SecondColumnLayout {
            UrlChooser {
                backendValue: backendValues.source
            }
        }

        Label {
            text: qsTr("Source Item")
            tooltip: qsTr("Defines an item to be used as the source of the texture.")
        }
        SecondColumnLayout {
            IdComboBox {
                typeFilter: "QtQuick.Item"
                Layout.fillWidth: true
                backendValue: backendValues.sourceItem
            }
        }

        Label {
            text: qsTr("U Scale")
            tooltip: qsTr("Defines how to scale the U texture coordinate when mapping to UV coordinates of a mesh.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: 0
                realDragRange: 10
                decimals: 0
                backendValue: backendValues.scaleU
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("V Scale")
            tooltip: qsTr("Defines how to scale the V texture coordinate when mapping to UV coordinates of a mesh.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: 0
                realDragRange: 10
                decimals: 0
                backendValue: backendValues.scaleV
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("Texture Mapping")
            tooltip: qsTr("Defines which method of mapping to use when sampling this texture.")
        }
        SecondColumnLayout {
            ComboBox {
                scope: "Texture"
                model: ["UV", "Environment", "LightProbe"]
                backendValue: backendValues.mappingMode
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("U Tiling")
            tooltip: qsTr("Controls how the texture is mapped when the U scaling value is greater than 1.")
        }
        SecondColumnLayout {
            ComboBox {
                scope: "Texture"
                model: ["Unknown", "ClampToEdge", "MirroredRepeat", "Repeat"]
                backendValue: backendValues.tilingModeHorizontal
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("V Tiling")
            tooltip: qsTr("Controls how the texture is mapped when the V scaling value is greater than 1.")
        }
        SecondColumnLayout {
            ComboBox {
                scope: "Texture"
                model: ["Unknown", "ClampToEdge", "MirroredRepeat", "Repeat"]
                backendValue: backendValues.tilingModeVertical
                Layout.fillWidth: true
            }
        }

        Label {
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
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("U Position")
            tooltip: qsTr("Offsets the U coordinate mapping from left to right.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.positionU
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("V Position")
            tooltip: qsTr("Offsets the V coordinate mapping from bottom to top.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.positionV
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("U Pivot")
            tooltip: qsTr("Sets the pivot U position.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.pivotU
                Layout.fillWidth: true
            }
        }

        Label {
            text: qsTr("V Pivot")
            tooltip: qsTr("Sets the pivot V position.")
        }
        SecondColumnLayout {
            SpinBox {
                maximumValue: 999999
                minimumValue: -999999
                realDragRange: 2
                decimals: 2
                stepSize: 0.1
                backendValue: backendValues.pivotV
                Layout.fillWidth: true
            }
        }

    }

}
