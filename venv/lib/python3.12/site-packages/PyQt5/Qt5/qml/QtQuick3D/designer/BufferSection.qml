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
import QtQuick.Layouts 1.15
import HelperWidgets 2.0
import StudioTheme 1.0 as StudioTheme

Column {
    width: parent.width

    Section {
        caption: qsTr("Buffer")
        width: parent.width

        SectionLayout {
            PropertyLabel {
                text: qsTr("Name")
                tooltip: qsTr("Buffer name.")
            }

            SecondColumnLayout {
                LineEdit {
                    backendValue: backendValues.name
                    showTranslateCheckBox: false
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                    width: implicitWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Format")
                tooltip: qsTr("Format of the buffer.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "Buffer"
                    model: ["Unknown", "R8", "R16", "R16F", "R32I", "R32UI", "R32F", "RG8", "RGBA8",
                            "RGB8", "SRGB8", "SRGB8A8", "RGB565", "RGBA16F", "RG16F", "RG32F",
                            "RGB32F", "RGBA32F", "R11G11B10", "RGB9E5", "Depth16", "Depth24",
                            "Depth32", "Depth24Stencil8"]
                    backendValue: backendValues.format
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Filter")
                tooltip: qsTr("Texture filter for the buffer.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "Buffer"
                    model: ["Unknown", "Nearest", "Linear"]
                    backendValue: backendValues.textureFilterOperation
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Coordinate Operation")
                tooltip: qsTr("Texture coordinate operation for the buffer.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "Buffer"
                    model: ["Unknown", "ClampToEdge", "MirroredRepeat", "Repeat"]
                    backendValue: backendValues.textureCoordOperation
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Allocation Flags")
                tooltip: qsTr("Allocation flags for the buffer.")
            }

            SecondColumnLayout {
                ComboBox {
                    scope: "Buffer"
                    model: ["None", "SceneLifetime"]
                    backendValue: backendValues.bufferFlags
                    implicitWidth: StudioTheme.Values.singleControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }

            PropertyLabel {
                text: qsTr("Size Multiplier")
                tooltip: qsTr("Defines the size multiplier for the buffer.")
            }

            SecondColumnLayout {
                SpinBox {
                    maximumValue: 10000
                    minimumValue: 0
                    decimals: 2
                    realDragRange: 30
                    backendValue: backendValues.sizeMultiplier
                    implicitWidth: StudioTheme.Values.twoControlColumnWidth
                                   + StudioTheme.Values.actionIndicatorWidth
                }

                ExpandingSpacer {}
            }
        }
    }
}
