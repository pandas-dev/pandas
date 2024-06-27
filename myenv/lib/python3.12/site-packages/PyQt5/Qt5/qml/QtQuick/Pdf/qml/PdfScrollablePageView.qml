/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtPDF module of the Qt Toolkit.
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
import QtQuick 2.14
import QtQuick.Controls 2.14
import QtQuick.Pdf 5.15
import QtQuick.Shapes 1.14
import Qt.labs.animation 1.0

Flickable {
    // public API
    // TODO 5.15: required property
    property var document: undefined
    property bool debug: false
    property alias status: image.status

    property alias selectedText: selection.text
    function selectAll() {
        selection.selectAll()
    }
    function copySelectionToClipboard() {
        selection.copyToClipboard()
    }

    // page navigation
    property alias currentPage: navigationStack.currentPage
    property alias backEnabled: navigationStack.backAvailable
    property alias forwardEnabled: navigationStack.forwardAvailable
    function back() { navigationStack.back() }
    function forward() { navigationStack.forward() }
    function goToPage(page) {
        if (page === navigationStack.currentPage)
            return
        goToLocation(page, Qt.point(0, 0), 0)
    }
    function goToLocation(page, location, zoom) {
        if (zoom > 0)
            root.renderScale = zoom
        navigationStack.push(page, location, zoom)
    }

    // page scaling
    property real renderScale: 1
    property real pageRotation: 0
    property alias sourceSize: image.sourceSize
    function resetScale() {
        paper.scale = 1
        root.renderScale = 1
    }
    function scaleToWidth(width, height) {
        var pagePointSize = document.pagePointSize(navigationStack.currentPage)
        root.renderScale = root.width / (paper.rot90 ? pagePointSize.height : pagePointSize.width)
        if (debug)
            console.log("scaling", pagePointSize, "to fit", root.width, "rotated?", paper.rot90, "scale", root.renderScale)
        root.contentX = 0
        root.contentY = 0
    }
    function scaleToPage(width, height) {
        var pagePointSize = document.pagePointSize(navigationStack.currentPage)
        root.renderScale = Math.min(
                    root.width / (paper.rot90 ? pagePointSize.height : pagePointSize.width),
                    root.height / (paper.rot90 ? pagePointSize.width : pagePointSize.height) )
        root.contentX = 0
        root.contentY = 0
    }

    // text search
    property alias searchModel: searchModel
    property alias searchString: searchModel.searchString
    function searchBack() { --searchModel.currentResult }
    function searchForward() { ++searchModel.currentResult }

    // implementation
    id: root
    PdfStyle { id: style }
    contentWidth: paper.width
    contentHeight: paper.height
    ScrollBar.vertical: ScrollBar {
        onActiveChanged:
            if (!active ) {
                var currentLocation = Qt.point((root.contentX + root.width / 2) / root.renderScale,
                                               (root.contentY + root.height / 2) / root.renderScale)
                navigationStack.update(navigationStack.currentPage, currentLocation, root.renderScale)
            }
    }
    ScrollBar.horizontal: ScrollBar {
        onActiveChanged:
            if (!active ) {
                var currentLocation = Qt.point((root.contentX + root.width / 2) / root.renderScale,
                                               (root.contentY + root.height / 2) / root.renderScale)
                navigationStack.update(navigationStack.currentPage, currentLocation, root.renderScale)
            }
    }

    onRenderScaleChanged: {
        image.sourceSize.width = document.pagePointSize(navigationStack.currentPage).width * renderScale
        image.sourceSize.height = 0
        paper.scale = 1
        var currentLocation = Qt.point((root.contentX + root.width / 2) / root.renderScale,
                                       (root.contentY + root.height / 2) / root.renderScale)
        navigationStack.update(navigationStack.currentPage, currentLocation, root.renderScale)
    }

    PdfSearchModel {
        id: searchModel
        document: root.document === undefined ? null : root.document
        // TODO maybe avoid jumping if the result is already fully visible in the viewport
        onCurrentResultBoundingRectChanged: root.goToLocation(currentPage,
            Qt.point(currentResultBoundingRect.x, currentResultBoundingRect.y), 0)
    }

    PdfNavigationStack {
        id: navigationStack
        onJumped: {
            root.renderScale = zoom
            var dx = Math.max(0, location.x * root.renderScale - root.width / 2) - root.contentX
            var dy = Math.max(0, location.y * root.renderScale - root.height / 2) - root.contentY
            // don't jump if location is in the viewport already, i.e. if the "error" between desired and actual contentX/Y is small
            if (Math.abs(dx) > root.width / 3)
                root.contentX += dx
            if (Math.abs(dy) > root.height / 3)
                root.contentY += dy
            if (root.debug) {
                console.log("going to zoom", zoom, "loc", location,
                            "on page", page, "ended up @", root.contentX + ", " + root.contentY)
            }
        }
        onCurrentPageChanged: searchModel.currentPage = currentPage
    }

    Rectangle {
        id: paper
        width: rot90 ? image.height : image.width
        height: rot90 ? image.width : image.height
        property real rotationModulus: Math.abs(root.pageRotation % 180)
        property bool rot90: rotationModulus > 45 && rotationModulus < 135

        Image {
            id: image
            currentFrame: navigationStack.currentPage
            source: document.status === PdfDocument.Ready ? document.source : ""
            asynchronous: true
            fillMode: Image.PreserveAspectFit
            rotation: root.pageRotation
            anchors.centerIn: parent
            property real pageScale: image.paintedWidth / document.pagePointSize(navigationStack.currentPage).width

            Shape {
                anchors.fill: parent
                visible: image.status === Image.Ready
                ShapePath {
                    strokeWidth: -1
                    fillColor: style.pageSearchResultsColor
                    scale: Qt.size(image.pageScale, image.pageScale)
                    PathMultiline {
                        paths: searchModel.currentPageBoundingPolygons
                    }
                }
                ShapePath {
                    strokeWidth: style.currentSearchResultStrokeWidth
                    strokeColor: style.currentSearchResultStrokeColor
                    fillColor: "transparent"
                    scale: Qt.size(image.pageScale, image.pageScale)
                    PathMultiline {
                        paths: searchModel.currentResultBoundingPolygons
                    }
                }
                ShapePath {
                    fillColor: style.selectionColor
                    scale: Qt.size(image.pageScale, image.pageScale)
                    PathMultiline {
                        paths: selection.geometry
                    }
                }
            }

            Repeater {
                model: PdfLinkModel {
                    id: linkModel
                    document: root.document
                    page: navigationStack.currentPage
                }
                delegate: Shape {
                    x: rect.x * image.pageScale
                    y: rect.y * image.pageScale
                    width: rect.width * image.pageScale
                    height: rect.height * image.pageScale
                    ShapePath {
                        strokeWidth: style.linkUnderscoreStrokeWidth
                        strokeColor: style.linkUnderscoreColor
                        strokeStyle: style.linkUnderscoreStrokeStyle
                        dashPattern: style.linkUnderscoreDashPattern
                        startX: 0; startY: height
                        PathLine { x: width; y: height }
                    }
                    MouseArea { // TODO switch to TapHandler / HoverHandler in 5.15
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                        onClicked: {
                            if (page >= 0)
                                navigationStack.push(page, Qt.point(0, 0), root.renderScale)
                            else
                                Qt.openUrlExternally(url)
                        }
                    }
                }
            }
            DragHandler {
                id: textSelectionDrag
                acceptedDevices: PointerDevice.Mouse | PointerDevice.Stylus
                target: null
            }
            TapHandler {
                id: mouseClickHandler
                acceptedDevices: PointerDevice.Mouse | PointerDevice.Stylus
            }
            TapHandler {
                id: touchTapHandler
                acceptedDevices: PointerDevice.TouchScreen
                onTapped: {
                    selection.clear()
                    selection.focus = true
                }
            }
        }

        PdfSelection {
            id: selection
            anchors.fill: parent
            document: root.document
            page: navigationStack.currentPage
            renderScale: image.pageScale
            fromPoint: textSelectionDrag.centroid.pressPosition
            toPoint: textSelectionDrag.centroid.position
            hold: !textSelectionDrag.active && !mouseClickHandler.pressed
            focus: true
        }

        PinchHandler {
            id: pinch
            minimumScale: 0.1
            maximumScale: root.renderScale < 4 ? 2 : 1
            minimumRotation: 0
            maximumRotation: 0
            enabled: image.sourceSize.width < 5000
            onActiveChanged:
                if (!active) {
                    var centroidInPoints = Qt.point(pinch.centroid.position.x / root.renderScale,
                                                    pinch.centroid.position.y / root.renderScale)
                    var centroidInFlickable = root.mapFromItem(paper, pinch.centroid.position.x, pinch.centroid.position.y)
                    var newSourceWidth = image.sourceSize.width * paper.scale
                    var ratio = newSourceWidth / image.sourceSize.width
                    if (root.debug)
                        console.log("pinch ended with centroid", pinch.centroid.position, centroidInPoints, "wrt flickable", centroidInFlickable,
                                    "page at", paper.x.toFixed(2), paper.y.toFixed(2),
                                    "contentX/Y were", root.contentX.toFixed(2), root.contentY.toFixed(2))
                    if (ratio > 1.1 || ratio < 0.9) {
                        var centroidOnPage = Qt.point(centroidInPoints.x * root.renderScale * ratio, centroidInPoints.y * root.renderScale * ratio)
                        paper.scale = 1
                        paper.x = 0
                        paper.y = 0
                        root.contentX = centroidOnPage.x - centroidInFlickable.x
                        root.contentY = centroidOnPage.y - centroidInFlickable.y
                        root.renderScale *= ratio // onRenderScaleChanged calls navigationStack.update() so we don't need to here
                        if (root.debug)
                            console.log("contentX/Y adjusted to", root.contentX.toFixed(2), root.contentY.toFixed(2))
                    } else {
                        paper.x = 0
                        paper.y = 0
                    }
                }
            grabPermissions: PointerHandler.CanTakeOverFromAnything
        }
    }
}
