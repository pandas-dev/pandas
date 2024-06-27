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

Rectangle {
    // public API
    // TODO 5.15: required property
    property var document: undefined
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
    function goToPage(page) { goToLocation(page, Qt.point(0, 0), 0) }
    function goToLocation(page, location, zoom) {
        if (zoom > 0)
            root.renderScale = zoom
        navigationStack.push(page, location, zoom)
    }

    // page scaling
    property real renderScale: 1
    property alias sourceSize: image.sourceSize
    function resetScale() {
        image.sourceSize.width = 0
        image.sourceSize.height = 0
        root.x = 0
        root.y = 0
        root.scale = 1
    }
    function scaleToWidth(width, height) {
        var halfRotation = Math.abs(root.rotation % 180)
        image.sourceSize = Qt.size((halfRotation > 45 && halfRotation < 135) ? height : width, 0)
        root.x = 0
        root.y = 0
        image.centerInSize = Qt.size(width, height)
        image.centerOnLoad = true
        image.vCenterOnLoad = (halfRotation > 45 && halfRotation < 135)
        root.scale = 1
    }
    function scaleToPage(width, height) {
        var windowAspect = width / height
        var halfRotation = Math.abs(root.rotation % 180)
        var pagePointSize = document.pagePointSize(navigationStack.currentPage)
        if (halfRotation > 45 && halfRotation < 135) {
            // rotated 90 or 270º
            var pageAspect = pagePointSize.height / pagePointSize.width
            if (windowAspect > pageAspect) {
                image.sourceSize = Qt.size(height, 0)
            } else {
                image.sourceSize = Qt.size(0, width)
            }
        } else {
            var pageAspect = pagePointSize.width / pagePointSize.height
            if (windowAspect > pageAspect) {
                image.sourceSize = Qt.size(0, height)
            } else {
                image.sourceSize = Qt.size(width, 0)
            }
        }
        image.centerInSize = Qt.size(width, height)
        image.centerOnLoad = true
        image.vCenterOnLoad = true
        root.scale = 1
    }

    // text search
    property alias searchModel: searchModel
    property alias searchString: searchModel.searchString
    function searchBack() { --searchModel.currentResult }
    function searchForward() { ++searchModel.currentResult }

    // implementation
    id: root
    width: image.width
    height: image.height

    PdfSelection {
        id: selection
        document: root.document
        page: navigationStack.currentPage
        fromPoint: Qt.point(textSelectionDrag.centroid.pressPosition.x / image.pageScale, textSelectionDrag.centroid.pressPosition.y / image.pageScale)
        toPoint: Qt.point(textSelectionDrag.centroid.position.x / image.pageScale, textSelectionDrag.centroid.position.y / image.pageScale)
        hold: !textSelectionDrag.active && !tapHandler.pressed
    }

    PdfSearchModel {
        id: searchModel
        document: root.document === undefined ? null : root.document
        onCurrentPageChanged: root.goToPage(currentPage)
    }

    PdfNavigationStack {
        id: navigationStack
        onCurrentPageChanged: searchModel.currentPage = currentPage
        // TODO onCurrentLocationChanged: position currentLocation.x and .y in middle // currentPageChanged() MUST occur first!
        onCurrentZoomChanged: root.renderScale = currentZoom
        // TODO deal with horizontal location (need WheelHandler or Flickable probably)
    }

    Image {
        id: image
        currentFrame: navigationStack.currentPage
        source: document.status === PdfDocument.Ready ? document.source : ""
        asynchronous: true
        fillMode: Image.PreserveAspectFit
        property bool centerOnLoad: false
        property bool vCenterOnLoad: false
        property size centerInSize
        property real pageScale: image.paintedWidth / document.pagePointSize(navigationStack.currentPage).width
        function reRenderIfNecessary() {
            var newSourceWidth = image.sourceSize.width * root.scale
            var ratio = newSourceWidth / image.sourceSize.width
            if (ratio > 1.1 || ratio < 0.9) {
                image.sourceSize.width = newSourceWidth
                image.sourceSize.height = 0
                root.scale = 1
            }
        }
        onStatusChanged:
            if (status == Image.Ready && centerOnLoad) {
                root.x = (centerInSize.width - image.implicitWidth) / 2
                root.y = vCenterOnLoad ? (centerInSize.height - image.implicitHeight) / 2 : 0
                centerOnLoad = false
                vCenterOnLoad = false
            }
    }
    onRenderScaleChanged: {
        image.sourceSize.width = document.pagePointSize(navigationStack.currentPage).width * renderScale
        image.sourceSize.height = 0
        root.scale = 1
    }

    Shape {
        anchors.fill: parent
        opacity: 0.25
        visible: image.status === Image.Ready
        ShapePath {
            strokeWidth: 1
            strokeColor: "cyan"
            fillColor: "steelblue"
            scale: Qt.size(image.pageScale, image.pageScale)
            PathMultiline {
                paths: searchModel.currentPageBoundingPolygons
            }
        }
        ShapePath {
            strokeWidth: 1
            strokeColor: "orange"
            fillColor: "cyan"
            scale: Qt.size(image.pageScale, image.pageScale)
            PathMultiline {
                paths: searchModel.currentResultBoundingPolygons
            }
        }
        ShapePath {
            fillColor: "orange"
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
        delegate: Rectangle {
            color: "transparent"
            border.color: "lightgrey"
            x: rect.x * image.pageScale
            y: rect.y * image.pageScale
            width: rect.width * image.pageScale
            height: rect.height * image.pageScale
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

    PinchHandler {
        id: pinch
        minimumScale: 0.1
        maximumScale: 10
        minimumRotation: 0
        maximumRotation: 0
        onActiveChanged: if (!active) image.reRenderIfNecessary()
        grabPermissions: PinchHandler.TakeOverForbidden // don't allow takeover if pinch has started
    }
    DragHandler {
        id: pageMovingTouchDrag
        acceptedDevices: PointerDevice.TouchScreen
    }
    DragHandler {
        id: pageMovingMiddleMouseDrag
        acceptedDevices: PointerDevice.Mouse | PointerDevice.Stylus
        acceptedButtons: Qt.MiddleButton
        snapMode: DragHandler.NoSnap
    }
    DragHandler {
        id: textSelectionDrag
        acceptedDevices: PointerDevice.Mouse | PointerDevice.Stylus
        target: null
    }
    TapHandler {
        id: tapHandler
        acceptedDevices: PointerDevice.Mouse | PointerDevice.Stylus
    }
    // prevent it from being scrolled out of view
    BoundaryRule on x {
        minimum: 100 - root.width
        maximum: root.parent.width - 100
    }
    BoundaryRule on y {
        minimum: 100 - root.height
        maximum: root.parent.height - 100
    }
}
