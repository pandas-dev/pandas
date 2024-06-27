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

ComboBox {
    id: comboBox

    property alias typeFilter: itemFilterModel.typeFilter

    manualMapping: true
    editable: true
    model: comboBox.addDefaultItem(itemFilterModel.itemModel)

    textInput.validator: RegExpValidator { regExp: /(^$|^[a-z_]\w*)/ }

    ItemFilterModel {
        id: itemFilterModel
        modelNodeBackendProperty: modelNodeBackend
    }

    property string defaultItem: qsTr("None")
    property string textValue: comboBox.backendValue.expression
    property bool block: false
    property bool dirty: true
    property var editRegExp: /^[a-z_]\w*/

    onTextValueChanged: {
        if (comboBox.block)
            return

        comboBox.setCurrentText(comboBox.textValue)
    }
    onModelChanged: comboBox.setCurrentText(comboBox.textValue)
    onCompressedActivated: function(index, reason) { comboBox.handleActivate(index) }
    Component.onCompleted: comboBox.setCurrentText(comboBox.textValue)

    onEditTextChanged: {
        comboBox.dirty = true
        colorLogic.errorState = !(editRegExp.exec(comboBox.editText) !== null
                                  || comboBox.editText === parenthesize(defaultItem))
    }
    onFocusChanged: {
        if (comboBox.dirty)
           comboBox.handleActivate(comboBox.currentIndex)
    }

    function handleActivate(index)
    {
        if (!comboBox.__isCompleted || comboBox.backendValue === undefined)
            return

        var cText = (index === -1) ? comboBox.editText : comboBox.textAt(index)
        comboBox.block = true
        comboBox.setCurrentText(cText)
        comboBox.block = false
    }

    function setCurrentText(text)
    {
        if (!comboBox.__isCompleted || comboBox.backendValue === undefined)
            return

        comboBox.currentIndex = comboBox.find(text)

        if (text === "") {
            comboBox.currentIndex = 0
            comboBox.editText = parenthesize(comboBox.defaultItem)
        } else {
            if (comboBox.currentIndex === -1)
                comboBox.editText = text
            else if (comboBox.currentIndex === 0)
                comboBox.editText = parenthesize(comboBox.defaultItem)
        }

        if (comboBox.currentIndex === 0) {
            comboBox.backendValue.resetValue()
        } else {
            if (comboBox.backendValue.expression !== comboBox.editText)
                comboBox.backendValue.expression = comboBox.editText
        }
        comboBox.dirty = false
    }

    function addDefaultItem(arr)
    {
        var copy = arr.slice()
        copy.unshift(parenthesize(comboBox.defaultItem))
        return copy
    }

    function parenthesize(value)
    {
        return "[" + value + "]"
    }
}
