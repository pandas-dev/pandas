#############################################################################
##
## Copyright (C) 2019 Riverbank Computing Limited.
## Copyright (C) 2006 Thorsten Marek.
## All right reserved.
##
## This file is part of PyQt.
##
## You may use this file under the terms of the GPL v2 or the revised BSD
## license as follows:
##
## "Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of the Riverbank Computing Limited nor the names
##     of its contributors may be used to endorse or promote products
##     derived from this software without specific prior written
##     permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
##
#############################################################################


import sys
import logging
import os
import re
from xml.etree.ElementTree import parse, SubElement

from .objcreator import QObjectCreator
from .properties import Properties


logger = logging.getLogger(__name__)
DEBUG = logger.debug

QtCore = None
QtWidgets = None


def _parse_alignment(alignment):
    """ Convert a C++ alignment to the corresponding flags. """

    align_flags = None
    for qt_align in alignment.split('|'):
        _, qt_align = qt_align.split('::')
        align = getattr(QtCore.Qt, qt_align)

        if align_flags is None:
            align_flags = align
        else:
            align_flags |= align

    return align_flags


def _layout_position(elem):
    """ Return either (), (0, alignment), (row, column, rowspan, colspan) or
    (row, column, rowspan, colspan, alignment) depending on the type of layout
    and its configuration.  The result will be suitable to use as arguments to
    the layout.
    """

    row = elem.attrib.get('row')
    column = elem.attrib.get('column')
    alignment = elem.attrib.get('alignment')

    # See if it is a box layout.
    if row is None or column is None:
        if alignment is None:
            return ()

        return (0, _parse_alignment(alignment))

    # It must be a grid or a form layout.
    row = int(row)
    column = int(column)

    rowspan = int(elem.attrib.get('rowspan', 1))
    colspan = int(elem.attrib.get('colspan', 1))

    if alignment is None:
        return (row, column, rowspan, colspan)

    return (row, column, rowspan, colspan, _parse_alignment(alignment))


class WidgetStack(list):
    topwidget = None
    def push(self, item):
        DEBUG("push %s %s" % (item.metaObject().className(),
                              item.objectName()))
        self.append(item)
        if isinstance(item, QtWidgets.QWidget):
            self.topwidget = item

    def popLayout(self):
        layout = list.pop(self)
        DEBUG("pop layout %s %s" % (layout.metaObject().className(),
                                    layout.objectName()))
        return layout

    def popWidget(self):
        widget = list.pop(self)
        DEBUG("pop widget %s %s" % (widget.metaObject().className(),
                                    widget.objectName()))
        for item in reversed(self):
            if isinstance(item, QtWidgets.QWidget):
                self.topwidget = item
                break
        else:
            self.topwidget = None
        DEBUG("new topwidget %s" % (self.topwidget,))
        return widget

    def peek(self):
        return self[-1]

    def topIsLayout(self):
        return isinstance(self[-1], QtWidgets.QLayout)

    def topIsLayoutWidget(self):
        # A plain QWidget is a layout widget unless it's parent is a
        # QMainWindow or a container widget.  Note that the corresponding uic
        # test is a little more complicated as it involves features not
        # supported by pyuic.

        if type(self[-1]) is not QtWidgets.QWidget:
            return False

        if len(self) < 2:
            return False

        parent = self[-2]

        return isinstance(parent, QtWidgets.QWidget) and type(parent) not in (
                QtWidgets.QMainWindow,
                QtWidgets.QStackedWidget,
                QtWidgets.QToolBox,
                QtWidgets.QTabWidget,
                QtWidgets.QScrollArea,
                QtWidgets.QMdiArea,
                QtWidgets.QWizard,
                QtWidgets.QDockWidget)


class ButtonGroup(object):
    """ Encapsulate the configuration of a button group and its implementation.
    """

    def __init__(self):
        """ Initialise the button group. """

        self.exclusive = True
        self.object = None


class UIParser(object):    
    def __init__(self, qtcore_module, qtgui_module, qtwidgets_module, creatorPolicy):
        self.factory = QObjectCreator(creatorPolicy)
        self.wprops = Properties(self.factory, qtcore_module, qtgui_module,
                qtwidgets_module)
        
        global QtCore, QtWidgets
        QtCore = qtcore_module
        QtWidgets = qtwidgets_module
        
        self.reset()

    def uniqueName(self, name):
        """UIParser.uniqueName(string) -> string

        Create a unique name from a string.
        >>> p = UIParser(QtCore, QtGui, QtWidgets)
        >>> p.uniqueName("foo")
        'foo'
        >>> p.uniqueName("foo")
        'foo1'
        """
        try:
            suffix = self.name_suffixes[name]
        except KeyError:
            self.name_suffixes[name] = 0
            return name

        suffix += 1
        self.name_suffixes[name] = suffix

        return "%s%i" % (name, suffix)

    def reset(self):
        try: self.wprops.reset()
        except AttributeError: pass
        self.toplevelWidget = None
        self.stack = WidgetStack()
        self.name_suffixes = {}
        self.defaults = {'spacing': -1, 'margin': -1}
        self.actions = []
        self.currentActionGroup = None
        self.resources = []
        self.button_groups = {}

    def setupObject(self, clsname, parent, branch, is_attribute=True):
        name = self.uniqueName(branch.attrib.get('name') or clsname[1:].lower())

        if parent is None:
            args = ()
        else:
            args = (parent, )

        obj = self.factory.createQObject(clsname, name, args, is_attribute)

        self.wprops.setProperties(obj, branch)
        obj.setObjectName(name)

        if is_attribute:
            setattr(self.toplevelWidget, name, obj)

        return obj

    def getProperty(self, elem, name):
        for prop in elem.findall('property'):
            if prop.attrib['name'] == name:
                return prop

        return None

    def createWidget(self, elem):
        self.column_counter = 0
        self.row_counter = 0
        self.item_nr = 0
        self.itemstack = []
        self.sorting_enabled = None

        widget_class = elem.attrib['class'].replace('::', '.')
        if widget_class == 'Line':
            widget_class = 'QFrame'
        
        # Ignore the parent if it is a container.
        parent = self.stack.topwidget
        if isinstance(parent, (QtWidgets.QDockWidget, QtWidgets.QMdiArea,
                               QtWidgets.QScrollArea, QtWidgets.QStackedWidget,
                               QtWidgets.QToolBox, QtWidgets.QTabWidget,
                               QtWidgets.QWizard)):
            parent = None

        self.stack.push(self.setupObject(widget_class, parent, elem))

        if isinstance(self.stack.topwidget, QtWidgets.QTableWidget):
            if self.getProperty(elem, 'columnCount') is None:
                self.stack.topwidget.setColumnCount(len(elem.findall("column")))

            if self.getProperty(elem, 'rowCount') is None:
                self.stack.topwidget.setRowCount(len(elem.findall("row")))

        self.traverseWidgetTree(elem)
        widget = self.stack.popWidget()

        if isinstance(widget, QtWidgets.QTreeView):
            self.handleHeaderView(elem, "header", widget.header())

        elif isinstance(widget, QtWidgets.QTableView):
            self.handleHeaderView(elem, "horizontalHeader",
                    widget.horizontalHeader())
            self.handleHeaderView(elem, "verticalHeader",
                    widget.verticalHeader())

        elif isinstance(widget, QtWidgets.QAbstractButton):
            bg_i18n = self.wprops.getAttribute(elem, "buttonGroup")
            if bg_i18n is not None:
                # This should be handled properly in case the problem arises
                # elsewhere as well.
                try:
                    # We are compiling the .ui file.
                    bg_name = bg_i18n.string
                except AttributeError:
                    # We are loading the .ui file.
                    bg_name = bg_i18n

				# Designer allows the creation of .ui files without explicit
				# button groups, even though uic then issues warnings.  We
				# handle it in two stages by first making sure it has a name
				# and then making sure one exists with that name.
                if not bg_name:
                    bg_name = 'buttonGroup'

                try:
                    bg = self.button_groups[bg_name]
                except KeyError:
                    bg = self.button_groups[bg_name] = ButtonGroup()

                if bg.object is None:
                    bg.object = self.factory.createQObject("QButtonGroup",
                            bg_name, (self.toplevelWidget, ))
                    setattr(self.toplevelWidget, bg_name, bg.object)

                    bg.object.setObjectName(bg_name)

                    if not bg.exclusive:
                        bg.object.setExclusive(False)

                bg.object.addButton(widget)

        if self.sorting_enabled is not None:
            widget.setSortingEnabled(self.sorting_enabled)
            self.sorting_enabled = None
        
        if self.stack.topIsLayout():
            lay = self.stack.peek()
            lp = elem.attrib['layout-position']

            if isinstance(lay, QtWidgets.QFormLayout):
                lay.setWidget(lp[0], self._form_layout_role(lp), widget)
            else:
                lay.addWidget(widget, *lp)

        topwidget = self.stack.topwidget

        if isinstance(topwidget, QtWidgets.QToolBox):
            icon = self.wprops.getAttribute(elem, "icon")
            if icon is not None:
                topwidget.addItem(widget, icon, self.wprops.getAttribute(elem, "label"))
            else:
                topwidget.addItem(widget, self.wprops.getAttribute(elem, "label"))

            tooltip = self.wprops.getAttribute(elem, "toolTip")
            if tooltip is not None:
                topwidget.setItemToolTip(topwidget.indexOf(widget), tooltip)
                
        elif isinstance(topwidget, QtWidgets.QTabWidget):
            icon = self.wprops.getAttribute(elem, "icon")
            if icon is not None:
                topwidget.addTab(widget, icon, self.wprops.getAttribute(elem, "title"))
            else:
                topwidget.addTab(widget, self.wprops.getAttribute(elem, "title"))

            tooltip = self.wprops.getAttribute(elem, "toolTip")
            if tooltip is not None:
                topwidget.setTabToolTip(topwidget.indexOf(widget), tooltip)
            
        elif isinstance(topwidget, QtWidgets.QWizard):
            topwidget.addPage(widget)
            
        elif isinstance(topwidget, QtWidgets.QStackedWidget):
            topwidget.addWidget(widget)
            
        elif isinstance(topwidget, (QtWidgets.QDockWidget, QtWidgets.QScrollArea)):
            topwidget.setWidget(widget)
            
        elif isinstance(topwidget, QtWidgets.QMainWindow):
            if type(widget) == QtWidgets.QWidget:
                topwidget.setCentralWidget(widget)
            elif isinstance(widget, QtWidgets.QToolBar):
                tbArea = self.wprops.getAttribute(elem, "toolBarArea")

                if tbArea is None:
                    topwidget.addToolBar(widget)
                else:
                    topwidget.addToolBar(tbArea, widget)

                tbBreak = self.wprops.getAttribute(elem, "toolBarBreak")

                if tbBreak:
                    topwidget.insertToolBarBreak(widget)

            elif isinstance(widget, QtWidgets.QMenuBar):
                topwidget.setMenuBar(widget)
            elif isinstance(widget, QtWidgets.QStatusBar):
                topwidget.setStatusBar(widget)
            elif isinstance(widget, QtWidgets.QDockWidget):
                dwArea = self.wprops.getAttribute(elem, "dockWidgetArea")
                topwidget.addDockWidget(QtCore.Qt.DockWidgetArea(dwArea),
                        widget)

    def handleHeaderView(self, elem, name, header):
        value = self.wprops.getAttribute(elem, name + "Visible")
        if value is not None:
            header.setVisible(value)

        value = self.wprops.getAttribute(elem, name + "CascadingSectionResizes")
        if value is not None:
            header.setCascadingSectionResizes(value)

        value = self.wprops.getAttribute(elem, name + "DefaultSectionSize")
        if value is not None:
            header.setDefaultSectionSize(value)

        value = self.wprops.getAttribute(elem, name + "HighlightSections")
        if value is not None:
            header.setHighlightSections(value)

        value = self.wprops.getAttribute(elem, name + "MinimumSectionSize")
        if value is not None:
            header.setMinimumSectionSize(value)

        value = self.wprops.getAttribute(elem, name + "ShowSortIndicator")
        if value is not None:
            header.setSortIndicatorShown(value)

        value = self.wprops.getAttribute(elem, name + "StretchLastSection")
        if value is not None:
            header.setStretchLastSection(value)

    def createSpacer(self, elem):
        width = elem.findtext("property/size/width")
        height = elem.findtext("property/size/height")

        if width is None or height is None:
            size_args = ()
        else:
            size_args = (int(width), int(height))

        sizeType = self.wprops.getProperty(elem, "sizeType",
                QtWidgets.QSizePolicy.Expanding)

        policy = (QtWidgets.QSizePolicy.Minimum, sizeType)

        if self.wprops.getProperty(elem, "orientation") == QtCore.Qt.Horizontal:
            policy = policy[1], policy[0]

        spacer = self.factory.createQObject("QSpacerItem",
                self.uniqueName("spacerItem"), size_args + policy,
                is_attribute=False)

        if self.stack.topIsLayout():
            lay = self.stack.peek()
            lp = elem.attrib['layout-position']

            if isinstance(lay, QtWidgets.QFormLayout):
                lay.setItem(lp[0], self._form_layout_role(lp), spacer)
            else:
                lay.addItem(spacer, *lp)

    def createLayout(self, elem):
        # We use an internal property to handle margins which will use separate
        # left, top, right and bottom margins if they are found to be
        # different.  The following will select, in order of preference,
        # separate margins, the same margin in all directions, and the default
        # margin.
        margin = -1 if self.stack.topIsLayout() else self.defaults['margin']
        margin = self.wprops.getProperty(elem, 'margin', margin)
        left = self.wprops.getProperty(elem, 'leftMargin', margin)
        top = self.wprops.getProperty(elem, 'topMargin', margin)
        right = self.wprops.getProperty(elem, 'rightMargin', margin)
        bottom = self.wprops.getProperty(elem, 'bottomMargin', margin)

        # A layout widget should, by default, have no margins.
        if self.stack.topIsLayoutWidget():
            if left < 0: left = 0
            if top < 0: top = 0
            if right < 0: right = 0
            if bottom < 0: bottom = 0

        if left >= 0 or top >= 0 or right >= 0 or bottom >= 0:
            # We inject the new internal property.
            cme = SubElement(elem, 'property', name='pyuicMargins')
            SubElement(cme, 'number').text = str(left)
            SubElement(cme, 'number').text = str(top)
            SubElement(cme, 'number').text = str(right)
            SubElement(cme, 'number').text = str(bottom)

        # We use an internal property to handle spacing which will use separate
        # horizontal and vertical spacing if they are found to be different.
        # The following will select, in order of preference, separate
        # horizontal and vertical spacing, the same spacing in both directions,
        # and the default spacing.
        spacing = self.wprops.getProperty(elem, 'spacing',
                self.defaults['spacing'])
        horiz = self.wprops.getProperty(elem, 'horizontalSpacing', spacing)
        vert = self.wprops.getProperty(elem, 'verticalSpacing', spacing)

        if horiz >= 0 or vert >= 0:
            # We inject the new internal property.
            cme = SubElement(elem, 'property', name='pyuicSpacing')
            SubElement(cme, 'number').text = str(horiz)
            SubElement(cme, 'number').text = str(vert)

        classname = elem.attrib["class"]
        if self.stack.topIsLayout():
            parent = None
        else:
            parent = self.stack.topwidget
        if "name" not in elem.attrib:
            elem.attrib["name"] = classname[1:].lower()
        self.stack.push(self.setupObject(classname, parent, elem))
        self.traverseWidgetTree(elem)

        layout = self.stack.popLayout()
        self.configureLayout(elem, layout)

        if self.stack.topIsLayout():
            top_layout = self.stack.peek()
            lp = elem.attrib['layout-position']

            if isinstance(top_layout, QtWidgets.QFormLayout):
                top_layout.setLayout(lp[0], self._form_layout_role(lp), layout)
            else:
                top_layout.addLayout(layout, *lp)

    def configureLayout(self, elem, layout):
        if isinstance(layout, QtWidgets.QGridLayout):
            self.setArray(elem, 'columnminimumwidth',
                    layout.setColumnMinimumWidth)
            self.setArray(elem, 'rowminimumheight',
                    layout.setRowMinimumHeight)
            self.setArray(elem, 'columnstretch', layout.setColumnStretch)
            self.setArray(elem, 'rowstretch', layout.setRowStretch)

        elif isinstance(layout, QtWidgets.QBoxLayout):
            self.setArray(elem, 'stretch', layout.setStretch)

    def setArray(self, elem, name, setter):
        array = elem.attrib.get(name)
        if array:
            for idx, value in enumerate(array.split(',')):
                value = int(value)
                if value > 0:
                    setter(idx, value)

    def disableSorting(self, w):
        if self.item_nr == 0:
            self.sorting_enabled = self.factory.invoke("__sortingEnabled",
                    w.isSortingEnabled)
            w.setSortingEnabled(False)

    def handleItem(self, elem):
        if self.stack.topIsLayout():
            elem[0].attrib['layout-position'] = _layout_position(elem)
            self.traverseWidgetTree(elem)
        else:
            w = self.stack.topwidget

            if isinstance(w, QtWidgets.QComboBox):
                text = self.wprops.getProperty(elem, "text")
                icon = self.wprops.getProperty(elem, "icon")

                if icon:
                    w.addItem(icon, '')
                else:
                    w.addItem('')

                w.setItemText(self.item_nr, text)

            elif isinstance(w, QtWidgets.QListWidget):
                self.disableSorting(w)
                item = self.createWidgetItem('QListWidgetItem', elem, w.item,
                        self.item_nr)
                w.addItem(item)

            elif isinstance(w, QtWidgets.QTreeWidget):
                if self.itemstack:
                    parent, _ = self.itemstack[-1]
                    _, nr_in_root = self.itemstack[0]
                else:
                    parent = w
                    nr_in_root = self.item_nr

                item = self.factory.createQObject("QTreeWidgetItem",
                        "item_%d" % len(self.itemstack), (parent, ), False)

                if self.item_nr == 0 and not self.itemstack:
                    self.sorting_enabled = self.factory.invoke("__sortingEnabled", w.isSortingEnabled)
                    w.setSortingEnabled(False)

                self.itemstack.append((item, self.item_nr))
                self.item_nr = 0

                # We have to access the item via the tree when setting the
                # text.
                titm = w.topLevelItem(nr_in_root)
                for child, nr_in_parent in self.itemstack[1:]:
                    titm = titm.child(nr_in_parent)

                column = -1
                for prop in elem.findall('property'):
                    c_prop = self.wprops.convert(prop)
                    c_prop_name = prop.attrib['name']

                    if c_prop_name == 'text':
                        column += 1
                        if c_prop:
                            titm.setText(column, c_prop)
                    elif c_prop_name == 'statusTip':
                        item.setStatusTip(column, c_prop)
                    elif c_prop_name == 'toolTip':
                        item.setToolTip(column, c_prop)
                    elif c_prop_name == 'whatsThis':
                        item.setWhatsThis(column, c_prop)
                    elif c_prop_name == 'font':
                        item.setFont(column, c_prop)
                    elif c_prop_name == 'icon':
                        item.setIcon(column, c_prop)
                    elif c_prop_name == 'background':
                        item.setBackground(column, c_prop)
                    elif c_prop_name == 'foreground':
                        item.setForeground(column, c_prop)
                    elif c_prop_name == 'flags':
                        item.setFlags(c_prop)
                    elif c_prop_name == 'checkState':
                        item.setCheckState(column, c_prop)

                self.traverseWidgetTree(elem)
                _, self.item_nr = self.itemstack.pop()

            elif isinstance(w, QtWidgets.QTableWidget):
                row = int(elem.attrib['row'])
                col = int(elem.attrib['column'])

                self.disableSorting(w)
                item = self.createWidgetItem('QTableWidgetItem', elem, w.item,
                        row, col)
                w.setItem(row, col, item)

            self.item_nr += 1

    def addAction(self, elem):
        self.actions.append((self.stack.topwidget, elem.attrib["name"]))

    @staticmethod
    def any_i18n(*args):
        """ Return True if any argument appears to be an i18n string. """

        for a in args:
            if a is not None and not isinstance(a, str):
                return True

        return False

    def createWidgetItem(self, item_type, elem, getter, *getter_args):
        """ Create a specific type of widget item. """

        item = self.factory.createQObject(item_type, "item", (), False)
        props = self.wprops

        # Note that not all types of widget items support the full set of
        # properties.

        text = props.getProperty(elem, 'text')
        status_tip = props.getProperty(elem, 'statusTip')
        tool_tip = props.getProperty(elem, 'toolTip')
        whats_this = props.getProperty(elem, 'whatsThis')

        if self.any_i18n(text, status_tip, tool_tip, whats_this):
            self.factory.invoke("item", getter, getter_args)

        if text:
            item.setText(text)

        if status_tip:
            item.setStatusTip(status_tip)

        if tool_tip:
            item.setToolTip(tool_tip)

        if whats_this:
            item.setWhatsThis(whats_this)

        text_alignment = props.getProperty(elem, 'textAlignment')
        if text_alignment:
            item.setTextAlignment(text_alignment)

        font = props.getProperty(elem, 'font')
        if font:
            item.setFont(font)

        icon = props.getProperty(elem, 'icon')
        if icon:
            item.setIcon(icon)

        background = props.getProperty(elem, 'background')
        if background:
            item.setBackground(background)

        foreground = props.getProperty(elem, 'foreground')
        if foreground:
            item.setForeground(foreground)

        flags = props.getProperty(elem, 'flags')
        if flags:
            item.setFlags(flags)

        check_state = props.getProperty(elem, 'checkState')
        if check_state:
            item.setCheckState(check_state)

        return item

    def addHeader(self, elem):
        w = self.stack.topwidget

        if isinstance(w, QtWidgets.QTreeWidget):
            props = self.wprops
            col = self.column_counter

            text = props.getProperty(elem, 'text')
            if text:
                w.headerItem().setText(col, text)

            status_tip = props.getProperty(elem, 'statusTip')
            if status_tip:
                w.headerItem().setStatusTip(col, status_tip)

            tool_tip = props.getProperty(elem, 'toolTip')
            if tool_tip:
                w.headerItem().setToolTip(col, tool_tip)

            whats_this = props.getProperty(elem, 'whatsThis')
            if whats_this:
                w.headerItem().setWhatsThis(col, whats_this)

            text_alignment = props.getProperty(elem, 'textAlignment')
            if text_alignment:
                w.headerItem().setTextAlignment(col, text_alignment)

            font = props.getProperty(elem, 'font')
            if font:
                w.headerItem().setFont(col, font)

            icon = props.getProperty(elem, 'icon')
            if icon:
                w.headerItem().setIcon(col, icon)

            background = props.getProperty(elem, 'background')
            if background:
                w.headerItem().setBackground(col, background)

            foreground = props.getProperty(elem, 'foreground')
            if foreground:
                w.headerItem().setForeground(col, foreground)

            self.column_counter += 1

        elif isinstance(w, QtWidgets.QTableWidget):
            if len(elem) != 0:
                if elem.tag == 'column':
                    item = self.createWidgetItem('QTableWidgetItem', elem,
                            w.horizontalHeaderItem, self.column_counter)
                    w.setHorizontalHeaderItem(self.column_counter, item)
                    self.column_counter += 1
                elif elem.tag == 'row':
                    item = self.createWidgetItem('QTableWidgetItem', elem,
                            w.verticalHeaderItem, self.row_counter)
                    w.setVerticalHeaderItem(self.row_counter, item)
                    self.row_counter += 1

    def setZOrder(self, elem):
        # Designer can generate empty zorder elements.
        if elem.text is None:
            return

        # Designer allows the z-order of spacer items to be specified even
        # though they can't be raised, so ignore any missing raise_() method.
        try:
            getattr(self.toplevelWidget, elem.text).raise_()
        except AttributeError:
            # Note that uic issues a warning message.
            pass

    def createAction(self, elem):
        self.setupObject("QAction", self.currentActionGroup or self.toplevelWidget,
                         elem)

    def createActionGroup(self, elem):
        action_group = self.setupObject("QActionGroup", self.toplevelWidget, elem)
        self.currentActionGroup = action_group
        self.traverseWidgetTree(elem)
        self.currentActionGroup = None

    widgetTreeItemHandlers = {
        "widget"    : createWidget,
        "addaction" : addAction,
        "layout"    : createLayout,
        "spacer"    : createSpacer,
        "item"      : handleItem,
        "action"    : createAction,
        "actiongroup": createActionGroup,
        "column"    : addHeader,
        "row"       : addHeader,
        "zorder"    : setZOrder,
        }

    def traverseWidgetTree(self, elem):
        for child in iter(elem):
            try:
                handler = self.widgetTreeItemHandlers[child.tag]
            except KeyError: 
                continue

            handler(self, child)

    def createUserInterface(self, elem):
        # Get the names of the class and widget.
        cname = elem.attrib["class"]
        wname = elem.attrib["name"]

        # If there was no widget name then derive it from the class name.
        if not wname:
            wname = cname

            if wname.startswith("Q"):
                wname = wname[1:]

            wname = wname[0].lower() + wname[1:]

        self.toplevelWidget = self.createToplevelWidget(cname, wname)
        self.toplevelWidget.setObjectName(wname)
        DEBUG("toplevel widget is %s",
              self.toplevelWidget.metaObject().className())
        self.wprops.setProperties(self.toplevelWidget, elem)
        self.stack.push(self.toplevelWidget)
        self.traverseWidgetTree(elem)
        self.stack.popWidget()
        self.addActions()
        self.setBuddies()
        self.setDelayedProps()
        
    def addActions(self):
        for widget, action_name in self.actions:
            if action_name == "separator":
                widget.addSeparator()
            else:
                DEBUG("add action %s to %s", action_name, widget.objectName())
                action_obj = getattr(self.toplevelWidget, action_name)
                if isinstance(action_obj, QtWidgets.QMenu):
                    widget.addAction(action_obj.menuAction())
                elif not isinstance(action_obj, QtWidgets.QActionGroup):
                    widget.addAction(action_obj)

    def setDelayedProps(self):
        for widget, layout, setter, args in self.wprops.delayed_props:
            if layout:
                widget = widget.layout()

            setter = getattr(widget, setter)
            setter(args)
            
    def setBuddies(self):
        for widget, buddy in self.wprops.buddies:
            DEBUG("%s is buddy of %s", buddy, widget.objectName())
            try:
                widget.setBuddy(getattr(self.toplevelWidget, buddy))
            except AttributeError:
                DEBUG("ERROR in ui spec: %s (buddy of %s) does not exist",
                      buddy, widget.objectName())

    def classname(self, elem):
        DEBUG("uiname is %s", elem.text)
        name = elem.text

        if name is None:
            name = ""

        self.uiname = name
        self.wprops.uiname = name
        self.setContext(name)

    def setContext(self, context):
        """
        Reimplemented by a sub-class if it needs to know the translation
        context.
        """
        pass

    def readDefaults(self, elem):
        self.defaults['margin'] = int(elem.attrib['margin'])
        self.defaults['spacing'] = int(elem.attrib['spacing'])

    def setTaborder(self, elem):
        lastwidget = None
        for widget_elem in elem:
            widget = getattr(self.toplevelWidget, widget_elem.text)

            if lastwidget is not None:
                self.toplevelWidget.setTabOrder(lastwidget, widget)

            lastwidget = widget

    def readResources(self, elem):
        """
        Read a "resources" tag and add the module to import to the parser's
        list of them.
        """
        try:
            iterator = getattr(elem, 'iter')
        except AttributeError:
            iterator = getattr(elem, 'getiterator')

        for include in iterator("include"):
            loc = include.attrib.get("location")

            # Apply the convention for naming the Python files generated by
            # pyrcc5.
            if loc and loc.endswith('.qrc'):
                mname = os.path.basename(loc[:-4] + self._resource_suffix)
                if mname not in self.resources:
                    self.resources.append(mname)

    def createConnections(self, elem):
        def name2object(obj):
            if obj == self.uiname:
                return self.toplevelWidget
            else:
                return getattr(self.toplevelWidget, obj)

        for conn in iter(elem):
            signal = conn.findtext('signal')
            signal_name, signal_args = signal.split('(')
            signal_args = signal_args[:-1].replace(' ', '')
            sender = name2object(conn.findtext('sender'))
            bound_signal = getattr(sender, signal_name)

            slot = self.factory.getSlot(name2object(conn.findtext('receiver')),
                    conn.findtext('slot').split('(')[0])

            if signal_args == '':
                bound_signal.connect(slot)
            else:
                signal_args = signal_args.split(',')

                if len(signal_args) == 1:
                    bound_signal[signal_args[0]].connect(slot)
                else:
                    bound_signal[tuple(signal_args)].connect(slot)

        QtCore.QMetaObject.connectSlotsByName(self.toplevelWidget)

    def customWidgets(self, elem):
        def header2module(header):
            """header2module(header) -> string

            Convert paths to C++ header files to according Python modules
            >>> header2module("foo/bar/baz.h")
            'foo.bar.baz'
            """
            if header.endswith(".h"):
                header = header[:-2]

            mpath = []
            for part in header.split('/'):
                # Ignore any empty parts or those that refer to the current
                # directory.
                if part not in ('', '.'):
                    if part == '..':
                        # We should allow this for Python3.
                        raise SyntaxError("custom widget header file name may not contain '..'.")

                    mpath.append(part)

            return '.'.join(mpath)
    
        for custom_widget in iter(elem):
            classname = custom_widget.findtext("class")
            self.factory.addCustomWidget(classname,
                                     custom_widget.findtext("extends") or "QWidget",
                                     header2module(custom_widget.findtext("header")))

    def createToplevelWidget(self, classname, widgetname):
        raise NotImplementedError

    def buttonGroups(self, elem):
        for button_group in iter(elem):
            if button_group.tag == 'buttongroup':
                bg_name = button_group.attrib['name']
                bg = ButtonGroup()
                self.button_groups[bg_name] = bg

                prop = self.getProperty(button_group, 'exclusive')
                if prop is not None:
                    if prop.findtext('bool') == 'false':
                        bg.exclusive = False

    # finalize will be called after the whole tree has been parsed and can be
    # overridden.
    def finalize(self):
        pass

    def parse(self, filename, resource_suffix):
        if hasattr(filename, 'read'):
            base_dir = ''
        else:
            # Allow the filename to be a QString.
            filename = str(filename)
            base_dir = os.path.dirname(filename)

        self.wprops.set_base_dir(base_dir)

        self._resource_suffix = resource_suffix

        # The order in which the different branches are handled is important.
        # The widget tree handler relies on all custom widgets being known, and
        # in order to create the connections, all widgets have to be populated.
        branchHandlers = (
            ("layoutdefault", self.readDefaults),
            ("class",         self.classname),
            ("buttongroups",  self.buttonGroups),
            ("customwidgets", self.customWidgets),
            ("widget",        self.createUserInterface),
            ("connections",   self.createConnections),
            ("tabstops",      self.setTaborder),
            ("resources",     self.readResources),
        )

        document = parse(filename)
        root = document.getroot()

        if root.tag != 'ui':
            raise SyntaxError("not created by Qt Designer")

        version = root.attrib.get('version')
        if version is None:
            raise SyntaxError("missing version number")

        # Right now, only version 4.0 is supported.
        if version != '4.0':
            raise SyntaxError("only Qt Designer files v4.0 are supported")

        for tagname, actor in branchHandlers:
            elem = document.find(tagname)
            if elem is not None:
                actor(elem)
        self.finalize()
        w = self.toplevelWidget
        self.reset()
        return w

    @staticmethod
    def _form_layout_role(layout_position):
        if layout_position[3] > 1:
            role = QtWidgets.QFormLayout.SpanningRole
        elif layout_position[1] == 1:
            role = QtWidgets.QFormLayout.FieldRole
        else:
            role = QtWidgets.QFormLayout.LabelRole

        return role
