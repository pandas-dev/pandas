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


import logging
import os.path
import sys

from .exceptions import NoSuchClassError, UnsupportedPropertyError
from .icon_cache import IconCache

if sys.hexversion >= 0x03000000:
    from .port_v3.ascii_upper import ascii_upper
else:
    from .port_v2.ascii_upper import ascii_upper


logger = logging.getLogger(__name__)
DEBUG = logger.debug


QtCore = None
QtGui = None
QtWidgets = None


def int_list(prop):
    return [int(child.text) for child in prop]

def float_list(prop):
    return [float(child.text) for child in prop]

bool_ = lambda v: v == "true"

def qfont_enum(v):
    return getattr(QtGui.QFont, v)

def needsWidget(func):
    func.needsWidget = True
    return func


class Properties(object):
    def __init__(self, factory, qtcore_module, qtgui_module, qtwidgets_module):
        self.factory = factory

        global QtCore, QtGui, QtWidgets
        QtCore = qtcore_module
        QtGui = qtgui_module
        QtWidgets = qtwidgets_module

        self._base_dir = ''

        self.reset()

    def set_base_dir(self, base_dir):
        """ Set the base directory to be used for all relative filenames. """

        self._base_dir = base_dir
        self.icon_cache.set_base_dir(base_dir)

    def reset(self):
        self.buddies = []
        self.delayed_props = []
        self.icon_cache = IconCache(self.factory, QtGui)

    def _pyEnumMember(self, cpp_name):
        try:
            prefix, membername = cpp_name.split("::")
        except ValueError:
            prefix = 'Qt'
            membername = cpp_name

        if prefix == 'Qt':
            return getattr(QtCore.Qt, membername)

        scope = self.factory.findQObjectType(prefix)
        if scope is None:
            raise NoSuchClassError(prefix)

        return getattr(scope, membername)

    def _set(self, prop):
        expr = [self._pyEnumMember(v) for v in prop.text.split('|')]

        value = expr[0]
        for v in expr[1:]:
            value |= v

        return value

    def _enum(self, prop):
        return self._pyEnumMember(prop.text)

    def _number(self, prop):
        return int(prop.text)

    _UInt = _uInt = _longLong = _uLongLong = _number

    def _double(self, prop):
        return float(prop.text)

    def _bool(self, prop):
        return prop.text == 'true'

    def _stringlist(self, prop):
        return [self._string(p, notr='true') for p in prop]

    def _string(self, prop, notr=None):
        text = prop.text

        if text is None:
            return ""

        if prop.get('notr', notr) == 'true':
            return text

        disambig = prop.get('comment')

        return QtWidgets.QApplication.translate(self.uiname, text, disambig)

    _char = _string

    def _cstring(self, prop):
        return str(prop.text)

    def _color(self, prop):
        args = int_list(prop)

        # Handle the optional alpha component.
        alpha = int(prop.get("alpha", "255"))

        if alpha != 255:
            args.append(alpha)

        return QtGui.QColor(*args)

    def _point(self, prop):
        return QtCore.QPoint(*int_list(prop))

    def _pointf(self, prop):
        return QtCore.QPointF(*float_list(prop))

    def _rect(self, prop):
        return QtCore.QRect(*int_list(prop))

    def _rectf(self, prop):
        return QtCore.QRectF(*float_list(prop))

    def _size(self, prop):
        return QtCore.QSize(*int_list(prop))

    def _sizef(self, prop):
        return QtCore.QSizeF(*float_list(prop))

    def _pixmap(self, prop):
        if prop.text:
            fname = prop.text.replace("\\", "\\\\")
            if self._base_dir != '' and fname[0] != ':' and not os.path.isabs(fname):
                fname = os.path.join(self._base_dir, fname)

            return QtGui.QPixmap(fname)

        # Don't bother to set the property if the pixmap is empty.
        return None

    def _iconset(self, prop):
        return self.icon_cache.get_icon(prop)

    def _url(self, prop):
        return QtCore.QUrl(prop[0].text)

    def _locale(self, prop):
        lang = getattr(QtCore.QLocale, prop.attrib['language'])
        country = getattr(QtCore.QLocale, prop.attrib['country'])
        return QtCore.QLocale(lang, country)

    def _date(self, prop):
        return QtCore.QDate(*int_list(prop))

    def _datetime(self, prop):
        args = int_list(prop)
        return QtCore.QDateTime(QtCore.QDate(*args[-3:]), QtCore.QTime(*args[:-3]))

    def _time(self, prop):
        return QtCore.QTime(*int_list(prop))

    def _gradient(self, prop):
        name = 'gradient'

        # Create the specific gradient.
        gtype = prop.get('type', '')

        if gtype == 'LinearGradient':
            startx = float(prop.get('startx'))
            starty = float(prop.get('starty'))
            endx = float(prop.get('endx'))
            endy = float(prop.get('endy'))
            gradient = self.factory.createQObject('QLinearGradient', name,
                    (startx, starty, endx, endy), is_attribute=False)

        elif gtype == 'ConicalGradient':
            centralx = float(prop.get('centralx'))
            centraly = float(prop.get('centraly'))
            angle = float(prop.get('angle'))
            gradient = self.factory.createQObject('QConicalGradient', name,
                    (centralx, centraly, angle), is_attribute=False)

        elif gtype == 'RadialGradient':
            centralx = float(prop.get('centralx'))
            centraly = float(prop.get('centraly'))
            radius = float(prop.get('radius'))
            focalx = float(prop.get('focalx'))
            focaly = float(prop.get('focaly'))
            gradient = self.factory.createQObject('QRadialGradient', name,
                    (centralx, centraly, radius, focalx, focaly),
                    is_attribute=False)

        else:
            raise UnsupportedPropertyError(prop.tag)

        # Set the common values.
        spread = prop.get('spread')
        if spread:
            gradient.setSpread(getattr(QtGui.QGradient, spread))

        cmode = prop.get('coordinatemode')
        if cmode:
            gradient.setCoordinateMode(getattr(QtGui.QGradient, cmode))

        # Get the gradient stops.
        for gstop in prop:
            if gstop.tag != 'gradientstop':
                raise UnsupportedPropertyError(gstop.tag)

            position = float(gstop.get('position'))
            color = self._color(gstop[0])

            gradient.setColorAt(position, color)

        return gradient

    def _palette(self, prop):
        palette = self.factory.createQObject("QPalette", "palette", (),
                is_attribute=False)

        for palette_elem in prop:
            sub_palette = getattr(QtGui.QPalette, palette_elem.tag.title())
            for role, color in enumerate(palette_elem):
                if color.tag == 'color':
                    # Handle simple colour descriptions where the role is
                    # implied by the colour's position.
                    palette.setColor(sub_palette,
                            QtGui.QPalette.ColorRole(role), self._color(color))
                elif color.tag == 'colorrole':
                    role = getattr(QtGui.QPalette, color.get('role'))
                    brush = self._brush(color[0])
                    palette.setBrush(sub_palette, role, brush)
                else:
                    raise UnsupportedPropertyError(color.tag)

        return palette

    def _brush(self, prop):
        brushstyle = prop.get('brushstyle')

        if brushstyle in ('LinearGradientPattern', 'ConicalGradientPattern', 'RadialGradientPattern'):
            gradient = self._gradient(prop[0])
            brush = self.factory.createQObject("QBrush", "brush", (gradient, ),
                    is_attribute=False)
        else:
            color = self._color(prop[0])
            brush = self.factory.createQObject("QBrush", "brush", (color, ),
                    is_attribute=False)

            brushstyle = getattr(QtCore.Qt, brushstyle)
            brush.setStyle(brushstyle)

        return brush

    #@needsWidget
    def _sizepolicy(self, prop, widget):
        values = [int(child.text) for child in prop]

        if len(values) == 2:
            # Qt v4.3.0 and later.
            horstretch, verstretch = values
            hsizetype = getattr(QtWidgets.QSizePolicy, prop.get('hsizetype'))
            vsizetype = getattr(QtWidgets.QSizePolicy, prop.get('vsizetype'))
        else:
            hsizetype, vsizetype, horstretch, verstretch = values
            hsizetype = QtWidgets.QSizePolicy.Policy(hsizetype)
            vsizetype = QtWidgets.QSizePolicy.Policy(vsizetype)

        sizePolicy = self.factory.createQObject('QSizePolicy', 'sizePolicy',
                (hsizetype, vsizetype), is_attribute=False)
        sizePolicy.setHorizontalStretch(horstretch)
        sizePolicy.setVerticalStretch(verstretch)
        sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
        return sizePolicy
    _sizepolicy = needsWidget(_sizepolicy)

    # font needs special handling/conversion of all child elements.
    _font_attributes = (("Family",          lambda s: s),
                        ("PointSize",       int),
                        ("Bold",            bool_),
                        ("Italic",          bool_),
                        ("Underline",       bool_),
                        ("Weight",          int),
                        ("StrikeOut",       bool_),
                        ("Kerning",         bool_),
                        ("StyleStrategy",   qfont_enum))

    def _font(self, prop):
        newfont = self.factory.createQObject("QFont", "font", (),
                                                     is_attribute = False)
        for attr, converter in self._font_attributes:
            v = prop.findtext("./%s" % (attr.lower(),))
            if v is None:
                continue

            getattr(newfont, "set%s" % (attr,))(converter(v))
        return newfont

    def _cursor(self, prop):
        return QtGui.QCursor(QtCore.Qt.CursorShape(int(prop.text)))

    def _cursorShape(self, prop):
        return QtGui.QCursor(getattr(QtCore.Qt, prop.text))

    def convert(self, prop, widget=None):
        try:
            func = getattr(self, "_" + prop[0].tag)
        except AttributeError:
            raise UnsupportedPropertyError(prop[0].tag)
        else:
            args = {}
            if getattr(func, "needsWidget", False):
                assert widget is not None
                args["widget"] = widget

            return func(prop[0], **args)


    def _getChild(self, elem_tag, elem, name, default=None):
        for prop in elem.findall(elem_tag):
            if prop.attrib["name"] == name:
                return self.convert(prop)
        else:
            return default

    def getProperty(self, elem, name, default=None):
        return self._getChild("property", elem, name, default)

    def getAttribute(self, elem, name, default=None):
        return self._getChild("attribute", elem, name, default)

    def setProperties(self, widget, elem):
        # Lines are sunken unless the frame shadow is explicitly set.
        set_sunken = (elem.attrib.get('class') == 'Line')

        for prop in elem.findall('property'):
            prop_name = prop.attrib['name']
            DEBUG("setting property %s" % (prop_name,))

            if prop_name == 'frameShadow':
                set_sunken = False

            try:
                stdset = bool(int(prop.attrib['stdset']))
            except KeyError:
                stdset = True

            if not stdset:
                self._setViaSetProperty(widget, prop)
            elif hasattr(self, prop_name):
                getattr(self, prop_name)(widget, prop)
            else:
                prop_value = self.convert(prop, widget)
                if prop_value is not None:
                    getattr(widget, 'set%s%s' % (ascii_upper(prop_name[0]), prop_name[1:]))(prop_value)

        if set_sunken:
            widget.setFrameShadow(QtWidgets.QFrame.Sunken)

    # SPECIAL PROPERTIES
    # If a property has a well-known value type but needs special,
    # context-dependent handling, the default behaviour can be overridden here.

    # Delayed properties will be set after the whole widget tree has been
    # populated.
    def _delayed_property(self, widget, prop):
        prop_value = self.convert(prop)
        if prop_value is not None:
            prop_name = prop.attrib["name"]
            self.delayed_props.append((widget, False,
                    'set%s%s' % (ascii_upper(prop_name[0]), prop_name[1:]),
                    prop_value))

    # These properties will be set with a widget.setProperty call rather than
    # calling the set<property> function.
    def _setViaSetProperty(self, widget, prop):
        prop_value = self.convert(prop, widget)
        if prop_value is not None:
            prop_name = prop.attrib['name']

            # This appears to be a Designer/uic hack where stdset=0 means that
            # the viewport should be used.
            if prop[0].tag == 'cursorShape':
                widget.viewport().setProperty(prop_name, prop_value)
            else:
                widget.setProperty(prop_name, prop_value)

    # Ignore the property.
    def _ignore(self, widget, prop):
        pass

    # Define properties that use the canned handlers.
    currentIndex = _delayed_property
    currentRow = _delayed_property

    showDropIndicator = _setViaSetProperty
    intValue = _setViaSetProperty
    value = _setViaSetProperty

    objectName = _ignore
    margin = _ignore
    leftMargin = _ignore
    topMargin = _ignore
    rightMargin = _ignore
    bottomMargin = _ignore
    spacing = _ignore
    horizontalSpacing = _ignore
    verticalSpacing = _ignore

    # tabSpacing is actually the spacing property of the widget's layout.
    def tabSpacing(self, widget, prop):
        prop_value = self.convert(prop)
        if prop_value is not None:
            self.delayed_props.append((widget, True, 'setSpacing', prop_value))

    # buddy setting has to be done after the whole widget tree has been
    # populated.  We can't use delay here because we cannot get the actual
    # buddy yet.
    def buddy(self, widget, prop):
        buddy_name = prop[0].text
        if buddy_name:
            self.buddies.append((widget, buddy_name))

    # geometry is handled specially if set on the toplevel widget.
    def geometry(self, widget, prop):
        if widget.objectName() == self.uiname:
            geom = int_list(prop[0])
            widget.resize(geom[2], geom[3])
        else:
            widget.setGeometry(self._rect(prop[0]))

    def orientation(self, widget, prop):
        # If the class is a QFrame, it's a line.
        if widget.metaObject().className() == 'QFrame':
            widget.setFrameShape(
                {'Qt::Horizontal': QtWidgets.QFrame.HLine,
                 'Qt::Vertical'  : QtWidgets.QFrame.VLine}[prop[0].text])
        else:
            widget.setOrientation(self._enum(prop[0]))

    # The isWrapping attribute of QListView is named inconsistently, it should
    # be wrapping.
    def isWrapping(self, widget, prop):
        widget.setWrapping(self.convert(prop))

    # This is a pseudo-property injected to deal with margins.
    def pyuicMargins(self, widget, prop):
        widget.setContentsMargins(*int_list(prop))

    # This is a pseudo-property injected to deal with spacing.
    def pyuicSpacing(self, widget, prop):
        horiz, vert = int_list(prop)

        if horiz == vert:
            widget.setSpacing(horiz)
        else:
            if horiz >= 0:
                widget.setHorizontalSpacing(horiz)

            if vert >= 0:
                widget.setVerticalSpacing(vert)
