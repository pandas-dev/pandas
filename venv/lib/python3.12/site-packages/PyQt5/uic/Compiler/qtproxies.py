#############################################################################
##
## Copyright (C) 2021 Riverbank Computing Limited.
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
import re

from .indenter import write_code
from .misc import Literal, moduleMember

if sys.hexversion >= 0x03000000:
    from ..port_v3.proxy_base import ProxyBase
    from ..port_v3.as_string import as_string
else:
    from ..port_v2.proxy_base import ProxyBase
    from ..port_v2.as_string import as_string


i18n_strings = []
i18n_context = ""

def i18n_print(string):
    i18n_strings.append(string)

def i18n_void_func(name):
    def _printer(self, *args):
        i18n_print("%s.%s(%s)" % (self, name, ", ".join(map(as_string, args))))
    return _printer

def i18n_func(name):
    def _printer(self, rname, *args):
        i18n_print("%s = %s.%s(%s)" % (rname, self, name, ", ".join(map(as_string, args))))
        return Literal(rname)

    return _printer

def strict_getattr(module, clsname):
    cls = getattr(module, clsname)
    if issubclass(cls, LiteralProxyClass):
        raise AttributeError(cls)
    else:
        return cls
    

class i18n_string(object):
    def __init__(self, string, disambig):
        self.string = string
        self.disambig = disambig

    def __str__(self):
        if self.disambig is None:
            return '_translate("%s", %s)' % (i18n_context, as_string(self.string))

        return '_translate("%s", %s, %s)' % (i18n_context, as_string(self.string), as_string(self.disambig))


# Classes with this flag will be handled as literal values. If functions are
# called on these classes, the literal value changes.
# Example:
# the code
# >>> QSize(9,10).expandedTo(...)
# will print just that code.
AS_ARGUMENT = 0x02

# Classes with this flag may have members that are signals which themselves
# will have a connect() member.
AS_SIGNAL = 0x01

# ATTENTION: currently, classes can either be literal or normal. If a class
# should need both kinds of behaviour, the code has to be changed.

class ProxyClassMember(object):
    def __init__(self, proxy, function_name, flags):
        self.proxy = proxy
        self.function_name = function_name
        self.flags = flags
        
    def __str__(self):
        return "%s.%s" % (self.proxy, self.function_name)
    
    def __call__(self, *args):
        if self.function_name == 'setProperty':
            str_args = (as_string(args[0]), as_string(args[1]))
        else:
            str_args = map(as_string, args)

        func_call = "%s.%s(%s)" % (self.proxy,
                                   self.function_name,
                                   ", ".join(str_args))
        if self.flags & AS_ARGUMENT:
            self.proxy._uic_name = func_call
            return self.proxy
        else:
            needs_translation = False
            for arg in args:
                if isinstance(arg, i18n_string):
                    needs_translation = True
            if needs_translation:
                i18n_print(func_call)
            else:
                if self.function_name == 'connect':
                    func_call += ' # type: ignore'

                write_code(func_call)                       

    def __getattribute__(self, attribute):
        """ Reimplemented to create a proxy connect() if requested and this
        might be a proxy for a signal.
        """

        try:
            return object.__getattribute__(self, attribute)
        except AttributeError:
            if attribute == 'connect' and self.flags & AS_SIGNAL:
                return ProxyClassMember(self, attribute, 0)

            raise

    def __getitem__(self, idx):
        """ Reimplemented to create a proxy member that should be a signal that
        passes arguments.  We handle signals without arguments before we get
        here and never apply the index notation to them.
        """

        return ProxySignalWithArguments(self.proxy, self.function_name, idx)


class ProxySignalWithArguments(object):
    """ This is a proxy for (what should be) a signal that passes arguments.
    """

    def __init__(self, sender, signal_name, signal_index):
        self._sender = sender
        self._signal_name = signal_name

        # Convert the signal index, which will be a single argument or a tuple
        # of arguments, to quoted strings.
        if isinstance(signal_index, tuple):
            self._signal_index = ','.join(["'%s'" % a for a in signal_index])
        else:
            self._signal_index = "'%s'" % signal_index

    def connect(self, slot):
        write_code("%s.%s[%s].connect(%s) # type: ignore" % (self._sender, self._signal_name, self._signal_index, slot))


class ProxyClass(ProxyBase):
    flags = 0

    def __init__(self, objectname, is_attribute, args=(), noInstantiation=False):
        if objectname:
            if is_attribute:
                objectname = "self." + objectname

            self._uic_name = objectname
        else:
            self._uic_name = "Unnamed"

        if not noInstantiation:
            funcall = "%s(%s)" % \
                    (moduleMember(self.module, self.__class__.__name__),
                    ", ".join(map(str, args)))

            if objectname:
                funcall = "%s = %s" % (objectname, funcall)

            write_code(funcall)
    
    def __str__(self):
        return self._uic_name

    def __getattribute__(self, attribute):
        try:
            return object.__getattribute__(self, attribute)
        except AttributeError:
            return ProxyClassMember(self, attribute, self.flags)


class LiteralProxyClass(ProxyClass):
    """LiteralObject(*args) -> new literal class

    a literal class can be used as argument in a function call

    >>> class Foo(LiteralProxyClass): pass
    >>> str(Foo(1,2,3)) == "Foo(1,2,3)"
    """
    flags = AS_ARGUMENT

    def __init__(self, *args):
        self._uic_name = "%s(%s)" % \
                     (moduleMember(self.module, self.__class__.__name__),
                      ", ".join(map(as_string, args)))
        

class ProxyNamespace(ProxyBase):
    pass


# These are all the Qt classes used by pyuic5 in their namespaces. If a class
# is missing, the compiler will fail, normally with an AttributeError.
#
# For adding new classes:
#     - utility classes used as literal values do not need to be listed
#       because they are created on the fly as subclasses of LiteralProxyClass
#     - classes which are *not* QWidgets inherit from ProxyClass and they
#       have to be listed explicitly in the correct namespace. These classes
#       are created via a ProxyQObjectCreator
#     - new QWidget-derived classes have to inherit from qtproxies.QWidget
#       If the widget does not need any special methods, it can be listed
#       in _qwidgets

class QtCore(ProxyNamespace):
    class Qt(ProxyNamespace):
        pass

    ## connectSlotsByName and connect have to be handled as class methods,
    ## otherwise they would be created as LiteralProxyClasses and never be
    ## printed
    class QMetaObject(ProxyClass):
        @classmethod
        def connectSlotsByName(cls, *args):
            ProxyClassMember(cls, "connectSlotsByName", 0)(*args)

    class QObject(ProxyClass):
        flags = AS_SIGNAL

        def metaObject(self):
            class _FakeMetaObject(object):
                def className(*args):
                    return self.__class__.__name__
            return _FakeMetaObject()

        def objectName(self):
            return self._uic_name.split(".")[-1]


class QtGui(ProxyNamespace):
    class QIcon(ProxyClass):
        class fromTheme(ProxyClass): pass

    class QConicalGradient(ProxyClass): pass
    class QLinearGradient(ProxyClass): pass
    class QRadialGradient(ProxyClass): pass
    class QBrush(ProxyClass): pass
    class QPainter(ProxyClass): pass
    class QPalette(ProxyClass): pass
    class QFont(ProxyClass): pass
    class QFontDatabase(ProxyClass): pass


# These sub-class QWidget but aren't themselves sub-classed.
_qwidgets = ("QCalendarWidget", "QDialogButtonBox", "QDockWidget", "QGroupBox",
        "QLineEdit", "QMainWindow", "QMenuBar", "QOpenGLWidget",
        "QProgressBar", "QStatusBar", "QToolBar", "QWizardPage")

class QtWidgets(ProxyNamespace):
    class QApplication(QtCore.QObject):
        @staticmethod
        def translate(uiname, text, disambig):
            return i18n_string(text or "", disambig)

    class QSpacerItem(ProxyClass): pass
    class QSizePolicy(ProxyClass): pass
    # QActions inherit from QObject for the meta-object stuff and the hierarchy
    # has to be correct since we have a isinstance(x, QtWidgets.QLayout) call
    # in the UI parser.
    class QAction(QtCore.QObject): pass
    class QActionGroup(QtCore.QObject): pass
    class QButtonGroup(QtCore.QObject): pass
    class QLayout(QtCore.QObject): pass
    class QGridLayout(QLayout): pass
    class QBoxLayout(QLayout): pass
    class QHBoxLayout(QBoxLayout): pass
    class QVBoxLayout(QBoxLayout): pass
    class QFormLayout(QLayout): pass
    
    class QWidget(QtCore.QObject):
        def font(self):
            return Literal("%s.font()" % self)

        def minimumSizeHint(self):
            return Literal("%s.minimumSizeHint()" % self)

        def sizePolicy(self):
            sp = LiteralProxyClass()
            sp._uic_name = "%s.sizePolicy()" % self
            return sp

    class QDialog(QWidget): pass
    class QColorDialog(QDialog): pass
    class QFileDialog(QDialog): pass
    class QFontDialog(QDialog): pass
    class QInputDialog(QDialog): pass
    class QMessageBox(QDialog): pass
    class QWizard(QDialog): pass

    class QAbstractSlider(QWidget): pass
    class QDial(QAbstractSlider): pass
    class QScrollBar(QAbstractSlider): pass
    class QSlider(QAbstractSlider): pass

    class QMenu(QWidget):
        def menuAction(self):
            return Literal("%s.menuAction()" % self)

    class QTabWidget(QWidget):
        def addTab(self, *args):
            text = args[-1]

            if isinstance(text, i18n_string):
                i18n_print("%s.setTabText(%s.indexOf(%s), %s)" % \
                        (self._uic_name, self._uic_name, args[0], text))
                args = args[:-1] + ("", )

            ProxyClassMember(self, "addTab", 0)(*args)

        def indexOf(self, page):
            return Literal("%s.indexOf(%s)" % (self, page))

    class QComboBox(QWidget): pass
    class QFontComboBox(QComboBox): pass

    class QAbstractSpinBox(QWidget): pass
    class QDoubleSpinBox(QAbstractSpinBox): pass
    class QSpinBox(QAbstractSpinBox): pass

    class QDateTimeEdit(QAbstractSpinBox): pass
    class QDateEdit(QDateTimeEdit): pass
    class QTimeEdit(QDateTimeEdit): pass

    class QFrame(QWidget): pass
    class QLabel(QFrame): pass
    class QLCDNumber(QFrame): pass
    class QSplitter(QFrame): pass
    class QStackedWidget(QFrame): pass

    class QToolBox(QFrame):
        def addItem(self, *args):
            text = args[-1]

            if isinstance(text, i18n_string):
                i18n_print("%s.setItemText(%s.indexOf(%s), %s)" % \
                        (self._uic_name, self._uic_name, args[0], text))
                args = args[:-1] + ("", )

            ProxyClassMember(self, "addItem", 0)(*args)

        def indexOf(self, page):
            return Literal("%s.indexOf(%s)" % (self, page))

        def layout(self):
            return QtWidgets.QLayout("%s.layout()" % self,
                    False, (), noInstantiation=True)

    class QAbstractScrollArea(QFrame):
        def viewport(self):
            return QtWidgets.QWidget("%s.viewport()" % self, False, (),
                    noInstantiation=True)

    class QGraphicsView(QAbstractScrollArea): pass
    class QMdiArea(QAbstractScrollArea): pass
    class QPlainTextEdit(QAbstractScrollArea): pass
    class QScrollArea(QAbstractScrollArea): pass

    class QTextEdit(QAbstractScrollArea): pass
    class QTextBrowser(QTextEdit): pass

    class QAbstractItemView(QAbstractScrollArea): pass
    class QColumnView(QAbstractItemView): pass
    class QHeaderView(QAbstractItemView): pass
    class QListView(QAbstractItemView): pass

    class QTableView(QAbstractItemView):
        def horizontalHeader(self):
            return QtWidgets.QHeaderView("%s.horizontalHeader()" % self,
                    False, (), noInstantiation=True)

        def verticalHeader(self):
            return QtWidgets.QHeaderView("%s.verticalHeader()" % self,
                    False, (), noInstantiation=True)

    class QTreeView(QAbstractItemView):
        def header(self):
            return QtWidgets.QHeaderView("%s.header()" % self,
                    False, (), noInstantiation=True)

    class QUndoView(QListView): pass

    class QListWidgetItem(ProxyClass): pass

    class QListWidget(QListView):
        setSortingEnabled = i18n_void_func("setSortingEnabled")
        isSortingEnabled = i18n_func("isSortingEnabled")
        item = i18n_func("item")

    class QTableWidgetItem(ProxyClass): pass

    class QTableWidget(QTableView):
        setSortingEnabled = i18n_void_func("setSortingEnabled")
        isSortingEnabled = i18n_func("isSortingEnabled")
        item = i18n_func("item")
        horizontalHeaderItem = i18n_func("horizontalHeaderItem")
        verticalHeaderItem = i18n_func("verticalHeaderItem")

    class QTreeWidgetItem(ProxyClass):
        def child(self, index):
            return QtWidgets.QTreeWidgetItem("%s.child(%i)" % (self, index),
                    False, (), noInstantiation=True)

    class QTreeWidget(QTreeView):
        setSortingEnabled = i18n_void_func("setSortingEnabled")
        isSortingEnabled = i18n_func("isSortingEnabled")

        def headerItem(self):
            return QtWidgets.QWidget("%s.headerItem()" % self, False, (),
                    noInstantiation=True)

        def topLevelItem(self, index):
            return QtWidgets.QTreeWidgetItem("%s.topLevelItem(%i)" % (self, index),
                    False, (), noInstantiation=True)

    class QAbstractButton(QWidget): pass
    class QCheckBox(QAbstractButton): pass
    class QRadioButton(QAbstractButton): pass
    class QToolButton(QAbstractButton): pass

    class QPushButton(QAbstractButton): pass
    class QCommandLinkButton(QPushButton): pass
    class QKeySequenceEdit(QWidget): pass

    # Add all remaining classes.
    for _class in _qwidgets:
        if _class not in locals():
            locals()[_class] = type(_class, (QWidget, ), {})
