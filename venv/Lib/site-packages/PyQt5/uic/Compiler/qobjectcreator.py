#############################################################################
##
## Copyright (C) 2018 Riverbank Computing Limited.
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
import sys

from .indenter import write_code
from .qtproxies import QtGui, QtWidgets, Literal, strict_getattr

if sys.hexversion >= 0x03000000:
    from ..port_v3.as_string import as_string
else:
    from ..port_v2.as_string import as_string


logger = logging.getLogger(__name__)
DEBUG = logger.debug


class _QtWrapper(object):
    @classmethod
    def search(cls, name):
        try:
            return strict_getattr(cls.module, name)
        except AttributeError:
            return None


class _QtGuiWrapper(_QtWrapper):
    module = QtGui


class _QtWidgetsWrapper(_QtWrapper):
    module = QtWidgets


class _ModuleWrapper(object):
    def __init__(self, name, classes):
        if "." in name:
            idx = name.rfind(".")
            self._package = name[:idx]
            self._module = name[idx + 1:]
        else:
            self._package = None
            self._module = name
            
        self._classes = classes
        self._used = False
    
    def search(self, cls):
        if cls in self._classes:
            self._used = True

            # Remove any C++ scope.
            cls = cls.split('.')[-1]

            return type(cls, (QtWidgets.QWidget,), {"module": self._module})
        else:
            return None

    def _writeImportCode(self):
        if self._used:
            if self._package is None:
                write_code("import %s" % self._module)
            else:
                write_code("from %s import %s" % (self._package, self._module))


class _CustomWidgetLoader(object):
    def __init__(self):
        self._widgets = {}
        self._usedWidgets = set()
        
    def addCustomWidget(self, widgetClass, baseClass, module):
        assert widgetClass not in self._widgets 
        self._widgets[widgetClass] = (baseClass, module)

    def _resolveBaseclass(self, baseClass):
        try:
            for x in range(0, 10):
                try: return strict_getattr(QtWidgets, baseClass)
                except AttributeError: pass
                
                baseClass = self._widgets[baseClass][0]
            else:
                raise ValueError("baseclass resolve took too long, check custom widgets")

        except KeyError:
            raise ValueError("unknown baseclass %s" % baseClass)
        
    def search(self, cls):
        try:
            baseClass = self._resolveBaseclass(self._widgets[cls][0])
            DEBUG("resolved baseclass of %s: %s" % (cls, baseClass))
        except KeyError:
            return None

        self._usedWidgets.add(cls)

        return type(cls, (baseClass, ), {"module" : ""})

    def _writeImportCode(self):
        imports = {}
        for widget in self._usedWidgets:
            _, module = self._widgets[widget]
            imports.setdefault(module, []).append(widget)

        for module, classes in sorted(imports.items()):
            write_code("from %s import %s" % (module, ", ".join(sorted(classes))))


class CompilerCreatorPolicy(object):
    def __init__(self):
        self._modules = []
        
    def createQtGuiWidgetsWrappers(self):
        return [_QtGuiWrapper, _QtWidgetsWrapper]

    def createModuleWrapper(self, name, classes):
        mw = _ModuleWrapper(name, classes)
        self._modules.append(mw)
        return mw

    def createCustomWidgetLoader(self):
        cw = _CustomWidgetLoader()
        self._modules.append(cw)
        return cw

    def instantiate(self, clsObject, objectname, ctor_args, is_attribute=True, no_instantiation=False):
        return clsObject(objectname, is_attribute, ctor_args, no_instantiation)

    def invoke(self, rname, method, args):
        return method(rname, *args)

    def getSlot(self, object, slotname):
        return Literal("%s.%s" % (object, slotname))

    def asString(self, s):
        return as_string(s)

    def _writeOutImports(self):
        for module in self._modules:
            module._writeImportCode()
