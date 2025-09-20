#############################################################################
##
## Copyright (C) 2015 Riverbank Computing Limited.
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

from PyQt5 import QtGui, QtWidgets


class _QtWrapper(object):
    @classmethod
    def search(cls, name):
        return getattr(cls.module, name, None)


class _QtGuiWrapper(_QtWrapper):
    module = QtGui


class _QtWidgetsWrapper(_QtWrapper):
    module = QtWidgets


class _ModuleWrapper(object):
    def __init__(self, moduleName, classes):
        self._moduleName = moduleName
        self._module = None
        self._classes = classes

    def search(self, cls):
        if cls in self._classes:
            if self._module is None:
                self._module = __import__(self._moduleName, {}, {}, self._classes)
            # Remove any C++ scope.
            cls = cls.split('.')[-1]

            return getattr(self._module, cls)

        return None


class _CustomWidgetLoader(object):
    def __init__(self, package):
        # should it stay this way?
        if '.' not in sys.path:
            sys.path.append('.')

        self._widgets = {}
        self._modules = {}
        self._package = package
        
    def addCustomWidget(self, widgetClass, baseClass, module):
        assert widgetClass not in self._widgets
        self._widgets[widgetClass] = module
    
    def search(self, cls):
        module_name = self._widgets.get(cls)
        if module_name is None:
            return None

        module = self._modules.get(module_name)
        if module is None:
            if module_name.startswith('.'):
                if self._package == '':
                    raise ImportError(
                            "relative import of %s without base package specified" % module_name)

                if self._package.startswith('.'):
                    raise ImportError(
                            "base package %s is relative" % self._package)

                mname = self._package + module_name
            else:
                mname = module_name

            try:
                module = __import__(mname, {}, {}, (cls,))
            except ValueError:
                # Raise a more helpful exception.
                raise ImportError("unable to import module %s" % mname)

            self._modules[module_name] = module

        return getattr(module, cls)


class LoaderCreatorPolicy(object):
    def __init__(self, package):
        self._package = package

    def createQtGuiWidgetsWrappers(self):
        return [_QtGuiWrapper, _QtWidgetsWrapper]
    
    def createModuleWrapper(self, moduleName, classes):
        return _ModuleWrapper(moduleName, classes)
    
    def createCustomWidgetLoader(self):
        return _CustomWidgetLoader(self._package)

    def instantiate(self, clsObject, objectName, ctor_args, is_attribute=True):
        return clsObject(*ctor_args)

    def invoke(self, rname, method, args):
        return method(*args)

    def getSlot(self, object, slotname):
        # Rename slots that correspond to Python keyword arguments.
        if slotname == 'raise':
            slotname += '_'

        return getattr(object, slotname)

    def asString(self, s):
        return s
