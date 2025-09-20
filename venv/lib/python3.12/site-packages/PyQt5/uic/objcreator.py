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


import os.path

from .exceptions import NoSuchWidgetError, WidgetPluginError


# The list of directories that are searched for widget plugins.  This is
# exposed as part of the API.
widgetPluginPath = [os.path.join(os.path.dirname(__file__), 'widget-plugins')]


MATCH = True
NO_MATCH = False
MODULE = 0
CW_FILTER = 1


class QObjectCreator(object):    
    def __init__(self, creatorPolicy):
        self._cpolicy = creatorPolicy

        self._cwFilters = []
        self._modules = self._cpolicy.createQtGuiWidgetsWrappers()

        # Get the optional plugins.
        for plugindir in widgetPluginPath:
            try:
                plugins = os.listdir(plugindir)
            except:
                plugins = []

            for filename in plugins:
                if not filename.endswith('.py'):
                    continue

                filename = os.path.join(plugindir, filename)

                plugin_globals = {
                    "MODULE": MODULE,
                    "CW_FILTER": CW_FILTER,
                    "MATCH": MATCH,
                    "NO_MATCH": NO_MATCH}

                plugin_locals = {}

                if self.load_plugin(filename, plugin_globals, plugin_locals):
                    pluginType = plugin_locals["pluginType"]
                    if pluginType == MODULE:
                        modinfo = plugin_locals["moduleInformation"]()
                        self._modules.append(self._cpolicy.createModuleWrapper(*modinfo))
                    elif pluginType == CW_FILTER:
                        self._cwFilters.append(plugin_locals["getFilter"]())
                    else:
                        raise WidgetPluginError("Unknown plugin type of %s" % filename)

        self._customWidgets = self._cpolicy.createCustomWidgetLoader()
        self._modules.append(self._customWidgets)

    def createQObject(self, classname, *args, **kwargs):
        # Handle regular and custom widgets.
        factory = self.findQObjectType(classname)

        if factory is None:
            # Handle scoped names, typically static factory methods.
            parts = classname.split('.')

            if len(parts) > 1:
                factory = self.findQObjectType(parts[0])

                if factory is not None:
                    for part in parts[1:]:
                        factory = getattr(factory, part, None)
                        if factory is None:
                            break

            if factory is None:
                raise NoSuchWidgetError(classname)

        return self._cpolicy.instantiate(factory, *args, **kwargs)

    def invoke(self, rname, method, args=()):
        return self._cpolicy.invoke(rname, method, args)

    def findQObjectType(self, classname):
        for module in self._modules:
            w = module.search(classname)
            if w is not None:
                return w
        return None

    def getSlot(self, obj, slotname):
        return self._cpolicy.getSlot(obj, slotname)

    def asString(self, s):
        return self._cpolicy.asString(s)

    def addCustomWidget(self, widgetClass, baseClass, module):
        for cwFilter in self._cwFilters:
            match, result = cwFilter(widgetClass, baseClass, module)
            if match:
                widgetClass, baseClass, module = result
                break

        self._customWidgets.addCustomWidget(widgetClass, baseClass, module)

    @staticmethod
    def load_plugin(filename, plugin_globals, plugin_locals):
        """ Load the plugin from the given file.  Return True if the plugin was
        loaded, or False if it wanted to be ignored.  Raise an exception if
        there was an error.
        """

        plugin = open(filename)

        try:
            exec(plugin.read(), plugin_globals, plugin_locals)
        except ImportError:
            return False
        except Exception as e:
            raise WidgetPluginError("%s: %s" % (e.__class__, str(e)))
        finally:
            plugin.close()

        return True
