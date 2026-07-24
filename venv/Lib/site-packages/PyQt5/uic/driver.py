#############################################################################
##
## Copyright (c) 2024 Riverbank Computing Limited <info@riverbankcomputing.com>
## 
## This file is part of PyQt5.
## 
## This file may be used under the terms of the GNU General Public License
## version 3.0 as published by the Free Software Foundation and appearing in
## the file LICENSE included in the packaging of this file.  Please review the
## following information to ensure the GNU General Public License version 3.0
## requirements will be met: http://www.gnu.org/copyleft/gpl.html.
## 
## If you do not wish to use this file under the terms of the GPL version 3.0
## then you may purchase a commercial license.  For more information contact
## info@riverbankcomputing.com.
## 
## This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
## WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
##
#############################################################################


import sys
import logging

from . import compileUi, loadUi


class Driver(object):
    """ This encapsulates access to the pyuic functionality so that it can be
    called by code that is Python v2/v3 specific.
    """

    LOGGER_NAME = 'PyQt5.uic'

    def __init__(self, opts, ui_file):
        """ Initialise the object.  opts is the parsed options.  ui_file is the
        name of the .ui file.
        """

        if opts.debug:
            logger = logging.getLogger(self.LOGGER_NAME)
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(name)s: %(message)s"))
            logger.addHandler(handler)
            logger.setLevel(logging.DEBUG)

        self._opts = opts
        self._ui_file = ui_file

    def invoke(self):
        """ Invoke the action as specified by the parsed options.  Returns 0 if
        there was no error.
        """

        if self._opts.preview:
            return self._preview()

        self._generate()

        return 0

    def _preview(self):
        """ Preview the .ui file.  Return the exit status to be passed back to
        the parent process.
        """

        from PyQt5 import QtWidgets

        app = QtWidgets.QApplication([self._ui_file])
        widget = loadUi(self._ui_file)
        widget.show()

        return app.exec_()

    def _generate(self):
        """ Generate the Python code. """

        needs_close = False

        if sys.hexversion >= 0x03000000:
            if self._opts.output == '-':
                from io import TextIOWrapper

                pyfile = TextIOWrapper(sys.stdout.buffer, encoding='utf8')
            else:
                pyfile = open(self._opts.output, 'wt', encoding='utf8')
                needs_close = True
        else:
            if self._opts.output == '-':
                pyfile = sys.stdout
            else:
                pyfile = open(self._opts.output, 'wt')
                needs_close = True

        import_from = self._opts.import_from

        if import_from:
            from_imports = True
        elif self._opts.from_imports:
            from_imports = True
            import_from = '.'
        else:
            from_imports = False

        compileUi(self._ui_file, pyfile, self._opts.execute, self._opts.indent,
                from_imports, self._opts.resource_suffix, import_from)

        if needs_close:
            pyfile.close()

    def on_IOError(self, e):
        """ Handle an IOError exception. """

        sys.stderr.write("Error: %s: \"%s\"\n" % (e.strerror, e.filename))

    def on_SyntaxError(self, e):
        """ Handle a SyntaxError exception. """

        sys.stderr.write("Error in input file: %s\n" % e)

    def on_NoSuchClassError(self, e):
        """ Handle a NoSuchClassError exception. """

        sys.stderr.write(str(e) + "\n")

    def on_NoSuchWidgetError(self, e):
        """ Handle a NoSuchWidgetError exception. """

        sys.stderr.write(str(e) + "\n")

    def on_Exception(self, e):
        """ Handle a generic exception. """

        if logging.getLogger(self.LOGGER_NAME).level == logging.DEBUG:
            import traceback

            traceback.print_exception(*sys.exc_info())
        else:
            from PyQt5 import QtCore

            sys.stderr.write("""An unexpected error occurred.
Check that you are using the latest version of PyQt5 and send an error report to
support@riverbankcomputing.com, including the following information:

  * your version of PyQt (%s)
  * the UI file that caused this error
  * the debug output of pyuic5 (use the -d flag when calling pyuic5)
""" % QtCore.PYQT_VERSION_STR)
