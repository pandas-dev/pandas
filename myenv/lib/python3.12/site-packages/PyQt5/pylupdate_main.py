# Copyright (c) 2023 Riverbank Computing Limited <info@riverbankcomputing.com>
# 
# This file is part of PyQt5.
# 
# This file may be used under the terms of the GNU General Public License
# version 3.0 as published by the Free Software Foundation and appearing in
# the file LICENSE included in the packaging of this file.  Please review the
# following information to ensure the GNU General Public License version 3.0
# requirements will be met: http://www.gnu.org/copyleft/gpl.html.
# 
# If you do not wish to use this file under the terms of the GPL version 3.0
# then you may purchase a commercial license.  For more information contact
# info@riverbankcomputing.com.
# 
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


import locale
import sys

from PyQt5.QtCore import (PYQT_VERSION_STR, QDir, QFile, QFileInfo, QIODevice,
        QTextStream)

from .pylupdate import *


def printUsage():
    sys.stderr.write(
"Usage:\n"
"    pylupdate5 [options] project-file\n"
"    pylupdate5 [options] source-files -ts ts-files\n"
"\n"
"Options:\n"
"    -help  Display this information and exit\n"
"    -version\n"
"           Display the version of pylupdate5 and exit\n"
"    -verbose\n"
"           Explain what is being done\n"
"    -noobsolete\n"
"           Drop all obsolete strings\n"
"    -tr-function name\n"
"           name() may be used instead of tr()\n"
"    -translate-function name\n"
"           name() may be used instead of translate()\n");


def updateTsFiles(fetchedTor, tsFileNames, codecForTr, noObsolete, verbose):
    dir = QDir()

    for t in tsFileNames:
        fn = dir.relativeFilePath(t)
        tor = MetaTranslator()
        out = MetaTranslator()

        tor.load(t)

        if codecForTr:
            tor.setCodec(codecForTr)

        merge(tor, fetchedTor, out, noObsolete, verbose, fn)

        if noObsolete:
            out.stripObsoleteMessages()

        out.stripEmptyContexts()

        if not out.save(t):
            sys.stderr.write("pylupdate5 error: Cannot save '%s'\n" % fn)


def _encoded_path(path):
    return path.encode(locale.getdefaultlocale()[1])


def main():
    # Initialise.

    defaultContext = "@default"
    fetchedTor = MetaTranslator()
    codecForTr = ''
    codecForSource = ''
    tsFileNames = []
    uiFileNames = []

    verbose = False
    noObsolete = False
    metSomething = False
    numFiles = 0
    standardSyntax = True
    metTsFlag = False
    tr_func = None
    translate_func = None

    # Parse the command line.

    for arg in sys.argv[1:]:
        if arg == "-ts":
            standardSyntax = False

    argc = len(sys.argv)
    i = 1

    while i < argc:
        arg = sys.argv[i]
        i += 1

        if arg == "-help":
            printUsage()
            sys.exit(0)

        if arg == "-version":
            sys.stderr.write("pylupdate5 v%s\n" % PYQT_VERSION_STR)
            sys.exit(0)

        if arg == "-noobsolete":
            noObsolete = True
            continue

        if arg == "-verbose":
            verbose = True
            continue

        if arg == "-ts":
            metTsFlag = True
            continue

        if arg == "-tr-function":
            if i >= argc:
                sys.stderr.write(
                        "pylupdate5 error: missing -tr-function name\n")
                sys.exit(2)

            tr_func = sys.argv[i]
            i += 1
            continue

        if arg == "-translate-function":
            if i >= argc:
                sys.stderr.write(
                        "pylupdate5 error: missing -translate-function name\n")
                sys.exit(2)

            translate_func = sys.argv[i]
            i += 1
            continue

        numFiles += 1

        fullText = ""

        if not metTsFlag:
            f = QFile(arg)

            if not f.open(QIODevice.ReadOnly):
                sys.stderr.write(
                        "pylupdate5 error: Cannot open file '%s'\n" % arg)
                sys.exit(1)

            t = QTextStream(f)
            fullText = t.readAll()
            f.close()

        if standardSyntax:
            oldDir = QDir.currentPath()
            QDir.setCurrent(QFileInfo(arg).path())

            fetchedTor = MetaTranslator()
            codecForTr = ''
            codecForSource = ''
            tsFileNames = []
            uiFileNames = []

            for key, value in proFileTagMap(fullText).items():
                for t in value.split(' '):
                    if key == "SOURCES":
                        fetchtr_py(
                                _encoded_path(
                                        QDir.current().absoluteFilePath(t)),
                                fetchedTor, defaultContext, True,
                                codecForSource, tr_func, translate_func)
                        metSomething = True

                    elif key == "TRANSLATIONS":
                        tsFileNames.append(QDir.current().absoluteFilePath(t))
                        metSomething = True

                    elif key in ("CODEC", "DEFAULTCODEC", "CODECFORTR"):
                        codecForTr = t
                        fetchedTor.setCodec(codecForTr)

                    elif key == "CODECFORSRC":
                        codecForSource = t

                    elif key == "FORMS":
                        fetchtr_ui(
                                _encoded_path(
                                        QDir.current().absoluteFilePath(t)),
                                fetchedTor, defaultContext, True)

            updateTsFiles(fetchedTor, tsFileNames, codecForTr, noObsolete,
                    verbose)

            if not metSomething:
                sys.stderr.write(
                        "pylupdate5 warning: File '%s' does not look like a "
                        "project file\n" % arg)
            elif len(tsFileNames) == 0:
                sys.stderr.write(
                        "pylupdate5 warning: Met no 'TRANSLATIONS' entry in "
                        "project file '%s'\n" % arg)

            QDir.setCurrent(oldDir)
        else:
            if metTsFlag:
                if arg.lower().endswith(".ts"):
                    fi = QFileInfo(arg)

                    if not fi.exists() or fi.isWritable():
                        tsFileNames.append(arg)
                    else:
                        sys.stderr.write(
                                "pylupdate5 warning: For some reason, I "
                                "cannot save '%s'\n" % arg)
                else:
                    sys.stderr.write(
                            "pylupdate5 error: File '%s' lacks .ts extension\n" % arg)
            else:
                fi = QFileInfo(arg)
                path = _encoded_path(fi.absoluteFilePath())

                if fi.suffix() in ("py", "pyw"):
                    fetchtr_py(path, fetchedTor, defaultContext, True,
                            codecForSource, tr_func, translate_func)
                else:
                    fetchtr_ui(path, fetchedTor, defaultContext, True)

    if not standardSyntax:
        updateTsFiles(fetchedTor, tsFileNames, codecForTr, noObsolete, verbose)

    if numFiles == 0:
        printUsage()
        sys.exit(1)


if __name__ == '__main__':
    main()
