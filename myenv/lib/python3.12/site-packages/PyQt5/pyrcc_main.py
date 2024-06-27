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


import sys

from PyQt5.QtCore import PYQT_VERSION_STR, QDir, QFile

from .pyrcc import *


# Initialise the globals.
verbose = False
compressLevel = CONSTANT_COMPRESSLEVEL_DEFAULT
compressThreshold = CONSTANT_COMPRESSTHRESHOLD_DEFAULT
resourceRoot = ''


def processResourceFile(filenamesIn, filenameOut, listFiles):
    if verbose:
        sys.stderr.write("PyQt5 resource compiler\n")

    # Setup.
    library = RCCResourceLibrary()
    library.setInputFiles(filenamesIn)
    library.setVerbose(verbose)
    library.setCompressLevel(compressLevel)
    library.setCompressThreshold(compressThreshold)
    library.setResourceRoot(resourceRoot)

    if not library.readFiles():
        return False

    if filenameOut == '-':
        filenameOut = ''

    if listFiles:
        # Open the output file or use stdout if not specified.
        if filenameOut:
            try:
                out_fd = open(filenameOut, 'w')
            except Exception:
                sys.stderr.write(
                        "Unable to open %s for writing\n" % filenameOut)
                return False
        else:
            out_fd = sys.stdout

        for df in library.dataFiles():
            out_fd.write("%s\n" % QDir.cleanPath(df))

        if out_fd is not sys.stdout:
            out_fd.close()

        return True

    return library.output(filenameOut)


def showHelp(error):
    sys.stderr.write("PyQt5 resource compiler\n")

    if error:
        sys.stderr.write("pyrcc5: %s\n" % error)

    sys.stderr.write(
"Usage: pyrcc5 [options] <inputs>\n"
"\n"
"Options:\n"
"    -o file           Write output to file rather than stdout\n"
"    -threshold level  Threshold to consider compressing files\n"
"    -compress level   Compress input files by level\n"
"    -root path        Prefix resource access path with root path\n"
"    -no-compress      Disable all compression\n"
"    -version          Display version\n"
"    -help             Display this information\n")


def main():
    # Parse the command line.  Note that this mimics the original C++ (warts
    # and all) in order to preserve backwards compatibility.
    global verbose
    global compressLevel
    global compressThreshold
    global resourceRoot

    outFilename = ''
    helpRequested = False
    listFiles = False
    files = []

    errorMsg = None
    argc = len(sys.argv)
    i = 1

    while i < argc:
        arg = sys.argv[i]
        i += 1

        if arg[0] == '-':
            opt = arg[1:]

            if opt == "o":
                if i >= argc:
                    errorMsg = "Missing output name"
                    break

                outFilename = sys.argv[i]
                i += 1

            elif opt == "root":
                if i >= argc:
                    errorMsg = "Missing root path"
                    break

                resourceRoot = QDir.cleanPath(sys.argv[i])
                i += 1

                if resourceRoot == '' or resourceRoot[0] != '/':
                    errorMsg = "Root must start with a /"
                    break

            elif opt == "compress":
                if i >= argc:
                    errorMsg = "Missing compression level"
                    break

                compressLevel = int(sys.argv[i])
                i += 1

            elif opt == "threshold":
                if i >= argc:
                    errorMsg = "Missing compression threshold"
                    break

                compressThreshold = int(sys.argv[i])
                i += 1

            elif opt == "verbose":
                verbose = True

            elif opt == "list":
                listFiles = True

            elif opt == "version":
                sys.stderr.write("pyrcc5 v%s\n" % PYQT_VERSION_STR)
                sys.exit(1)

            elif opt == "help" or opt == "h":
                helpRequested = True

            elif opt == "no-compress":
                compressLevel = -2

            else:
                errorMsg = "Unknown option: '%s'" % arg
                break
        else:
            if not QFile.exists(arg):
                sys.stderr.write(
                        "%s: File does not exist '%s'\n" % (sys.argv[0], arg))
                sys.exit(1)

            files.append(arg)

    # Handle any errors or a request for help.
    if len(files) == 0 or errorMsg is not None or helpRequested:
        showHelp(errorMsg)
        sys.exit(1)

    if not processResourceFile(files, outFilename, listFiles):
        sys.exit(1)


if __name__ == '__main__':
    main()
