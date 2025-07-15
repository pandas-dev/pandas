# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

from __future__ import unicode_literals
import os.path as op
from send2trash.compat import text_type
from send2trash.util import preprocess_paths
from platform import version
import pythoncom
import pywintypes
from win32com.shell import shell, shellcon
from send2trash.win.IFileOperationProgressSink import create_sink


def send2trash(paths):
    paths = preprocess_paths(paths)
    if not paths:
        return
    # convert data type
    paths = [text_type(path, "mbcs") if not isinstance(path, text_type) else path for path in paths]
    # convert to full paths
    paths = [op.abspath(path) if not op.isabs(path) else path for path in paths]
    # remove the leading \\?\ if present
    paths = [path[4:] if path.startswith("\\\\?\\") else path for path in paths]
    # Need to initialize the com before using
    pythoncom.CoInitialize()
    # create instance of file operation object
    fileop = pythoncom.CoCreateInstance(
        shell.CLSID_FileOperation,
        None,
        pythoncom.CLSCTX_ALL,
        shell.IID_IFileOperation,
    )
    # default flags to use
    flags = shellcon.FOF_NOCONFIRMATION | shellcon.FOF_NOERRORUI | shellcon.FOF_SILENT | shellcon.FOFX_EARLYFAILURE
    # determine rest of the flags based on OS version
    # use newer recommended flags if available
    if int(version().split(".", 1)[0]) >= 8:
        flags |= 0x20000000 | 0x00080000  # FOFX_ADDUNDORECORD win 8+  # FOFX_RECYCLEONDELETE win 8+
    else:
        flags |= shellcon.FOF_ALLOWUNDO
    # set the flags
    fileop.SetOperationFlags(flags)
    # actually try to perform the operation, this section may throw a
    # pywintypes.com_error which does not seem to create as nice of an
    # error as OSError so wrapping with try to convert
    sink = create_sink()
    try:
        for path in paths:
            item = shell.SHCreateItemFromParsingName(path, None, shell.IID_IShellItem)
            fileop.DeleteItem(item, sink)
        result = fileop.PerformOperations()
        aborted = fileop.GetAnyOperationsAborted()
        # if non-zero result or aborted throw an exception
        if result or aborted:
            raise OSError(None, None, paths, result)
    except pywintypes.com_error as error:
        # convert to standard OS error, allows other code to get a
        # normal errno
        raise OSError(None, error.strerror, path, error.hresult)
    finally:
        # Need to make sure we call this once fore every init
        pythoncom.CoUninitialize()
