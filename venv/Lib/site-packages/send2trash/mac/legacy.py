# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

from __future__ import unicode_literals

from ctypes import cdll, byref, Structure, c_char, c_char_p
from ctypes.util import find_library

from send2trash.util import preprocess_paths

Foundation = cdll.LoadLibrary(find_library("Foundation"))
CoreServices = cdll.LoadLibrary(find_library("CoreServices"))

GetMacOSStatusCommentString = Foundation.GetMacOSStatusCommentString
GetMacOSStatusCommentString.restype = c_char_p
FSPathMakeRefWithOptions = CoreServices.FSPathMakeRefWithOptions
FSMoveObjectToTrashSync = CoreServices.FSMoveObjectToTrashSync

kFSPathMakeRefDefaultOptions = 0
kFSPathMakeRefDoNotFollowLeafSymlink = 0x01

kFSFileOperationDefaultOptions = 0
kFSFileOperationOverwrite = 0x01
kFSFileOperationSkipSourcePermissionErrors = 0x02
kFSFileOperationDoNotMoveAcrossVolumes = 0x04
kFSFileOperationSkipPreflight = 0x08


class FSRef(Structure):
    _fields_ = [("hidden", c_char * 80)]


def check_op_result(op_result):
    if op_result:
        msg = GetMacOSStatusCommentString(op_result).decode("utf-8")
        raise OSError(msg)


def send2trash(paths):
    paths = preprocess_paths(paths)
    paths = [path.encode("utf-8") if not isinstance(path, bytes) else path for path in paths]
    for path in paths:
        fp = FSRef()
        opts = kFSPathMakeRefDoNotFollowLeafSymlink
        op_result = FSPathMakeRefWithOptions(path, opts, byref(fp), None)
        check_op_result(op_result)
        opts = kFSFileOperationDefaultOptions
        op_result = FSMoveObjectToTrashSync(byref(fp), None, opts)
        check_op_result(op_result)
