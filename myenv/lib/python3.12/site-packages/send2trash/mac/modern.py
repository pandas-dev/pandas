# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

from Foundation import NSFileManager, NSURL
from send2trash.compat import text_type
from send2trash.util import preprocess_paths


def check_op_result(op_result):
    # First value will be false on failure
    if not op_result[0]:
        # Error is in third value, localized failure reason matchs ctypes version
        raise OSError(op_result[2].localizedFailureReason())


def send2trash(paths):
    paths = preprocess_paths(paths)
    paths = [path.decode("utf-8") if not isinstance(path, text_type) else path for path in paths]
    for path in paths:
        file_url = NSURL.fileURLWithPath_(path)
        fm = NSFileManager.defaultManager()
        op_result = fm.trashItemAtURL_resultingItemURL_error_(file_url, None, None)
        check_op_result(op_result)
