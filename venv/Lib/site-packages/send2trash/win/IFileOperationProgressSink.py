# Sample implementation of IFileOperationProgressSink that just prints
# some basic info

import pythoncom
from win32com.shell import shell, shellcon
from win32com.server.policy import DesignatedWrapPolicy


class FileOperationProgressSink(DesignatedWrapPolicy):
    _com_interfaces_ = [shell.IID_IFileOperationProgressSink]
    _public_methods_ = [
        "StartOperations",
        "FinishOperations",
        "PreRenameItem",
        "PostRenameItem",
        "PreMoveItem",
        "PostMoveItem",
        "PreCopyItem",
        "PostCopyItem",
        "PreDeleteItem",
        "PostDeleteItem",
        "PreNewItem",
        "PostNewItem",
        "UpdateProgress",
        "ResetTimer",
        "PauseTimer",
        "ResumeTimer",
    ]

    def __init__(self):
        self._wrap_(self)
        self.newItem = None

    def PreDeleteItem(self, flags, item):
        # Can detect cases where to stop via flags and condition below, however the operation
        # does not actual stop, we can resort to raising an exception as that does stop things
        # but that may need some additional considerations before implementing.
        return 0 if flags & shellcon.TSF_DELETE_RECYCLE_IF_POSSIBLE else 0x80004005  # S_OK, or E_FAIL

    def PostDeleteItem(self, flags, item, hr_delete, newly_created):
        if newly_created:
            self.newItem = newly_created.GetDisplayName(shellcon.SHGDN_FORPARSING)


def create_sink():
    return pythoncom.WrapObject(FileOperationProgressSink(), shell.IID_IFileOperationProgressSink)
