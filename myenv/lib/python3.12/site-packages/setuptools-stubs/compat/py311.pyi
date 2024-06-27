from _typeshed import StrOrBytesPath
from shutil import _OnExcCallback

def shutil_rmtree(path: StrOrBytesPath, ignore_errors: bool = False, onexc: _OnExcCallback | None = None) -> None: ...
