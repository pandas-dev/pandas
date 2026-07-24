# Magic utility that "redirects" to pythoncomXX.dll
from typing import Any

import pywintypes

pywintypes.__import_pywin32_system_module__("pythoncom", globals())


# This module dynamically re-exports from a C-Extension.
# `__getattr__() -> Any` prevents checkers attribute access errors
# without the use of external type stubs.
# So keep an empty stub behind `if TYPE_CHECKING` if pythoncom.frozen support is removed.
def __getattr__(name: str) -> Any:
    if name == "frozen":
        import sys
        import warnings

        warnings.warn(
            f"`pythoncom.frozen` used to expose `Py_FrozenFlag` from the C API. "
            + "`Py_FrozenFlag` is deprecated since Python 3.12. "
            + "Ever since pywin32 b200, loading the `win32com` module has silently "
            + "been replacing `pythoncom.frozen` with `sys.frozen`. "
            + 'Use `getattr(sys, "frozen", False)` directly instead.',
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(sys, "frozen", False)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
