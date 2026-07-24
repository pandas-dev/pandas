"""\
The `winxpgui` module is obsolete and has been completely replaced \
by `win32gui` and `win32console.GetConsoleWindow`. Use those instead. \
"""

import warnings

from win32console import (
    GetConsoleWindow as GetConsoleWindow,  # noqa: PLC0414 # Explicit re-export
)
from win32gui import *

warnings.warn(str(__doc__), category=DeprecationWarning, stacklevel=2)
