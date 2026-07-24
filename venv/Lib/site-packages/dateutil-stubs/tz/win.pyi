import sys
from ctypes import _NameTypes, _Pointer, c_wchar
from datetime import datetime, timedelta
from typing import Any, ClassVar, Final

from dateutil.tz import tzwin as tzwin, tzwinlocal as tzwinlocal

if sys.platform == "win32":
    from winreg import _KeyType

    __all__ = ["tzwin", "tzwinlocal", "tzres"]

    ONEWEEK: timedelta
    TZKEYNAMENT: Final[str]
    TZKEYNAME9X: Final[str]
    TZLOCALKEYNAME: Final[str]
    TZKEYNAME: Final[str]

    class tzres:
        p_wchar: ClassVar[type[_Pointer[c_wchar]]]
        tzres_loc: _NameTypes
        def __init__(self, tzres_loc: _NameTypes = "tzres.dll") -> None: ...
        def load_name(self, offset: int) -> str: ...
        def name_from_string(self, tzname_str: str) -> str: ...

    def picknthweekday(year: int, month: int, dayofweek: int, hour: int, minute: int, whichweek: int) -> datetime: ...
    def valuestodict(key: _KeyType) -> dict[str, Any]: ...  # keys and values in dict are results of winreg.EnumValue() function
