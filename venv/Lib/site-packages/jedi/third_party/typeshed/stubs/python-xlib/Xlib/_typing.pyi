from collections.abc import Callable
from typing import TypeVar
from typing_extensions import TypeAlias

from Xlib.error import XError
from Xlib.protocol.rq import Request

_T = TypeVar("_T")
ErrorHandler: TypeAlias = Callable[[XError, Request | None], _T]
