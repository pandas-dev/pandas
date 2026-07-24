from logging import Logger
from traceback import StackSummary
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

log: Logger

@type_check_only
class _StackInfo(TypedDict):
    path: str
    line: int
    label: str

@type_check_only
class _ThrowableAttrs(TypedDict):
    id: str
    message: NotRequired[str]
    type: str
    remote: bool
    stack: NotRequired[list[_StackInfo]]

class Throwable:
    id: str
    message: str
    type: str
    remote: bool
    stack: list[_StackInfo] | None
    def __init__(self, exception: Exception, stack: StackSummary, remote: bool = False) -> None: ...
    def to_dict(self) -> _ThrowableAttrs: ...
