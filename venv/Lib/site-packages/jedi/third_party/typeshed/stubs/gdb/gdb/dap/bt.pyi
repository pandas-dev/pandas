from _typeshed import Unused
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

from .sources import Source
from .varref import _ValueFormat

@type_check_only
class _StackFrameFormat(_ValueFormat, total=False):
    parameters: bool
    parameterTypes: bool
    parameterNames: bool
    parameterValues: bool
    line: bool
    module: bool
    includeAll: bool

@type_check_only
class _StackFrame(TypedDict):
    id: int
    name: str
    line: int
    column: int
    instructionPointerReference: str
    moduleId: NotRequired[str | None]
    source: NotRequired[Source]

@type_check_only
class _StackTraceResult(TypedDict):
    stackFrames: list[_StackFrame]

def check_stack_frame(
    *,
    # From source:
    #   Note that StackFrameFormat extends ValueFormat, which is why
    #   "hex" appears here.
    hex: bool = False,
    parameters: bool = False,
    parameterTypes: bool = False,
    parameterNames: bool = False,
    parameterValues: bool = False,
    line: bool = False,
    module: bool = False,
    includeAll: bool = False,
    **rest: Unused,
) -> _StackFrameFormat: ...
def stacktrace(
    *, levels: int = 0, startFrame: int = 0, threadId: int, format: _StackFrameFormat | None = None, **extra: Unused
) -> _StackTraceResult: ...
