from _typeshed import Unused
from typing import Literal, TypedDict, type_check_only
from typing_extensions import TypeAlias

@type_check_only
class _ContinueRequestResult(TypedDict):
    allThreadsContinued: bool

_Granularity: TypeAlias = Literal["statement", "instruction"]

def next(*, threadId: int, singleThread: bool = False, granularity: _Granularity = "statement", **args: Unused) -> None: ...
def step_in(*, threadId: int, singleThread: bool = False, granularity: _Granularity = "statement", **args: Unused) -> None: ...
def step_out(*, threadId: int, singleThread: bool = False, **args: Unused): ...
def continue_request(*, threadId: int, singleThread: bool = False, **args: Unused) -> _ContinueRequestResult: ...
