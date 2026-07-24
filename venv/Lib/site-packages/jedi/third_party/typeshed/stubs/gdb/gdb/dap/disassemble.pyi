from _typeshed import Unused
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

from .sources import Source

@type_check_only
class _Instruction(TypedDict):
    address: str
    instruction: str
    instructionBytes: str
    symbol: NotRequired[str]  # only set if there's a corresponding label
    line: NotRequired[int]  # only set if source is available
    location: NotRequired[Source]  # only set if source is available

@type_check_only
class _DisassembleResult(TypedDict):
    instructions: list[_Instruction]

def disassemble(
    *, memoryReference: str, offset: int = 0, instructionOffset: int = 0, instructionCount: int, **extra: Unused
) -> _DisassembleResult: ...
