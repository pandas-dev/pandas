import types
from typing import Final

PYTHON_VERSION: Final[tuple[int, int]]
CONDITIONAL_JUMPS: Final[list[str]]
UNCONDITIONAL_JUMPS: Final[list[str]]
ENDS_FUNCTION: Final[list[str]]
HAVE_REL_REFERENCE: Final[list[str]]
HAVE_ABS_REFERENCE: Final[list[str]]
REL_REFERENCE_IS_INVERTED: Final[list[str]]

def rel_reference_scale(opname: str) -> int: ...

REVERSE_CMP_OP: Final[list[int]]

def jump_arg_bytes(arg: int) -> int: ...
def add_bytes_to_jump_arg(arg: int, size: int) -> int: ...

class ExceptionTableEntry:
    def __init__(self, start_offset: int, end_offset: int, target: int, depth: int, lasti: bool) -> None: ...
    def __eq__(self, other: object) -> bool: ...

class ExceptionTable:
    def __init__(self, entries: list[ExceptionTableEntry]) -> None: ...
    def __eq__(self, other: object) -> bool: ...

def generate_exceptiontable(original_code: types.CodeType, exception_table_entries: list[ExceptionTableEntry]) -> bytes: ...
