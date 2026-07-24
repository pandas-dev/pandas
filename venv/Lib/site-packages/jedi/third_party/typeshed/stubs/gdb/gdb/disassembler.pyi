from collections.abc import Sequence
from typing import Final, final
from typing_extensions import disjoint_base

import gdb
from gdb import Architecture, Progspace

class Disassembler:
    def __init__(self, name: str) -> None: ...
    def __call__(self, info): ...

@disjoint_base
class DisassembleInfo:
    address: int
    architecture: Architecture
    progspace: Progspace
    def __init__(self, info: DisassembleInfo) -> None: ...
    def address_part(self, address: int) -> DisassemblerAddressPart: ...
    def is_valid(self) -> bool: ...
    def read_memory(self, len: int, offset: int = 0): ...
    def text_part(self, style: int, string: str) -> DisassemblerTextPart: ...

class DisassemblerPart:
    def __init__(self, /, *args, **kwargs) -> None: ...

@final
class DisassemblerAddressPart(DisassemblerPart):
    address: int
    string: str

@final
class DisassemblerTextPart(DisassemblerPart):
    string: str
    style: int

@final
class DisassemblerResult:
    def __init__(self, length: int, string: str | None = None, parts: Sequence[DisassemblerPart] | None = None) -> None: ...
    length: int
    parts: Sequence[DisassemblerPart]
    string: str

STYLE_TEXT: Final = 0
STYLE_MNEMONIC: Final = 1
STYLE_SUB_MNEMONIC: Final = 2
STYLE_ASSEMBLER_DIRECTIVE: Final = 3
STYLE_REGISTER: Final = 4
STYLE_IMMEDIATE: Final = 5
STYLE_ADDRESS: Final = 6
STYLE_ADDRESS_OFFSET: Final = 7
STYLE_SYMBOL: Final = 8
STYLE_COMMENT_START: Final = 9

def builtin_disassemble(info: DisassembleInfo) -> None: ...

class maint_info_py_disassemblers_cmd(gdb.Command):
    def __init__(self) -> None: ...
    def invoke(self, args, from_tty): ...

def register_disassembler(disassembler: type[Disassembler], architecture: str | None = None): ...
