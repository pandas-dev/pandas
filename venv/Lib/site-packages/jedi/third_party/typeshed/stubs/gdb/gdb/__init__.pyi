# The GDB Python API is implemented in C, so the type hints below were made
# reading the documentation
# (https://sourceware.org/gdb/onlinedocs/gdb/Python-API.html).

import _typeshed
import threading
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import AbstractContextManager
from typing import Any, Final, Generic, Literal, Protocol, TypedDict, TypeVar, final, overload, type_check_only
from typing_extensions import TypeAlias, deprecated, disjoint_base

import gdb.FrameDecorator
import gdb.types
from gdb.missing_debug import MissingDebugHandler
from gdb.missing_files import MissingFileHandler

# The following submodules are automatically imported
from . import events as events, printing as printing, prompt as prompt, types as types

# Basic

VERSION: str
HOST_CONFIG: str
TARGET_CONFIG: str

PYTHONDIR: str

STDOUT: int
STDERR: int
STDLOG: int

@overload
def execute(command: str, from_tty: bool = False, to_string: Literal[False] = False) -> None: ...
@overload
def execute(command: str, *, to_string: Literal[True]) -> str: ...
@overload
def execute(command: str, from_tty: bool, to_string: Literal[True]) -> str: ...
@overload
def execute(command: str, from_tty: bool = False, to_string: bool = False) -> str | None: ...
def breakpoints() -> Sequence[Breakpoint]: ...
def rbreak(regex: str, minsyms: bool = ..., throttle: int = ..., symtabs: Iterator[Symtab] = ...) -> list[Breakpoint]: ...
def parameter(parameter: str, /) -> bool | int | str | None: ...
def set_parameter(name: str, value: bool | int | str | None) -> None: ...
def with_parameter(name: str, value: bool | int | str | None) -> AbstractContextManager[None]: ...
def history(number: int, /) -> Value: ...
def add_history(value: Value, /) -> int: ...
def history_count() -> int: ...
def convenience_variable(name: str, /) -> Value | None: ...
def set_convenience_variable(name: str, value: _ValueOrNative | None, /) -> None: ...
def parse_and_eval(expression: str, global_context: bool = False) -> Value: ...
def format_address(address: int, progspace: Progspace = ..., architecture: Architecture = ...): ...
def find_pc_line(pc: int | Value) -> Symtab_and_line: ...
def post_event(event: Callable[[], object], /) -> None: ...
def write(string: str, stream: int = ...) -> None: ...
def flush(stream: int = ...) -> None: ...
def target_charset() -> str: ...
def target_wide_charset() -> str: ...
def host_charset() -> str: ...
def solib_name(addr: int) -> str | None: ...
def decode_line(expression: str = ..., /) -> tuple[str | None, tuple[Symtab_and_line, ...] | None]: ...
def architecture_names() -> list[str]: ...
def connections() -> list[TargetConnection]: ...

prompt_hook: Callable[[str], str | None]

# Exceptions

class error(RuntimeError): ...
class MemoryError(error): ...
class GdbError(Exception): ...

# Values

_ValueOrNative: TypeAlias = bool | int | float | str | Value | LazyString
_ValueOrInt: TypeAlias = Value | int

@disjoint_base
class Value:
    address: Value | None
    is_optimized_out: bool
    type: Type
    dynamic_type: Type
    is_lazy: bool
    bytes: bytes

    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __float__(self) -> float: ...
    def __add__(self, other: _ValueOrNative, /) -> Value: ...
    def __radd__(self, other: _ValueOrNative, /) -> Value: ...
    def __sub__(self, other: _ValueOrNative, /) -> Value: ...
    def __rsub__(self, other: _ValueOrNative, /) -> Value: ...
    def __mul__(self, other: _ValueOrNative, /) -> Value: ...
    def __rmul__(self, other: _ValueOrNative, /) -> Value: ...
    def __truediv__(self, other: _ValueOrNative, /) -> Value: ...
    def __rtruediv__(self, other: _ValueOrNative, /) -> Value: ...
    def __mod__(self, other: _ValueOrNative, /) -> Value: ...
    def __rmod__(self, other: _ValueOrNative, /) -> Value: ...
    def __pow__(self, other: _ValueOrNative, mod: None = None, /) -> Value: ...
    def __and__(self, other: _ValueOrNative, /) -> Value: ...
    def __or__(self, other: _ValueOrNative, /) -> Value: ...
    def __xor__(self, other: _ValueOrNative, /) -> Value: ...
    def __lshift__(self, other: _ValueOrNative, /) -> Value: ...
    def __rshift__(self, other: _ValueOrNative, /) -> Value: ...
    def __eq__(self, other: _ValueOrNative, /) -> bool: ...  # type: ignore[override]
    def __ne__(self, other: _ValueOrNative, /) -> bool: ...  # type: ignore[override]
    def __lt__(self, other: _ValueOrNative, /) -> bool: ...
    def __le__(self, other: _ValueOrNative, /) -> bool: ...
    def __gt__(self, other: _ValueOrNative, /) -> bool: ...
    def __ge__(self, other: _ValueOrNative, /) -> bool: ...
    def __getitem__(self, key: int | str | Field, /) -> Value: ...
    def __call__(self, *args: _ValueOrNative) -> Value: ...
    def __init__(self, val: _ValueOrNative, type: Type | None = None) -> None: ...
    def cast(self, type: Type) -> Value: ...
    def dereference(self) -> Value: ...
    def referenced_value(self) -> Value: ...
    def reference_value(self) -> Value: ...
    def rvalue_reference_value(self) -> Value: ...
    def const_value(self) -> Value: ...
    def dynamic_cast(self, type: Type) -> Value: ...
    def reinterpret_cast(self, type: Type) -> Value: ...
    def format_string(
        self,
        raw: bool = ...,
        pretty_arrays: bool = ...,
        pretty_structs: bool = ...,
        array_indexes: bool = ...,
        symbols: bool = ...,
        unions: bool = ...,
        address: bool = ...,
        styling: bool = ...,
        nibbles: bool = ...,
        deref_refs: bool = ...,
        actual_objects: bool = ...,
        static_members: bool = ...,
        max_characters: int = ...,
        max_elements: int = ...,
        max_depth: int = ...,
        repeat_threshold: int = ...,
        format: str = ...,
    ) -> str: ...
    def string(self, encoding: str = ..., errors: str = ..., length: int = ...) -> str: ...
    def lazy_string(self, encoding: str = ..., length: int = ...) -> LazyString: ...
    def fetch_lazy(self) -> None: ...
    def assign(self, val): ...
    def to_array(self): ...

# Types

def lookup_type(name: str, block: Block = ...) -> Type: ...
@final
class Type(Mapping[str, Field]):
    alignof: int
    code: int
    dynamic: bool
    name: str
    sizeof: int
    tag: str | None
    objfile: Objfile | None
    is_scalar: bool
    is_signed: bool
    is_array_like: bool
    is_string_like: bool

    def fields(self) -> list[Field]: ...
    def array(self, n1: int | Value, n2: int | Value = ...) -> Type: ...
    def vector(self, n1: int, n2: int = ...) -> Type: ...
    def iteritems(self) -> TypeIterator[tuple[str, Field]]: ...
    def iterkeys(self) -> TypeIterator[str]: ...
    def itervalues(self) -> TypeIterator[Field]: ...
    def const(self) -> Type: ...
    def volatile(self) -> Type: ...
    def unqualified(self) -> Type: ...
    def range(self) -> tuple[int, int]: ...
    def reference(self) -> Type: ...
    def pointer(self) -> Type: ...
    def strip_typedefs(self) -> Type: ...
    def target(self) -> Type: ...
    def template_argument(self, n: int, block: Block = ...) -> Type: ...
    def optimized_out(self) -> Value: ...
    def get(self, key: str, default: Any = ...) -> Field | Any: ...
    def has_key(self, key: str) -> bool: ...
    def __len__(self) -> int: ...
    def __getitem__(self, key: str, /) -> Field: ...
    def __iter__(self) -> TypeIterator[str]: ...

_T = TypeVar("_T")

@final
class TypeIterator(Iterator[_T]):
    def __next__(self) -> _T: ...

@final
class Field:
    bitpos: int | None
    enumval: int
    name: str | None
    artificial: bool
    is_base_class: bool
    bitsize: int
    type: Type | None
    parent_type: Type

TYPE_CODE_BITSTRING: Final = -1
TYPE_CODE_PTR: Final = 1
TYPE_CODE_ARRAY: Final = 2
TYPE_CODE_STRUCT: Final = 3
TYPE_CODE_UNION: Final = 4
TYPE_CODE_ENUM: Final = 5
TYPE_CODE_FLAGS: Final = 6
TYPE_CODE_FUNC: Final = 7
TYPE_CODE_INT: Final = 8
TYPE_CODE_FLT: Final = 9
TYPE_CODE_VOID: Final = 10
TYPE_CODE_SET: Final = 11
TYPE_CODE_RANGE: Final = 12
TYPE_CODE_STRING: Final = 13
TYPE_CODE_ERROR: Final = 14
TYPE_CODE_METHOD: Final = 15
TYPE_CODE_METHODPTR: Final = 16
TYPE_CODE_MEMBERPTR: Final = 17
TYPE_CODE_REF: Final = 18
TYPE_CODE_RVALUE_REF: Final = 19
TYPE_CODE_CHAR: Final = 20
TYPE_CODE_BOOL: Final = 21
TYPE_CODE_COMPLEX: Final = 22
TYPE_CODE_TYPEDEF: Final = 23
TYPE_CODE_NAMESPACE: Final = 24
TYPE_CODE_DECFLOAT: Final = 25
TYPE_CODE_MODULE: Final = 26
TYPE_CODE_INTERNAL_FUNCTION: Final = 27
TYPE_CODE_XMETHOD: Final = 28
TYPE_CODE_FIXED_POINT: Final = 29
TYPE_CODE_NAMELIST: Final = 30

SEARCH_UNDEF_DOMAIN: Final[int]
SEARCH_VAR_DOMAIN: Final[int]
SEARCH_STRUCT_DOMAIN: Final[int]
SEARCH_MODULE_DOMAIN: Final[int]
SEARCH_LABEL_DOMAIN: Final[int]
SEARCH_COMMON_BLOCK_DOMAIN: Final[int]
SEARCH_TYPE_DOMAIN: Final[int]
SEARCH_FUNCTION_DOMAIN: Final[int]

# Pretty Printing

@type_check_only
class _PrettyPrinter(Protocol):
    # TODO: The "children" and "display_hint" methods are optional for
    # pretty-printers. Unfortunately, there is no such thing as an optional
    # method in the type system at the moment.
    #
    # def children(self) -> Iterator[tuple[str, _ValueOrNative]]: ...
    # def display_hint(self) -> str | None: ...
    def to_string(self) -> str | LazyString: ...

_PrettyPrinterLookupFunction: TypeAlias = Callable[[Value], _PrettyPrinter | None]

def default_visualizer(value: Value, /) -> _PrettyPrinter | None: ...

# Selecting Pretty-Printers

pretty_printers: list[_PrettyPrinterLookupFunction]
type_printers: list[gdb.types._TypePrinter]

# Filtering Frames

@type_check_only
class _FrameFilter(Protocol):
    name: str
    enabled: bool
    priority: int

    def filter(
        self, iterator: Iterator[gdb.FrameDecorator.FrameDecorator | gdb.FrameDecorator.DAPFrameDecorator]
    ) -> Iterator[gdb.FrameDecorator.FrameDecorator | gdb.FrameDecorator.DAPFrameDecorator]: ...

frame_filters: dict[str, _FrameFilter]

# Unwinding Frames

@final
class PendingFrame:
    def read_register(self, register: str | RegisterDescriptor | int) -> Value: ...
    def create_unwind_info(self, frame_id: object, /) -> UnwindInfo: ...
    def architecture(self) -> Architecture: ...
    def language(self): ...
    def level(self) -> int: ...
    def name(self) -> str: ...
    def pc(self) -> int: ...
    def block(self) -> Block: ...
    def find_sal(self) -> Symtab_and_line: ...
    def function(self) -> Symbol: ...
    def is_valid(self) -> bool: ...

@final
class UnwindInfo:
    def add_saved_register(self, register: str | RegisterDescriptor | int, value: Value) -> None: ...

@type_check_only
class _Unwinder(Protocol):
    @property
    def name(self) -> str: ...
    enabled: bool

    def __call__(self, pending_frame: PendingFrame) -> UnwindInfo | None: ...

frame_unwinders: list[_Unwinder]

# Inferiors

def inferiors() -> tuple[Inferior, ...]: ...
def selected_inferior() -> Inferior: ...

_BufferType: TypeAlias = _typeshed.ReadableBuffer

@final
class Inferior:
    num: int
    connection: TargetConnection | None
    connection_num: int | None
    pid: int
    was_attached: bool
    progspace: Progspace
    main_name: str | None
    @property
    def arguments(self) -> str | None: ...
    @arguments.setter
    def arguments(self, args: str | Sequence[str]) -> None: ...
    def is_valid(self) -> bool: ...
    def threads(self) -> tuple[InferiorThread, ...]: ...
    def architecture(self) -> Architecture: ...
    def read_memory(self, address: _ValueOrInt, length: int) -> memoryview: ...
    def write_memory(self, address: _ValueOrInt, buffer: _BufferType, length: int = ...) -> None: ...
    def search_memory(self, address: _ValueOrInt, length: int, pattern: _BufferType) -> int | None: ...
    def thread_from_handle(self, handle: Value) -> InferiorThread: ...
    @deprecated("Use gdb.thread_from_handle() instead.")
    def thread_from_thread_handle(self, handle: Value) -> InferiorThread: ...
    def set_env(self, name: str, value: str) -> None: ...
    def unset_env(self, name: str) -> None: ...
    def clear_env(self) -> None: ...

# Threads

class Thread(threading.Thread): ...

def selected_thread() -> InferiorThread: ...
@final
class InferiorThread:
    name: str | None
    details: str | None
    num: int
    global_num: int
    ptid: tuple[int, int, int]
    ptid_string: str
    inferior: Inferior

    def is_valid(self) -> bool: ...
    def switch(self) -> None: ...
    def is_stopped(self) -> bool: ...
    def is_running(self) -> bool: ...
    def is_exited(self) -> bool: ...
    def handle(self) -> bytes: ...

# Recordings

def start_recording(method: str = ..., format: str = ..., /) -> Record: ...
def current_recording() -> Record | None: ...
def stop_recording() -> None: ...

class Record:
    method: str
    format: str | None
    begin: Instruction
    end: Instruction
    replay_position: Instruction | None
    instruction_history: list[Instruction]
    function_call_history: list[RecordFunctionSegment]

    def goto(self, instruction: Instruction, /) -> None: ...
    def clear(self) -> None: ...

class Instruction:
    pc: int
    data: memoryview
    decoded: str
    size: int

class RecordInstruction(Instruction):
    number: int
    sal: Symtab_and_line | None
    is_speculative: bool

class RecordGap(Instruction):
    number: int
    error_code: int
    error_string: str

class RecordFunctionSegment:
    number: int
    symbol: Symbol | None
    level: int | None
    instructions: list[RecordInstruction | RecordGap]
    up: RecordFunctionSegment | None
    prev: RecordFunctionSegment | None
    next: RecordFunctionSegment | None

# CLI Commands

@disjoint_base
class Command:
    def __init__(self, name: str, command_class: int, completer_class: int = ..., prefix: bool = ...) -> None: ...
    def dont_repeat(self) -> None: ...
    def invoke(self, argument: str, from_tty: bool) -> None: ...
    def complete(self, text: str, word: str) -> object: ...

def string_to_argv(argv: str, /) -> list[str]: ...

COMMAND_NONE: int
COMMAND_RUNNING: int
COMMAND_DATA: int
COMMAND_STACK: int
COMMAND_FILES: int
COMMAND_SUPPORT: int
COMMAND_STATUS: int
COMMAND_BREAKPOINTS: int
COMMAND_TRACEPOINTS: int
COMMAND_TUI: int
COMMAND_USER: int
COMMAND_OBSCURE: int
COMMAND_MAINTENANCE: int

COMPLETE_NONE: int
COMPLETE_FILENAME: int
COMPLETE_LOCATION: int
COMPLETE_COMMAND: int
COMPLETE_SYMBOL: int
COMPLETE_EXPRESSION: int

# GDB/MI Commands

@disjoint_base
class MICommand:
    name: str
    installed: bool

    def __init__(self, name: str) -> None: ...
    def invoke(self, arguments: list[str]) -> dict[str, object] | None: ...

# Parameters

@disjoint_base
class Parameter:
    set_doc: str
    show_doc: str
    value: object

    def __init__(self, name: str, command_class: int, parameter_class: int, enum_sequence: Sequence[str] = ...) -> None: ...
    def get_set_string(self) -> str: ...
    def get_show_string(self, svalue: str) -> str: ...

PARAM_BOOLEAN: int
PARAM_AUTO_BOOLEAN: int
PARAM_UINTEGER: int
PARAM_INTEGER: int
PARAM_STRING: int
PARAM_STRING_NOESCAPE: int
PARAM_OPTIONAL_FILENAME: int
PARAM_FILENAME: int
PARAM_ZINTEGER: int
PARAM_ZUINTEGER: int
PARAM_ZUINTEGER_UNLIMITED: int
PARAM_ENUM: int

# Convenience functions

class Function:
    def __init__(self, name: str) -> None: ...
    def invoke(self, *args: Value) -> _ValueOrNative: ...

# Progspaces

def current_progspace() -> Progspace | None: ...
def progspaces() -> Sequence[Progspace]: ...
@final
class Progspace:
    executable_filename: str | None
    filename: str | None
    symbol_file: Objfile | None
    pretty_printers: list[_PrettyPrinterLookupFunction]
    type_printers: list[gdb.types._TypePrinter]
    frame_filters: dict[str, _FrameFilter]
    frame_unwinders: list[_Unwinder]
    missing_file_handlers: Sequence[tuple[Literal["debug"], MissingDebugHandler] | tuple[Literal["file"], MissingFileHandler]]

    def block_for_pc(self, pc: int, /) -> Block | None: ...
    def find_pc_line(self, pc: int, /) -> Symtab_and_line: ...
    def is_valid(self) -> bool: ...
    def objfile_for_address(self, address: int, /) -> Objfile | None: ...
    def objfiles(self) -> Sequence[Objfile]: ...
    def solib_name(self, address: int, /) -> str | None: ...

# Objfiles

def current_objfile() -> Objfile | None: ...
def objfiles() -> list[Objfile]: ...
def lookup_objfile(name: str, by_build_id: bool = ...) -> Objfile | None: ...
@final
class Objfile:
    filename: str | None
    username: str | None
    owner: Objfile | None
    build_id: str | None
    progspace: Progspace | None
    pretty_printers: list[_PrettyPrinterLookupFunction]
    type_printers: list[gdb.types._TypePrinter]
    frame_filters: dict[str, _FrameFilter]
    frame_unwinders: list[_Unwinder]
    is_file: bool

    def is_valid(self) -> bool: ...
    def add_separate_debug_file(self, file: str) -> None: ...
    def lookup_global_symbol(self, name: str, domain: int = ...) -> Symbol | None: ...
    def lookup_static_symbol(self, name: str, domain: int = ...) -> Symbol | None: ...

# Frames

def selected_frame() -> Frame: ...
def newest_frame() -> Frame: ...
def frame_stop_reason_string(code: int, /) -> str: ...
def invalidate_cached_frames() -> None: ...

NORMAL_FRAME: int
DUMMY_FRAME: int
INLINE_FRAME: int
TAILCALL_FRAME: int
SIGTRAMP_FRAME: int
ARCH_FRAME: int
SENTINEL_FRAME: int

FRAME_UNWIND_NO_REASON: int
FRAME_UNWIND_NULL_ID: int
FRAME_UNWIND_OUTERMOST: int
FRAME_UNWIND_UNAVAILABLE: int
FRAME_UNWIND_INNER_ID: int
FRAME_UNWIND_SAME_ID: int
FRAME_UNWIND_NO_SAVED_PC: int
FRAME_UNWIND_MEMORY_ERROR: int

@final
class Frame:
    def is_valid(self) -> bool: ...
    def name(self) -> str | None: ...
    def architecture(self) -> Architecture: ...
    def type(self) -> int: ...
    def unwind_stop_reason(self) -> int: ...
    def pc(self) -> int: ...
    def block(self) -> Block: ...
    def function(self) -> Symbol: ...
    def older(self) -> Frame | None: ...
    def newer(self) -> Frame | None: ...
    def find_sal(self) -> Symtab_and_line: ...
    def read_register(self, register: str | RegisterDescriptor | int) -> Value: ...
    def read_var(self, variable: str | Symbol, block: Block | None = ...) -> Value: ...
    def select(self) -> None: ...
    def level(self) -> int: ...
    def static_link(self) -> Frame | None: ...
    def language(self): ...

# Blocks

def block_for_pc(pc: int) -> Block | None: ...
@final
class Block:
    start: int
    end: int
    function: Symbol | None
    superblock: Block | None
    global_block: Block
    static_block: Block | None
    is_global: bool
    is_static: bool

    def is_valid(self) -> bool: ...
    def __iter__(self) -> BlockIterator: ...

@final
class BlockIterator:
    def is_valid(self) -> bool: ...
    def __iter__(self: _typeshed.Self) -> _typeshed.Self: ...
    def __next__(self) -> Symbol: ...

# Symbols

def lookup_symbol(name: str, block: Block | None = ..., domain: int = ...) -> tuple[Symbol | None, bool]: ...
def lookup_global_symbol(name: str, domain: int = ...) -> Symbol | None: ...
def lookup_static_symbol(name: str, domain: int = ...) -> Symbol | None: ...
def lookup_static_symbols(name: str, domain: int = ...) -> list[Symbol]: ...
@final
class Symbol:
    type: Type | None
    symtab: Symtab
    line: int
    name: str
    linkage_name: str
    print_name: str
    addr_class: int
    needs_frame: bool
    is_argument: bool
    is_constant: bool
    is_function: bool
    is_variable: bool
    is_artificial: bool

    def is_valid(self) -> bool: ...
    def value(self, frame: Frame = ..., /) -> Value: ...

SYMBOL_UNDEF_DOMAIN: Final = 0
SYMBOL_VAR_DOMAIN: Final = 1
SYMBOL_STRUCT_DOMAIN: Final = 2
SYMBOL_MODULE_DOMAIN: Final = 3
SYMBOL_LABEL_DOMAIN: Final = 4
SYMBOL_COMMON_BLOCK_DOMAIN: Final = 5
SYMBOL_TYPE_DOMAIN: Final = 6
SYMBOL_FUNCTION_DOMAIN: Final = 7

SYMBOL_LOC_UNDEF: Final = 0
SYMBOL_LOC_CONST: Final = 1
SYMBOL_LOC_STATIC: Final = 2
SYMBOL_LOC_REGISTER: Final = 3
SYMBOL_LOC_ARG: Final = 4
SYMBOL_LOC_REF_ARG: Final = 5
SYMBOL_LOC_REGPARM_ADDR: Final = 6
SYMBOL_LOC_LOCAL: Final = 7
SYMBOL_LOC_TYPEDEF: Final = 8
SYMBOL_LOC_LABEL: Final = 9
SYMBOL_LOC_BLOCK: Final = 10
SYMBOL_LOC_CONST_BYTES: Final = 11
SYMBOL_LOC_UNRESOLVED: Final = 12
SYMBOL_LOC_OPTIMIZED_OUT: Final = 13
SYMBOL_LOC_COMPUTED: Final = 14
SYMBOL_LOC_COMMON_BLOCK: Final = 15

# Symbol tables

@final
class Symtab_and_line:
    symtab: Symtab
    pc: int
    last: int
    line: int

    def is_valid(self) -> bool: ...

@final
class Symtab:
    filename: str
    objfile: Objfile
    producer: str

    def is_valid(self) -> bool: ...
    def fullname(self) -> str: ...
    def global_block(self) -> Block: ...
    def static_block(self) -> Block: ...
    def linetable(self) -> LineTable: ...

# Line Tables

@final
class LineTableEntry:
    line: int
    pc: int

@final
class LineTableIterator(Iterator[LineTableEntry]):
    def __next__(self) -> LineTableEntry: ...
    def is_valid(self) -> bool: ...

@final
class LineTable:
    def __iter__(self) -> LineTableIterator: ...
    def line(self, line: int, /) -> tuple[LineTableEntry, ...]: ...
    def has_line(self, line: int, /) -> bool: ...
    def source_lines(self) -> list[int]: ...
    def is_valid(self) -> bool: ...

# Breakpoints

@disjoint_base
class Breakpoint:
    # The where="spec" form of __init__().  See py-breakpoints.c:bppy_init():keywords for the positional order.
    @overload
    def __init__(
        self,
        # where
        spec: str = ...,
        # options
        type: int = ...,
        wp_class: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...

    # The where="location" form of __init__().  A watchpoint (`type=BP_WATCHPOINT`) cannot be created with this form.
    #
    # We exclude the `wp_class` (watchpoint class) option here, even though py-breakpoints.c accepts it.  It doesn't make sense
    # unless type==BP_WATCHPOINT, and is silently ignored in those cases; allowing it in those cases is likely an oversight, not
    # an intentional allowance.
    #
    # We repeat this 7 times because the type system doesn't have simple a way for us to say "at least one of `function`, `label`,
    # or `line`", so we must repeat it for each combination of the 3.
    #
    # The len=3 combination.
    @overload
    def __init__(
        self,
        *,
        # where
        source: str = ...,
        function: str,
        label: str,
        line: int | str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...
    # The 3 len=2 combinations.
    @overload
    def __init__(
        self,
        *,
        source: str = ...,
        # where
        label: str,
        line: int | str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        source: str = ...,
        # where
        function: str,
        line: int | str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        source: str = ...,
        # where
        function: str,
        label: str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...
    # The 3 len=1 combinations.
    @overload
    def __init__(
        self,
        *,
        source: str = ...,
        # where
        function: str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        source: str = ...,
        # where
        label: str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...
    @overload
    def __init__(
        self,
        *,
        source: str = ...,
        # where
        line: int | str,
        # options
        type: int = ...,
        internal: bool = ...,
        temporary: bool = ...,
        qualified: bool = ...,
    ) -> None: ...

    # Methods.
    def stop(self) -> bool: ...
    def is_valid(self) -> bool: ...
    def delete(self) -> None: ...

    enabled: bool
    silent: bool
    pending: bool
    thread: int | None
    task: str | None
    ignore_count: int
    number: int
    type: int
    visible: bool
    temporary: bool
    hit_count: int
    location: str | None
    locations: Sequence[BreakpointLocation]
    inferior: int | None
    expression: str | None
    condition: str | None
    commands: str | None

@final
class BreakpointLocation:
    address: int
    enabled: bool
    fullname: str
    function: str | None
    owner: Breakpoint
    source: tuple[str, int]
    thread_groups: Sequence[int]

BP_NONE: int
BP_BREAKPOINT: int
BP_HARDWARE_BREAKPOINT: int
BP_WATCHPOINT: int
BP_HARDWARE_WATCHPOINT: int
BP_READ_WATCHPOINT: int
BP_ACCESS_WATCHPOINT: int
BP_CATCHPOINT: int

WP_READ: int
WP_WRITE: int
WP_ACCESS: int

# Finish Breakpoints

@disjoint_base
class FinishBreakpoint(Breakpoint):
    return_value: Value | None

    def __init__(self, frame: Frame = ..., internal: bool = ...) -> None: ...
    def out_of_scope(self) -> None: ...

# Lazy strings

class LazyString:
    def value(self) -> Value: ...

    address: Value
    length: int
    encoding: str
    type: Type

# Architectures

@type_check_only
class _Instruction(TypedDict):
    addr: int
    asm: str
    length: int

@final
class Architecture:
    def name(self) -> str: ...
    def disassemble(self, start_pc: int, end_pc: int = ..., count: int = ...) -> list[_Instruction]: ...
    def integer_type(self, size: int, signed: bool = ...) -> Type: ...
    def registers(self, reggroup: str = ...) -> RegisterDescriptorIterator: ...
    def register_groups(self) -> RegisterGroupsIterator: ...

# Registers

@final
class RegisterDescriptor:
    name: str

@final
class RegisterDescriptorIterator(Iterator[RegisterDescriptor]):
    def __next__(self) -> RegisterDescriptor: ...
    def find(self, name: str) -> RegisterDescriptor | None: ...

@final
class RegisterGroup:
    name: str

@final
class RegisterGroupsIterator(Iterator[RegisterGroup]):
    def __next__(self) -> RegisterGroup: ...

# Connections

@disjoint_base
class TargetConnection:
    def is_valid(self) -> bool: ...

    num: int
    type: str
    description: str
    details: str | None

@final
class RemoteTargetConnection(TargetConnection):
    def send_packet(self, packet: str | bytes) -> bytes: ...

# TUI Windows

def register_window_type(name: str, factory: Callable[[TuiWindow], _Window]) -> None: ...

class TuiWindow:
    width: int
    height: int
    title: str

    def is_valid(self) -> bool: ...
    def erase(self) -> None: ...
    def write(self, string: str, full_window: bool = ...) -> None: ...

@type_check_only
class _Window(Protocol):
    def close(self) -> None: ...
    def render(self) -> None: ...
    def hscroll(self, num: int) -> None: ...
    def vscroll(self, num: int) -> None: ...
    def click(self, x: int, y: int, button: int) -> None: ...

# Events
@disjoint_base
class Event: ...

class ThreadEvent(Event):
    inferior_thread: InferiorThread

class ContinueEvent(ThreadEvent): ...

class ExitedEvent(Event):
    exit_code: int
    inferior: Inferior

class ThreadExitedEvent(ThreadEvent): ...

class StopEvent(ThreadEvent):
    details: dict[str, object]

class BreakpointEvent(StopEvent):
    breakpoints: Sequence[Breakpoint]
    breakpoint: Breakpoint

class NewObjFileEvent(Event):
    new_objfile: Objfile

class FreeObjFileEvent(Event):
    objfile: Objfile

class ClearObjFilesEvent(Event):
    progspace: Progspace

class NewProgspaceEvent(Event): ...
class FreeProgspaceEvent(Event): ...

class SignalEvent(StopEvent):
    stop_signal: str

@type_check_only
class _InferiorCallEvent(Event): ...

class InferiorCallPreEvent(_InferiorCallEvent):
    ptid: InferiorThread
    address: Value

class InferiorCallPostEvent(_InferiorCallEvent):
    ptid: InferiorThread
    address: Value

class MemoryChangedEvent(Event):
    address: Value
    length: int

class RegisterChangedEvent(Event):
    frame: Frame
    regnum: str

class NewInferiorEvent(Event):
    inferior: Inferior

class InferiorDeletedEvent(Event):
    inferior: Inferior

class NewThreadEvent(ThreadEvent): ...

class GdbExitingEvent(Event):
    exit_code: int

class ConnectionEvent(Event):
    connection: TargetConnection

class ExecutableChangedEvent(Event): ...

class TuiEnabledEvent(Event):
    enabled: bool

_ET = TypeVar("_ET", bound=Event | Breakpoint | None)

@final
class EventRegistry(Generic[_ET]):
    def connect(self, object: Callable[[_ET], object], /) -> None: ...
    def disconnect(self, object: Callable[[_ET], object], /) -> None: ...

class ValuePrinter: ...

def blocked_signals(): ...
def notify_mi(name: str, data: dict[str, object] | None = None): ...
def interrupt(): ...
def execute_mi(command: str, *args: str) -> dict[str, object]: ...
