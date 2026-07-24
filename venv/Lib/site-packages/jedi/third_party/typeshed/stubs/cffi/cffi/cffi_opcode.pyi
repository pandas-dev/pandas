from typing import Final

class CffiOp:
    op: int | None
    arg: str | None
    def __init__(self, op: int | None, arg: str | None) -> None: ...
    def as_c_expr(self) -> str: ...
    def as_python_bytes(self) -> str: ...

def format_four_bytes(num: int) -> str: ...

OP_PRIMITIVE: Final = 1
OP_POINTER: Final = 3
OP_ARRAY: Final = 5
OP_OPEN_ARRAY: Final = 7
OP_STRUCT_UNION: Final = 9
OP_ENUM: Final = 11
OP_FUNCTION: Final = 13
OP_FUNCTION_END: Final = 15
OP_NOOP: Final = 17
OP_BITFIELD: Final = 19
OP_TYPENAME: Final = 21
OP_CPYTHON_BLTN_V: Final = 23
OP_CPYTHON_BLTN_N: Final = 25
OP_CPYTHON_BLTN_O: Final = 27
OP_CONSTANT: Final = 29
OP_CONSTANT_INT: Final = 31
OP_GLOBAL_VAR: Final = 33
OP_DLOPEN_FUNC: Final = 35
OP_DLOPEN_CONST: Final = 37
OP_GLOBAL_VAR_F: Final = 39
OP_EXTERN_PYTHON: Final = 41
PRIM_VOID: Final = 0
PRIM_BOOL: Final = 1
PRIM_CHAR: Final = 2
PRIM_SCHAR: Final = 3
PRIM_UCHAR: Final = 4
PRIM_SHORT: Final = 5
PRIM_USHORT: Final = 6
PRIM_INT: Final = 7
PRIM_UINT: Final = 8
PRIM_LONG: Final = 9
PRIM_ULONG: Final = 10
PRIM_LONGLONG: Final = 11
PRIM_ULONGLONG: Final = 12
PRIM_FLOAT: Final = 13
PRIM_DOUBLE: Final = 14
PRIM_LONGDOUBLE: Final = 15
PRIM_WCHAR: Final = 16
PRIM_INT8: Final = 17
PRIM_UINT8: Final = 18
PRIM_INT16: Final = 19
PRIM_UINT16: Final = 20
PRIM_INT32: Final = 21
PRIM_UINT32: Final = 22
PRIM_INT64: Final = 23
PRIM_UINT64: Final = 24
PRIM_INTPTR: Final = 25
PRIM_UINTPTR: Final = 26
PRIM_PTRDIFF: Final = 27
PRIM_SIZE: Final = 28
PRIM_SSIZE: Final = 29
PRIM_INT_LEAST8: Final = 30
PRIM_UINT_LEAST8: Final = 31
PRIM_INT_LEAST16: Final = 32
PRIM_UINT_LEAST16: Final = 33
PRIM_INT_LEAST32: Final = 34
PRIM_UINT_LEAST32: Final = 35
PRIM_INT_LEAST64: Final = 36
PRIM_UINT_LEAST64: Final = 37
PRIM_INT_FAST8: Final = 38
PRIM_UINT_FAST8: Final = 39
PRIM_INT_FAST16: Final = 40
PRIM_UINT_FAST16: Final = 41
PRIM_INT_FAST32: Final = 42
PRIM_UINT_FAST32: Final = 43
PRIM_INT_FAST64: Final = 44
PRIM_UINT_FAST64: Final = 45
PRIM_INTMAX: Final = 46
PRIM_UINTMAX: Final = 47
PRIM_FLOATCOMPLEX: Final = 48
PRIM_DOUBLECOMPLEX: Final = 49
PRIM_CHAR16: Final = 50
PRIM_CHAR32: Final = 51
PRIMITIVE_TO_INDEX: Final[dict[str, int]]
F_UNION: Final = 1
F_CHECK_FIELDS: Final = 2
F_PACKED: Final = 4
F_EXTERNAL: Final = 8
F_OPAQUE: Final = 16
G_FLAGS: Final[dict[bytes, bytes]]
CLASS_NAME: Final[dict[int, str]]
