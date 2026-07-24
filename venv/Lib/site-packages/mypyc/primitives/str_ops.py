"""Primitive str ops."""

from __future__ import annotations

from mypyc.ir.deps import STR_EXTRA_OPS
from mypyc.ir.ops import ERR_MAGIC, ERR_NEVER
from mypyc.ir.rtypes import (
    RType,
    bit_rprimitive,
    bool_rprimitive,
    bytes_rprimitive,
    c_int_rprimitive,
    c_pyssize_t_rprimitive,
    int64_rprimitive,
    int_rprimitive,
    list_rprimitive,
    object_rprimitive,
    pointer_rprimitive,
    short_int_rprimitive,
    str_rprimitive,
    tuple_rprimitive,
)
from mypyc.primitives.registry import (
    ERR_NEG_INT,
    binary_op,
    custom_op,
    custom_primitive_op,
    function_op,
    load_address_op,
    method_op,
)

# Get the 'str' type object.
load_address_op(name="builtins.str", type=object_rprimitive, src="PyUnicode_Type")

# str(obj)
str_op = function_op(
    name="builtins.str",
    arg_types=[object_rprimitive],
    return_type=str_rprimitive,
    c_function_name="PyObject_Str",
    error_kind=ERR_MAGIC,
)

# repr(obj)
function_op(
    name="builtins.repr",
    arg_types=[object_rprimitive],
    return_type=str_rprimitive,
    c_function_name="PyObject_Repr",
    error_kind=ERR_MAGIC,
)

# translate isinstance(obj, str)
isinstance_str = function_op(
    name="builtins.isinstance",
    arg_types=[object_rprimitive],
    return_type=bit_rprimitive,
    c_function_name="PyUnicode_Check",
    error_kind=ERR_NEVER,
)

# str1 + str2
binary_op(
    name="+",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="PyUnicode_Concat",
    error_kind=ERR_MAGIC,
)

# str1 += str2
#
# PyUnicode_Append makes an effort to reuse the LHS when the refcount
# is 1. This is super dodgy but oh well, the interpreter does it.
binary_op(
    name="+=",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Append",
    error_kind=ERR_MAGIC,
    steals=[True, False],
)

# str * int
binary_op(
    name="*",
    arg_types=[str_rprimitive, int_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Multiply",
    error_kind=ERR_MAGIC,
)

# int * str
binary_op(
    name="*",
    arg_types=[int_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Multiply",
    error_kind=ERR_MAGIC,
    ordering=[1, 0],
)

# str1 == str2 (very common operation, so we provide our own)
str_eq = custom_primitive_op(
    name="str_eq",
    c_function_name="CPyStr_Equal",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=bool_rprimitive,
    error_kind=ERR_NEVER,
)

str_eq_literal = custom_primitive_op(
    name="str_eq_literal",
    c_function_name="CPyStr_EqualLiteral",
    arg_types=[str_rprimitive, str_rprimitive, c_pyssize_t_rprimitive],
    return_type=bool_rprimitive,
    error_kind=ERR_NEVER,
)

unicode_compare = custom_op(
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="PyUnicode_Compare",
    error_kind=ERR_NEVER,
)

# str[index] (for an int index)
method_op(
    name="__getitem__",
    arg_types=[str_rprimitive, int_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_GetItem",
    error_kind=ERR_MAGIC,
)

# This is unsafe since it assumes that the index is within reasonable bounds.
# In the future this might do no bounds checking at all.
str_get_item_unsafe_op = custom_op(
    arg_types=[str_rprimitive, c_pyssize_t_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_GetItemUnsafe",
    error_kind=ERR_MAGIC,
)

# str[begin:end]
str_slice_op = custom_op(
    arg_types=[str_rprimitive, int_rprimitive, int_rprimitive],
    return_type=object_rprimitive,
    c_function_name="CPyStr_GetSlice",
    error_kind=ERR_MAGIC,
)

# item in str
binary_op(
    name="in",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="PyUnicode_Contains",
    error_kind=ERR_NEG_INT,
    truncated_type=bool_rprimitive,
    ordering=[1, 0],
)

# str.find(...) and str.rfind(...)
str_find_types: list[RType] = [str_rprimitive, str_rprimitive, int_rprimitive, int_rprimitive]
str_find_functions = ["CPyStr_Find", "CPyStr_Find", "CPyStr_FindWithEnd"]
str_find_constants: list[list[tuple[int, RType]]] = [[(0, c_int_rprimitive)], [], []]
str_rfind_constants: list[list[tuple[int, RType]]] = [[(0, c_int_rprimitive)], [], []]
for i in range(len(str_find_types) - 1):
    method_op(
        name="find",
        arg_types=str_find_types[0 : i + 2],
        return_type=int_rprimitive,
        c_function_name=str_find_functions[i],
        extra_int_constants=str_find_constants[i] + [(1, c_int_rprimitive)],
        error_kind=ERR_MAGIC,
    )
    method_op(
        name="rfind",
        arg_types=str_find_types[0 : i + 2],
        return_type=int_rprimitive,
        c_function_name=str_find_functions[i],
        extra_int_constants=str_rfind_constants[i] + [(-1, c_int_rprimitive)],
        error_kind=ERR_MAGIC,
    )

# str.join(obj)
method_op(
    name="join",
    arg_types=[str_rprimitive, object_rprimitive],
    return_type=str_rprimitive,
    c_function_name="PyUnicode_Join",
    error_kind=ERR_MAGIC,
)

str_build_op = custom_op(
    arg_types=[c_pyssize_t_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Build",
    error_kind=ERR_MAGIC,
    var_arg_type=str_rprimitive,
)

# str.strip, str.lstrip, str.rstrip
for strip_prefix in ["l", "r", ""]:
    method_op(
        name=f"{strip_prefix}strip",
        arg_types=[str_rprimitive, str_rprimitive],
        return_type=str_rprimitive,
        c_function_name=f"CPyStr_{strip_prefix.upper()}Strip",
        error_kind=ERR_NEVER,
    )
    method_op(
        name=f"{strip_prefix}strip",
        arg_types=[str_rprimitive],
        return_type=str_rprimitive,
        c_function_name=f"CPyStr_{strip_prefix.upper()}Strip",
        # This 0 below is implicitly treated as NULL in C.
        extra_int_constants=[(0, c_int_rprimitive)],
        error_kind=ERR_NEVER,
    )

# str.startswith(str)
method_op(
    name="startswith",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="CPyStr_Startswith",
    truncated_type=bool_rprimitive,
    error_kind=ERR_NEVER,
)

# str.startswith(tuple)
method_op(
    name="startswith",
    arg_types=[str_rprimitive, tuple_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyStr_Startswith",
    error_kind=ERR_MAGIC,
)

# str.endswith(str)
method_op(
    name="endswith",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="CPyStr_Endswith",
    truncated_type=bool_rprimitive,
    error_kind=ERR_NEVER,
)

# str.endswith(tuple)
method_op(
    name="endswith",
    arg_types=[str_rprimitive, tuple_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyStr_Endswith",
    error_kind=ERR_MAGIC,
)

# str.removeprefix(str)
method_op(
    name="removeprefix",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Removeprefix",
    error_kind=ERR_NEVER,
)

# str.removesuffix(str)
method_op(
    name="removesuffix",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Removesuffix",
    error_kind=ERR_NEVER,
)

# str.split(...) and str.rsplit(...)
str_split_types: list[RType] = [str_rprimitive, str_rprimitive, int_rprimitive]
str_split_functions = ["PyUnicode_Split", "PyUnicode_Split", "CPyStr_Split"]
str_rsplit_functions = ["PyUnicode_RSplit", "PyUnicode_RSplit", "CPyStr_RSplit"]
str_split_constants: list[list[tuple[int, RType]]] = [
    [(0, pointer_rprimitive), (-1, c_int_rprimitive)],
    [(-1, c_int_rprimitive)],
    [],
]
for i in range(len(str_split_types)):
    method_op(
        name="split",
        arg_types=str_split_types[0 : i + 1],
        return_type=list_rprimitive,
        c_function_name=str_split_functions[i],
        extra_int_constants=str_split_constants[i],
        error_kind=ERR_MAGIC,
    )
    method_op(
        name="rsplit",
        arg_types=str_split_types[0 : i + 1],
        return_type=list_rprimitive,
        c_function_name=str_rsplit_functions[i],
        extra_int_constants=str_split_constants[i],
        error_kind=ERR_MAGIC,
    )

# str.splitlines(...)
str_splitlines_types: list[RType] = [str_rprimitive, bool_rprimitive]
str_splitlines_constants: list[list[tuple[int, RType]]] = [[(0, c_int_rprimitive)], []]
for i in range(2):
    method_op(
        name="splitlines",
        arg_types=str_splitlines_types[0 : i + 1],
        return_type=list_rprimitive,
        c_function_name="PyUnicode_Splitlines",
        extra_int_constants=str_splitlines_constants[i],
        error_kind=ERR_NEVER,
    )

# str.partition(str)
method_op(
    name="partition",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=tuple_rprimitive,
    c_function_name="PyUnicode_Partition",
    error_kind=ERR_MAGIC,
)

# str.rpartition(str)
method_op(
    name="rpartition",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=tuple_rprimitive,
    c_function_name="PyUnicode_RPartition",
    error_kind=ERR_MAGIC,
)

# str.count(substring)
method_op(
    name="count",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=c_pyssize_t_rprimitive,
    c_function_name="CPyStr_Count",
    error_kind=ERR_NEG_INT,
    extra_int_constants=[(0, c_pyssize_t_rprimitive)],
)

# str.count(substring, start)
method_op(
    name="count",
    arg_types=[str_rprimitive, str_rprimitive, int_rprimitive],
    return_type=c_pyssize_t_rprimitive,
    c_function_name="CPyStr_Count",
    error_kind=ERR_NEG_INT,
)

# str.count(substring, start, end)
method_op(
    name="count",
    arg_types=[str_rprimitive, str_rprimitive, int_rprimitive, int_rprimitive],
    return_type=c_pyssize_t_rprimitive,
    c_function_name="CPyStr_CountFull",
    error_kind=ERR_NEG_INT,
)

# str.replace(old, new)
method_op(
    name="replace",
    arg_types=[str_rprimitive, str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="PyUnicode_Replace",
    error_kind=ERR_MAGIC,
    extra_int_constants=[(-1, c_int_rprimitive)],
)

# str.replace(old, new, count)
method_op(
    name="replace",
    arg_types=[str_rprimitive, str_rprimitive, str_rprimitive, int_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Replace",
    error_kind=ERR_MAGIC,
)

# check if a string is true (isn't an empty string)
str_check_if_true = custom_op(
    arg_types=[str_rprimitive],
    return_type=bit_rprimitive,
    c_function_name="CPyStr_IsTrue",
    error_kind=ERR_NEVER,
)

str_ssize_t_size_op = custom_op(
    arg_types=[str_rprimitive],
    return_type=c_pyssize_t_rprimitive,
    c_function_name="CPyStr_Size_size_t",
    error_kind=ERR_NEG_INT,
)

# str.lower()
method_op(
    name="lower",
    arg_types=[str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Lower",
    error_kind=ERR_MAGIC,
)

# str.upper()
method_op(
    name="upper",
    arg_types=[str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Upper",
    error_kind=ERR_MAGIC,
)

method_op(
    name="isspace",
    arg_types=[str_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyStr_IsSpace",
    error_kind=ERR_NEVER,
)

method_op(
    name="isalnum",
    arg_types=[str_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyStr_IsAlnum",
    error_kind=ERR_NEVER,
)

method_op(
    name="isdigit",
    arg_types=[str_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyStr_IsDigit",
    error_kind=ERR_NEVER,
)


# obj.decode()
method_op(
    name="decode",
    arg_types=[bytes_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPy_Decode",
    error_kind=ERR_MAGIC,
    extra_int_constants=[(0, pointer_rprimitive), (0, pointer_rprimitive)],
)

# obj.decode(encoding)
method_op(
    name="decode",
    arg_types=[bytes_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPy_Decode",
    error_kind=ERR_MAGIC,
    extra_int_constants=[(0, pointer_rprimitive)],
)

# bytes.decode(encoding, errors)
method_op(
    name="decode",
    arg_types=[bytes_rprimitive, str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPy_Decode",
    error_kind=ERR_MAGIC,
)

# bytes.decode(encoding) - utf8 strict specialization
bytes_decode_utf8_strict = custom_op(
    arg_types=[bytes_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPy_DecodeUTF8",
    error_kind=ERR_MAGIC,
)

# bytes.decode(encoding) - ascii strict specialization
bytes_decode_ascii_strict = custom_op(
    arg_types=[bytes_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPy_DecodeASCII",
    error_kind=ERR_MAGIC,
)

# bytes.decode(encoding) - latin1 strict specialization
bytes_decode_latin1_strict = custom_op(
    arg_types=[bytes_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPy_DecodeLatin1",
    error_kind=ERR_MAGIC,
)

# str.encode()
method_op(
    name="encode",
    arg_types=[str_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPy_Encode",
    error_kind=ERR_MAGIC,
    extra_int_constants=[(0, pointer_rprimitive), (0, pointer_rprimitive)],
)

# str.encode(encoding)
method_op(
    name="encode",
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPy_Encode",
    error_kind=ERR_MAGIC,
    extra_int_constants=[(0, pointer_rprimitive)],
)

# str.encode(encoding) - utf8 strict specialization
str_encode_utf8_strict = custom_op(
    arg_types=[str_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="PyUnicode_AsUTF8String",
    error_kind=ERR_MAGIC,
)

# str.encode(encoding) - ascii strict specialization
str_encode_ascii_strict = custom_op(
    arg_types=[str_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="PyUnicode_AsASCIIString",
    error_kind=ERR_MAGIC,
)

# str.encode(encoding) - latin1 strict specialization
str_encode_latin1_strict = custom_op(
    arg_types=[str_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="PyUnicode_AsLatin1String",
    error_kind=ERR_MAGIC,
)

# str.encode(encoding, errors)
method_op(
    name="encode",
    arg_types=[str_rprimitive, str_rprimitive, str_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPy_Encode",
    error_kind=ERR_MAGIC,
)

function_op(
    name="builtins.ord",
    arg_types=[str_rprimitive],
    return_type=int_rprimitive,
    c_function_name="CPyStr_Ord",
    error_kind=ERR_MAGIC,
)

# Optimized str indexing for ord(s[i])

# str index adjustment - convert negative index to positive
str_adjust_index_op = custom_primitive_op(
    name="str_adjust_index",
    arg_types=[str_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="CPyStr_AdjustIndex",
    error_kind=ERR_NEVER,
    dependencies=[STR_EXTRA_OPS],
)

# str range check - check if index is in valid range
str_range_check_op = custom_primitive_op(
    name="str_range_check",
    arg_types=[str_rprimitive, int64_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyStr_RangeCheck",
    error_kind=ERR_NEVER,
    dependencies=[STR_EXTRA_OPS],
)

# str.__getitem__() as int - get character at index as int (ord value) - no bounds checking
str_get_item_unsafe_as_int_op = custom_primitive_op(
    name="str_get_item_unsafe_as_int",
    arg_types=[str_rprimitive, int64_rprimitive],
    return_type=short_int_rprimitive,
    c_function_name="CPyStr_GetItemUnsafeAsInt",
    error_kind=ERR_NEVER,
    dependencies=[STR_EXTRA_OPS],
)
