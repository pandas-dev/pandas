"""Primitive bytes ops."""

from __future__ import annotations

from mypyc.ir.deps import BYTES_EXTRA_OPS
from mypyc.ir.ops import ERR_MAGIC, ERR_NEVER
from mypyc.ir.rtypes import (
    RUnion,
    bit_rprimitive,
    bool_rprimitive,
    bytes_rprimitive,
    c_int_rprimitive,
    c_pyssize_t_rprimitive,
    dict_rprimitive,
    int64_rprimitive,
    int_rprimitive,
    list_rprimitive,
    object_rprimitive,
    str_rprimitive,
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

# Get the 'bytes' type object.
load_address_op(name="builtins.bytes", type=object_rprimitive, src="PyBytes_Type")

# bytes(obj)
function_op(
    name="builtins.bytes",
    arg_types=[RUnion([list_rprimitive, dict_rprimitive, str_rprimitive])],
    return_type=bytes_rprimitive,
    c_function_name="PyBytes_FromObject",
    error_kind=ERR_MAGIC,
)

# translate isinstance(obj, bytes)
isinstance_bytes = function_op(
    name="builtins.isinstance",
    arg_types=[object_rprimitive],
    return_type=bit_rprimitive,
    c_function_name="PyBytes_Check",
    error_kind=ERR_NEVER,
)

# bytes ==/!= (return -1/0/1)
bytes_compare = custom_op(
    arg_types=[bytes_rprimitive, bytes_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="CPyBytes_Compare",
    error_kind=ERR_NEG_INT,
)

# bytes + bytes
binary_op(
    name="+",
    arg_types=[bytes_rprimitive, bytes_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_Concat",
    error_kind=ERR_MAGIC,
    steals=[True, False],
)

# bytes * int
binary_op(
    name="*",
    arg_types=[bytes_rprimitive, int_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_Multiply",
    error_kind=ERR_MAGIC,
)

# int * bytes
binary_op(
    name="*",
    arg_types=[int_rprimitive, bytes_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_Multiply",
    error_kind=ERR_MAGIC,
    ordering=[1, 0],
)

# bytes[begin:end]
bytes_slice_op = custom_op(
    arg_types=[bytes_rprimitive, int_rprimitive, int_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_GetSlice",
    error_kind=ERR_MAGIC,
)

# bytes[index]
# bytearray[index]
method_op(
    name="__getitem__",
    arg_types=[bytes_rprimitive, int_rprimitive],
    return_type=int_rprimitive,
    c_function_name="CPyBytes_GetItem",
    error_kind=ERR_MAGIC,
)

# bytes.join(obj)
method_op(
    name="join",
    arg_types=[bytes_rprimitive, object_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_Join",
    error_kind=ERR_MAGIC,
)

# bytes.translate(table)
method_op(
    name="translate",
    arg_types=[bytes_rprimitive, object_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_Translate",
    error_kind=ERR_MAGIC,
    dependencies=[BYTES_EXTRA_OPS],
)

# bytes.startswith(bytes)
method_op(
    name="startswith",
    arg_types=[bytes_rprimitive, bytes_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="CPyBytes_Startswith",
    truncated_type=bool_rprimitive,
    error_kind=ERR_MAGIC,
)

# bytes.endswith(bytes)
method_op(
    name="endswith",
    arg_types=[bytes_rprimitive, bytes_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name="CPyBytes_Endswith",
    truncated_type=bool_rprimitive,
    error_kind=ERR_MAGIC,
)

# Join bytes objects and return a new bytes.
# The first argument is the total number of the following bytes.
bytes_build_op = custom_op(
    arg_types=[c_pyssize_t_rprimitive],
    return_type=bytes_rprimitive,
    c_function_name="CPyBytes_Build",
    error_kind=ERR_MAGIC,
    var_arg_type=bytes_rprimitive,
)

function_op(
    name="builtins.ord",
    arg_types=[bytes_rprimitive],
    return_type=int_rprimitive,
    c_function_name="CPyBytes_Ord",
    error_kind=ERR_MAGIC,
)

# Optimized bytes.__getitem__ operations

# bytes index adjustment - convert negative index to positive
bytes_adjust_index_op = custom_primitive_op(
    name="bytes_adjust_index",
    arg_types=[bytes_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="CPyBytes_AdjustIndex",
    error_kind=ERR_NEVER,
    dependencies=[BYTES_EXTRA_OPS],
)

# bytes range check - check if index is in valid range
bytes_range_check_op = custom_primitive_op(
    name="bytes_range_check",
    arg_types=[bytes_rprimitive, int64_rprimitive],
    return_type=bool_rprimitive,
    c_function_name="CPyBytes_RangeCheck",
    error_kind=ERR_NEVER,
    dependencies=[BYTES_EXTRA_OPS],
)

# bytes.__getitem__() - get byte at index (no bounds checking)
bytes_get_item_unsafe_op = custom_primitive_op(
    name="bytes_get_item_unsafe",
    arg_types=[bytes_rprimitive, int64_rprimitive],
    return_type=int_rprimitive,
    c_function_name="CPyBytes_GetItemUnsafe",
    error_kind=ERR_NEVER,
    dependencies=[BYTES_EXTRA_OPS],
)
