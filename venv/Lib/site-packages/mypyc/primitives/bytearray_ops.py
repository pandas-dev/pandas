"""Primitive bytearray ops.

NOTE: Most of these should be added to bytearray_extra_ops.c, which requires the
      BYTEARRAY_EXTRA_OPS primitive dependency, since these are used relatively rarely and we
      don't want to compile them unless needed.
"""

from __future__ import annotations

from mypyc.ir.deps import BYTEARRAY_EXTRA_OPS
from mypyc.ir.ops import ERR_MAGIC, ERR_NEVER
from mypyc.ir.rtypes import bit_rprimitive, bytearray_rprimitive, object_rprimitive
from mypyc.primitives.registry import custom_primitive_op, function_op, load_address_op

# Get the 'bytearray' type object.
load_address_op(name="builtins.bytearray", type=object_rprimitive, src="PyByteArray_Type")

# bytearray(obj)
function_op(
    name="builtins.bytearray",
    arg_types=[object_rprimitive],
    return_type=bytearray_rprimitive,
    c_function_name="PyByteArray_FromObject",
    error_kind=ERR_MAGIC,
)

# bytearray() -- construct empty bytearray
function_op(
    name="builtins.bytearray",
    arg_types=[],
    return_type=bytearray_rprimitive,
    c_function_name="CPyByteArray_New",
    error_kind=ERR_MAGIC,
    dependencies=[BYTEARRAY_EXTRA_OPS],
)

# isinstance(obj, bytearray)
isinstance_bytearray = custom_primitive_op(
    name="builtins.isinstance",
    arg_types=[object_rprimitive],
    return_type=bit_rprimitive,
    c_function_name="PyByteArray_Check",
    error_kind=ERR_NEVER,
)
