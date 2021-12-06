"""Primitive str ops."""

from typing import List, Tuple

from mypyc.ir.ops import ERR_MAGIC, ERR_NEVER
from mypyc.ir.rtypes import (
    RType, object_rprimitive, str_rprimitive, int_rprimitive, list_rprimitive,
    c_int_rprimitive, pointer_rprimitive, bool_rprimitive, bit_rprimitive
)
from mypyc.primitives.registry import (
    method_op, binary_op, function_op,
    load_address_op, custom_op
)


# Get the 'str' type object.
load_address_op(
    name='builtins.str',
    type=object_rprimitive,
    src='PyUnicode_Type')

# str(obj)
function_op(
    name='builtins.str',
    arg_types=[object_rprimitive],
    return_type=str_rprimitive,
    c_function_name='PyObject_Str',
    error_kind=ERR_MAGIC)

# str1 + str2
binary_op(name='+',
          arg_types=[str_rprimitive, str_rprimitive],
          return_type=str_rprimitive,
          c_function_name='PyUnicode_Concat',
          error_kind=ERR_MAGIC)

# str.join(obj)
method_op(
    name='join',
    arg_types=[str_rprimitive, object_rprimitive],
    return_type=str_rprimitive,
    c_function_name='PyUnicode_Join',
    error_kind=ERR_MAGIC
)

# str.startswith(str)
method_op(
    name='startswith',
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=bool_rprimitive,
    c_function_name='CPyStr_Startswith',
    error_kind=ERR_NEVER
)

# str.endswith(str)
method_op(
    name='endswith',
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=bool_rprimitive,
    c_function_name='CPyStr_Endswith',
    error_kind=ERR_NEVER
)

# str[index] (for an int index)
method_op(
    name='__getitem__',
    arg_types=[str_rprimitive, int_rprimitive],
    return_type=str_rprimitive,
    c_function_name='CPyStr_GetItem',
    error_kind=ERR_MAGIC
)

# str.split(...)
str_split_types = [str_rprimitive, str_rprimitive, int_rprimitive]  # type: List[RType]
str_split_functions = ['PyUnicode_Split', 'PyUnicode_Split', 'CPyStr_Split']
str_split_constants = [[(0, pointer_rprimitive), (-1, c_int_rprimitive)],
                       [(-1, c_int_rprimitive)],
                       []] \
                       # type: List[List[Tuple[int, RType]]]
for i in range(len(str_split_types)):
    method_op(
        name='split',
        arg_types=str_split_types[0:i+1],
        return_type=list_rprimitive,
        c_function_name=str_split_functions[i],
        extra_int_constants=str_split_constants[i],
        error_kind=ERR_MAGIC)

# str1 += str2
#
# PyUnicodeAppend makes an effort to reuse the LHS when the refcount
# is 1. This is super dodgy but oh well, the interpreter does it.
binary_op(name='+=',
          arg_types=[str_rprimitive, str_rprimitive],
          return_type=str_rprimitive,
          c_function_name='CPyStr_Append',
          error_kind=ERR_MAGIC,
          steals=[True, False])

unicode_compare = custom_op(
    arg_types=[str_rprimitive, str_rprimitive],
    return_type=c_int_rprimitive,
    c_function_name='PyUnicode_Compare',
    error_kind=ERR_NEVER)

# str[begin:end]
str_slice_op = custom_op(
    arg_types=[str_rprimitive, int_rprimitive, int_rprimitive],
    return_type=object_rprimitive,
    c_function_name='CPyStr_GetSlice',
    error_kind=ERR_MAGIC)

# str.replace(old, new)
method_op(
    name='replace',
    arg_types=[str_rprimitive, str_rprimitive, str_rprimitive],
    return_type=str_rprimitive,
    c_function_name="PyUnicode_Replace",
    error_kind=ERR_MAGIC,
    extra_int_constants=[(-1, c_int_rprimitive)])

# str.replace(old, new, count)
method_op(
    name='replace',
    arg_types=[str_rprimitive, str_rprimitive, str_rprimitive, int_rprimitive],
    return_type=str_rprimitive,
    c_function_name="CPyStr_Replace",
    error_kind=ERR_MAGIC)

# check if a string is true (isn't an empty string)
str_check_if_true = custom_op(
    arg_types=[str_rprimitive],
    return_type=bit_rprimitive,
    c_function_name="CPyStr_IsTrue",
    error_kind=ERR_NEVER,
)
