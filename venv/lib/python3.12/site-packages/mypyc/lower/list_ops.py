from __future__ import annotations

from mypyc.common import PLATFORM_SIZE
from mypyc.ir.ops import GetElementPtr, IncRef, Integer, IntOp, LoadMem, SetMem, Value
from mypyc.ir.rtypes import (
    PyListObject,
    c_pyssize_t_rprimitive,
    object_rprimitive,
    pointer_rprimitive,
)
from mypyc.irbuild.ll_builder import LowLevelIRBuilder
from mypyc.lower.registry import lower_primitive_op


@lower_primitive_op("buf_init_item")
def buf_init_item(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    """Initialize an item in a buffer of "PyObject *" values at given index.

    This can be used to initialize the data buffer of a freshly allocated list
    object.
    """
    base = args[0]
    index_value = args[1]
    value = args[2]
    assert isinstance(index_value, Integer)
    index = index_value.numeric_value()
    if index == 0:
        ptr = base
    else:
        ptr = builder.add(
            IntOp(
                pointer_rprimitive,
                base,
                Integer(index * PLATFORM_SIZE, c_pyssize_t_rprimitive),
                IntOp.ADD,
                line,
            )
        )
    return builder.add(SetMem(object_rprimitive, ptr, value, line))


@lower_primitive_op("list_items")
def list_items(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    ob_item_ptr = builder.add(GetElementPtr(args[0], PyListObject, "ob_item", line))
    return builder.add(LoadMem(pointer_rprimitive, ob_item_ptr, line))


def list_item_ptr(builder: LowLevelIRBuilder, obj: Value, index: Value, line: int) -> Value:
    """Get a pointer to a list item (index must be valid and non-negative).

    Type of index must be c_pyssize_t_rprimitive, and obj must refer to a list object.
    """
    # List items are represented as an array of pointers. Pointer to the item obj[index] is
    # <pointer to first item> + index * <pointer size>.
    items = list_items(builder, [obj], line)
    delta = builder.add(
        IntOp(
            c_pyssize_t_rprimitive,
            index,
            Integer(PLATFORM_SIZE, c_pyssize_t_rprimitive),
            IntOp.MUL,
        )
    )
    return builder.add(IntOp(pointer_rprimitive, items, delta, IntOp.ADD))


@lower_primitive_op("list_get_item_unsafe")
def list_get_item_unsafe(builder: LowLevelIRBuilder, args: list[Value], line: int) -> Value:
    index = builder.coerce(args[1], c_pyssize_t_rprimitive, line)
    item_ptr = list_item_ptr(builder, args[0], index, line)
    value = builder.add(LoadMem(object_rprimitive, item_ptr, line))
    builder.add(IncRef(value))
    return value
