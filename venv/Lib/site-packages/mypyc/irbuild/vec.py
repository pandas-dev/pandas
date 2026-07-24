"""Generate IR for librt.vecs.vec operations"""

from __future__ import annotations

from typing import TYPE_CHECKING, Final

from mypyc.common import IS_32_BIT_PLATFORM, PLATFORM_SIZE
from mypyc.ir.ops import (
    ERR_MAGIC,
    Assign,
    BasicBlock,
    Branch,
    CallC,
    ComparisonOp,
    DecRef,
    GetElement,
    Integer,
    IntOp,
    RaiseStandardError,
    Register,
    SetElement,
    TupleGet,
    TupleSet,
    Unborrow,
    Undef,
    Unreachable,
    Value,
)
from mypyc.ir.rtypes import (
    RInstance,
    RPrimitive,
    RTuple,
    RType,
    RUnion,
    RVec,
    VecNestedBufItem,
    bool_rprimitive,
    c_pyssize_t_rprimitive,
    c_size_t_rprimitive,
    int32_rprimitive,
    int64_rprimitive,
    is_c_py_ssize_t_rprimitive,
    is_int32_rprimitive,
    is_int64_rprimitive,
    is_int_rprimitive,
    is_short_int_rprimitive,
    list_rprimitive,
    object_rprimitive,
    optional_value_type,
    pointer_rprimitive,
    tuple_rprimitive,
    vec_api_by_item_type,
    vec_item_type_tags,
)

if TYPE_CHECKING:
    from mypyc.irbuild.ll_builder import LowLevelIRBuilder


def as_platform_int(builder: LowLevelIRBuilder, v: Value, line: int) -> Value:
    rtype = v.type
    if is_c_py_ssize_t_rprimitive(rtype):
        return v
    if isinstance(v, Integer):
        if is_short_int_rprimitive(rtype) or is_int_rprimitive(rtype):
            return Integer(v.value // 2, c_pyssize_t_rprimitive)
        return Integer(v.value, c_pyssize_t_rprimitive)
    if isinstance(rtype, RPrimitive):
        if PLATFORM_SIZE == 8 and is_int64_rprimitive(rtype):
            return v
        if PLATFORM_SIZE == 4 and is_int32_rprimitive(rtype):
            return v
    return builder.coerce(v, c_pyssize_t_rprimitive, line)


def vec_create(
    builder: LowLevelIRBuilder,
    vtype: RVec,
    length: int | Value,
    line: int,
    *,
    capacity: Value | None = None,
) -> Value:
    if isinstance(length, int):
        length = Integer(length, c_pyssize_t_rprimitive)
    length = as_platform_int(builder, length, line)
    if capacity is not None:
        capacity = as_platform_int(builder, capacity, line)
    else:
        capacity = length

    item_type = vtype.item_type
    api_name = vec_api_by_item_type.get(item_type)
    if api_name is not None:
        call = CallC(
            f"{api_name}.alloc",
            [length, capacity],
            vtype,
            False,
            False,
            error_kind=ERR_MAGIC,
            line=line,
        )
        return builder.add(call)

    typeobj, optional, depth = vec_item_type_info(builder, item_type, line)
    if typeobj is not None:
        typeval: Value
        if isinstance(typeobj, Integer):
            typeval = typeobj
        else:
            # Create an integer which will hold the type object * as an integral value.
            # Assign implicitly coerces between pointer/integer types.
            typeval = Register(pointer_rprimitive)
            builder.add(Assign(typeval, typeobj))
            if optional:
                typeval = builder.add(
                    IntOp(pointer_rprimitive, typeval, Integer(1, pointer_rprimitive), IntOp.OR)
                )
        if depth == 0:
            call = CallC(
                "VecTApi.alloc",
                [length, capacity, typeval],
                vtype,
                False,
                False,
                error_kind=ERR_MAGIC,
                line=line,
            )
            return builder.add(call)
        else:
            call = CallC(
                "VecNestedApi.alloc",
                [length, capacity, typeval, Integer(depth, int32_rprimitive)],
                vtype,
                False,
                False,
                error_kind=ERR_MAGIC,
                line=line,
            )
            return builder.add(call)

    assert False, "unsupported: %s" % vtype


def vec_create_initialized(
    builder: LowLevelIRBuilder,
    vtype: RVec,
    length: int | Value,
    init: Value,
    line: int,
    *,
    capacity: Value | None = None,
) -> Value:
    """Create vec with items initialized to the given value."""
    if isinstance(length, int):
        length = Integer(length, c_pyssize_t_rprimitive)
    length = as_platform_int(builder, length, line)

    item_type = vtype.item_type
    init = builder.coerce(init, item_type, line)
    vec = vec_create(builder, vtype, length, line, capacity=capacity)

    items_start = vec_items(builder, vec)
    step = step_size(item_type)
    items_end = builder.int_add(items_start, builder.int_mul(length, step))

    for_loop = builder.begin_for(
        items_start, items_end, Integer(step, c_pyssize_t_rprimitive), signed=False
    )
    vec_set_mem_item(builder, for_loop.index, item_type, init)
    for_loop.finish()

    builder.keep_alive([vec], line)
    return vec


def vec_create_from_values(
    builder: LowLevelIRBuilder,
    vtype: RVec,
    values: list[Value],
    line: int,
    *,
    capacity: Value | None = None,
) -> Value:
    vec = vec_create(builder, vtype, len(values), line, capacity=capacity)
    ptr = vec_items(builder, vec)
    item_type = vtype.item_type
    step = step_size(item_type)
    for value in values:
        vec_set_mem_item(builder, ptr, item_type, value)
        ptr = builder.int_add(ptr, step)
    builder.keep_alive([vec], line)
    return vec


def step_size(item_type: RType) -> int:
    if isinstance(item_type, RPrimitive):
        return item_type.size
    elif isinstance(item_type, RVec):
        return PLATFORM_SIZE * 2
    else:
        return PLATFORM_SIZE


VEC_TYPE_INFO_I64: Final = 2


def vec_item_type_info(
    builder: LowLevelIRBuilder, typ: RType, line: int
) -> tuple[Value | None, bool, int]:
    if isinstance(typ, RPrimitive) and typ.is_refcounted:
        return builder.load_builtin(typ.name, line), False, 0
    elif isinstance(typ, RInstance):
        return builder.load_native_type_object(typ.name), False, 0
    elif typ in vec_item_type_tags:
        return Integer(vec_item_type_tags[typ], c_size_t_rprimitive), False, 0
    elif isinstance(typ, RUnion):
        non_opt = optional_value_type(typ)
        assert non_opt is not None
        typeval, _, _ = vec_item_type_info(builder, non_opt, line)
        if typeval is not None:
            return typeval, True, 0
    elif isinstance(typ, RVec):
        typeval, optional, depth = vec_item_type_info(builder, typ.item_type, line)
        if typeval is not None:
            return typeval, optional, depth + 1
    return None, False, 0


def vec_len(builder: LowLevelIRBuilder, val: Value) -> Value:
    """Return len(<vec>) as i64."""
    len_val = vec_len_native(builder, val)
    if IS_32_BIT_PLATFORM:
        return builder.coerce(len_val, int64_rprimitive, -1)
    return len_val


def vec_len_native(builder: LowLevelIRBuilder, val: Value) -> Value:
    """Return len(<vec>) as platform integer type (32-bit/64-bit)."""
    return builder.get_element(val, "len")


def vec_items(builder: LowLevelIRBuilder, vecobj: Value) -> Value:
    """Return pointer to first item in vec.

    The items field points directly to the first element in the buffer.
    """
    return builder.get_element(vecobj, "items")


def vec_item_ptr(builder: LowLevelIRBuilder, vecobj: Value, index: Value) -> Value:
    items_addr = vec_items(builder, vecobj)
    assert isinstance(vecobj.type, RVec)
    # TODO: Do we need to care about alignment?
    item_type = vecobj.type.item_type
    if isinstance(item_type, RPrimitive):
        item_size = item_type.size
    elif isinstance(item_type, RVec):
        item_size = 2 * PLATFORM_SIZE
    else:
        item_size = object_rprimitive.size
    delta = builder.int_mul(index, item_size)
    return builder.int_add(items_addr, delta)


def vec_load_mem_item(
    builder: LowLevelIRBuilder, ptr: Value, item_type: RType, *, can_borrow: bool = False
) -> Value:
    """Load a vec item from storage, converting nested vec slots to RVec values."""
    return builder.load_mem(ptr, item_type, borrow=can_borrow)


def vec_set_mem_item(
    builder: LowLevelIRBuilder, ptr: Value, item_type: RType, item: Value
) -> None:
    """Store a vec item, converting RVec values to nested storage items."""
    builder.set_mem(ptr, item_type, item)


def vec_check_and_adjust_index(
    builder: LowLevelIRBuilder, lenv: Value, index: Value, line: int
) -> Value:
    r = Register(int64_rprimitive)
    index = builder.coerce(index, int64_rprimitive, line)
    lenv = builder.coerce(lenv, int64_rprimitive, line)
    ok, ok2, ok3 = BasicBlock(), BasicBlock(), BasicBlock()
    fail, fail2 = BasicBlock(), BasicBlock()
    is_less = builder.comparison_op(index, lenv, ComparisonOp.ULT, line)
    builder.add_bool_branch(is_less, ok2, fail)
    builder.activate_block(fail)

    x = builder.int_add(index, lenv)
    is_less2 = builder.comparison_op(x, lenv, ComparisonOp.ULT, line)
    builder.add_bool_branch(is_less2, ok, fail2)

    builder.activate_block(fail2)
    # TODO: Include index in exception
    builder.add(RaiseStandardError(RaiseStandardError.INDEX_ERROR, None, line))
    builder.add(Unreachable())

    builder.activate_block(ok)
    builder.assign(r, x)
    builder.goto(ok3)

    builder.activate_block(ok2)
    builder.assign(r, index)
    builder.goto(ok3)

    builder.activate_block(ok3)
    return r


def vec_get_item(
    builder: LowLevelIRBuilder, base: Value, index: Value, line: int, *, can_borrow: bool = False
) -> Value:
    """Generate inlined vec __getitem__ call.

    We inline this, since it's simple but performance-critical.
    """
    # TODO: Support more item types
    # TODO: Support more index types
    len_val = vec_len(builder, base)
    index = vec_check_and_adjust_index(builder, len_val, index, line)
    return vec_get_item_unsafe(builder, base, index, line, can_borrow=can_borrow)


def vec_get_item_unsafe(
    builder: LowLevelIRBuilder, base: Value, index: Value, line: int, *, can_borrow: bool = False
) -> Value:
    """Get vec item, assuming index is non-negative and within bounds."""
    assert isinstance(base.type, RVec)
    index = as_platform_int(builder, index, line)
    vtype = base.type
    item_addr = vec_item_ptr(builder, base, index)
    result = vec_load_mem_item(builder, item_addr, vtype.item_type, can_borrow=can_borrow)
    builder.keep_alives.append(base)
    return result


def vec_set_item(
    builder: LowLevelIRBuilder, base: Value, index: Value, item: Value, line: int
) -> None:
    assert isinstance(base.type, RVec)
    vtype = base.type
    len_val = vec_len(builder, base)
    index = vec_check_and_adjust_index(builder, len_val, index, line)
    index = builder.coerce(index, c_pyssize_t_rprimitive, line)
    item_addr = vec_item_ptr(builder, base, index)
    item_type = vtype.item_type
    item = builder.coerce(item, item_type, line)
    if item_type.is_refcounted:
        # Read an unborrowed reference to cause a decref to be
        # generated for the old item.
        old_item = vec_load_mem_item(builder, item_addr, item_type, can_borrow=True)
        builder.add(DecRef(old_item))
    vec_set_mem_item(builder, item_addr, item_type, item)
    builder.keep_alive([base], line)


def vec_init_item_unsafe(
    builder: LowLevelIRBuilder, base: Value, index: Value, item: Value, line: int
) -> None:
    assert isinstance(base.type, RVec)
    index = as_platform_int(builder, index, line)
    vtype = base.type
    item_addr = vec_item_ptr(builder, base, index)
    item_type = vtype.item_type
    item = builder.coerce(item, item_type, line)
    vec_set_mem_item(builder, item_addr, item_type, item)
    builder.keep_alive([base], line)


def convert_to_t_ext_item(builder: LowLevelIRBuilder, item: Value) -> Value:
    vec_len = builder.add(GetElement(item, "len"))
    vec_items = builder.add(GetElement(item, "items"))
    temp = builder.add(SetElement(Undef(VecNestedBufItem), "len", vec_len))
    return builder.add(SetElement(temp, "items", vec_items))


def convert_from_t_ext_item(builder: LowLevelIRBuilder, item: Value, vec_type: RVec) -> Value:
    """Convert an owned VecNestedBufItem to the corresponding RVec value."""
    vec_len = builder.add(GetElement(item, "len"))
    vec_items = builder.add(GetElement(item, "items"))
    temp = builder.add(SetElement(Undef(vec_type), "len", vec_len))
    return builder.add(SetElement(temp, "items", vec_items))


def vec_item_type(builder: LowLevelIRBuilder, item_type: RType, line: int) -> Value:
    typeobj, optional, depth = vec_item_type_info(builder, item_type, line)
    assert typeobj is not None
    if isinstance(typeobj, Integer):
        return typeobj
    else:
        # Create an integer which will hold the type object * as an integral value.
        # Assign implicitly coerces between pointer/integer types.
        typeval: Value
        typeval = Register(pointer_rprimitive)
        builder.add(Assign(typeval, typeobj))
        if optional:
            typeval = builder.add(
                IntOp(pointer_rprimitive, typeval, Integer(1, pointer_rprimitive), IntOp.OR)
            )
        return typeval


def vec_append(builder: LowLevelIRBuilder, vec: Value, item: Value, line: int) -> Value:
    vec_type = vec.type
    assert isinstance(vec_type, RVec)
    item_type = vec_type.item_type
    coerced_item = builder.coerce(item, item_type, line)
    item_type_arg = []
    api_name = vec_api_by_item_type.get(item_type)
    if api_name is not None:
        name = f"{api_name}.append"
    elif vec_type.depth() == 0:
        name = "VecTApi.append"
        item_type_arg = [vec_item_type(builder, item_type, line)]
    else:
        coerced_item = convert_to_t_ext_item(builder, coerced_item)
        name = "VecNestedApi.append"
    call = builder.add(
        CallC(
            name,
            [vec, coerced_item] + item_type_arg,
            vec_type,
            steals=[True, False] + ([False] if item_type_arg else []),
            is_borrowed=False,
            error_kind=ERR_MAGIC,
            line=line,
        )
    )
    if vec_type.depth() > 0:
        builder.keep_alive([item], line)
    return call


def vec_extend(builder: LowLevelIRBuilder, vec: Value, iterable: Value, line: int) -> Value:
    vec_type = vec.type
    assert isinstance(vec_type, RVec)
    item_type = vec_type.item_type
    if isinstance(iterable.type, RVec) and iterable.type == vec_type:
        suffix = "_vec"
        src = iterable
    else:
        suffix = ""
        src = builder.coerce(iterable, object_rprimitive, line)
    item_type_arg: list[Value] = []
    api_name = vec_api_by_item_type.get(item_type)
    if api_name is not None:
        name = f"{api_name}.extend{suffix}"
    elif vec_type.depth() == 0:
        name = f"VecTApi.extend{suffix}"
        item_type_arg = [vec_item_type(builder, item_type, line)]
    else:
        name = f"VecNestedApi.extend{suffix}"
    return builder.add(
        CallC(
            name,
            [vec, src] + item_type_arg,
            vec_type,
            steals=[True, False] + ([False] if item_type_arg else []),
            is_borrowed=False,
            error_kind=ERR_MAGIC,
            line=line,
        )
    )


def vec_pop(builder: LowLevelIRBuilder, base: Value, index: Value, line: int) -> Value:
    assert isinstance(base.type, RVec)
    vec_type = base.type
    item_type = vec_type.item_type
    index = as_platform_int(builder, index, line)

    api_name = vec_api_by_item_type.get(item_type)
    if api_name is not None:
        name = f"{api_name}.pop"
    elif vec_type.depth() == 0:
        name = "VecTApi.pop"
    else:
        name = "VecNestedApi.pop"
        # Nested vecs return a generic vec struct.
        item_type = VecNestedBufItem
    result = builder.add(
        CallC(
            name,
            [base, index],
            RTuple([vec_type, item_type]),
            steals=[True, False],
            is_borrowed=False,
            error_kind=ERR_MAGIC,
            line=line,
        )
    )
    if vec_type.depth() > 0:
        orig = result
        x = builder.add(TupleGet(result, 0, borrow=True))
        x = builder.add(Unborrow(x))
        y = builder.add(TupleGet(result, 1, borrow=True))
        y = builder.add(Unborrow(y))
        assert isinstance(vec_type.item_type, RVec)
        z = convert_from_t_ext_item(builder, y, vec_type.item_type)
        result = builder.add(TupleSet([x, z], line))
        builder.keep_alive([orig], line, steal=True)
    return result


def vec_remove(builder: LowLevelIRBuilder, vec: Value, item: Value, line: int) -> Value:
    assert isinstance(vec.type, RVec)
    vec_type = vec.type
    item_type = vec_type.item_type
    coerced_item = builder.coerce(item, item_type, line)

    if item_type in vec_api_by_item_type:
        name = f"{vec_api_by_item_type[item_type]}.remove"
    elif vec_type.depth() == 0:
        name = "VecTApi.remove"
    else:
        coerced_item = convert_to_t_ext_item(builder, coerced_item)
        name = "VecNestedApi.remove"
    call = builder.add(
        CallC(
            name,
            [vec, coerced_item],
            vec_type,
            steals=[True, False],
            is_borrowed=False,
            error_kind=ERR_MAGIC,
            line=line,
        )
    )
    if vec_type.depth() > 0:
        builder.keep_alive([item], line)
    return call


def vec_contains(builder: LowLevelIRBuilder, vec: Value, target: Value, line: int) -> Value:
    assert isinstance(vec.type, RVec)
    vec_type = vec.type
    item_type = vec_type.item_type
    target = builder.coerce(target, item_type, line)

    step = step_size(item_type)
    len_val = vec_len_native(builder, vec)
    items_start = vec_items(builder, vec)
    items_end = builder.int_add(items_start, builder.int_mul(len_val, step))

    true, end = BasicBlock(), BasicBlock()

    for_loop = builder.begin_for(
        items_start, items_end, Integer(step, c_pyssize_t_rprimitive), signed=False
    )
    item = vec_load_mem_item(builder, for_loop.index, item_type, can_borrow=True)
    comp = builder.binary_op(item, target, "==", line)
    false = BasicBlock()
    builder.add(Branch(comp, true, false, Branch.BOOL))
    builder.activate_block(false)
    for_loop.finish()

    builder.keep_alive([vec], line)

    res = Register(bool_rprimitive)
    builder.assign(res, Integer(0, bool_rprimitive))
    builder.goto(end)
    builder.activate_block(true)
    builder.assign(res, Integer(1, bool_rprimitive))
    builder.goto_and_activate(end)
    return res


def vec_slice(
    builder: LowLevelIRBuilder, vec: Value, begin: Value, end: Value, line: int
) -> Value:
    assert isinstance(vec.type, RVec)
    vec_type = vec.type
    item_type = vec_type.item_type
    begin = builder.coerce(begin, int64_rprimitive, line)
    end = builder.coerce(end, int64_rprimitive, line)
    api_name = vec_api_by_item_type.get(item_type)
    if api_name is not None:
        name = f"{api_name}.slice"
    elif vec_type.depth() == 0:
        name = "VecTApi.slice"
    else:
        name = "VecNestedApi.slice"
    call = CallC(
        name,
        [vec, begin, end],
        vec_type,
        steals=[False, False, False],
        is_borrowed=False,
        error_kind=ERR_MAGIC,
        line=line,
    )
    return builder.add(call)


def vec_to_list(builder: LowLevelIRBuilder, vec: Value, line: int) -> Value | None:
    return _vec_to_sequence(builder, vec, line, "to_list", list_rprimitive)


def vec_to_tuple(builder: LowLevelIRBuilder, vec: Value, line: int) -> Value | None:
    return _vec_to_sequence(builder, vec, line, "to_tuple", tuple_rprimitive)


def supports_vec_to_sequence(vec_type: RVec) -> bool:
    return vec_api_by_item_type.get(vec_type.item_type) is not None or vec_type.depth() == 0


def _vec_to_sequence(
    builder: LowLevelIRBuilder, vec: Value, line: int, method: str, result_type: RType
) -> Value | None:
    vec_type = vec.type
    assert isinstance(vec_type, RVec)
    item_type = vec_type.item_type
    api_name = vec_api_by_item_type.get(item_type)
    if api_name is not None:
        name = f"{api_name}.{method}"
    elif supports_vec_to_sequence(vec_type):
        name = f"VecTApi.{method}"
    else:
        return None
    return builder.add(
        CallC(
            name,
            [vec],
            result_type,
            steals=[True],
            is_borrowed=False,
            error_kind=ERR_MAGIC,
            line=line,
        )
    )
