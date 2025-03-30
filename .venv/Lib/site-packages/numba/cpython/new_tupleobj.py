"""
Implementation of tuple objects
"""

import operator

from numba.core.imputils import (lower_builtin, lower_getattr_generic,
                                    lower_cast, lower_constant, iternext_impl,
                                    impl_ret_borrowed, impl_ret_untracked,
                                    RefType)
from numba.core import typing, types, cgutils
from numba.core.extending import overload_method, overload, intrinsic


@lower_builtin(types.NamedTupleClass, types.VarArg(types.Any))
def namedtuple_constructor(context, builder, sig, args):
    # A namedtuple has the same representation as a regular tuple
    # the arguments need casting (lower_cast) from the types in the ctor args
    # to those in the ctor return type, this is to handle cases such as a
    # literal present in the args, but a type present in the return type.
    newargs = []
    for i, arg in enumerate(args):
        casted = context.cast(builder, arg, sig.args[i], sig.return_type[i])
        newargs.append(casted)
    res = context.make_tuple(builder, sig.return_type, tuple(newargs))
    # The tuple's contents are borrowed
    return impl_ret_borrowed(context, builder, sig.return_type, res)

@lower_builtin(operator.add, types.BaseTuple, types.BaseTuple)
def tuple_add(context, builder, sig, args):
    left, right = [cgutils.unpack_tuple(builder, x) for x in args]
    res = context.make_tuple(builder, sig.return_type, left + right)
    # The tuple's contents are borrowed
    return impl_ret_borrowed(context, builder, sig.return_type, res)

def tuple_cmp_ordered(context, builder, op, sig, args):
    tu, tv = sig.args
    u, v = args
    res = cgutils.alloca_once_value(builder, cgutils.true_bit)
    bbend = builder.append_basic_block("cmp_end")
    for i, (ta, tb) in enumerate(zip(tu.types, tv.types)):
        a = builder.extract_value(u, i)
        b = builder.extract_value(v, i)
        not_equal = context.generic_compare(builder, operator.ne, (ta, tb), (a, b))
        with builder.if_then(not_equal):
            pred = context.generic_compare(builder, op, (ta, tb), (a, b))
            builder.store(pred, res)
            builder.branch(bbend)
    # Everything matched equal => compare lengths
    len_compare = op(len(tu.types), len(tv.types))
    pred = context.get_constant(types.boolean, len_compare)
    builder.store(pred, res)
    builder.branch(bbend)
    builder.position_at_end(bbend)
    return builder.load(res)


@lower_builtin(operator.eq, types.BaseTuple, types.BaseTuple)
def tuple_eq(context, builder, sig, args):
    tu, tv = sig.args
    u, v = args
    if len(tu.types) != len(tv.types):
        res = context.get_constant(types.boolean, False)
        return impl_ret_untracked(context, builder, sig.return_type, res)
    res = context.get_constant(types.boolean, True)
    for i, (ta, tb) in enumerate(zip(tu.types, tv.types)):
        a = builder.extract_value(u, i)
        b = builder.extract_value(v, i)
        pred = context.generic_compare(builder, operator.eq, (ta, tb), (a, b))
        res = builder.and_(res, pred)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(operator.ne, types.BaseTuple, types.BaseTuple)
def tuple_ne(context, builder, sig, args):
    res = builder.not_(tuple_eq(context, builder, sig, args))
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(operator.lt, types.BaseTuple, types.BaseTuple)
def tuple_lt(context, builder, sig, args):
    res = tuple_cmp_ordered(context, builder, operator.lt, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(operator.le, types.BaseTuple, types.BaseTuple)
def tuple_le(context, builder, sig, args):
    res = tuple_cmp_ordered(context, builder, operator.le, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(operator.gt, types.BaseTuple, types.BaseTuple)
def tuple_gt(context, builder, sig, args):
    res = tuple_cmp_ordered(context, builder, operator.gt, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)

@lower_builtin(operator.ge, types.BaseTuple, types.BaseTuple)
def tuple_ge(context, builder, sig, args):
    res = tuple_cmp_ordered(context, builder, operator.ge, sig, args)
    return impl_ret_untracked(context, builder, sig.return_type, res)

# for hashing see hashing.py

@lower_getattr_generic(types.BaseNamedTuple)
def namedtuple_getattr(context, builder, typ, value, attr):
    """
    Fetch a namedtuple's field.
    """
    index = typ.fields.index(attr)
    res = builder.extract_value(value, index)
    return impl_ret_borrowed(context, builder, typ[index], res)


@lower_constant(types.UniTuple)
@lower_constant(types.NamedUniTuple)
def unituple_constant(context, builder, ty, pyval):
    """
    Create a homogeneous tuple constant.
    """
    consts = [context.get_constant_generic(builder, ty.dtype, v)
              for v in pyval]
    return impl_ret_borrowed(
        context, builder, ty, cgutils.pack_array(builder, consts),
    )

@lower_constant(types.Tuple)
@lower_constant(types.NamedTuple)
def unituple_constant(context, builder, ty, pyval):
    """
    Create a heterogeneous tuple constant.
    """
    consts = [context.get_constant_generic(builder, ty.types[i], v)
              for i, v in enumerate(pyval)]
    return impl_ret_borrowed(
        context, builder, ty, cgutils.pack_struct(builder, consts),
    )


#------------------------------------------------------------------------------
# Tuple iterators

@lower_builtin('getiter', types.UniTuple)
@lower_builtin('getiter', types.NamedUniTuple)
def getiter_unituple(context, builder, sig, args):
    [tupty] = sig.args
    [tup] = args

    iterval = context.make_helper(builder, types.UniTupleIter(tupty))

    index0 = context.get_constant(types.intp, 0)
    indexptr = cgutils.alloca_once(builder, index0.type)
    builder.store(index0, indexptr)

    iterval.index = indexptr
    iterval.tuple = tup

    res = iterval._getvalue()
    return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin('iternext', types.UniTupleIter)
@iternext_impl(RefType.BORROWED)
def iternext_unituple(context, builder, sig, args, result):
    [tupiterty] = sig.args
    [tupiter] = args

    iterval = context.make_helper(builder, tupiterty, value=tupiter)

    tup = iterval.tuple
    idxptr = iterval.index
    idx = builder.load(idxptr)
    count = context.get_constant(types.intp, tupiterty.container.count)

    is_valid = builder.icmp_signed('<', idx, count)
    result.set_valid(is_valid)

    with builder.if_then(is_valid):
        getitem_sig = typing.signature(tupiterty.container.dtype,
                                       tupiterty.container,
                                       types.intp)
        getitem_out = getitem_unituple(context, builder, getitem_sig,
                                       [tup, idx])
        # As a iternext_impl function, this will incref the yieled value.
        # We need to release the new reference from getitem_unituple.
        if context.enable_nrt:
            context.nrt.decref(builder, tupiterty.container.dtype, getitem_out)
        result.yield_(getitem_out)
        nidx = builder.add(idx, context.get_constant(types.intp, 1))
        builder.store(nidx, iterval.index)


@overload(operator.getitem)
def getitem_literal_idx(tup, idx):
    """
    Overloads BaseTuple getitem to cover cases where constant
    inference and RewriteConstGetitems cannot replace it
    with a static_getitem.
    """
    if not (isinstance(tup, types.BaseTuple)
            and isinstance(idx, types.IntegerLiteral)):
        return None

    idx_val = idx.literal_value
    def getitem_literal_idx_impl(tup, idx):
        return tup[idx_val]

    return getitem_literal_idx_impl


@lower_builtin('typed_getitem', types.BaseTuple, types.Any)
def getitem_typed(context, builder, sig, args):
    tupty, _ = sig.args
    tup, idx = args
    errmsg_oob = ("tuple index out of range",)

    if len(tupty) == 0:
        # Empty tuple.

        # Always branch and raise IndexError
        with builder.if_then(cgutils.true_bit):
            context.call_conv.return_user_exc(builder, IndexError,
                                              errmsg_oob)
        # This is unreachable in runtime,
        # but it exists to not terminate the current basicblock.
        res = context.get_constant_null(sig.return_type)
        return impl_ret_untracked(context, builder,
                                  sig.return_type, res)
    else:
        # The tuple is not empty

        bbelse = builder.append_basic_block("typed_switch.else")
        bbend = builder.append_basic_block("typed_switch.end")
        switch = builder.switch(idx, bbelse)

        with builder.goto_block(bbelse):
            context.call_conv.return_user_exc(builder, IndexError,
                                            errmsg_oob)

        lrtty = context.get_value_type(sig.return_type)
        voidptrty = context.get_value_type(types.voidptr)
        with builder.goto_block(bbend):
            phinode = builder.phi(voidptrty)

        for i in range(tupty.count):
            ki = context.get_constant(types.intp, i)
            bbi = builder.append_basic_block("typed_switch.%d" % i)
            switch.add_case(ki, bbi)
            # handle negative indexing, create case (-tuple.count + i) to
            # reference same block as i
            kin = context.get_constant(types.intp, -tupty.count + i)
            switch.add_case(kin, bbi)
            with builder.goto_block(bbi):
                value = builder.extract_value(tup, i)
                # Dragon warning...
                # The fact the code has made it this far suggests that type
                # inference decided whatever was being done with the item pulled
                # from the tuple was legitimate, it is not the job of lowering
                # to argue about that. However, here lies a problem, the tuple
                # lowering is implemented as a switch table with each case
                # writing to a phi node slot that is returned. The type of this
                # phi node slot needs to be "correct" for the current type but
                # it also needs to survive stores being made to it from the
                # other cases that will in effect never run. To do this a stack
                # slot is made for each case for the specific type and then cast
                # to a void pointer type, this is then added as an incoming on
                # the phi node, at the end of the switch the phi node is then
                # cast back to the required return type for this typed_getitem.
                # The only further complication is that if the value is not a
                # pointer then the void* juggle won't work so a cast is made
                # prior to store, again, that type inference has permitted it
                # suggests this is safe.
                # End Dragon warning...
                DOCAST = context.typing_context.unify_types(sig.args[0][i],
                                        sig.return_type) == sig.return_type
                if DOCAST:
                    value_slot = builder.alloca(lrtty,
                                                name="TYPED_VALUE_SLOT%s" % i)
                    casted = context.cast(builder, value, sig.args[0][i],
                                        sig.return_type)
                    builder.store(casted, value_slot)
                else:
                    value_slot = builder.alloca(value.type,
                                                name="TYPED_VALUE_SLOT%s" % i)
                    builder.store(value, value_slot)
                phinode.add_incoming(builder.bitcast(value_slot, voidptrty),
                                     bbi)
                builder.branch(bbend)

        builder.position_at_end(bbend)
        res = builder.bitcast(phinode, lrtty.as_pointer())
        res = builder.load(res)
        return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin(operator.getitem, types.UniTuple, types.py_int)
@lower_builtin(operator.getitem, types.UniTuple, types.np_intp)
@lower_builtin(operator.getitem, types.UniTuple, types.np_uintp)
@lower_builtin(operator.getitem, types.NamedUniTuple, types.py_int)
@lower_builtin(operator.getitem, types.NamedUniTuple, types.np_intp)
@lower_builtin(operator.getitem, types.NamedUniTuple, types.np_uintp)
def getitem_unituple(context, builder, sig, args):
    tupty, _ = sig.args
    tup, idx = args

    errmsg_oob = ("tuple index out of range",)

    if len(tupty) == 0:
        # Empty tuple.

        # Always branch and raise IndexError
        with builder.if_then(cgutils.true_bit):
            context.call_conv.return_user_exc(builder, IndexError,
                                              errmsg_oob)
        # This is unreachable in runtime,
        # but it exists to not terminate the current basicblock.
        res = context.get_constant_null(sig.return_type)
        return impl_ret_untracked(context, builder,
                                  sig.return_type, res)
    else:
        # The tuple is not empty
        bbelse = builder.append_basic_block("switch.else")
        bbend = builder.append_basic_block("switch.end")
        switch = builder.switch(idx, bbelse)

        with builder.goto_block(bbelse):
            context.call_conv.return_user_exc(builder, IndexError,
                                              errmsg_oob)

        lrtty = context.get_value_type(tupty.dtype)
        with builder.goto_block(bbend):
            phinode = builder.phi(lrtty)

        for i in range(tupty.count):
            ki = context.get_constant(types.intp, i)
            bbi = builder.append_basic_block("switch.%d" % i)
            switch.add_case(ki, bbi)
            # handle negative indexing, create case (-tuple.count + i) to
            # reference same block as i
            kin = context.get_constant(types.intp, -tupty.count + i)
            switch.add_case(kin, bbi)
            with builder.goto_block(bbi):
                value = builder.extract_value(tup, i)
                builder.branch(bbend)
                phinode.add_incoming(value, bbi)

        builder.position_at_end(bbend)
        res = phinode
        assert sig.return_type == tupty.dtype
        return impl_ret_borrowed(context, builder, sig.return_type, res)


@lower_builtin('static_getitem', types.LiteralStrKeyDict, types.StringLiteral)
@lower_builtin('static_getitem', types.LiteralList, types.IntegerLiteral)
@lower_builtin('static_getitem', types.LiteralList, types.SliceLiteral)
@lower_builtin('static_getitem', types.BaseTuple, types.IntegerLiteral)
@lower_builtin('static_getitem', types.BaseTuple, types.SliceLiteral)
def static_getitem_tuple(context, builder, sig, args):
    tupty, idxty = sig.args
    tup, idx = args
    if isinstance(idx, int):
        if idx < 0:
            idx += len(tupty)
        if not 0 <= idx < len(tupty):
            raise IndexError("cannot index at %d in %s" % (idx, tupty))
        res = builder.extract_value(tup, idx)
    elif isinstance(idx, slice):
        items = cgutils.unpack_tuple(builder, tup)[idx]
        res = context.make_tuple(builder, sig.return_type, items)
    elif isinstance(tupty, types.LiteralStrKeyDict):
        # pretend to be a dictionary
        idx_val = idxty.literal_value
        idx_offset = tupty.fields.index(idx_val)
        res = builder.extract_value(tup, idx_offset)
    else:
        raise NotImplementedError("unexpected index %r for %s"
                                  % (idx, sig.args[0]))
    return impl_ret_borrowed(context, builder, sig.return_type, res)


#------------------------------------------------------------------------------
# Implicit conversion

@lower_cast(types.BaseTuple, types.BaseTuple)
def tuple_to_tuple(context, builder, fromty, toty, val):
    if (isinstance(fromty, types.BaseNamedTuple)
        or isinstance(toty, types.BaseNamedTuple)):
        # Disallowed by typing layer
        raise NotImplementedError

    if len(fromty) != len(toty):
        # Disallowed by typing layer
        raise NotImplementedError

    olditems = cgutils.unpack_tuple(builder, val, len(fromty))
    items = [context.cast(builder, v, f, t)
             for v, f, t in zip(olditems, fromty, toty)]
    return context.make_tuple(builder, toty, items)


#------------------------------------------------------------------------------
# Methods

@overload_method(types.BaseTuple, 'index')
def tuple_index(tup, value):

    def tuple_index_impl(tup, value):
        for i in range(len(tup)):
            if tup[i] == value:
                return i
        raise ValueError("tuple.index(x): x not in tuple")

    return tuple_index_impl


@overload(operator.contains)
def in_seq_empty_tuple(x, y):
    if isinstance(x, types.Tuple) and not x.types:
        return lambda x, y: False
