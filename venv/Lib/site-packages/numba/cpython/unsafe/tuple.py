"""
This file provides internal compiler utilities that support certain special
operations with tuple and workarounds for limitations enforced in userland.
"""

from numba.core import types, typing, errors
from numba.core.cgutils import alloca_once
from numba.core.extending import intrinsic


@intrinsic
def tuple_setitem(typingctx, tup, idx, val):
    """Return a copy of the tuple with item at *idx* replaced with *val*.

    Operation: ``out = tup[:idx] + (val,) + tup[idx + 1:]

    **Warning**

    - No boundchecking.
    - The dtype of the tuple cannot be changed.
      *val* is always cast to the existing dtype of the tuple.
    """
    def codegen(context, builder, signature, args):
        tup, idx, val = args
        stack = alloca_once(builder, tup.type)
        builder.store(tup, stack)
        # Unsafe load on unchecked bounds.  Poison value maybe returned.
        offptr = builder.gep(stack, [idx.type(0), idx], inbounds=True)
        builder.store(val, offptr)
        return builder.load(stack)

    sig = tup(tup, idx, tup.dtype)
    return sig, codegen


@intrinsic
def build_full_slice_tuple(tyctx, sz):
    """Creates a sz-tuple of full slices."""
    if not isinstance(sz, types.IntegerLiteral):
        raise errors.RequireLiteralValue(sz)

    size = int(sz.literal_value)
    tuple_type = types.UniTuple(dtype=types.slice2_type, count=size)
    sig = tuple_type(sz)

    def codegen(context, builder, signature, args):
        def impl(length, empty_tuple):
            out = empty_tuple
            for i in range(length):
                out = tuple_setitem(out, i, slice(None, None))
            return out

        inner_argtypes = [types.intp, tuple_type]
        inner_sig = typing.signature(tuple_type, *inner_argtypes)
        ll_idx_type = context.get_value_type(types.intp)
        # Allocate an empty tuple
        empty_tuple = context.get_constant_undef(tuple_type)
        inner_args = [ll_idx_type(size), empty_tuple]

        res = context.compile_internal(builder, impl, inner_sig, inner_args)
        return res

    return sig, codegen


@intrinsic
def unpack_single_tuple(tyctx, tup):
    """This exists to handle the situation y = (*x,), the interpreter injects a
    call to it in the case of a single value unpack. It's not possible at
    interpreting time to differentiate between an unpack on a variable sized
    container e.g. list and a fixed one, e.g. tuple. This function handles the
    situation should it arise.
    """
    # See issue #6534
    if not isinstance(tup, types.BaseTuple):
        msg = (f"Only tuples are supported when unpacking a single item, "
               f"got type: {tup}")
        raise errors.UnsupportedError(msg)

    sig = tup(tup)

    def codegen(context, builder, signature, args):
        return args[0] # there's only one tuple and it's a simple pass through
    return sig, codegen
