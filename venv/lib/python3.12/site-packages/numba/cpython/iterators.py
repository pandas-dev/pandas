"""
Implementation of various iterable and iterator types.
"""

from numba.core import types, cgutils
from numba.core.imputils import (
    lower_builtin, iternext_impl, call_iternext, call_getiter,
    impl_ret_borrowed, impl_ret_new_ref, RefType)



@lower_builtin('getiter', types.IteratorType)
def iterator_getiter(context, builder, sig, args):
    [it] = args
    return impl_ret_borrowed(context, builder, sig.return_type, it)

#-------------------------------------------------------------------------------
# builtin `enumerate` implementation

@lower_builtin(enumerate, types.IterableType)
@lower_builtin(enumerate, types.IterableType, types.Integer)
def make_enumerate_object(context, builder, sig, args):
    assert len(args) == 1 or len(args) == 2 # enumerate(it) or enumerate(it, start)
    srcty = sig.args[0]

    if len(args) == 1:
        src = args[0]
        start_val = context.get_constant(types.intp, 0)
    elif len(args) == 2:
        src = args[0]
        start_val = context.cast(builder, args[1], sig.args[1], types.intp)

    iterobj = call_getiter(context, builder, srcty, src)

    enum = context.make_helper(builder, sig.return_type)

    countptr = cgutils.alloca_once(builder, start_val.type)
    builder.store(start_val, countptr)

    enum.count = countptr
    enum.iter = iterobj

    res = enum._getvalue()
    return impl_ret_new_ref(context, builder, sig.return_type, res)

@lower_builtin('iternext', types.EnumerateType)
@iternext_impl(RefType.NEW)
def iternext_enumerate(context, builder, sig, args, result):
    [enumty] = sig.args
    [enum] = args

    enum = context.make_helper(builder, enumty, value=enum)

    count = builder.load(enum.count)
    ncount = builder.add(count, context.get_constant(types.intp, 1))
    builder.store(ncount, enum.count)

    srcres = call_iternext(context, builder, enumty.source_type, enum.iter)
    is_valid = srcres.is_valid()
    result.set_valid(is_valid)

    with builder.if_then(is_valid):
        srcval = srcres.yielded_value()
        result.yield_(context.make_tuple(builder, enumty.yield_type,
                                         [count, srcval]))


#-------------------------------------------------------------------------------
# builtin `zip` implementation

@lower_builtin(zip, types.VarArg(types.Any))
def make_zip_object(context, builder, sig, args):
    zip_type = sig.return_type

    assert len(args) == len(zip_type.source_types)

    zipobj = context.make_helper(builder, zip_type)

    for i, (arg, srcty) in enumerate(zip(args, sig.args)):
        zipobj[i] = call_getiter(context, builder, srcty, arg)

    res = zipobj._getvalue()
    return impl_ret_new_ref(context, builder, sig.return_type, res)

@lower_builtin('iternext', types.ZipType)
@iternext_impl(RefType.NEW)
def iternext_zip(context, builder, sig, args, result):
    [zip_type] = sig.args
    [zipobj] = args

    zipobj = context.make_helper(builder, zip_type, value=zipobj)

    if len(zipobj) == 0:
        # zip() is an empty iterator
        result.set_exhausted()
        return

    p_ret_tup = cgutils.alloca_once(builder,
                                    context.get_value_type(zip_type.yield_type))
    p_is_valid = cgutils.alloca_once_value(builder, value=cgutils.true_bit)

    for i, (iterobj, srcty) in enumerate(zip(zipobj, zip_type.source_types)):
        is_valid = builder.load(p_is_valid)
        # Avoid calling the remaining iternext if a iterator has been exhausted
        with builder.if_then(is_valid):
            srcres = call_iternext(context, builder, srcty, iterobj)
            is_valid = builder.and_(is_valid, srcres.is_valid())
            builder.store(is_valid, p_is_valid)
            val = srcres.yielded_value()
            ptr = cgutils.gep_inbounds(builder, p_ret_tup, 0, i)
            builder.store(val, ptr)

    is_valid = builder.load(p_is_valid)
    result.set_valid(is_valid)

    with builder.if_then(is_valid):
        result.yield_(builder.load(p_ret_tup))


#-------------------------------------------------------------------------------
# generator implementation

@lower_builtin('iternext', types.Generator)
@iternext_impl(RefType.BORROWED)
def iternext_zip(context, builder, sig, args, result):
    genty, = sig.args
    gen, = args
    impl = context.get_generator_impl(genty)
    status, retval = impl(context, builder, sig, args)
    context.add_linking_libs(getattr(impl, 'libs', ()))

    with cgutils.if_likely(builder, status.is_ok):
        result.set_valid(True)
        result.yield_(retval)
    with cgutils.if_unlikely(builder, status.is_stop_iteration):
        result.set_exhausted()
    with cgutils.if_unlikely(builder,
                             builder.and_(status.is_error,
                                          builder.not_(status.is_stop_iteration))):
        context.call_conv.return_status_propagate(builder, status)
