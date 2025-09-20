import operator

from numba.core import types, typing, cgutils

from numba.core.imputils import (lower_cast, lower_builtin,
                                 lower_getattr_generic, impl_ret_untracked,
                                 lower_setattr_generic)


def always_return_true_impl(context, builder, sig, args):
    return cgutils.true_bit


def always_return_false_impl(context, builder, sig, args):
    return cgutils.false_bit


def optional_is_none(context, builder, sig, args):
    """
    Check if an Optional value is invalid
    """
    [lty, rty] = sig.args
    [lval, rval] = args

    # Make sure None is on the right
    if lty == types.none:
        lty, rty = rty, lty
        lval, rval = rval, lval

    opt_type = lty
    opt_val = lval

    opt = context.make_helper(builder, opt_type, opt_val)
    res = builder.not_(cgutils.as_bool_bit(builder, opt.valid))
    return impl_ret_untracked(context, builder, sig.return_type, res)


# None is/not None
lower_builtin(operator.is_, types.none, types.none)(always_return_true_impl)

# Optional is None
lower_builtin(operator.is_, types.Optional, types.none)(optional_is_none)
lower_builtin(operator.is_, types.none, types.Optional)(optional_is_none)


@lower_getattr_generic(types.Optional)
def optional_getattr(context, builder, typ, value, attr):
    """
    Optional.__getattr__ => redirect to the wrapped type.
    """
    inner_type = typ.type
    val = context.cast(builder, value, typ, inner_type)
    imp = context.get_getattr(inner_type, attr)
    return imp(context, builder, inner_type, val, attr)


@lower_setattr_generic(types.Optional)
def optional_setattr(context, builder, sig, args, attr):
    """
    Optional.__setattr__ => redirect to the wrapped type.
    """
    basety, valty = sig.args
    target, val = args
    target_type = basety.type
    target = context.cast(builder, target, basety, target_type)

    newsig = typing.signature(sig.return_type, target_type, valty)
    imp = context.get_setattr(attr, newsig)
    return imp(builder, (target, val))


@lower_cast(types.Optional, types.Optional)
def optional_to_optional(context, builder, fromty, toty, val):
    """
    The handling of optional->optional cast must be special cased for
    correct propagation of None value.  Given type T and U. casting of
    T? to U? (? denotes optional) should always succeed.   If the from-value
    is None, the None value the casted value (U?) should be None; otherwise,
    the from-value is casted to U. This is different from casting T? to U,
    which requires the from-value must not be None.
    """
    optval = context.make_helper(builder, fromty, value=val)
    validbit = cgutils.as_bool_bit(builder, optval.valid)
    # Create uninitialized optional value
    outoptval = context.make_helper(builder, toty)

    with builder.if_else(validbit) as (is_valid, is_not_valid):
        with is_valid:
            # Cast internal value
            outoptval.valid = cgutils.true_bit
            outoptval.data = context.cast(builder, optval.data,
                                          fromty.type, toty.type)

        with is_not_valid:
            # Store None to result
            outoptval.valid = cgutils.false_bit
            outoptval.data = cgutils.get_null_value(
                outoptval.data.type)

    return outoptval._getvalue()


@lower_cast(types.Any, types.Optional)
def any_to_optional(context, builder, fromty, toty, val):
    if fromty == types.none:
        return context.make_optional_none(builder, toty.type)
    else:
        val = context.cast(builder, val, fromty, toty.type)
        return context.make_optional_value(builder, toty.type, val)


@lower_cast(types.Optional, types.Any)
@lower_cast(types.Optional, types.Boolean)
def optional_to_any(context, builder, fromty, toty, val):
    optval = context.make_helper(builder, fromty, value=val)
    validbit = cgutils.as_bool_bit(builder, optval.valid)
    with builder.if_then(builder.not_(validbit), likely=False):
        msg = "expected %s, got None" % (fromty.type,)
        context.call_conv.return_user_exc(builder, TypeError, (msg,))

    return context.cast(builder, optval.data, fromty.type, toty)
