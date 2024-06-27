"""
This file implements the lowering for `dict()`
"""
from numba.core import types
from numba.core.imputils import lower_builtin


_message_dict_support = """
Unsupported use of `dict()` with keyword argument(s). \
The only supported uses are `dict()` or `dict(*iterable)`.
""".strip()


@lower_builtin(dict, types.IterableType)
def dict_constructor(context, builder, sig, args):
    from numba.typed import Dict

    dicttype = sig.return_type
    kt, vt = dicttype.key_type, dicttype.value_type

    def dict_impl(iterable):
        res = Dict.empty(kt, vt)
        for k, v in iterable:
            res[k] = v
        return res

    return context.compile_internal(builder, dict_impl, sig, args)


@lower_builtin(dict)
def impl_dict(context, builder, sig, args):
    """
    The `dict()` implementation simply forwards the work to `Dict.empty()`.
    """
    from numba.typed import Dict

    dicttype = sig.return_type
    kt, vt = dicttype.key_type, dicttype.value_type

    def call_ctor():
        return Dict.empty(kt, vt)

    return context.compile_internal(builder, call_ctor, sig, args)
