"""
Exception handling intrinsics.
"""

from numba.core import types, errors, cgutils
from numba.core.extending import intrinsic


@intrinsic
def exception_check(typingctx):
    """An intrinsic to check if an exception is raised
    """
    def codegen(context, builder, signature, args):
        nrt = context.nrt
        return nrt.eh_check(builder)

    restype = types.boolean
    return restype(), codegen


@intrinsic
def mark_try_block(typingctx):
    """An intrinsic to mark the start of a *try* block.
    """
    def codegen(context, builder, signature, args):
        nrt = context.nrt
        nrt.eh_try(builder)
        return context.get_dummy_value()

    restype = types.none
    return restype(), codegen


@intrinsic
def end_try_block(typingctx):
    """An intrinsic to mark the end of a *try* block.
    """
    def codegen(context, builder, signature, args):
        nrt = context.nrt
        nrt.eh_end_try(builder)
        return context.get_dummy_value()

    restype = types.none
    return restype(), codegen


@intrinsic
def exception_match(typingctx, exc_value, exc_class):
    """Basically do ``isinstance(exc_value, exc_class)`` for exception objects.
    Used in ``except Exception:`` syntax.
    """
    # Check for our limitation
    if exc_class.exc_class is not Exception:
        msg = "Exception matching is limited to {}"
        raise errors.UnsupportedError(msg.format(Exception))

    def codegen(context, builder, signature, args):
        # Intentionally always True.
        return cgutils.true_bit

    restype = types.boolean
    return restype(exc_value, exc_class), codegen
