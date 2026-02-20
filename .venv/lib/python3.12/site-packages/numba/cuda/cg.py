from numba.core import types
from numba.core.extending import overload, overload_method
from numba.core.typing import signature
from numba.cuda import nvvmutils
from numba.cuda.extending import intrinsic
from numba.cuda.types import grid_group, GridGroup as GridGroupClass


class GridGroup:
    """A cooperative group representing the entire grid"""

    def sync() -> None:
        """Synchronize this grid group"""


def this_grid() -> GridGroup:
    """Get the current grid group."""
    return GridGroup()


@intrinsic
def _this_grid(typingctx):
    sig = signature(grid_group)

    def codegen(context, builder, sig, args):
        one = context.get_constant(types.int32, 1)
        mod = builder.module
        return builder.call(
            nvvmutils.declare_cudaCGGetIntrinsicHandle(mod),
            (one,))

    return sig, codegen


@overload(this_grid, target='cuda')
def _ol_this_grid():
    def impl():
        return _this_grid()

    return impl


@intrinsic
def _grid_group_sync(typingctx, group):
    sig = signature(types.int32, group)

    def codegen(context, builder, sig, args):
        flags = context.get_constant(types.int32, 0)
        mod = builder.module
        return builder.call(
            nvvmutils.declare_cudaCGSynchronize(mod),
            (*args, flags))

    return sig, codegen


@overload_method(GridGroupClass, 'sync', target='cuda')
def _ol_grid_group_sync(group):
    def impl(group):
        return _grid_group_sync(group)

    return impl
