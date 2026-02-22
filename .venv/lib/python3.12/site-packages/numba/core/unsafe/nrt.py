"""
Contains unsafe intrinsic that calls NRT C API
"""

from numba.core import types
from numba.core.typing import signature
from numba.core.extending import intrinsic


@intrinsic
def NRT_get_api(tyctx):
    """NRT_get_api()

    Calls NRT_get_api() from the NRT C API
    Returns LLVM Type i8* (void pointer)
    """
    def codegen(cgctx, builder, sig, args):
        return cgctx.nrt.get_nrt_api(builder)
    sig = signature(types.voidptr)
    return sig, codegen
