"""
Implementation of some CFFI functions
"""


from numba.core.imputils import Registry
from numba.core import types
from numba.np import arrayobj

registry = Registry('cffiimpl')

@registry.lower('ffi.from_buffer', types.Buffer)
def from_buffer(context, builder, sig, args):
    assert len(sig.args) == 1
    assert len(args) == 1
    [fromty] = sig.args
    [val] = args
    # Type inference should have prevented passing a buffer from an
    # array to a pointer of the wrong type
    assert fromty.dtype == sig.return_type.dtype
    ary = arrayobj.make_array(fromty)(context, builder, val)
    return ary.data
