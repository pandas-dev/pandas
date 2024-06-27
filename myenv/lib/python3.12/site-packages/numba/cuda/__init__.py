from numba import runtests
from numba.core import config

if config.ENABLE_CUDASIM:
    from .simulator_init import *
else:
    from .device_init import *
    from .device_init import _auto_device

from numba.cuda.compiler import (compile, compile_for_current_device,
                                 compile_ptx, compile_ptx_for_current_device)

def test(*args, **kwargs):
    if not is_available():
        raise cuda_error()

    return runtests.main("numba.cuda.tests", *args, **kwargs)
