import sys

from .api import *
from .vector_types import vector_types
from .reduction import Reduce
from .cudadrv.devicearray import (device_array, device_array_like, pinned,
                                  pinned_array, pinned_array_like,
                                  mapped_array, to_device, auto_device)
from .cudadrv import devicearray
from .cudadrv.devices import require_context, gpus
from .cudadrv.devices import get_context as current_context
from .cudadrv.runtime import runtime
from numba.core import config
reduce = Reduce

# Register simulated vector types as module level variables
for name, svty in vector_types.items():
    setattr(sys.modules[__name__], name, svty)
    for alias in svty.aliases:
        setattr(sys.modules[__name__], alias, svty)
del vector_types, name, svty, alias

# Ensure that any user code attempting to import cudadrv etc. gets the
# simulator's version and not the real version if the simulator is enabled.
if config.ENABLE_CUDASIM:
    import sys
    from numba.cuda.simulator import cudadrv
    sys.modules['numba.cuda.cudadrv'] = cudadrv
    sys.modules['numba.cuda.cudadrv.devicearray'] = cudadrv.devicearray
    sys.modules['numba.cuda.cudadrv.devices'] = cudadrv.devices
    sys.modules['numba.cuda.cudadrv.driver'] = cudadrv.driver
    sys.modules['numba.cuda.cudadrv.runtime'] = cudadrv.runtime
    sys.modules['numba.cuda.cudadrv.drvapi'] = cudadrv.drvapi
    sys.modules['numba.cuda.cudadrv.error'] = cudadrv.error
    sys.modules['numba.cuda.cudadrv.nvvm'] = cudadrv.nvvm

    from . import compiler
    sys.modules['numba.cuda.compiler'] = compiler
