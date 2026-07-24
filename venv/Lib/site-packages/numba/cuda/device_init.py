# Re export
import sys
from numba.cuda import cg
from .stubs import (threadIdx, blockIdx, blockDim, gridDim, laneid, warpsize,
                    syncwarp, shared, local, const, atomic,
                    shfl_sync_intrinsic, vote_sync_intrinsic, match_any_sync,
                    match_all_sync, threadfence_block, threadfence_system,
                    threadfence, selp, popc, brev, clz, ffs, fma, cbrt,
                    activemask, lanemask_lt, nanosleep, fp16,
                    _vector_type_stubs)
from .intrinsics import (grid, gridsize, syncthreads, syncthreads_and,
                         syncthreads_count, syncthreads_or)
from .cudadrv.error import CudaSupportError
from numba.cuda.cudadrv.driver import (BaseCUDAMemoryManager,
                                       HostOnlyCUDAMemoryManager,
                                       GetIpcHandleMixin, MemoryPointer,
                                       MappedMemory, PinnedMemory, MemoryInfo,
                                       IpcHandle, set_memory_manager)
from numba.cuda.cudadrv.runtime import runtime
from .cudadrv import nvvm
from numba.cuda import initialize
from .errors import KernelRuntimeError

from .decorators import jit, declare_device
from .api import *
from .api import _auto_device
from .args import In, Out, InOut

from .intrinsic_wrapper import (all_sync, any_sync, eq_sync, ballot_sync,
                                shfl_sync, shfl_up_sync, shfl_down_sync,
                                shfl_xor_sync)

from .kernels import reduction

reduce = Reduce = reduction.Reduce

# Expose vector type constructors and aliases as module level attributes.
for vector_type_stub in _vector_type_stubs:
    setattr(sys.modules[__name__], vector_type_stub.__name__, vector_type_stub)
    for alias in vector_type_stub.aliases:
        setattr(sys.modules[__name__], alias, vector_type_stub)
del vector_type_stub, _vector_type_stubs


def is_available():
    """Returns a boolean to indicate the availability of a CUDA GPU.

    This will initialize the driver if it hasn't been initialized.
    """
    # whilst `driver.is_available` will init the driver itself,
    # the driver initialization may raise and as a result break
    # test discovery/orchestration as `cuda.is_available` is often
    # used as a guard for whether to run a CUDA test, the try/except
    # below is to handle this case.
    driver_is_available = False
    try:
        driver_is_available = driver.driver.is_available
    except CudaSupportError:
        pass

    return driver_is_available and nvvm.is_available()


def is_supported_version():
    """Returns True if the CUDA Runtime is a supported version.

    Unsupported versions (e.g. newer versions than those known to Numba)
    may still work; this function provides a facility to check whether the
    current Numba version is tested and known to work with the current
    runtime version. If the current version is unsupported, the caller can
    decide how to act. Options include:

    - Continuing silently,
    - Emitting a warning,
    - Generating an error or otherwise preventing the use of CUDA.
    """

    return runtime.is_supported_version()


def cuda_error():
    """Returns None if there was no error initializing the CUDA driver.
    If there was an error initializing the driver, a string describing the
    error is returned.
    """
    return driver.driver.initialization_error


initialize.initialize_all()
