import os
import platform
import shutil

from numba.tests.support import SerialMixin
from numba.cuda.cuda_paths import get_conda_ctk
from numba.cuda.cudadrv import driver, devices, libs
from numba.core import config
from numba.tests.support import TestCase
from pathlib import Path
import unittest

numba_cuda_dir = Path(__file__).parent
test_data_dir = numba_cuda_dir / 'tests' / 'data'


class CUDATestCase(SerialMixin, TestCase):
    """
    For tests that use a CUDA device. Test methods in a CUDATestCase must not
    be run out of module order, because the ContextResettingTestCase may reset
    the context and destroy resources used by a normal CUDATestCase if any of
    its tests are run between tests from a CUDATestCase.
    """

    def setUp(self):
        self._low_occupancy_warnings = config.CUDA_LOW_OCCUPANCY_WARNINGS
        self._warn_on_implicit_copy = config.CUDA_WARN_ON_IMPLICIT_COPY

        # Disable warnings about low gpu utilization in the test suite
        config.CUDA_LOW_OCCUPANCY_WARNINGS = 0
        # Disable warnings about host arrays in the test suite
        config.CUDA_WARN_ON_IMPLICIT_COPY = 0

    def tearDown(self):
        config.CUDA_LOW_OCCUPANCY_WARNINGS = self._low_occupancy_warnings
        config.CUDA_WARN_ON_IMPLICIT_COPY = self._warn_on_implicit_copy

    def skip_if_lto(self, reason):
        # Some linkers need the compute capability to be specified, so we
        # always specify it here.
        cc = devices.get_context().device.compute_capability
        linker = driver.Linker.new(cc=cc)
        if linker.lto:
            self.skipTest(reason)


class ContextResettingTestCase(CUDATestCase):
    """
    For tests where the context needs to be reset after each test. Typically
    these inspect or modify parts of the context that would usually be expected
    to be internal implementation details (such as the state of allocations and
    deallocations, etc.).
    """

    def tearDown(self):
        super().tearDown()
        from numba.cuda.cudadrv.devices import reset
        reset()


def ensure_supported_ccs_initialized():
    from numba.cuda import is_available as cuda_is_available
    from numba.cuda.cudadrv import nvvm

    if cuda_is_available():
        # Ensure that cudart.so is loaded and the list of supported compute
        # capabilities in the nvvm module is populated before a fork. This is
        # needed because some compilation tests don't require a CUDA context,
        # but do use NVVM, and it is required that libcudart.so should be
        # loaded before a fork (note that the requirement is not explicitly
        # documented).
        nvvm.get_supported_ccs()


def skip_on_cudasim(reason):
    """Skip this test if running on the CUDA simulator"""
    return unittest.skipIf(config.ENABLE_CUDASIM, reason)


def skip_unless_cudasim(reason):
    """Skip this test if running on CUDA hardware"""
    return unittest.skipUnless(config.ENABLE_CUDASIM, reason)


def skip_unless_conda_cudatoolkit(reason):
    """Skip test if the CUDA toolkit was not installed by Conda"""
    return unittest.skipUnless(get_conda_ctk() is not None, reason)


def skip_if_external_memmgr(reason):
    """Skip test if an EMM Plugin is in use"""
    return unittest.skipIf(config.CUDA_MEMORY_MANAGER != 'default', reason)


def skip_under_cuda_memcheck(reason):
    return unittest.skipIf(os.environ.get('CUDA_MEMCHECK') is not None, reason)


def skip_without_nvdisasm(reason):
    nvdisasm_path = shutil.which('nvdisasm')
    return unittest.skipIf(nvdisasm_path is None, reason)


def skip_with_nvdisasm(reason):
    nvdisasm_path = shutil.which('nvdisasm')
    return unittest.skipIf(nvdisasm_path is not None, reason)


def skip_on_arm(reason):
    cpu = platform.processor()
    is_arm = cpu.startswith('arm') or cpu.startswith('aarch')
    return unittest.skipIf(is_arm, reason)


def skip_if_cuda_includes_missing(fn):
    # Skip when cuda.h is not available - generally this should indicate
    # whether the CUDA includes are available or not
    cuda_h = os.path.join(config.CUDA_INCLUDE_PATH, 'cuda.h')
    cuda_h_file = (os.path.exists(cuda_h) and os.path.isfile(cuda_h))
    reason = 'CUDA include dir not available on this system'
    return unittest.skipUnless(cuda_h_file, reason)(fn)


def skip_if_mvc_enabled(reason):
    """Skip a test if Minor Version Compatibility is enabled"""
    return unittest.skipIf(config.CUDA_ENABLE_MINOR_VERSION_COMPATIBILITY,
                           reason)


def skip_if_mvc_libraries_unavailable(fn):
    libs_available = False
    try:
        import cubinlinker  # noqa: F401
        import ptxcompiler  # noqa: F401
        libs_available = True
    except ImportError:
        pass

    return unittest.skipUnless(libs_available,
                               "Requires cubinlinker and ptxcompiler")(fn)


def cc_X_or_above(major, minor):
    if not config.ENABLE_CUDASIM:
        cc = devices.get_context().device.compute_capability
        return cc >= (major, minor)
    else:
        return True


def skip_unless_cc_50(fn):
    return unittest.skipUnless(cc_X_or_above(5, 0), "requires cc >= 5.0")(fn)


def skip_unless_cc_53(fn):
    return unittest.skipUnless(cc_X_or_above(5, 3), "requires cc >= 5.3")(fn)


def skip_unless_cc_60(fn):
    return unittest.skipUnless(cc_X_or_above(6, 0), "requires cc >= 6.0")(fn)


def skip_unless_cc_75(fn):
    return unittest.skipUnless(cc_X_or_above(7, 5), "requires cc >= 7.5")(fn)


def xfail_unless_cudasim(fn):
    if config.ENABLE_CUDASIM:
        return fn
    else:
        return unittest.expectedFailure(fn)


def skip_with_cuda_python(reason):
    return unittest.skipIf(driver.USE_NV_BINDING, reason)


def cudadevrt_missing():
    if config.ENABLE_CUDASIM:
        return False
    try:
        path = libs.get_cudalib('cudadevrt', static=True)
        libs.check_static_lib(path)
    except FileNotFoundError:
        return True
    return False


def skip_if_cudadevrt_missing(fn):
    return unittest.skipIf(cudadevrt_missing(), 'cudadevrt missing')(fn)


class ForeignArray(object):
    """
    Class for emulating an array coming from another library through the CUDA
    Array interface. This just hides a DeviceNDArray so that it doesn't look
    like a DeviceNDArray.
    """

    def __init__(self, arr):
        self._arr = arr
        self.__cuda_array_interface__ = arr.__cuda_array_interface__
