"""
CUDA Runtime wrapper.

This provides a very minimal set of bindings, since the Runtime API is not
really used in Numba except for querying the Runtime version.
"""

import ctypes
import functools
import sys

from numba.core import config
from numba.cuda.cudadrv.driver import ERROR_MAP, make_logger
from numba.cuda.cudadrv.error import CudaSupportError, CudaRuntimeError
from numba.cuda.cudadrv.libs import open_cudalib
from numba.cuda.cudadrv.rtapi import API_PROTOTYPES
from numba.cuda.cudadrv import enums


class CudaRuntimeAPIError(CudaRuntimeError):
    """
    Raised when there is an error accessing a C API from the CUDA Runtime.
    """
    def __init__(self, code, msg):
        self.code = code
        self.msg = msg
        super().__init__(code, msg)

    def __str__(self):
        return "[%s] %s" % (self.code, self.msg)


class Runtime:
    """
    Runtime object that lazily binds runtime API functions.
    """

    def __init__(self):
        self.is_initialized = False

    def _initialize(self):
        # lazily initialize logger
        global _logger
        _logger = make_logger()

        if config.DISABLE_CUDA:
            msg = ("CUDA is disabled due to setting NUMBA_DISABLE_CUDA=1 "
                   "in the environment, or because CUDA is unsupported on "
                   "32-bit systems.")
            raise CudaSupportError(msg)
        self.lib = open_cudalib('cudart')

        self.is_initialized = True

    def __getattr__(self, fname):
        # First request of a runtime API function
        try:
            proto = API_PROTOTYPES[fname]
        except KeyError:
            raise AttributeError(fname)
        restype = proto[0]
        argtypes = proto[1:]

        if not self.is_initialized:
            self._initialize()

        # Find function in runtime library
        libfn = self._find_api(fname)
        libfn.restype = restype
        libfn.argtypes = argtypes

        safe_call = self._wrap_api_call(fname, libfn)
        setattr(self, fname, safe_call)
        return safe_call

    def _wrap_api_call(self, fname, libfn):
        @functools.wraps(libfn)
        def safe_cuda_api_call(*args):
            _logger.debug('call runtime api: %s', libfn.__name__)
            retcode = libfn(*args)
            self._check_error(fname, retcode)
        return safe_cuda_api_call

    def _check_error(self, fname, retcode):
        if retcode != enums.CUDA_SUCCESS:
            errname = ERROR_MAP.get(retcode, "cudaErrorUnknown")
            msg = "Call to %s results in %s" % (fname, errname)
            _logger.error(msg)
            raise CudaRuntimeAPIError(retcode, msg)

    def _find_api(self, fname):
        try:
            return getattr(self.lib, fname)
        except AttributeError:
            pass

        # Not found.
        # Delay missing function error to use
        def absent_function(*args, **kws):
            msg = "runtime missing function: %s."
            raise CudaRuntimeError(msg % fname)

        setattr(self, fname, absent_function)
        return absent_function

    def get_version(self):
        """
        Returns the CUDA Runtime version as a tuple (major, minor).
        """
        rtver = ctypes.c_int()
        self.cudaRuntimeGetVersion(ctypes.byref(rtver))
        # The version is encoded as (1000 * major) + (10 * minor)
        major = rtver.value // 1000
        minor = (rtver.value - (major * 1000)) // 10
        return (major, minor)

    def is_supported_version(self):
        """
        Returns True if the CUDA Runtime is a supported version.
        """

        return self.get_version() in self.supported_versions

    @property
    def supported_versions(self):
        """A tuple of all supported CUDA toolkit versions. Versions are given in
        the form ``(major_version, minor_version)``."""
        if sys.platform not in ('linux', 'win32') or config.MACHINE_BITS != 64:
            # Only 64-bit Linux and Windows are supported
            return ()
        return ((11, 0), (11, 1), (11, 2), (11, 3), (11, 4), (11, 5), (11, 6),
                (11, 7))


runtime = Runtime()


def get_version():
    """
    Return the runtime version as a tuple of (major, minor)
    """
    return runtime.get_version()
