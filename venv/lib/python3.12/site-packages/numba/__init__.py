"""
Expose top-level symbols that are safe for import *
"""

import platform
import re
import sys
import warnings


# ---------------------- WARNING WARNING WARNING ----------------------------
# THIS MUST RUN FIRST, DO NOT MOVE... SEE DOCSTRING IN _ensure_critical_deps
def _ensure_critical_deps():
    """
    Make sure the Python, NumPy and SciPy present are supported versions.
    This has to be done _before_ importing anything from Numba such that
    incompatible versions can be reported to the user. If this occurs _after_
    importing things from Numba and there's an issue in e.g. a Numba c-ext, a
    SystemError might have occurred which prevents reporting the likely cause of
    the problem (incompatible versions of critical dependencies).
    """
    #NOTE THIS CODE SHOULD NOT IMPORT ANYTHING FROM NUMBA!

    def extract_version(mod):
        return tuple(map(int, mod.__version__.split('.')[:2]))

    PYVERSION = sys.version_info[:2]

    if PYVERSION < (3, 10):
        msg = ("Numba needs Python 3.10 or greater. Got Python "
               f"{PYVERSION[0]}.{PYVERSION[1]}.")
        raise ImportError(msg)

    import numpy as np
    numpy_version = extract_version(np)

    if numpy_version < (1, 22):
        msg = (f"Numba needs NumPy 1.22 or greater. Got NumPy "
               f"{numpy_version[0]}.{numpy_version[1]}.")
        raise ImportError(msg)

    if numpy_version > (2, 3):
        msg = (f"Numba needs NumPy 2.3 or less. Got NumPy "
                f"{numpy_version[0]}.{numpy_version[1]}.")
        raise ImportError(msg)

    try:
        import scipy
    except ImportError:
        pass
    else:
        sp_version = extract_version(scipy)
        if sp_version < (1, 0):
            msg = ("Numba requires SciPy version 1.0 or greater. Got SciPy "
                   f"{scipy.__version__}.")
            raise ImportError(msg)


_ensure_critical_deps()
# END DO NOT MOVE
# ---------------------- WARNING WARNING WARNING ----------------------------


from ._version import get_versions
from numba.misc.init_utils import generate_version_info

__version__ = get_versions()['version']
version_info = generate_version_info(__version__)
del get_versions
del generate_version_info


from numba.core import config
from numba.core import types, errors

# Re-export typeof
from numba.misc.special import (
    typeof, prange, pndindex, gdb, gdb_breakpoint, gdb_init,
    literally, literal_unroll,
)

# Re-export error classes
from numba.core.errors import *

# Re-export types itself
import numba.core.types as types

# Re-export all type names
from numba.core.types import *

# Re-export decorators
from numba.core.decorators import (cfunc, jit, njit, stencil,
                                   jit_module)

# Re-export vectorize decorators and the thread layer querying function
from numba.np.ufunc import (vectorize, guvectorize, threading_layer,
                            get_num_threads, set_num_threads,
                            set_parallel_chunksize, get_parallel_chunksize,
                            get_thread_id)

# Re-export Numpy helpers
from numba.np.numpy_support import carray, farray, from_dtype

# Re-export experimental
from numba import experimental

# Initialize withcontexts
import numba.core.withcontexts
from numba.core.withcontexts import objmode_context as objmode
from numba.core.withcontexts import parallel_chunksize

# Initialize target extensions
import numba.core.target_extension

# Initialize typed containers
import numba.typed

# Keep this for backward compatibility.
def test(argv, **kwds):
    # To speed up the import time, avoid importing `unittest` and other test
    # dependencies unless the user is actually trying to run tests.
    from numba.testing import _runtests as runtests
    return runtests.main(argv, **kwds)

__all__ = """
    cfunc
    from_dtype
    guvectorize
    jit
    experimental
    njit
    stencil
    jit_module
    typeof
    prange
    gdb
    gdb_breakpoint
    gdb_init
    vectorize
    objmode
    literal_unroll
    get_num_threads
    set_num_threads
    set_parallel_chunksize
    get_parallel_chunksize
    parallel_chunksize
    """.split() + types.__all__ + errors.__all__


_min_llvmlite_version = (0, 46, 0)
_min_llvm_version = (14, 0, 0)

def _ensure_llvm():
    """
    Make sure llvmlite is operational.
    """
    import warnings
    import llvmlite

    # Only look at the major, minor and bugfix version numbers.
    # Ignore other stuffs
    regex = re.compile(r'(\d+)\.(\d+).(\d+)')
    m = regex.match(llvmlite.__version__)
    if m:
        ver = tuple(map(int, m.groups()))
        if ver < _min_llvmlite_version:
            msg = ("Numba requires at least version %d.%d.%d of llvmlite.\n"
                   "Installed version is %s.\n"
                   "Please update llvmlite." %
                   (_min_llvmlite_version + (llvmlite.__version__,)))
            raise ImportError(msg)
    else:
        # Not matching?
        warnings.warn("llvmlite version format not recognized!")

    from llvmlite.binding import llvm_version_info, check_jit_execution

    if llvm_version_info < _min_llvm_version:
        msg = ("Numba requires at least version %d.%d.%d of LLVM.\n"
               "Installed llvmlite is built against version %d.%d.%d.\n"
               "Please update llvmlite." %
               (_min_llvm_version + llvm_version_info))
        raise ImportError(msg)

    check_jit_execution()


def _try_enable_svml():
    """
    Tries to enable SVML if configuration permits use and the library is found.
    """
    if not config.DISABLE_INTEL_SVML:
        try:
            if sys.platform.startswith('linux'):
                llvmlite.binding.load_library_permanently("libsvml.so")
            elif sys.platform.startswith('darwin'):
                llvmlite.binding.load_library_permanently("libsvml.dylib")
            elif sys.platform.startswith('win'):
                llvmlite.binding.load_library_permanently("svml_dispmd")
            else:
                return False
            # The SVML library is loaded, therefore SVML *could* be supported.
            # Now see if LLVM has been compiled with the SVML support patch.
            # If llvmlite has the checking function `has_svml` and it returns
            # True, then LLVM was compiled with SVML support and the setup
            # for SVML can proceed. We err on the side of caution and if the
            # checking function is missing, regardless of that being fine for
            # most 0.23.{0,1} llvmlite instances (i.e. conda or pip installed),
            # we assume that SVML was not compiled in. llvmlite 0.23.2 is a
            # bugfix release with the checking function present that will always
            # produce correct behaviour. For context see: #3006.
            try:
                if not getattr(llvmlite.binding.targets, "has_svml")():
                    # has detection function, but no svml compiled in, therefore
                    # disable SVML
                    return False
            except AttributeError:
                if platform.machine() == 'x86_64' and config.DEBUG:
                    msg = ("SVML was found but llvmlite >= 0.23.2 is "
                           "needed to support it.")
                    warnings.warn(msg)
                # does not have detection function, cannot detect reliably,
                # disable SVML.
                return False

            # All is well, detection function present and reports SVML is
            # compiled in, set the vector library to SVML.
            llvmlite.binding.set_option('SVML', '-vector-library=SVML')
            return True
        except:
            if platform.machine() == 'x86_64' and config.DEBUG:
                warnings.warn("SVML was not found/could not be loaded.")
    return False

_ensure_llvm()

# we know llvmlite is working as the above tests passed, import it now as SVML
# needs to mutate runtime options (sets the `-vector-library`).
import llvmlite

"""
Is set to True if Intel SVML is in use.
"""
config.USING_SVML = _try_enable_svml()


# ---------------------- WARNING WARNING WARNING ----------------------------
# The following imports occur below here (SVML init) because somewhere in their
# import sequence they have a `@njit` wrapped function. This triggers too early
# a bind to the underlying LLVM libraries which then irretrievably sets the LLVM
# SVML state to "no SVML". See https://github.com/numba/numba/issues/4689 for
# context.
# ---------------------- WARNING WARNING WARNING ----------------------------
