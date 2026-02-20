"""This module is used to switch between C and Python implementations of the
wrappers.
"""

import os

from .wrappers import BoundFunctionWrapper, CallableObjectProxy, FunctionWrapper
from .wrappers import ObjectProxy as BaseObjectProxy
from .wrappers import PartialCallableObjectProxy, _FunctionWrapperBase

# Try to use C extensions if not disabled.

_using_c_extension = False

_use_extensions = not os.environ.get("WRAPT_DISABLE_EXTENSIONS")

if _use_extensions:
    try:
        from ._wrappers import (  # type: ignore[no-redef,import-not-found]
            BoundFunctionWrapper,
            CallableObjectProxy,
            FunctionWrapper,
        )
        from ._wrappers import ObjectProxy as BaseObjectProxy  # type: ignore[no-redef]
        from ._wrappers import (  # type: ignore[no-redef,import-not-found]
            PartialCallableObjectProxy,
            _FunctionWrapperBase,
        )

        _using_c_extension = True
    except ImportError:
        # C extensions not available, using Python implementations
        pass


def partial(*args, **kwargs):
    """Create a callable object proxy with partial application of the given
    arguments and keywords. This behaves the same as `functools.partial`, but
    implemented using the `ObjectProxy` class to provide better support for
    introspection.
    """
    return PartialCallableObjectProxy(*args, **kwargs)
