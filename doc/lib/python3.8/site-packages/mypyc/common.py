import sys
from typing import Dict, Any, Tuple
import sys

from typing_extensions import Final

PREFIX = 'CPyPy_'  # type: Final # Python wrappers
NATIVE_PREFIX = 'CPyDef_'  # type: Final # Native functions etc.
DUNDER_PREFIX = 'CPyDunder_'  # type: Final # Wrappers for exposing dunder methods to the API
REG_PREFIX = 'cpy_r_'  # type: Final # Registers
STATIC_PREFIX = 'CPyStatic_'  # type: Final # Static variables (for literals etc.)
TYPE_PREFIX = 'CPyType_'  # type: Final # Type object struct
MODULE_PREFIX = 'CPyModule_'  # type: Final # Cached modules
ATTR_PREFIX = '_'  # type: Final # Attributes

ENV_ATTR_NAME = '__mypyc_env__'  # type: Final
NEXT_LABEL_ATTR_NAME = '__mypyc_next_label__'  # type: Final
TEMP_ATTR_NAME = '__mypyc_temp__'  # type: Final
LAMBDA_NAME = '__mypyc_lambda__'  # type: Final
PROPSET_PREFIX = '__mypyc_setter__'  # type: Final
SELF_NAME = '__mypyc_self__'  # type: Final

# Max short int we accept as a literal is based on 32-bit platforms,
# so that we can just always emit the same code.

TOP_LEVEL_NAME = '__top_level__'  # type: Final # Special function representing module top level

# Maximal number of subclasses for a class to trigger fast path in isinstance() checks.
FAST_ISINSTANCE_MAX_SUBCLASSES = 2  # type: Final

IS_32_BIT_PLATFORM = sys.maxsize < (1 << 31)  # type: Final

PLATFORM_SIZE = 4 if IS_32_BIT_PLATFORM else 8

# Python 3.5 on macOS uses a hybrid 32/64-bit build that requires some workarounds.
# The same generated C will be compiled in both 32 and 64 bit modes when building mypy
# wheels (for an unknown reason).
#
# Note that we use "in ['darwin']" because of https://github.com/mypyc/mypyc/issues/761.
IS_MIXED_32_64_BIT_BUILD = sys.platform in ['darwin'] and sys.version_info < (3, 6)  # type: Final

# Maximum value for a short tagged integer.
MAX_SHORT_INT = sys.maxsize >> 1  # type: Final

# Maximum value for a short tagged integer represented as a C integer literal.
#
# Note: Assume that the compiled code uses the same bit width as mypyc, except for
#       Python 3.5 on macOS.
MAX_LITERAL_SHORT_INT = (sys.maxsize >> 1 if not IS_MIXED_32_64_BIT_BUILD
                         else 2**30 - 1)  # type: Final

# Runtime C library files
RUNTIME_C_FILES = [
    'init.c',
    'getargs.c',
    'getargsfast.c',
    'int_ops.c',
    'list_ops.c',
    'dict_ops.c',
    'str_ops.c',
    'set_ops.c',
    'tuple_ops.c',
    'exc_ops.c',
    'misc_ops.c',
    'generic_ops.c',
]  # type: Final


JsonDict = Dict[str, Any]


def decorator_helper_name(func_name: str) -> str:
    return '__mypyc_{}_decorator_helper__'.format(func_name)


def shared_lib_name(group_name: str) -> str:
    """Given a group name, return the actual name of its extension module.

    (This just adds a suffix to the final component.)
    """
    return '{}__mypyc'.format(group_name)


def short_name(name: str) -> str:
    if name.startswith('builtins.'):
        return name[9:]
    return name


def use_fastcall(capi_version: Tuple[int, int]) -> bool:
    # We can use METH_FASTCALL for faster wrapper functions on Python 3.7+.
    return capi_version >= (3, 7)


def use_vectorcall(capi_version: Tuple[int, int]) -> bool:
    # We can use vectorcalls to make calls on Python 3.8+ (PEP 590).
    return capi_version >= (3, 8)


def use_method_vectorcall(capi_version: Tuple[int, int]) -> bool:
    # We can use a dedicated vectorcall API to call methods on Python 3.9+.
    return capi_version >= (3, 9)
