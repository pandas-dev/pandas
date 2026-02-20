#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

# Hey Ruff, please ignore the next violations
# ruff: noqa: E402 - Module level import not at top of file
# ruff: noqa: F401 - `var` imported but unused

import contextlib
import os
import platform
from enum import Enum

import numpy as np

_HAS_NUMBA = False
try:
    import numba

    _HAS_NUMBA = True
except ImportError:
    pass
# Do the platform check once at module level
IS_WASM = platform.machine() == "wasm32"
# IS_WASM = True  # for testing (comment this line out for production)
"""
Flag for WebAssembly platform.
"""

if not IS_WASM:
    import numexpr

from .version import __array_api_version__, __version__

__version__ = __version__
__array_api_version__ = __array_api_version__
"""
Python-Blosc2 version.
"""


class Codec(Enum):
    """
    Available codecs.
    """

    BLOSCLZ = 0
    LZ4 = 1
    LZ4HC = 2
    ZLIB = 4
    ZSTD = 5
    NDLZ = 32
    ZFP_ACC = 33
    ZFP_PREC = 34
    ZFP_RATE = 35
    #: Needs to be installed with ``pip install blosc2-openhtj2k``
    OPENHTJ2K = 36
    #: Needs to be installed with ``pip install blosc2-grok``
    GROK = 37
    #: Needs to be installed with ``pip install blosc2-openzl``
    OPENZL = 38


class Filter(Enum):
    """
    Available filters.
    """

    NOFILTER = 0
    SHUFFLE = 1
    BITSHUFFLE = 2
    DELTA = 3
    TRUNC_PREC = 4
    NDCELL = 32
    NDMEAN = 33
    BYTEDELTA = 35
    INT_TRUNC = 36


class SplitMode(Enum):
    """
    Available split modes.
    """

    ALWAYS_SPLIT = 1
    NEVER_SPLIT = 2
    AUTO_SPLIT = 3
    FORWARD_COMPAT_SPLIT = 4


class SpecialValue(Enum):
    """
    Possible special values in a chunk.
    """

    NOT_SPECIAL = 0
    ZERO = 1
    NAN = 2
    VALUE = 3
    UNINIT = 4


class Tuner(Enum):
    """
    Available tuners.
    """

    #: A 'simple' tuner. This is the default in the Blosc2 library
    STUNE = 0
    #: A more sophisticated tuner that can select different codecs/filters for different chunks
    #: (more info `here <https://github.com/Blosc/blosc2_btune/>`_); Needs to be installed with
    #: ``pip install blosc2-btune``
    BTUNE = 32


class FPAccuracy(Enum):
    """
    Floating point accuracy modes for Blosc2 computing with lazy expressions.

    This is only relevant when using floating point dtypes with miniexpr.
    """

    #: Use 1.0 ULPs (Units in the Last Place) for floating point functions
    HIGH = 1
    #: Use 3.5 ULPs (Units in the Last Place) for floating point functions
    MEDIUM = 2
    #: Use default accuracy. This is MEDIUM, which should be enough for most applications.
    DEFAULT = MEDIUM


from .blosc2_ext import (
    DEFINED_CODECS_STOP,
    EXTENDED_HEADER_LENGTH,
    GLOBAL_REGISTERED_CODECS_STOP,
    MAX_BLOCKSIZE,
    MAX_BUFFERSIZE,
    MAX_DIM,
    MAX_OVERHEAD,
    MAX_TYPESIZE,
    MIN_HEADER_LENGTH,
    USER_REGISTERED_CODECS_STOP,
    VERSION_DATE,
    VERSION_STRING,
)

DEFINED_CODECS_STOP = DEFINED_CODECS_STOP
"""
Maximum possible Blosc2-defined codec id."""

GLOBAL_REGISTERED_CODECS_STOP = GLOBAL_REGISTERED_CODECS_STOP
"""
Maximum possible Blosc2 global registered codec id."""

USER_REGISTERED_CODECS_STOP = USER_REGISTERED_CODECS_STOP
"""
Maximum possible Blosc2 user registered codec id."""

EXTENDED_HEADER_LENGTH = EXTENDED_HEADER_LENGTH
"""
Blosc2 extended header length in bytes."""

MAX_BUFFERSIZE = MAX_BUFFERSIZE
"""
Maximum buffer size in bytes for a Blosc2 chunk."""

MAX_FAST_PATH_SIZE = 2**30
"""
Maximum size in bytes for a fast path evaluation.
"""

MAX_OVERHEAD = MAX_OVERHEAD
"""
Maximum overhead during compression (in bytes). This is
equal to :py:obj:`blosc2.EXTENDED_HEADER_LENGTH <EXTENDED_HEADER_LENGTH>`."""

MAX_TYPESIZE = MAX_TYPESIZE
"""
Blosc2 maximum type size (in bytes)."""

MIN_HEADER_LENGTH = MIN_HEADER_LENGTH
"""
Blosc2 minimum header length (in bytes)."""

VERSION_DATE = VERSION_DATE
"""
The C-Blosc2 version's date."""

VERSION_STRING = VERSION_STRING
"""
The C-Blosc2 version's string."""


# For array-api compatibility
iinfo = np.iinfo
finfo = np.finfo


def isdtype(a_dtype: np.dtype, kind: str | np.dtype | tuple):
    """
    Returns a boolean indicating whether a provided dtype is of a specified data type "kind".

    Parameters
    ----------
    dtype: dtype
        The input dtype.

    kind: str | dtype | Tuple[str, dtype]
        Data type kind.

        If kind is a dtype, return boolean indicating whether the input dtype is equal to the dtype specified by kind.

        If kind is a string, return boolean indicating whether the input dtype is of a specified data type kind.
        The following dtype kinds are supporte:

            * 'bool': boolean data types (e.g., bool).

            * 'signed integer': signed integer data types (e.g., int8, int16, int32, int64).

            * 'unsigned integer': unsigned integer data types (e.g., uint8, uint16, uint32, uint64).

            * 'integral': integer data types. Shorthand for ('signed integer', 'unsigned integer').

            * 'real floating': real-valued floating-point data types (e.g., float32, float64).

            * 'complex floating': complex floating-point data types (e.g., complex64, complex128).

            * 'numeric': numeric data types. Shorthand for ('integral', 'real floating', 'complex floating').

    Returns
    -------
    out: bool
        Boolean indicating whether a provided dtype is of a specified data type kind.
    """
    kind = (kind,) if not isinstance(kind, tuple) else kind
    for _ in kind:
        if a_dtype == kind:
            return True

    _complex, _signedint, _uint, _rfloat = False, False, False, False
    if a_dtype in (complex64, complex128):
        _complex = True
        if "complex floating" in kind:
            return True
    if a_dtype == bool_ and "bool" in kind:
        return True
    if a_dtype in (int8, int16, int32, int64):
        _signedint = True
        if "signed integer" in kind:
            return True
    if a_dtype in (uint8, uint16, uint32, uint64):
        _uint = True
        if "unsigned integer" in kind:
            return True
    if a_dtype in (float16, float32, float64):
        _rfloat = True
        if "real floating" in kind:
            return True
    if "integral" in kind and (_signedint or _uint):
        return True
    return "numeric" in kind and (
        _signedint or _uint or _rfloat or _complex
    )  # checked everything, otherwise False


# dtypes for array-api
str_ = np.str_
bytes_ = np.bytes_
object_ = np.object_

from numpy import (
    bool_,
    complex64,
    complex128,
    e,
    euler_gamma,
    float16,
    float32,
    float64,
    inf,
    int8,
    int16,
    int32,
    int64,
    nan,
    newaxis,
    pi,
    uint8,
    uint16,
    uint32,
    uint64,
)

bool = bool

DEFAULT_COMPLEX = complex128
"""
Default complex floating dtype."""

DEFAULT_FLOAT = float64
"""
Default real floating dtype."""

DEFAULT_INT = int64
"""
Default integer dtype."""

DEFAULT_INDEX = int64
"""
Default indexing dtype."""


class Info:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


def __array_namespace_info__() -> Info:
    """
    Return information about the array namespace following the Array API specification.
    """

    def _raise(exc):
        raise exc

    return Info(
        capabilities=lambda: {
            "boolean indexing": True,
            "data-dependent shapes": False,
            "max dimensions": MAX_DIM,
        },
        default_device=lambda: "cpu",
        default_dtypes=lambda device=None: {
            "real floating": DEFAULT_FLOAT,
            "complex floating": DEFAULT_COMPLEX,
            "integral": DEFAULT_INT,
            "indexing": DEFAULT_INDEX,
        }
        if (device == "cpu" or device is None)
        else _raise(ValueError("Only cpu devices allowed")),
        dtypes=lambda device=None, kind=None: np.__array_namespace_info__().dtypes(kind=kind, device=device)
        if (device == "cpu" or device is None)
        else _raise(ValueError("Only cpu devices allowed")),
        devices=lambda: ["cpu"],
        name="blosc2",
        version=__version__,
    )


# Public API for container module
from .core import (
    clib_info,
    compress,
    compress2,
    compressor_list,
    compute_chunks_blocks,
    decompress,
    decompress2,
    detect_number_of_cores,
    free_resources,
    from_cframe,
    get_blocksize,
    get_cbuffer_sizes,
    get_clib,
    get_compressor,
    get_cpu_info,
    load_array,
    load_tensor,
    ndarray_from_cframe,
    pack,
    pack_array,
    pack_array2,
    pack_tensor,
    print_versions,
    register_codec,
    register_filter,
    remove_urlpath,
    save_array,
    save_tensor,
    schunk_from_cframe,
    set_blocksize,
    set_compressor,
    set_nthreads,
    set_releasegil,
    unpack,
    unpack_array,
    unpack_array2,
    unpack_tensor,
)

# Internal Blosc threading
# Get CPU info
cpu_info = get_cpu_info()
nthreads = ncores = cpu_info.get("count", 1)
"""Number of threads to be used in compression/decompression.
"""
# Protection against too many threads
nthreads = min(nthreads, 64)
# Experiments say that, when using a large number of threads, it is better to not use them all
if nthreads > 16:
    nthreads -= nthreads // 8
if not IS_WASM:
    # WASM does not support threading
    # Only call set_num_threads if within NUMEXPR_MAX_THREADS limit to avoid warning
    numexpr_max_env = os.environ.get("NUMEXPR_MAX_THREADS")
    numexpr_max: int | None = None
    if numexpr_max_env is not None:
        with contextlib.suppress(ValueError):
            numexpr_max = int(numexpr_max_env)
    if numexpr_max is None or nthreads <= numexpr_max:
        numexpr.set_num_threads(nthreads)

# This import must be before ndarray and schunk
from .storage import (  # noqa: I001
    CParams,
    cparams_dflts,
    DParams,
    dparams_dflts,
    Storage,
    storage_dflts,
)

from .ndarray import (
    Array,
    NDArray,
    NDField,
    Operand,
    are_partitions_aligned,
    are_partitions_behaved,
    arange,
    broadcast_to,
    linspace,
    eye,
    asarray,
    astype,
    indices,
    sort,
    reshape,
    copy,
    concat,
    expand_dims,
    empty,
    empty_like,
    frombuffer,
    fromiter,
    get_slice_nchunks,
    meshgrid,
    nans,
    uninit,
    zeros,
    zeros_like,
    ones,
    ones_like,
    full,
    full_like,
    save,
    stack,
)
from .embed_store import EmbedStore, estore_from_cframe
from .dict_store import DictStore
from .tree_store import TreeStore

from .c2array import c2context, C2Array, URLPath

from .lazyexpr import (
    LazyExpr,
    lazyudf,
    lazyexpr,
    LazyArray,
    LazyUDF,
    _open_lazyarray,
    get_expr_operands,
    validate_expr,
    evaluate,
    result_type,
    can_cast,
)
from .proxy import Proxy, ProxySource, ProxyNDSource, ProxyNDField, SimpleProxy, jit, as_simpleproxy

from .schunk import SChunk, open
from . import linalg
from .linalg import tensordot, vecdot, permute_dims, matrix_transpose, matmul, transpose, diagonal, outer
from .utils import linalg_funcs as linalg_funcs_list
from . import fft

# Registry for postfilters
postfilter_funcs = {}
"""
Registry for postfilter functions. For more info see
 :func:`SChunk.postfilter <blosc2.schunk.SChunk.postfilter>`"""
# Registry for prefilters
prefilter_funcs = {}
"""
Registry for prefilter functions. For more info see
 :func:`SChunk.prefilter <blosc2.schunk.SChunk.prefilter>`"""

# Registry for user-defined codecs
ucodecs_registry = {}
"""
Registry for user-defined codecs. For more info see
 :func:`blosc2.register_codec <blosc2.register_codec>`"""
# Registry for user-defined filters
ufilters_registry = {}
"""
Registry for user-defined filters. For more info see
 :func:`blosc2.register_filter <blosc2.register_filter>`"""

blosclib_version = f"{VERSION_STRING} ({VERSION_DATE})"
"""
The blosc2 version + date.
"""

# Private global variables
_disable_overloaded_equal = False
"""
Disable the overloaded equal operator.
"""

# Delayed imports for avoiding overwriting of python builtins
from .ndarray import (
    abs,
    acos,
    acosh,
    add,
    all,
    any,
    arccos,
    arccosh,
    arcsin,
    arcsinh,
    arctan,
    arctan2,
    arctanh,
    argmax,
    argmin,
    array_from_ffi_ptr,
    asin,
    asinh,
    atan,
    atan2,
    atanh,
    bitwise_and,
    bitwise_invert,
    bitwise_left_shift,
    bitwise_or,
    bitwise_right_shift,
    bitwise_xor,
    ceil,
    clip,
    conj,
    contains,
    copysign,
    cos,
    cosh,
    count_nonzero,
    divide,
    equal,
    exp,
    expm1,
    floor,
    floor_divide,
    greater,
    greater_equal,
    hypot,
    imag,
    isfinite,
    isinf,
    isnan,
    lazywhere,
    less,
    less_equal,
    log,
    log1p,
    log2,
    log10,
    logaddexp,
    logical_and,
    logical_not,
    logical_or,
    logical_xor,
    max,
    maximum,
    mean,
    min,
    minimum,
    multiply,
    negative,
    nextafter,
    not_equal,
    positive,
    pow,
    prod,
    real,
    reciprocal,
    remainder,
    round,
    sign,
    signbit,
    sin,
    sinh,
    sqrt,
    square,
    squeeze,
    std,
    subtract,
    sum,
    take,
    take_along_axis,
    tan,
    tanh,
    trunc,
    var,
    where,
)

__all__ = [  # noqa : RUF022
    # Constants
    "EXTENDED_HEADER_LENGTH",
    "MAX_BUFFERSIZE",
    "MAX_TYPESIZE",
    "MIN_HEADER_LENGTH",
    "VERSION_DATE",
    "VERSION_STRING",
    # Default dtypes
    "DEFAULT_COMPLEX",
    "DEFAULT_FLOAT",
    "DEFAULT_INDEX",
    "DEFAULT_INT",
    # Mathematical constants
    "e",
    "pi",
    "inf",
    "nan",
    "newaxis",
    # Classes
    "C2Array",
    "CParams",
    # Enums
    "Codec",
    "DParams",
    "DictStore",
    "EmbedStore",
    "Filter",
    "LazyArray",
    "LazyExpr",
    "LazyUDF",
    "NDArray",
    "NDField",
    "Operand",
    "Proxy",
    "ProxyNDField",
    "ProxyNDSource",
    "ProxySource",
    "SChunk",
    "SimpleProxy",
    "SpecialValue",
    "SplitMode",
    "Storage",
    "TreeStore",
    "Tuner",
    "URLPath",
    # Version
    "__version__",
    # Utils
    "linalg_funcs_list",
    # Functions
    "abs",
    "acos",
    "acosh",
    "add",
    "all",
    "any",
    "arange",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctan2",
    "arctanh",
    "are_partitions_aligned",
    "are_partitions_behaved",
    "argmax",
    "argmin",
    "array_from_ffi_ptr",
    "asarray",
    "asin",
    "asinh",
    "as_simpleproxy",
    "astype",
    "atan",
    "atan2",
    "atanh",
    "bitwise_and",
    "bitwise_invert",
    "bitwise_left_shift",
    "bitwise_or",
    "bitwise_right_shift",
    "bitwise_xor",
    "broadcast_to",
    "can_cast",
    "ceil",
    "clib_info",
    "clip",
    "compress",
    "compress2",
    "compressor_list",
    "compute_chunks_blocks",
    "concat",
    "conj",
    "contains",
    "copy",
    "copysign",
    "cos",
    "cosh",
    "count_nonzero",
    "cparams_dflts",
    "cpu_info",
    "decompress",
    "decompress2",
    "detect_number_of_cores",
    "divide",
    "dparams_dflts",
    "empty",
    "empty_like",
    "equal",
    "estore_from_cframe",
    "exp",
    "expand_dims",
    "expm1",
    "eye",
    "finfo",
    "floor",
    "floor_divide",
    "free_resources",
    "from_cframe",
    "frombuffer",
    "fromiter",
    "full",
    "full_like",
    "get_blocksize",
    "get_cbuffer_sizes",
    "get_clib",
    "get_compressor",
    "get_cpu_info",
    "get_expr_operands",
    "get_slice_nchunks",
    "greater",
    "greater_equal",
    "hypot",
    "imag",
    "iinfo",
    "indices",
    "isdtype",
    "isfinite",
    "isinf",
    "isnan",
    "jit",
    "lazyexpr",
    "lazyudf",
    "lazywhere",
    "less",
    "less_equal",
    "linspace",
    "load_array",
    "load_tensor",
    "log",
    "log1p",
    "log2",
    "log10",
    "logaddexp",
    "logical_and",
    "logical_not",
    "logical_or",
    "logical_xor",
    "matmul",
    "matrix_transpose",
    "max",
    "maximum",
    "mean",
    "meshgrid",
    "min",
    "minimum",
    "multiply",
    "nans",
    "ndarray_from_cframe",
    "negative",
    "nextafter",
    "not_equal",
    "ones",
    "ones_like",
    "open",
    "pack",
    "pack_array",
    "pack_array2",
    "pack_tensor",
    "permute_dims",
    "positive",
    "postfilter_funcs",
    "pow",
    "prefilter_funcs",
    "print_versions",
    "prod",
    "real",
    "reciprocal",
    "register_codec",
    "register_filter",
    "remainder",
    "remove_urlpath",
    "reshape",
    "result_type",
    "round",
    "save",
    "save_array",
    "save_tensor",
    "schunk_from_cframe",
    "set_blocksize",
    "set_compressor",
    "set_nthreads",
    "set_releasegil",
    "sign",
    "signbit",
    "sin",
    "sinh",
    "sort",
    "sqrt",
    "square",
    "squeeze",
    "stack",
    "std",
    "storage_dflts",
    "subtract",
    "sum",
    "take",
    "take_along_axis",
    "tan",
    "tanh",
    "tensordot",
    "transpose",
    "trunc",
    "uninit",
    "unpack",
    "unpack_array",
    "unpack_array2",
    "unpack_tensor",
    "validate_expr",
    "var",
    "vecdot",
    "where",
    "zeros",
    "zeros_like",
]
