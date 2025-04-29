#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

# Hey Ruff, please ignore the next violations
# ruff: noqa: E402 - Module level import not at top of file
# ruff: noqa: F401 - `var` imported but unused

import platform
from enum import Enum

# Do the platform check once at module level
IS_WASM = platform.machine() == "wasm32"
# IS_WASM = True  # for testing (comment this line out for production)
"""
Flag for WebAssembly platform.
"""

if not IS_WASM:
    import numexpr

from .version import __version__

__version__ = __version__
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


from .blosc2_ext import (
    DEFINED_CODECS_STOP,
    EXTENDED_HEADER_LENGTH,
    GLOBAL_REGISTERED_CODECS_STOP,
    MAX_BLOCKSIZE,
    MAX_BUFFERSIZE,
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
    NDArray,
    NDField,
    Operand,
    are_partitions_aligned,
    are_partitions_behaved,
    arange,
    linspace,
    eye,
    asarray,
    indices,
    sort,
    reshape,
    copy,
    empty,
    frombuffer,
    fromiter,
    get_slice_nchunks,
    nans,
    uninit,
    zeros,
    ones,
    full,
    save,
    matmul,
)

from .c2array import c2context, C2Array, URLPath

from .lazyexpr import (
    LazyExpr,
    lazyudf,
    lazyexpr,
    LazyArray,
    _open_lazyarray,
    get_expr_operands,
    validate_expr,
    evaluate,
)
from .proxy import Proxy, ProxySource, ProxyNDSource, ProxyNDField, SimpleProxy, jit

from .schunk import SChunk, open


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
    all,
    any,
    arccos,
    arccosh,
    arcsin,
    arcsinh,
    arctan,
    arctan2,
    arctanh,
    conj,
    contains,
    cos,
    cosh,
    exp,
    expm1,
    imag,
    lazywhere,
    log,
    log1p,
    log10,
    max,
    mean,
    min,
    prod,
    real,
    sin,
    sinh,
    sqrt,
    std,
    sum,
    tan,
    tanh,
    var,
    where,
)

__all__ = [
    "EXTENDED_HEADER_LENGTH",
    "MAX_BUFFERSIZE",
    "MAX_TYPESIZE",
    "MIN_HEADER_LENGTH",
    "VERSION_DATE",
    "VERSION_STRING",
    "CParams",
    "DParams",
    "SChunk",
    "Storage",
    "__version__",
    "clib_info",
    "compress",
    "compress2",
    "compressor_list",
    "compute_chunks_blocks",
    "cparams_dflts",
    "cpu_info",
    "decompress",
    "decompress2",
    "detect_number_of_cores",
    "dparams_dflts",
    "free_resources",
    "get_blocksize",
    "get_clib",
    "get_compressor",
    "get_slice_nchunks",
    "lazyexpr",
    "lazyudf",
    "lazywhere",
    "load_array",
    "nthreads",
    "open",
    "pack",
    "pack_array",
    "pack_array2",
    "print_versions",
    "remove_urlpath",
    "save_array",
    "set_blocksize",
    "set_compressor",
    "set_nthreads",
    "set_releasegil",
    "storage_dflts",
    "unpack",
    "unpack_array",
    "unpack_array2",
]
