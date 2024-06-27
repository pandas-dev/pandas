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

from enum import Enum

import cpuinfo


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
    OPENHTJ2K = 36
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
    #: (more info `here <https://github.com/Blosc/blosc2_btune/>`_).
    BTUNE = 32


from .blosc2_ext import (
    DEFINED_CODECS_STOP,
    EXTENDED_HEADER_LENGTH,
    MAX_BUFFERSIZE,
    MAX_OVERHEAD,
    MAX_TYPESIZE,
    MIN_HEADER_LENGTH,
    VERSION_DATE,
    VERSION_STRING,
)

DEFINED_CODECS_STOP = DEFINED_CODECS_STOP
"""
Maximum possible Blosc2-defined codec id."""

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

cpu_info = cpuinfo.get_cpu_info()

# Public API for container module
from .core import (  # noqa: I001
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

from .ndarray import (  # noqa: I001
    NDArray,
    asarray,
    copy,
    empty,
    frombuffer,
    full,
    uninit,
    zeros,
    get_slice_nchunks,
    sin,
    cos,
    tan,
    sinh,
    cosh,
    tanh,
    arcsin,
    arccos,
    arctan,
    arctan2,
    arcsinh,
    arccosh,
    arctanh,
    exp,
    expm1,
    log,
    log10,
    log1p,
    sqrt,
    conj,
    real,
    imag,
    contains,
    abs,
)

from .lazyexpr import LazyExpr

from .schunk import SChunk, open
from .version import __version__

__version__ = __version__
"""
Python-Blosc2 version.
"""

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
# Internal Blosc threading
nthreads = ncores = detect_number_of_cores()
if ncores > 1:
    nthreads = ncores // 2
"""Number of threads to be used in compression/decompression.
"""
# Protection against too many threads
nthreads = min(nthreads, 16)
set_nthreads(nthreads)

# Defaults for compression params
cparams_dflts = {
    "codec": Codec.ZSTD,
    "codec_meta": 0,
    "clevel": 1,
    "use_dict": False,
    "typesize": 8,
    "nthreads": nthreads,
    "blocksize": 0,
    "splitmode": SplitMode.ALWAYS_SPLIT,
    "schunk": None,
    "filters": [
        Filter.NOFILTER,
        Filter.NOFILTER,
        Filter.NOFILTER,
        Filter.NOFILTER,
        Filter.NOFILTER,
        Filter.SHUFFLE,
    ],
    "filters_meta": [0, 0, 0, 0, 0, 0],
    "prefilter": None,
    "preparams": None,
    "tuner": Tuner.STUNE,
    "instr_codec": False,
}
"""
Compression params defaults.
"""

# Defaults for decompression params
dparams_dflts = {"nthreads": nthreads, "schunk": None, "postfilter": None, "postparams": None}
"""
Decompression params defaults.
"""
# Default for storage
storage_dflts = {"contiguous": False, "urlpath": None, "cparams": None, "dparams": None, "io": None}
"""
Storage params defaults. This is meant only for :ref:`SChunk <SChunk>` or :ref:`NDArray <NDArray>`.
"""

_disable_overloaded_equal = False

__all__ = [
    "__version__",
    "compress",
    "decompress",
    "cparams_dflts",
    "dparams_dflts",
    "storage_dflts",
    "set_compressor",
    "free_resources",
    "set_nthreads",
    "clib_info",
    "get_clib",
    "compressor_list",
    "set_blocksize",
    "pack",
    "unpack",
    "pack_array",
    "pack_array2",
    "save_array",
    "unpack_array",
    "unpack_array2",
    "load_array",
    "get_compressor",
    "set_releasegil",
    "detect_number_of_cores",
    "print_versions",
    "get_blocksize",
    "MAX_TYPESIZE",
    "MAX_BUFFERSIZE",
    "VERSION_STRING",
    "VERSION_DATE",
    "MIN_HEADER_LENGTH",
    "EXTENDED_HEADER_LENGTH",
    "compress2",
    "decompress2",
    "SChunk",
    "open",
    "remove_urlpath",
    "nthreads",
    "compute_chunks_blocks",
    "cpu_info",
    "get_slice_nchunks",
]
