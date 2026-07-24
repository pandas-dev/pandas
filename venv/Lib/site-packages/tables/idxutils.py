"""Utilities to be used mainly by the Index class."""

from __future__ import annotations

import math
from typing import Literal, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .index import Index


# Hints for chunk/slice/block/superblock computations:
# - The slicesize should not exceed 2**32 elements (because of
# implementation reasons).  Such an extreme case would make the
# sorting algorithms to consume up to 64 GB of memory.
# - In general, one should favor a small chunksize ( < 128 KB) if one
# wants to reduce the latency for indexed queries. However, keep in
# mind that a very low value of chunksize for big datasets may hurt
# the performance by requiring the HDF5 to use a lot of memory and CPU
# for its internal B-Tree.


def csformula(nrows: int) -> float:
    """Return the fitted chunksize (a float value) for nrows."""
    # This formula has been computed using two points:
    # 2**12 = m * 2**(n + log10(10**6))
    # 2**15 = m * 2**(n + log10(10**9))
    # where 2**12 and 2**15 are reasonable values for chunksizes for indexes
    # with 10**6 and 10**9 elements respectively.
    # Yes, return a floating point number!
    return 64 * 2 ** math.log10(nrows)


def limit_er(expectedrows: int) -> int:
    """Protection against creating too small or too large chunks or slices."""
    if expectedrows < 10**5:
        expectedrows = 10**5
    elif expectedrows > 10**12:
        expectedrows = 10**12
    return expectedrows


def computechunksize(expectedrows: int) -> int:
    """Get the optimum chunksize based on expectedrows."""
    expectedrows = limit_er(expectedrows)
    zone = int(math.log10(expectedrows))
    nrows = 10**zone
    return int(csformula(nrows))


def computeslicesize(expectedrows: int, memlevel: int) -> int:
    """Get the optimum slicesize based on expectedrows and memorylevel."""
    expectedrows = limit_er(expectedrows)
    # First, the optimum chunksize
    cs = csformula(expectedrows)
    # Now, the actual chunksize
    chunksize = computechunksize(expectedrows)
    # The optimal slicesize
    ss = int(cs * memlevel**2)
    # We *need* slicesize to be an exact multiple of the actual chunksize
    ss = (ss // chunksize) * chunksize
    ss *= 4  # slicesize should be at least divisible by 4
    # ss cannot be bigger than 2**31 - 1 elements because of fundamental
    # reasons (this limitation comes mainly from the way of compute
    # indices for indexes, but also because C keysort is not implemented
    # yet for the string type).  Besides, it cannot be larger than
    # 2**30, because limitations of the optimized binary search code
    # (in idx-opt.c, the line ``mid = lo + (hi-lo)/2;`` will overflow
    # for values of ``lo`` and ``hi`` >= 2**30).  Finally, ss must be a
    # multiple of 4, so 2**30 must definitely be an upper limit.
    if ss > 2**30:
        ss = 2**30
    return ss


def computeblocksize(
    expectedrows: int, compoundsize: int, lowercompoundsize: int
) -> int:
    """Calculate the optimum number of superblocks made from compounds blocks.

    This is useful for computing the sizes of both blocks and
    superblocks (using the PyTables terminology for blocks in indexes).

    """
    nlowerblocks = (expectedrows // lowercompoundsize) + 1
    if nlowerblocks > 2**20:
        # Protection against too large number of compound blocks
        nlowerblocks = 2**20
    size = int(lowercompoundsize * nlowerblocks)
    # We *need* superblocksize to be an exact multiple of the actual
    # compoundblock size (a ceil must be performed here!)
    size = ((size // compoundsize) + 1) * compoundsize
    return size


def calc_chunksize(
    expectedrows: int,
    optlevel: int = 6,
    indsize: int = 4,
    memlevel: int = 4,
    node: Index | None = None,
) -> tuple[int, int, int, int]:
    """Calculate the HDF5 chunk size for index and sorted arrays.

    The logic to do that is based purely in experiments playing with
    different chunksizes and compression flag. It is obvious that using
    big chunks optimizes the I/O speed, but if they are too large, the
    uncompressor takes too much time. This might (should) be further
    optimized by doing more experiments.

    """
    chunksize = computechunksize(expectedrows)
    slicesize = computeslicesize(expectedrows, memlevel)

    # Avoid excessive slicesize in Indexes,
    # see https://github.com/PyTables/PyTables/issues/879
    if node is not None:
        maxsize = (
            node._v_file.params["BUFFER_TIMES"]
            * node._v_file.params["IO_BUFFER_SIZE"]
        )
        while (slicesize * node.dtype.itemsize) > maxsize:
            slicesize = slicesize // 2

    # Correct the slicesize and the chunksize based on optlevel
    if indsize == 1:  # ultralight
        chunksize, slicesize = ccs_ultralight(optlevel, chunksize, slicesize)
    elif indsize == 2:  # light
        chunksize, slicesize = ccs_light(optlevel, chunksize, slicesize)
    elif indsize == 4:  # medium
        chunksize, slicesize = ccs_medium(optlevel, chunksize, slicesize)
    elif indsize == 8:  # full
        chunksize, slicesize = ccs_full(optlevel, chunksize, slicesize)

    # Finally, compute blocksize and superblocksize
    blocksize = computeblocksize(expectedrows, slicesize, chunksize)
    superblocksize = computeblocksize(expectedrows, blocksize, slicesize)
    # The size for different blocks information
    sizes = (superblocksize, blocksize, slicesize, chunksize)
    return sizes


def ccs_ultralight(
    optlevel: int, chunksize: int, slicesize: int
) -> tuple[int, int]:
    """Correct the slicesize and the chunksize based on optlevel."""
    if optlevel in (0, 1, 2):
        slicesize //= 2
        slicesize += optlevel * slicesize
    elif optlevel in (3, 4, 5):
        slicesize *= optlevel - 1
    elif optlevel in (6, 7, 8):
        slicesize *= optlevel - 1
    elif optlevel == 9:
        slicesize *= optlevel - 1
    return chunksize, slicesize


def ccs_light(
    optlevel: int, chunksize: int, slicesize: int
) -> tuple[int, int]:
    """Correct the slicesize and the chunksize based on optlevel."""
    if optlevel in (0, 1, 2):
        slicesize //= 2
    elif optlevel in (3, 4, 5):
        pass
    elif optlevel in (6, 7, 8):
        chunksize //= 2
    elif optlevel == 9:
        # Reducing the chunksize and enlarging the slicesize is the
        # best way to reduce the entropy with the current algorithm.
        chunksize //= 2
        slicesize *= 2
    return chunksize, slicesize


def ccs_medium(
    optlevel: int, chunksize: int, slicesize: int
) -> tuple[int, int]:
    """Correct the slicesize and the chunksize based on optlevel."""
    if optlevel in (0, 1, 2):
        slicesize //= 2
    elif optlevel in (3, 4, 5):
        pass
    elif optlevel in (6, 7, 8):
        chunksize //= 2
    elif optlevel == 9:
        # Reducing the chunksize and enlarging the slicesize is the
        # best way to reduce the entropy with the current algorithm.
        chunksize //= 2
        slicesize *= 2
    return chunksize, slicesize


def ccs_full(optlevel: int, chunksize: int, slicesize: int) -> tuple[int, int]:
    """Correct the slicesize and the chunksize based on optlevel."""
    if optlevel in (0, 1, 2):
        slicesize //= 2
    elif optlevel in (3, 4, 5):
        pass
    elif optlevel in (6, 7, 8):
        chunksize //= 2
    elif optlevel == 9:
        # Reducing the chunksize and enlarging the slicesize is the
        # best way to reduce the entropy with the current algorithm.
        chunksize //= 2
        slicesize *= 2
    return chunksize, slicesize


def calcoptlevels(
    nblocks: int, optlevel: int, indsize: int
) -> tuple[bool, bool, bool, bool]:
    """Compute the optimizations to be done.

    The calculation is based on the number of blocks, optlevel and
    indexing mode.

    """
    if indsize == 2:  # light
        return col_light(nblocks, optlevel)
    elif indsize == 4:  # medium
        return col_medium(nblocks, optlevel)
    elif indsize == 8:  # full
        return col_full(nblocks, optlevel)


def col_light(nblocks: int, optlevel: int) -> tuple[bool, bool, bool, bool]:
    """Compute the optimizations to be done for light indexes."""
    optmedian, optstarts, optstops, optfull = (False,) * 4

    if 0 < optlevel <= 3:
        optmedian = True
    elif 3 < optlevel <= 6:
        optmedian, optstarts = (True, True)
    elif 6 < optlevel <= 9:
        optmedian, optstarts, optstops = (True, True, True)

    return optmedian, optstarts, optstops, optfull


def col_medium(nblocks: int, optlevel: int) -> tuple[bool, bool, bool, bool]:
    """Compute the optimizations to be done for medium indexes."""
    optmedian, optstarts, optstops, optfull = (False,) * 4

    # Medium case
    if nblocks <= 1:
        if 0 < optlevel <= 3:
            optmedian = True
        elif 3 < optlevel <= 6:
            optmedian, optstarts = (True, True)
        elif 6 < optlevel <= 9:
            optfull = 1
    else:  # More than a block
        if 0 < optlevel <= 3:
            optfull = 1
        elif 3 < optlevel <= 6:
            optfull = 2
        elif 6 < optlevel <= 9:
            optfull = 3

    return optmedian, optstarts, optstops, optfull


def col_full(nblocks: int, optlevel: int) -> tuple[bool, bool, bool, bool]:
    """Compute the optimizations to be done for full indexes."""
    optmedian, optstarts, optstops, optfull = (False,) * 4

    # Full case
    if nblocks <= 1:
        if 0 < optlevel <= 3:
            optmedian = True
        elif 3 < optlevel <= 6:
            optmedian, optstarts = (True, True)
        elif 6 < optlevel <= 9:
            optfull = 1
    else:  # More than a block
        if 0 < optlevel <= 3:
            optfull = 1
        elif 3 < optlevel <= 6:
            optfull = 2
        elif 6 < optlevel <= 9:
            optfull = 3

    return optmedian, optstarts, optstops, optfull


def get_reduction_level(
    indsize: int, optlevel: int, slicesize: int, chunksize: int
) -> int:
    """Compute the reduction level based on indsize and optlevel."""
    rlevels = [
        [8, 8, 8, 8, 4, 4, 4, 2, 2, 1],  # 8-bit indices (ultralight)
        [4, 4, 4, 4, 2, 2, 2, 1, 1, 1],  # 16-bit indices (light)
        [2, 2, 2, 2, 1, 1, 1, 1, 1, 1],  # 32-bit indices (medium)
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],  # 64-bit indices (full)
    ]
    isizes = {1: 0, 2: 1, 4: 2, 8: 3}
    rlevel = rlevels[isizes[indsize]][optlevel]
    # The next cases should only happen in tests
    if rlevel >= slicesize:
        rlevel = 1
    if slicesize <= chunksize * rlevel:
        rlevel = 1
    if indsize == 8:
        # Ensure that, for full indexes we will never perform a reduction.
        # This is required because of implementation assumptions.
        assert rlevel == 1
    return rlevel


# Python implementations of NextAfter and NextAfterF
#
# These implementations exist because the standard function
# nextafterf is not available on Microsoft platforms.
#
# These implementations are based on the IEEE representation of
# floats and doubles.
# Author:  Shack Toms - shack@livedata.com
#
# Thanks to Shack Toms shack@livedata.com for NextAfter and NextAfterF
# implementations in Python. 2004-10-01
# epsilon  = math.ldexp(1.0, -53) # smallest double such that
#                                 # 0.5 + epsilon != 0.5
# epsilonF = math.ldexp(1.0, -24) # smallest float such that 0.5 + epsilonF
# != 0.5
# maxFloat = float(2**1024 - 2**971)  # From the IEEE 754 standard
# maxFloatF = float(2**128 - 2**104)  # From the IEEE 754 standard
# minFloat  = math.ldexp(1.0, -1022) # min positive normalized double
# minFloatF = math.ldexp(1.0, -126)  # min positive normalized float
# smallEpsilon  = math.ldexp(1.0, -1074) # smallest increment for
#                                        # doubles < minFloat
# smallEpsilonF = math.ldexp(1.0, -149)  # smallest increment for
#                                        # floats < minFloatF
infinity = math.ldexp(1.0, 1023) * 2
infinityf = math.ldexp(1.0, 128)
# Finf = float("inf")  # Infinite in the IEEE 754 standard (not avail in Win)

# A portable representation of NaN
# if sys.byteorder == "little":
#     testNaN = struct.unpack("d", '\x01\x00\x00\x00\x00\x00\xf0\x7f')[0]
# elif sys.byteorder == "big":
#     testNaN = struct.unpack("d", '\x7f\xf0\x00\x00\x00\x00\x00\x01')[0]
# else:
#     raise ValueError("Byteorder '%s' not supported!" % sys.byteorder)
# This one seems better
# testNaN = infinity - infinity

# "infinity" for several types
infinitymap = {
    "bool": [0, 1],
    "int8": [-(2**7), 2**7 - 1],
    "uint8": [0, 2**8 - 1],
    "int16": [-(2**15), 2**15 - 1],
    "uint16": [0, 2**16 - 1],
    "int32": [-(2**31), 2**31 - 1],
    "uint32": [0, 2**32 - 1],
    "int64": [-(2**63), 2**63 - 1],
    "uint64": [0, 2**64 - 1],
    "float32": [-infinityf, infinityf],
    "float64": [-infinity, infinity],
}

if hasattr(np, "float16"):
    infinitymap["float16"] = [-np.float16(np.inf), np.float16(np.inf)]
if hasattr(np, "float96"):
    infinitymap["float96"] = [-np.float96(np.inf), np.float96(np.inf)]
if hasattr(np, "float128"):
    infinitymap["float128"] = [-np.float128(np.inf), np.float128(np.inf)]

# Utility functions


def inftype(
    dtype: np.dtype, itemsize: int, sign: Literal[-1, 1] = 1
) -> bytes | float | int:
    """Return a superior limit for maximum representable data type."""
    assert sign in [-1, +1]

    if dtype.kind == "S":
        if sign < 0:
            return b"\x00" * itemsize
        else:
            return b"\xff" * itemsize
    try:
        return infinitymap[dtype.name][sign >= 0]
    except KeyError:
        raise TypeError("Type %s is not supported" % dtype.name)


def string_next_after(
    x: bytes, direction: Literal[-1, 1], itemsize: int
) -> bytes:
    """Return the next neighbor of x in the specified direction."""
    assert direction in [-1, +1]

    # Pad the string with \x00 chars until itemsize completion
    padsize = itemsize - len(x)
    if padsize > 0:
        x += b"\x00" * padsize
    # int.to_bytes is not available in Python < 3.2
    # xlist = [i.to_bytes(1, sys.byteorder) for i in x]
    xlist = [bytes([i]) for i in x]
    xlist.reverse()
    i = 0
    if direction > 0:
        if xlist == b"\xff" * itemsize:
            # Maximum value, return this
            return b"".join(xlist)
        for xchar in xlist:
            if ord(xchar) < 0xFF:
                xlist[i] = chr(ord(xchar) + 1).encode("ascii")
                break
            else:
                xlist[i] = b"\x00"
            i += 1
    else:
        if xlist == b"\x00" * itemsize:
            # Minimum value, return this
            return b"".join(xlist)
        for xchar in xlist:
            if ord(xchar) > 0x00:
                xlist[i] = chr(ord(xchar) - 1).encode("ascii")
                break
            else:
                xlist[i] = b"\xff"
            i += 1
    xlist.reverse()
    return b"".join(xlist)


def int_type_next_after(
    x: float | int, direction: Literal[-1, 1], itemsize: int
) -> int:
    """Return the next neighbor of x in the specified direction."""
    assert direction in [-1, +1]

    # x is guaranteed to be either an int or a float
    if direction < 0:
        if isinstance(x, int):
            return x - 1
        else:
            # return int(PyNextAfter(x, x - 1))
            return int(np.nextafter(x, x - 1))
    else:
        if isinstance(x, int):
            return x + 1
        else:
            # return int(PyNextAfter(x,x + 1)) + 1
            return int(np.nextafter(x, x + 1)) + 1


def bool_type_next_after(
    x: bool, direction: Literal[-1, 1], itemsize: int
) -> bool:
    """Return the next representable neighbor of x in the specified direction."""
    assert direction in [-1, +1]

    # x is guaranteed to be either a boolean
    if direction < 0:
        return False
    else:
        return True


def nextafter(
    x: bool | bytes | float | int,
    direction: Literal[-1, 0, 1],
    dtype: np.dtype,
    itemsize: int,
) -> bool | bytes | int | float:
    """Return the next representable neighbor of x in the specified direction."""
    assert direction in [-1, 0, +1]
    assert dtype.kind == "S" or type(x) in (bool, float, int)

    if direction == 0:
        return x

    if dtype.kind == "S":
        return string_next_after(x, direction, itemsize)

    if dtype.kind in ["b"]:
        return bool_type_next_after(x, direction, itemsize)
    elif dtype.kind in ["i", "u"]:
        return int_type_next_after(x, direction, itemsize)
    elif dtype.kind == "f":
        if direction < 0:
            return np.nextafter(x, x - 1)
        else:
            return np.nextafter(x, x + 1)

    # elif dtype.name == "float32":
    #    if direction < 0:
    #        return PyNextAfterF(x,x-1)
    #    else:
    #        return PyNextAfterF(x,x + 1)
    # elif dtype.name == "float64":
    #    if direction < 0:
    #        return PyNextAfter(x,x-1)
    #    else:
    #        return PyNextAfter(x,x + 1)

    raise TypeError("data type ``%s`` is not supported" % dtype)
