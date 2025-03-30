"""
Hash implementations for Numba types
"""

import math
import numpy as np
import sys
import ctypes
import warnings
from collections import namedtuple

import llvmlite.binding as ll
from llvmlite import ir

from numba import literal_unroll
from numba.core.extending import (
    overload, overload_method, intrinsic, register_jitable)
from numba.core import errors
from numba.core import types
from numba.core.unsafe.bytes import grab_byte, grab_uint64_t
from numba.cpython.randomimpl import (const_int, get_next_int, get_next_int32,
                                      get_state_ptr)

# This is Py_hash_t, which is a Py_ssize_t, which has sizeof(size_t):
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Include/pyport.h#L91-L96    # noqa: E501
_hash_width = sys.hash_info.width
_Py_hash_t = getattr(types, 'py_int')
_Py_uhash_t = getattr(types, 'py_int')

# Constants from CPython source, obtained by various means:
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Include/pyhash.h    # noqa: E501
_PyHASH_INF = sys.hash_info.inf
_PyHASH_NAN = sys.hash_info.nan
_PyHASH_MODULUS = _Py_uhash_t(sys.hash_info.modulus)
_PyHASH_BITS = 31 if types.py_int.bitwidth == 32 else 61  # mersenne primes
_PyHASH_MULTIPLIER = 0xf4243  # 1000003UL
_PyHASH_IMAG = _PyHASH_MULTIPLIER
_PyLong_SHIFT = sys.int_info.bits_per_digit
_Py_HASH_CUTOFF = sys.hash_info.cutoff
_Py_hashfunc_name = sys.hash_info.algorithm


# This stub/overload pair are used to force branch pruning to remove the dead
# branch based on the potential `None` type of the hash_func which works better
# if the predicate for the prune in an ir.Arg. The obj is an arg to allow for
# a custom error message.
def _defer_hash(hash_func):
    pass


@overload(_defer_hash)
def ol_defer_hash(obj, hash_func):
    err_msg = f"unhashable type: '{obj}'"

    def impl(obj, hash_func):
        if hash_func is None:
            raise TypeError(err_msg)
        else:
            return hash_func()
    return impl


# hash(obj) is implemented by calling obj.__hash__()
@overload(hash)
def hash_overload(obj):
    attempt_generic_msg = ("No __hash__ is defined for object of type "
                           f"'{obj}' and a generic hash() cannot be "
                           "performed as there is no suitable object "
                           "represention in Numba compiled code!")

    def impl(obj):
        if hasattr(obj, '__hash__'):
            return _defer_hash(obj, getattr(obj, '__hash__'))
        else:
            raise TypeError(attempt_generic_msg)
    return impl


@register_jitable
def process_return(val):
    asint = _Py_hash_t(val)
    if (asint == int(-1)):
        asint = int(-2)
    return asint


# This is a translation of CPython's _Py_HashDouble:
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Python/pyhash.c#L34-L129   # noqa: E501
# NOTE: In Python 3.10 hash of nan is now hash of the pointer to the PyObject
# containing said nan. Numba cannot replicate this as there is no object, so it
# elects to replicate the behaviour i.e. hash of nan is something "unique" which
# satisfies https://bugs.python.org/issue43475.


@register_jitable(locals={'x': _Py_uhash_t,
                          'y': _Py_uhash_t,
                          'm': types.np_double,
                          'e': types.np_intc,
                          'sign': types.np_intc,
                          '_PyHASH_MODULUS': _Py_uhash_t,
                          '_PyHASH_BITS': types.np_intc
                          })
def _Py_HashDouble(v):
    if not np.isfinite(v):
        if (np.isinf(v)):
            if (v > 0):
                return _PyHASH_INF
            else:
                return -_PyHASH_INF
        else:
            # Python 3.10 does not use `_PyHASH_NAN`.
            # https://github.com/python/cpython/blob/2c4792264f9218692a1bd87398a60591f756b171/Python/pyhash.c#L102   # noqa: E501
            # Numba returns a pseudo-random number to reflect the spirit of the
            # change.
            x = _prng_random_hash()
            return process_return(x)

    m, e = math.frexp(v)

    sign = 1
    if (m < 0):
        sign = -1
        m = -m

    # process 28 bits at a time;  this should work well both for binary
    #  and hexadecimal floating point.
    x = 0
    while (m):
        x = ((x << 28) & _PyHASH_MODULUS) | x >> (_PyHASH_BITS - 28)
        m *= 268435456.0  # /* 2**28 */
        e -= 28
        y = int(m)  # /* pull out integer part */
        m -= y
        x += y
        if x >= _PyHASH_MODULUS:
            x -= _PyHASH_MODULUS
    # /* adjust for the exponent;  first reduce it modulo _PyHASH_BITS */
    if e >= 0:
        e = e % _PyHASH_BITS
    else:
        e = _PyHASH_BITS - 1 - ((-1 - e) % _PyHASH_BITS)

    x = ((x << e) & _PyHASH_MODULUS) | x >> (_PyHASH_BITS - e)

    x = x * sign
    return process_return(x)


@intrinsic
def _fpext(tyctx, val):
    def impl(cgctx, builder, signature, args):
        val = args[0]
        return builder.fpext(val, ir.DoubleType())
    sig = types.c_float64(types.c_float32)
    return sig, impl


@intrinsic
def _prng_random_hash(tyctx):

    def impl(cgctx, builder, signature, args):
        state_ptr = get_state_ptr(cgctx, builder, "internal")
        bits = const_int(_hash_width)

        # Why not just use get_next_int() with the correct bitwidth?
        # get_next_int() always returns an i64, because the bitwidth it is
        # passed may not be a compile-time constant, so it needs to allocate
        # the largest unit of storage that may be required. Therefore, if the
        # hash width is 32, then we need to use get_next_int32() to ensure we
        # don't return a wider-than-expected hash, even if everything above
        # the low 32 bits would have been zero.
        if _hash_width == 32:
            value = get_next_int32(cgctx, builder, state_ptr)
        else:
            value = get_next_int(cgctx, builder, state_ptr, bits, False)

        return value

    sig = _Py_hash_t()
    return sig, impl


# This is a translation of CPython's long_hash, but restricted to the numerical
# domain reachable by int64/uint64 (i.e. no BigInt like support):
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Objects/longobject.c#L2934-L2989    # noqa: E501
# obdigit is a uint32_t which is typedef'd to digit
# int32_t is typedef'd to sdigit


@register_jitable(locals={'x': _Py_uhash_t,
                          'p1': _Py_uhash_t,
                          'p2': _Py_uhash_t,
                          'p3': _Py_uhash_t,
                          'p4': _Py_uhash_t,
                          '_PyHASH_MODULUS': _Py_uhash_t,
                          '_PyHASH_BITS': types.c_int32,
                          '_PyLong_SHIFT': types.c_int32,})
def _long_impl(val):
    # This function assumes val came from a long int repr with val being a
    # uint64_t this means having to split the input into PyLong_SHIFT size
    # chunks in an unsigned hash wide type, max numba can handle is a 64bit int

    # mask to select low _PyLong_SHIFT bits
    _tmp_shift = 32 - _PyLong_SHIFT
    mask_shift = (~types.c_uint32(0x0)) >> _tmp_shift

    # a 64bit wide max means Numba only needs 3 x 30 bit values max,
    # or 5 x 15 bit values max on 32bit platforms
    i = (64 // _PyLong_SHIFT) + 1

    # alg as per hash_long
    x = 0
    p3 = (_PyHASH_BITS - _PyLong_SHIFT)
    for idx in range(i - 1, -1, -1):
        p1 = x << _PyLong_SHIFT
        p2 = p1 & _PyHASH_MODULUS
        p4 = x >> p3
        x = p2 | p4
        # the shift and mask splits out the `ob_digit` parts of a Long repr
        x += types.c_uint32((val >> idx * _PyLong_SHIFT) & mask_shift)
        if x >= _PyHASH_MODULUS:
            x -= _PyHASH_MODULUS
    return _Py_hash_t(x)


# This has no CPython equivalent, CPython uses long_hash.
@overload_method(types.Integer, '__hash__')
@overload_method(types.Boolean, '__hash__')
def int_hash(val):

    _HASH_I64_MIN = -2 if sys.maxsize <= 2 ** 32 else -4
    _SIGNED_MIN = types.c_int64(-0x8000000000000000)

    # Find a suitable type to hold a "big" value, i.e. iinfo(ty).min/max
    # this is to ensure e.g. int32.min is handled ok as it's abs() is its value
    _BIG = types.c_int64 if getattr(val, 'signed', False) else types.c_uint64

    # this is a bit involved due to the CPython repr of ints
    def impl(val):
        # If the magnitude is under PyHASH_MODULUS, just return the
        # value val as the hash, couple of special cases if val == val:
        # 1. it's 0, in which case return 0
        # 2. it's signed int minimum value, return the value CPython computes
        # but Numba cannot as there's no type wide enough to hold the shifts.
        #
        # If the magnitude is greater than PyHASH_MODULUS then... if the value
        # is negative then negate it switch the sign on the hash once computed
        # and use the standard wide unsigned hash implementation
        val = _BIG(val)
        mag = abs(val)
        if mag < _PyHASH_MODULUS:
            if val == 0:
                ret = 0
            elif val == _SIGNED_MIN:  # e.g. int64 min, -0x8000000000000000
                ret = _Py_hash_t(_HASH_I64_MIN)
            else:
                ret = _Py_hash_t(val)
        else:
            needs_negate = False
            if val < 0:
                val = -val
                needs_negate = True
            ret = _long_impl(val)
            if needs_negate:
                ret = -ret
        return process_return(ret)
    return impl

# This is a translation of CPython's float_hash:
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Objects/floatobject.c#L528-L532    # noqa: E501


@overload_method(types.Float, '__hash__')
def float_hash(val):
    if val.bitwidth == 64:
        def impl(val):
            hashed = _Py_HashDouble(val)
            return hashed
    else:
        def impl(val):
            # widen the 32bit float to 64bit
            fpextended = np.float64(_fpext(val))
            hashed = _Py_HashDouble(fpextended)
            return hashed
    return impl

# This is a translation of CPython's complex_hash:
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Objects/complexobject.c#L408-L428    # noqa: E501


@overload_method(types.Complex, '__hash__')
def complex_hash(val):
    def impl(val):
        hashreal = hash(val.real)
        hashimag = hash(val.imag)
        # Note:  if the imaginary part is 0, hashimag is 0 now,
        # so the following returns hashreal unchanged.  This is
        # important because numbers of different types that
        # compare equal must have the same hash value, so that
        # hash(x + 0*j) must equal hash(x).
        combined = hashreal + _PyHASH_IMAG * hashimag
        return process_return(combined)
    return impl


# Python 3.8 strengthened its hash alg for tuples.
# This is a translation of CPython's tuplehash for Python >=3.8
# https://github.com/python/cpython/blob/b738237d6792acba85b1f6e6c8993a812c7fd815/Objects/tupleobject.c#L338-L391    # noqa: E501

# These consts are needed for this alg variant, they are from:
# https://github.com/python/cpython/blob/b738237d6792acba85b1f6e6c8993a812c7fd815/Objects/tupleobject.c#L353-L363    # noqa: E501

if _Py_uhash_t.bitwidth // 8 > 4:
    _PyHASH_XXPRIME_1 = _Py_uhash_t(11400714785074694791)
    _PyHASH_XXPRIME_2 = _Py_uhash_t(14029467366897019727)
    _PyHASH_XXPRIME_5 = _Py_uhash_t(2870177450012600261)

    @register_jitable(locals={'x': types.np_uint64})
    def _PyHASH_XXROTATE(x):
        # Rotate left 31 bits
        return ((x << types.np_uint64(31)) | (x >> types.np_uint64(33)))
else:
    _PyHASH_XXPRIME_1 = _Py_uhash_t(2654435761)
    _PyHASH_XXPRIME_2 = _Py_uhash_t(2246822519)
    _PyHASH_XXPRIME_5 = _Py_uhash_t(374761393)

    @register_jitable(locals={'x': types.np_uint64})
    def _PyHASH_XXROTATE(x):
        # Rotate left 13 bits
        return ((x << types.np_uint64(13)) | (x >> types.np_uint64(19)))


@register_jitable(locals={'acc': _Py_uhash_t, 'lane': _Py_uhash_t,
                          '_PyHASH_XXPRIME_5': _Py_uhash_t,
                          '_PyHASH_XXPRIME_1': _Py_uhash_t,
                          'tl': _Py_uhash_t})
def _tuple_hash(tup):
    tl = len(tup)
    acc = _PyHASH_XXPRIME_5
    for x in literal_unroll(tup):
        lane = hash(x)
        if lane == _Py_uhash_t(-1):
            return -1
        acc += lane * _PyHASH_XXPRIME_2
        acc = _PyHASH_XXROTATE(acc)
        acc *= _PyHASH_XXPRIME_1

    acc += tl ^ (_PyHASH_XXPRIME_5 ^ _Py_uhash_t(3527539))

    if acc == _Py_uhash_t(-1):
        return process_return(1546275796)

    return process_return(acc)


@overload_method(types.BaseTuple, '__hash__')
def tuple_hash(val):
    def impl(val):
        return _tuple_hash(val)
    return impl


# ------------------------------------------------------------------------------
# String/bytes hashing needs hashseed info, this is from:
# https://stackoverflow.com/a/41088757
# with thanks to Martijn Pieters
#
# Developer note:
# CPython makes use of an internal "hashsecret" which is essentially a struct
# containing some state that is set on CPython initialization and contains magic
# numbers used particularly in unicode/string hashing. This code binds to the
# Python runtime libraries in use by the current process and reads the
# "hashsecret" state so that it can be used by Numba. As this is done at runtime
# the behaviour and influence of the PYTHONHASHSEED environment variable is
# accommodated.

from ctypes import (  # noqa
    c_size_t,
    c_ubyte,
    c_uint64,
    pythonapi,
    Structure,
    Union,
)  # noqa


class FNV(Structure):
    _fields_ = [
        ('prefix', c_size_t),
        ('suffix', c_size_t)
    ]


class SIPHASH(Structure):
    _fields_ = [
        ('k0', c_uint64),
        ('k1', c_uint64),
    ]


class DJBX33A(Structure):
    _fields_ = [
        ('padding', c_ubyte * 16),
        ('suffix', c_size_t),
    ]


class EXPAT(Structure):
    _fields_ = [
        ('padding', c_ubyte * 16),
        ('hashsalt', c_size_t),
    ]


class _Py_HashSecret_t(Union):
    _fields_ = [
        # ensure 24 bytes
        ('uc', c_ubyte * 24),
        # two Py_hash_t for FNV
        ('fnv', FNV),
        # two uint64 for SipHash24
        ('siphash', SIPHASH),
        # a different (!) Py_hash_t for small string optimization
        ('djbx33a', DJBX33A),
        ('expat', EXPAT),
    ]


_hashsecret_entry = namedtuple('_hashsecret_entry', ['symbol', 'value'])


# Only a few members are needed at present
def _build_hashsecret():
    """Read hash secret from the Python process

    Returns
    -------
    info : dict
        - keys are "djbx33a_suffix", "siphash_k0", siphash_k1".
        - values are the namedtuple[symbol:str, value:int]
    """
    # Read hashsecret and inject it into the LLVM symbol map under the
    # prefix `_numba_hashsecret_`.
    pyhashsecret = _Py_HashSecret_t.in_dll(pythonapi, '_Py_HashSecret')
    info = {}

    def inject(name, val):
        symbol_name = "_numba_hashsecret_{}".format(name)
        val = ctypes.c_uint64(val)
        addr = ctypes.addressof(val)
        ll.add_symbol(symbol_name, addr)
        info[name] = _hashsecret_entry(symbol=symbol_name, value=val)

    inject('djbx33a_suffix', pyhashsecret.djbx33a.suffix)
    inject('siphash_k0', pyhashsecret.siphash.k0)
    inject('siphash_k1', pyhashsecret.siphash.k1)
    return info


_hashsecret = _build_hashsecret()


# ------------------------------------------------------------------------------


if _Py_hashfunc_name in ('siphash13', 'siphash24', 'fnv'):

    # Check for use of the FNV hashing alg, warn users that it's not implemented
    # and functionality relying of properties derived from hashing will be fine
    # but hash values themselves are likely to be different.
    if _Py_hashfunc_name == 'fnv':
        msg = ("FNV hashing is not implemented in Numba. See PEP 456 "
               "https://www.python.org/dev/peps/pep-0456/ "
               "for rationale over not using FNV. Numba will continue to work, "
               "but hashes for built in types will be computed using "
               "siphash24. This will permit e.g. dictionaries to continue to "
               "behave as expected, however anything relying on the value of "
               "the hash opposed to hash as a derived property is likely to "
               "not work as expected.")
        warnings.warn(msg)

    # This is a translation of CPython's siphash24 function:
    # https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Python/pyhash.c#L287-L413    # noqa: E501
    # and also, since Py 3.11, a translation of CPython's siphash13 function:
    # https://github.com/python/cpython/blob/9dda9020abcf0d51d59b283a89c58c8e1fb0f574/Python/pyhash.c#L376-L424
    # the only differences are in the use of SINGLE_ROUND in siphash13 vs.
    # DOUBLE_ROUND in siphash24, and that siphash13 has an extra "ROUND" applied
    # just before the final XORing of components to create the return value.

    # /* *********************************************************************
    # <MIT License>
    # Copyright (c) 2013  Marek Majkowski <marek@popcount.org>

    # Permission is hereby granted, free of charge, to any person obtaining a
    # copy of this software and associated documentation files (the "Software"),
    # to deal in the Software without restriction, including without limitation
    # the rights to use, copy, modify, merge, publish, distribute, sublicense,
    # and/or sell copies of the Software, and to permit persons to whom the
    # Software is furnished to do so, subject to the following conditions:

    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.

    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    # THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    # FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    # DEALINGS IN THE SOFTWARE.
    # </MIT License>

    # Original location:
    # https://github.com/majek/csiphash/

    # Solution inspired by code from:
    # Samuel Neves (supercop/crypto_auth/siphash24/little)
    #djb (supercop/crypto_auth/siphash24/little2)
    # Jean-Philippe Aumasson (https://131002.net/siphash/siphash24.c)

    # Modified for Python by Christian Heimes:
    # - C89 / MSVC compatibility
    # - _rotl64() on Windows
    # - letoh64() fallback
    # */

    @register_jitable(locals={'x': types.np_uint64,
                              'b': types.np_uint64, })
    def _ROTATE(x, b):
        return types.np_uint64(((x) << (b)) |
                               ((x) >> (types.np_uint64(64) - (b))))

    @register_jitable(locals={'a': types.np_uint64,
                              'b': types.np_uint64,
                              'c': types.np_uint64,
                              'd': types.np_uint64,
                              's': types.np_uint64,
                              't': types.np_uint64, })
    def _HALF_ROUND(a, b, c, d, s, t):
        a += b
        c += d
        b = _ROTATE(b, s) ^ a
        d = _ROTATE(d, t) ^ c
        a = _ROTATE(a, 32)
        return a, b, c, d

    @register_jitable(locals={'v0': types.np_uint64,
                              'v1': types.np_uint64,
                              'v2': types.np_uint64,
                              'v3': types.np_uint64, })
    def _SINGLE_ROUND(v0, v1, v2, v3):
        v0, v1, v2, v3 = _HALF_ROUND(v0, v1, v2, v3, 13, 16)
        v2, v1, v0, v3 = _HALF_ROUND(v2, v1, v0, v3, 17, 21)
        return v0, v1, v2, v3

    @register_jitable(locals={'v0': types.np_uint64,
                              'v1': types.np_uint64,
                              'v2': types.np_uint64,
                              'v3': types.np_uint64, })
    def _DOUBLE_ROUND(v0, v1, v2, v3):
        v0, v1, v2, v3 = _SINGLE_ROUND(v0, v1, v2, v3)
        v0, v1, v2, v3 = _SINGLE_ROUND(v0, v1, v2, v3)
        return v0, v1, v2, v3

    def _gen_siphash(alg):
        if alg == 'siphash13':
            _ROUNDER = _SINGLE_ROUND
            _EXTRA_ROUND = True
        elif alg == 'siphash24':
            _ROUNDER = _DOUBLE_ROUND
            _EXTRA_ROUND = False
        else:
            assert 0, 'unreachable'

        @register_jitable(locals={'v0': types.np_uint64,
                                  'v1': types.np_uint64,
                                  'v2': types.np_uint64,
                                  'v3': types.np_uint64,
                                  'b': types.np_uint64,
                                  'mi': types.np_uint64,
                                  't': types.np_uint64,
                                  'mask': types.np_uint64,
                                  'jmp': types.np_uint64,
                                  'ohexefef': types.np_uint64})
        def _siphash(k0, k1, src, src_sz):
            b = types.np_uint64(src_sz) << 56
            v0 = k0 ^ types.np_uint64(0x736f6d6570736575)
            v1 = k1 ^ types.np_uint64(0x646f72616e646f6d)
            v2 = k0 ^ types.np_uint64(0x6c7967656e657261)
            v3 = k1 ^ types.np_uint64(0x7465646279746573)

            idx = 0
            while (src_sz >= 8):
                mi = grab_uint64_t(src, idx)
                idx += 1
                src_sz -= 8
                v3 ^= mi
                v0, v1, v2, v3 = _ROUNDER(v0, v1, v2, v3)
                v0 ^= mi

            # this is the switch fallthrough:
            # https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Python/pyhash.c#L390-L400    # noqa: E501
            t = types.np_uint64(0x0)
            boffset = idx * 8
            ohexefef = types.np_uint64(0xff)
            if src_sz >= 7:
                jmp = (6 * 8)
                mask = ~types.np_uint64(ohexefef << jmp)
                t = (t & mask) | (types.np_uint64(grab_byte(src, boffset + 6))
                                  << jmp)
            if src_sz >= 6:
                jmp = (5 * 8)
                mask = ~types.np_uint64(ohexefef << jmp)
                t = (t & mask) | (types.np_uint64(grab_byte(src, boffset + 5))
                                  << jmp)
            if src_sz >= 5:
                jmp = (4 * 8)
                mask = ~types.np_uint64(ohexefef << jmp)
                t = (t & mask) | (types.np_uint64(grab_byte(src, boffset + 4))
                                  << jmp)
            if src_sz >= 4:
                t &= types.np_uint64(0xffffffff00000000)
                for i in range(4):
                    jmp = i * 8
                    mask = ~types.np_uint64(ohexefef << jmp)
                    t = (t & mask) | \
                        (types.np_uint64(grab_byte(src, boffset + i)) << jmp)
            if src_sz >= 3:
                jmp = (2 * 8)
                mask = ~types.np_uint64(ohexefef << jmp)
                t = (t & mask) | (types.np_uint64(grab_byte(src, boffset + 2))
                                  << jmp)
            if src_sz >= 2:
                jmp = (1 * 8)
                mask = ~types.np_uint64(ohexefef << jmp)
                t = (t & mask) | (types.np_uint64(grab_byte(src, boffset + 1))
                                  << jmp)
            if src_sz >= 1:
                mask = ~(ohexefef)
                t = (t & mask) | (types.np_uint64(grab_byte(src, boffset + 0)))

            b |= t
            v3 ^= b
            v0, v1, v2, v3 = _ROUNDER(v0, v1, v2, v3)
            v0 ^= b
            v2 ^= ohexefef
            v0, v1, v2, v3 = _ROUNDER(v0, v1, v2, v3)
            v0, v1, v2, v3 = _ROUNDER(v0, v1, v2, v3)
            if _EXTRA_ROUND:
                v0, v1, v2, v3 = _ROUNDER(v0, v1, v2, v3)
            t = (v0 ^ v1) ^ (v2 ^ v3)
            return t

        return _siphash

    _siphash13 = _gen_siphash('siphash13')
    _siphash24 = _gen_siphash('siphash24')

    _siphasher = _siphash13 if _Py_hashfunc_name == 'siphash13' else _siphash24

else:
    msg = "Unsupported hashing algorithm in use %s" % _Py_hashfunc_name
    raise ValueError(msg)


@intrinsic
def _inject_hashsecret_read(tyctx, name):
    """Emit code to load the hashsecret.
    """
    if not isinstance(name, types.StringLiteral):
        raise errors.TypingError("requires literal string")

    sym = _hashsecret[name.literal_value].symbol
    resty = types.np_uint64
    sig = resty(name)

    def impl(cgctx, builder, sig, args):
        mod = builder.module
        try:
            # Search for existing global
            gv = mod.get_global(sym)
        except KeyError:
            # Inject the symbol if not already exist.
            gv = ir.GlobalVariable(mod, ir.IntType(64), name=sym)
        v = builder.load(gv)
        return v

    return sig, impl


def _load_hashsecret(name):
    return _hashsecret[name].value


@overload(_load_hashsecret)
def _impl_load_hashsecret(name):
    def imp(name):
        return _inject_hashsecret_read(name)
    return imp


# This is a translation of CPythons's _Py_HashBytes:
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Python/pyhash.c#L145-L191    # noqa: E501


@register_jitable(locals={'_hash': _Py_uhash_t})
def _Py_HashBytes(val, _len):
    if (_len == 0):
        return process_return(0)

    if (_len < _Py_HASH_CUTOFF):
        # TODO: this branch needs testing, needs a CPython setup for it!
        # /* Optimize hashing of very small strings with inline DJBX33A. */
        _hash = _Py_uhash_t(5381)  # /* DJBX33A starts with 5381 */
        for idx in range(_len):
            _hash = ((_hash << 5) + _hash) + np.uint8(grab_byte(val, idx))

        _hash ^= _len
        _hash ^= _load_hashsecret('djbx33a_suffix')
    else:
        tmp = _siphasher(types.np_uint64(_load_hashsecret('siphash_k0')),
                         types.np_uint64(_load_hashsecret('siphash_k1')),
                         val, _len)
        _hash = process_return(tmp)
    return process_return(_hash)

# This is an approximate translation of CPython's unicode_hash:
# https://github.com/python/cpython/blob/d1dd6be613381b996b9071443ef081de8e5f3aff/Objects/unicodeobject.c#L11635-L11663    # noqa: E501


@overload_method(types.UnicodeType, '__hash__')
def unicode_hash(val):
    from numba.cpython.unicode import _kind_to_byte_width

    def impl(val):
        kindwidth = _kind_to_byte_width(val._kind)
        _len = len(val)
        # use the cache if possible
        current_hash = val._hash
        if current_hash != -1:
            return current_hash
        else:
            # cannot write hash value to cache in the unicode struct due to
            # pass by value on the struct making the struct member immutable
            return _Py_HashBytes(val._data, kindwidth * _len)

    return impl
