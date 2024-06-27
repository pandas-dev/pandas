import sys
import operator

import numpy as np
from llvmlite.ir import IntType, Constant

from numba.core.cgutils import is_nonelike
from numba.core.extending import (
    models,
    register_model,
    make_attribute_wrapper,
    unbox,
    box,
    NativeValue,
    overload,
    overload_method,
    intrinsic,
    register_jitable,
)
from numba.core.imputils import (lower_constant, lower_cast, lower_builtin,
                                 iternext_impl, impl_ret_new_ref, RefType)
from numba.core.datamodel import register_default, StructModel
from numba.core import types, cgutils
from numba.core.utils import PYVERSION
from numba.core.pythonapi import (
    PY_UNICODE_1BYTE_KIND,
    PY_UNICODE_2BYTE_KIND,
    PY_UNICODE_4BYTE_KIND,
)
from numba._helperlib import c_helpers
from numba.cpython.hashing import _Py_hash_t
from numba.core.unsafe.bytes import memcpy_region
from numba.core.errors import TypingError
from numba.cpython.unicode_support import (_Py_TOUPPER, _Py_TOLOWER, _Py_UCS4,
                                           _Py_ISALNUM,
                                           _PyUnicode_ToUpperFull,
                                           _PyUnicode_ToLowerFull,
                                           _PyUnicode_ToFoldedFull,
                                           _PyUnicode_ToTitleFull,
                                           _PyUnicode_IsPrintable,
                                           _PyUnicode_IsSpace,
                                           _Py_ISSPACE,
                                           _PyUnicode_IsXidStart,
                                           _PyUnicode_IsXidContinue,
                                           _PyUnicode_IsCased,
                                           _PyUnicode_IsCaseIgnorable,
                                           _PyUnicode_IsUppercase,
                                           _PyUnicode_IsLowercase,
                                           _PyUnicode_IsLineBreak,
                                           _Py_ISLINEBREAK,
                                           _Py_ISLINEFEED,
                                           _Py_ISCARRIAGERETURN,
                                           _PyUnicode_IsTitlecase,
                                           _Py_ISLOWER,
                                           _Py_ISUPPER,
                                           _Py_TAB,
                                           _Py_LINEFEED,
                                           _Py_CARRIAGE_RETURN,
                                           _Py_SPACE,
                                           _PyUnicode_IsAlpha,
                                           _PyUnicode_IsNumeric,
                                           _Py_ISALPHA,
                                           _PyUnicode_IsDigit,
                                           _PyUnicode_IsDecimalDigit)
from numba.cpython import slicing

if PYVERSION in ((3, 9), (3, 10), (3, 11)):
    from numba.core.pythonapi import PY_UNICODE_WCHAR_KIND

# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L84-L85    # noqa: E501
_MAX_UNICODE = 0x10ffff

# https://github.com/python/cpython/blob/1960eb005e04b7ad8a91018088cfdb0646bc1ca0/Objects/stringlib/fastsearch.h#L31    # noqa: E501
_BLOOM_WIDTH = types.intp.bitwidth

# DATA MODEL


@register_model(types.UnicodeType)
class UnicodeModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('data', types.voidptr),
            ('length', types.intp),
            ('kind', types.int32),
            ('is_ascii', types.uint32),
            ('hash', _Py_hash_t),
            ('meminfo', types.MemInfoPointer(types.voidptr)),
            # A pointer to the owner python str/unicode object
            ('parent', types.pyobject),
        ]
        models.StructModel.__init__(self, dmm, fe_type, members)


make_attribute_wrapper(types.UnicodeType, 'data', '_data')
make_attribute_wrapper(types.UnicodeType, 'length', '_length')
make_attribute_wrapper(types.UnicodeType, 'kind', '_kind')
make_attribute_wrapper(types.UnicodeType, 'is_ascii', '_is_ascii')
make_attribute_wrapper(types.UnicodeType, 'hash', '_hash')


@register_default(types.UnicodeIteratorType)
class UnicodeIteratorModel(StructModel):
    def __init__(self, dmm, fe_type):
        members = [('index', types.EphemeralPointer(types.uintp)),
                   ('data', fe_type.data)]
        super(UnicodeIteratorModel, self).__init__(dmm, fe_type, members)

# CAST


def compile_time_get_string_data(obj):
    """Get string data from a python string for use at compile-time to embed
    the string data into the LLVM module.
    """
    from ctypes import (
        CFUNCTYPE, c_void_p, c_int, c_uint, c_ssize_t, c_ubyte, py_object,
        POINTER, byref,
    )

    extract_unicode_fn = c_helpers['extract_unicode']
    proto = CFUNCTYPE(c_void_p, py_object, POINTER(c_ssize_t), POINTER(c_int),
                      POINTER(c_uint), POINTER(c_ssize_t))
    fn = proto(extract_unicode_fn)
    length = c_ssize_t()
    kind = c_int()
    is_ascii = c_uint()
    hashv = c_ssize_t()
    data = fn(obj, byref(length), byref(kind), byref(is_ascii), byref(hashv))
    if data is None:
        raise ValueError("cannot extract unicode data from the given string")
    length = length.value
    kind = kind.value
    is_ascii = is_ascii.value
    nbytes = (length + 1) * _kind_to_byte_width(kind)
    out = (c_ubyte * nbytes).from_address(data)
    return bytes(out), length, kind, is_ascii, hashv.value


def make_string_from_constant(context, builder, typ, literal_string):
    """
    Get string data by `compile_time_get_string_data()` and return a
    unicode_type LLVM value
    """
    databytes, length, kind, is_ascii, hashv = \
        compile_time_get_string_data(literal_string)
    mod = builder.module
    gv = context.insert_const_bytes(mod, databytes)
    uni_str = cgutils.create_struct_proxy(typ)(context, builder)
    uni_str.data = gv
    uni_str.length = uni_str.length.type(length)
    uni_str.kind = uni_str.kind.type(kind)
    uni_str.is_ascii = uni_str.is_ascii.type(is_ascii)
    # Set hash to -1 to indicate that it should be computed.
    # We cannot bake in the hash value because of hashseed randomization.
    uni_str.hash = uni_str.hash.type(-1)
    return uni_str._getvalue()


@lower_cast(types.StringLiteral, types.unicode_type)
def cast_from_literal(context, builder, fromty, toty, val):
    return make_string_from_constant(
        context, builder, toty, fromty.literal_value,
    )


# CONSTANT

@lower_constant(types.unicode_type)
def constant_unicode(context, builder, typ, pyval):
    return make_string_from_constant(context, builder, typ, pyval)


# BOXING


@unbox(types.UnicodeType)
def unbox_unicode_str(typ, obj, c):
    """
    Convert a unicode str object to a native unicode structure.
    """
    ok, data, length, kind, is_ascii, hashv = \
        c.pyapi.string_as_string_size_and_kind(obj)
    uni_str = cgutils.create_struct_proxy(typ)(c.context, c.builder)
    uni_str.data = data
    uni_str.length = length
    uni_str.kind = kind
    uni_str.is_ascii = is_ascii
    uni_str.hash = hashv
    uni_str.meminfo = c.pyapi.nrt_meminfo_new_from_pyobject(
        data,  # the borrowed data pointer
        obj,   # the owner pyobject; the call will incref it.
    )
    uni_str.parent = obj

    is_error = cgutils.is_not_null(c.builder, c.pyapi.err_occurred())
    return NativeValue(uni_str._getvalue(), is_error=is_error)


@box(types.UnicodeType)
def box_unicode_str(typ, val, c):
    """
    Convert a native unicode structure to a unicode string
    """
    uni_str = cgutils.create_struct_proxy(typ)(c.context, c.builder, value=val)
    res = c.pyapi.string_from_kind_and_data(
        uni_str.kind, uni_str.data, uni_str.length)
    # hash isn't needed now, just compute it so it ends up in the unicodeobject
    # hash cache, cpython doesn't always do this, depends how a string was
    # created it's safe, just burns the cycles required to hash on @box
    c.pyapi.object_hash(res)
    c.context.nrt.decref(c.builder, typ, val)
    return res


# HELPER FUNCTIONS


def make_deref_codegen(bitsize):
    def codegen(context, builder, signature, args):
        data, idx = args
        ptr = builder.bitcast(data, IntType(bitsize).as_pointer())
        ch = builder.load(builder.gep(ptr, [idx]))
        return builder.zext(ch, IntType(32))

    return codegen


@intrinsic
def deref_uint8(typingctx, data, offset):
    sig = types.uint32(types.voidptr, types.intp)
    return sig, make_deref_codegen(8)


@intrinsic
def deref_uint16(typingctx, data, offset):
    sig = types.uint32(types.voidptr, types.intp)
    return sig, make_deref_codegen(16)


@intrinsic
def deref_uint32(typingctx, data, offset):
    sig = types.uint32(types.voidptr, types.intp)
    return sig, make_deref_codegen(32)


@intrinsic
def _malloc_string(typingctx, kind, char_bytes, length, is_ascii):
    """make empty string with data buffer of size alloc_bytes.

    Must set length and kind values for string after it is returned
    """
    def details(context, builder, signature, args):
        [kind_val, char_bytes_val, length_val, is_ascii_val] = args

        # fill the struct
        uni_str_ctor = cgutils.create_struct_proxy(types.unicode_type)
        uni_str = uni_str_ctor(context, builder)
        # add null padding character
        nbytes_val = builder.mul(char_bytes_val,
                                 builder.add(length_val,
                                             Constant(length_val.type, 1)))
        uni_str.meminfo = context.nrt.meminfo_alloc(builder, nbytes_val)
        uni_str.kind = kind_val
        uni_str.is_ascii = is_ascii_val
        uni_str.length = length_val
        # empty string has hash value -1 to indicate "need to compute hash"
        uni_str.hash = context.get_constant(_Py_hash_t, -1)
        uni_str.data = context.nrt.meminfo_data(builder, uni_str.meminfo)
        # Set parent to NULL
        uni_str.parent = cgutils.get_null_value(uni_str.parent.type)
        return uni_str._getvalue()

    sig = types.unicode_type(types.int32, types.intp, types.intp, types.uint32)
    return sig, details


@register_jitable
def _empty_string(kind, length, is_ascii=0):
    char_width = _kind_to_byte_width(kind)
    s = _malloc_string(kind, char_width, length, is_ascii)
    _set_code_point(s, length, np.uint32(0))    # Write NULL character
    return s


# Disable RefCt for performance.
@register_jitable(_nrt=False)
def _get_code_point(a, i):
    if a._kind == PY_UNICODE_1BYTE_KIND:
        return deref_uint8(a._data, i)
    elif a._kind == PY_UNICODE_2BYTE_KIND:
        return deref_uint16(a._data, i)
    elif a._kind == PY_UNICODE_4BYTE_KIND:
        return deref_uint32(a._data, i)
    else:
        # there's also a wchar kind, but that's one of the above,
        # so skipping for this example
        return 0

####


def make_set_codegen(bitsize):
    def codegen(context, builder, signature, args):
        data, idx, ch = args
        if bitsize < 32:
            ch = builder.trunc(ch, IntType(bitsize))
        ptr = builder.bitcast(data, IntType(bitsize).as_pointer())
        builder.store(ch, builder.gep(ptr, [idx]))
        return context.get_dummy_value()

    return codegen


@intrinsic
def set_uint8(typingctx, data, idx, ch):
    sig = types.void(types.voidptr, types.int64, types.uint32)
    return sig, make_set_codegen(8)


@intrinsic
def set_uint16(typingctx, data, idx, ch):
    sig = types.void(types.voidptr, types.int64, types.uint32)
    return sig, make_set_codegen(16)


@intrinsic
def set_uint32(typingctx, data, idx, ch):
    sig = types.void(types.voidptr, types.int64, types.uint32)
    return sig, make_set_codegen(32)


@register_jitable(_nrt=False)
def _set_code_point(a, i, ch):
    # WARNING: This method is very dangerous:
    #   * Assumes that data contents can be changed (only allowed for new
    #     strings)
    #   * Assumes that the kind of unicode string is sufficiently wide to
    #     accept ch.  Will truncate ch to make it fit.
    #   * Assumes that i is within the valid boundaries of the function
    if a._kind == PY_UNICODE_1BYTE_KIND:
        set_uint8(a._data, i, ch)
    elif a._kind == PY_UNICODE_2BYTE_KIND:
        set_uint16(a._data, i, ch)
    elif a._kind == PY_UNICODE_4BYTE_KIND:
        set_uint32(a._data, i, ch)
    else:
        raise AssertionError(
            "Unexpected unicode representation in _set_code_point")


if PYVERSION in ((3, 12),):
    @register_jitable
    def _pick_kind(kind1, kind2):
        if kind1 == PY_UNICODE_1BYTE_KIND:
            return kind2
        elif kind1 == PY_UNICODE_2BYTE_KIND:
            if kind2 == PY_UNICODE_4BYTE_KIND:
                return kind2
            else:
                return kind1
        elif kind1 == PY_UNICODE_4BYTE_KIND:
            return kind1
        else:
            raise AssertionError(
                "Unexpected unicode representation in _pick_kind")
elif PYVERSION in ((3, 9), (3, 10), (3, 11)):
    @register_jitable
    def _pick_kind(kind1, kind2):
        if (kind1 == PY_UNICODE_WCHAR_KIND or kind2 == PY_UNICODE_WCHAR_KIND):
            raise AssertionError("PY_UNICODE_WCHAR_KIND unsupported")

        if kind1 == PY_UNICODE_1BYTE_KIND:
            return kind2
        elif kind1 == PY_UNICODE_2BYTE_KIND:
            if kind2 == PY_UNICODE_4BYTE_KIND:
                return kind2
            else:
                return kind1
        elif kind1 == PY_UNICODE_4BYTE_KIND:
            return kind1
        else:
            raise AssertionError(
                "Unexpected unicode representation in _pick_kind")
else:
    raise NotImplementedError(PYVERSION)


@register_jitable
def _pick_ascii(is_ascii1, is_ascii2):
    if is_ascii1 == 1 and is_ascii2 == 1:
        return types.uint32(1)
    return types.uint32(0)


if PYVERSION in ((3, 12),):
    @register_jitable
    def _kind_to_byte_width(kind):
        if kind == PY_UNICODE_1BYTE_KIND:
            return 1
        elif kind == PY_UNICODE_2BYTE_KIND:
            return 2
        elif kind == PY_UNICODE_4BYTE_KIND:
            return 4
        else:
            raise AssertionError("Unexpected unicode encoding encountered")
elif PYVERSION in ((3, 9), (3, 10), (3, 11)):
    @register_jitable
    def _kind_to_byte_width(kind):
        if kind == PY_UNICODE_1BYTE_KIND:
            return 1
        elif kind == PY_UNICODE_2BYTE_KIND:
            return 2
        elif kind == PY_UNICODE_4BYTE_KIND:
            return 4
        elif kind == PY_UNICODE_WCHAR_KIND:
            raise AssertionError("PY_UNICODE_WCHAR_KIND unsupported")
        else:
            raise AssertionError("Unexpected unicode encoding encountered")
else:
    raise NotImplementedError(PYVERSION)


@register_jitable(_nrt=False)
def _cmp_region(a, a_offset, b, b_offset, n):
    if n == 0:
        return 0
    elif a_offset + n > a._length:
        return -1
    elif b_offset + n > b._length:
        return 1

    for i in range(n):
        a_chr = _get_code_point(a, a_offset + i)
        b_chr = _get_code_point(b, b_offset + i)
        if a_chr < b_chr:
            return -1
        elif a_chr > b_chr:
            return 1

    return 0


@register_jitable
def _codepoint_to_kind(cp):
    """
    Compute the minimum unicode kind needed to hold a given codepoint
    """
    if cp < 256:
        return PY_UNICODE_1BYTE_KIND
    elif cp < 65536:
        return PY_UNICODE_2BYTE_KIND
    else:
        # Maximum code point of Unicode 6.0: 0x10ffff (1,114,111)
        MAX_UNICODE = 0x10ffff
        if cp > MAX_UNICODE:
            msg = "Invalid codepoint. Found value greater than Unicode maximum"
            raise ValueError(msg)
        return PY_UNICODE_4BYTE_KIND


@register_jitable
def _codepoint_is_ascii(ch):
    """
    Returns true if a codepoint is in the ASCII range
    """
    return ch < 128


# PUBLIC API


@overload(len)
def unicode_len(s):
    if isinstance(s, types.UnicodeType):
        def len_impl(s):
            return s._length
        return len_impl


@overload(operator.eq)
def unicode_eq(a, b):
    if not (a.is_internal and b.is_internal):
        return
    if isinstance(a, types.Optional):
        check_a = a.type
    else:
        check_a = a
    if isinstance(b, types.Optional):
        check_b = b.type
    else:
        check_b = b
    accept = (types.UnicodeType, types.StringLiteral, types.UnicodeCharSeq)
    a_unicode = isinstance(check_a, accept)
    b_unicode = isinstance(check_b, accept)
    if a_unicode and b_unicode:
        def eq_impl(a, b):
            # handle Optionals at runtime
            a_none = a is None
            b_none = b is None
            if a_none or b_none:
                if a_none and b_none:
                    return True
                else:
                    return False
            # the str() is for UnicodeCharSeq, it's a nop else
            a = str(a)
            b = str(b)
            if len(a) != len(b):
                return False
            return _cmp_region(a, 0, b, 0, len(a)) == 0
        return eq_impl
    elif a_unicode ^ b_unicode:
        # one of the things is unicode, everything compares False
        def eq_impl(a, b):
            return False
        return eq_impl


@overload(operator.ne)
def unicode_ne(a, b):
    if not (a.is_internal and b.is_internal):
        return
    accept = (types.UnicodeType, types.StringLiteral, types.UnicodeCharSeq)
    a_unicode = isinstance(a, accept)
    b_unicode = isinstance(b, accept)
    if a_unicode and b_unicode:
        def ne_impl(a, b):
            return not (a == b)
        return ne_impl
    elif a_unicode ^ b_unicode:
        # one of the things is unicode, everything compares True
        def eq_impl(a, b):
            return True
        return eq_impl


@overload(operator.lt)
def unicode_lt(a, b):
    a_unicode = isinstance(a, (types.UnicodeType, types.StringLiteral))
    b_unicode = isinstance(b, (types.UnicodeType, types.StringLiteral))
    if a_unicode and b_unicode:
        def lt_impl(a, b):
            minlen = min(len(a), len(b))
            eqcode = _cmp_region(a, 0, b, 0, minlen)
            if eqcode == -1:
                return True
            elif eqcode == 0:
                return len(a) < len(b)
            return False
        return lt_impl


@overload(operator.gt)
def unicode_gt(a, b):
    a_unicode = isinstance(a, (types.UnicodeType, types.StringLiteral))
    b_unicode = isinstance(b, (types.UnicodeType, types.StringLiteral))
    if a_unicode and b_unicode:
        def gt_impl(a, b):
            minlen = min(len(a), len(b))
            eqcode = _cmp_region(a, 0, b, 0, minlen)
            if eqcode == 1:
                return True
            elif eqcode == 0:
                return len(a) > len(b)
            return False
        return gt_impl


@overload(operator.le)
def unicode_le(a, b):
    a_unicode = isinstance(a, (types.UnicodeType, types.StringLiteral))
    b_unicode = isinstance(b, (types.UnicodeType, types.StringLiteral))
    if a_unicode and b_unicode:
        def le_impl(a, b):
            return not (a > b)
        return le_impl


@overload(operator.ge)
def unicode_ge(a, b):
    a_unicode = isinstance(a, (types.UnicodeType, types.StringLiteral))
    b_unicode = isinstance(b, (types.UnicodeType, types.StringLiteral))
    if a_unicode and b_unicode:
        def ge_impl(a, b):
            return not (a < b)
        return ge_impl


@overload(operator.contains)
def unicode_contains(a, b):
    if isinstance(a, types.UnicodeType) and isinstance(b, types.UnicodeType):
        def contains_impl(a, b):
            # note parameter swap: contains(a, b) == b in a
            return _find(a, b) > -1
        return contains_impl


def unicode_idx_check_type(ty, name):
    """Check object belongs to one of specific types
    ty: type
        Type of the object
    name: str
        Name of the object
    """
    thety = ty
    # if the type is omitted, the concrete type is the value
    if isinstance(ty, types.Omitted):
        thety = ty.value
    # if the type is optional, the concrete type is the captured type
    elif isinstance(ty, types.Optional):
        thety = ty.type

    accepted = (types.Integer, types.NoneType)
    if thety is not None and not isinstance(thety, accepted):
        raise TypingError('"{}" must be {}, not {}'.format(name, accepted, ty))


def unicode_sub_check_type(ty, name):
    """Check object belongs to unicode type"""
    if not isinstance(ty, types.UnicodeType):
        msg = '"{}" must be {}, not {}'.format(name, types.UnicodeType, ty)
        raise TypingError(msg)


# FAST SEARCH algorithm implementation from cpython

@register_jitable
def _bloom_add(mask, ch):
    mask |= (1 << (ch & (_BLOOM_WIDTH - 1)))
    return mask


@register_jitable
def _bloom_check(mask, ch):
    return mask & (1 << (ch & (_BLOOM_WIDTH - 1)))


# https://github.com/python/cpython/blob/1960eb005e04b7ad8a91018088cfdb0646bc1ca0/Objects/stringlib/fastsearch.h#L550    # noqa: E501
@register_jitable
def _default_find(data, substr, start, end):
    """Left finder."""
    m = len(substr)
    if m == 0:
        return start

    gap = mlast = m - 1
    last = _get_code_point(substr, mlast)

    zero = types.intp(0)
    mask = _bloom_add(zero, last)
    for i in range(mlast):
        ch = _get_code_point(substr, i)
        mask = _bloom_add(mask, ch)
        if ch == last:
            gap = mlast - i - 1

    i = start
    while i <= end - m:
        ch = _get_code_point(data, mlast + i)
        if ch == last:
            j = 0
            while j < mlast:
                haystack_ch = _get_code_point(data, i + j)
                needle_ch = _get_code_point(substr, j)
                if haystack_ch != needle_ch:
                    break
                j += 1
            if j == mlast:
                # got a match
                return i

            ch = _get_code_point(data, mlast + i + 1)
            if _bloom_check(mask, ch) == 0:
                i += m
            else:
                i += gap
        else:
            ch = _get_code_point(data, mlast + i + 1)
            if _bloom_check(mask, ch) == 0:
                i += m
        i += 1

    return -1


@register_jitable
def _default_rfind(data, substr, start, end):
    """Right finder."""
    m = len(substr)
    if m == 0:
        return end

    skip = mlast = m - 1
    mfirst = _get_code_point(substr, 0)
    mask = _bloom_add(0, mfirst)
    i = mlast
    while i > 0:
        ch = _get_code_point(substr, i)
        mask = _bloom_add(mask, ch)
        if ch == mfirst:
            skip = i - 1
        i -= 1

    i = end - m
    while i >= start:
        ch = _get_code_point(data, i)
        if ch == mfirst:
            j = mlast
            while j > 0:
                haystack_ch = _get_code_point(data, i + j)
                needle_ch = _get_code_point(substr, j)
                if haystack_ch != needle_ch:
                    break
                j -= 1

            if j == 0:
                # got a match
                return i

            ch = _get_code_point(data, i - 1)
            if i > start and _bloom_check(mask, ch) == 0:
                i -= m
            else:
                i -= skip

        else:
            ch = _get_code_point(data, i - 1)
            if i > start and _bloom_check(mask, ch) == 0:
                i -= m
        i -= 1

    return -1


def generate_finder(find_func):
    """Generate finder either left or right."""
    def impl(data, substr, start=None, end=None):
        length = len(data)
        sub_length = len(substr)
        if start is None:
            start = 0
        if end is None:
            end = length

        start, end = _adjust_indices(length, start, end)
        if end - start < sub_length:
            return -1

        return find_func(data, substr, start, end)

    return impl


_find = register_jitable(generate_finder(_default_find))
_rfind = register_jitable(generate_finder(_default_rfind))


@overload_method(types.UnicodeType, 'find')
def unicode_find(data, substr, start=None, end=None):
    """Implements str.find()"""
    if isinstance(substr, types.UnicodeCharSeq):
        def find_impl(data, substr, start=None, end=None):
            return data.find(str(substr))
        return find_impl

    unicode_idx_check_type(start, 'start')
    unicode_idx_check_type(end, 'end')
    unicode_sub_check_type(substr, 'substr')

    return _find


@overload_method(types.UnicodeType, 'rfind')
def unicode_rfind(data, substr, start=None, end=None):
    """Implements str.rfind()"""
    if isinstance(substr, types.UnicodeCharSeq):
        def rfind_impl(data, substr, start=None, end=None):
            return data.rfind(str(substr))
        return rfind_impl

    unicode_idx_check_type(start, 'start')
    unicode_idx_check_type(end, 'end')
    unicode_sub_check_type(substr, 'substr')

    return _rfind


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12831-L12857    # noqa: E501
@overload_method(types.UnicodeType, 'rindex')
def unicode_rindex(s, sub, start=None, end=None):
    """Implements str.rindex()"""
    unicode_idx_check_type(start, 'start')
    unicode_idx_check_type(end, 'end')
    unicode_sub_check_type(sub, 'sub')

    def rindex_impl(s, sub, start=None, end=None):
        result = s.rfind(sub, start, end)
        if result < 0:
            raise ValueError('substring not found')

        return result

    return rindex_impl


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11692-L11718    # noqa: E501
@overload_method(types.UnicodeType, 'index')
def unicode_index(s, sub, start=None, end=None):
    """Implements str.index()"""
    unicode_idx_check_type(start, 'start')
    unicode_idx_check_type(end, 'end')
    unicode_sub_check_type(sub, 'sub')

    def index_impl(s, sub, start=None, end=None):
        result = s.find(sub, start, end)
        if result < 0:
            raise ValueError('substring not found')

        return result

    return index_impl


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12922-L12976    # noqa: E501
@overload_method(types.UnicodeType, 'partition')
def unicode_partition(data, sep):
    """Implements str.partition()"""
    thety = sep
    # if the type is omitted, the concrete type is the value
    if isinstance(sep, types.Omitted):
        thety = sep.value
    # if the type is optional, the concrete type is the captured type
    elif isinstance(sep, types.Optional):
        thety = sep.type

    accepted = (types.UnicodeType, types.UnicodeCharSeq)
    if thety is not None and not isinstance(thety, accepted):
        msg = '"{}" must be {}, not {}'.format('sep', accepted, sep)
        raise TypingError(msg)

    def impl(data, sep):
        # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/stringlib/partition.h#L7-L60    # noqa: E501
        sep = str(sep)
        empty_str = _empty_string(data._kind, 0, data._is_ascii)
        sep_length = len(sep)
        if data._kind < sep._kind or len(data) < sep_length:
            return data, empty_str, empty_str

        if sep_length == 0:
            raise ValueError('empty separator')

        pos = data.find(sep)
        if pos < 0:
            return data, empty_str, empty_str

        return data[0:pos], sep, data[pos + sep_length:len(data)]

    return impl


@overload_method(types.UnicodeType, 'count')
def unicode_count(src, sub, start=None, end=None):

    _count_args_types_check(start)
    _count_args_types_check(end)

    if isinstance(sub, types.UnicodeType):
        def count_impl(src, sub, start=None, end=None):
            count = 0
            src_len = len(src)
            sub_len = len(sub)

            start = _normalize_slice_idx_count(start, src_len, 0)
            end = _normalize_slice_idx_count(end, src_len, src_len)

            if end - start < 0 or start > src_len:
                return 0

            src = src[start : end]
            src_len = len(src)
            start, end = 0, src_len
            if sub_len == 0:
                return src_len + 1

            while (start + sub_len <= src_len):
                if src[start : start + sub_len] == sub:
                    count += 1
                    start += sub_len
                else:
                    start += 1
            return count
        return count_impl
    error_msg = "The substring must be a UnicodeType, not {}"
    raise TypingError(error_msg.format(type(sub)))


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12979-L13033    # noqa: E501
@overload_method(types.UnicodeType, 'rpartition')
def unicode_rpartition(data, sep):
    """Implements str.rpartition()"""
    thety = sep
    # if the type is omitted, the concrete type is the value
    if isinstance(sep, types.Omitted):
        thety = sep.value
    # if the type is optional, the concrete type is the captured type
    elif isinstance(sep, types.Optional):
        thety = sep.type

    accepted = (types.UnicodeType, types.UnicodeCharSeq)
    if thety is not None and not isinstance(thety, accepted):
        msg = '"{}" must be {}, not {}'.format('sep', accepted, sep)
        raise TypingError(msg)

    def impl(data, sep):
        # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/stringlib/partition.h#L62-L115    # noqa: E501
        sep = str(sep)
        empty_str = _empty_string(data._kind, 0, data._is_ascii)
        sep_length = len(sep)
        if data._kind < sep._kind or len(data) < sep_length:
            return empty_str, empty_str, data

        if sep_length == 0:
            raise ValueError('empty separator')

        pos = data.rfind(sep)
        if pos < 0:
            return empty_str, empty_str, data

        return data[0:pos], sep, data[pos + sep_length:len(data)]

    return impl


# https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L9342-L9354    # noqa: E501
@register_jitable
def _adjust_indices(length, start, end):
    if end > length:
        end = length
    if end < 0:
        end += length
        if end < 0:
            end = 0
    if start < 0:
        start += length
        if start < 0:
            start = 0

    return start, end


@overload_method(types.UnicodeType, 'startswith')
def unicode_startswith(s, prefix, start=None, end=None):
    if not is_nonelike(start) and not isinstance(start, types.Integer):
        raise TypingError(
            "When specified, the arg 'start' must be an Integer or None")

    if not is_nonelike(end) and not isinstance(end, types.Integer):
        raise TypingError(
            "When specified, the arg 'end' must be an Integer or None")

    if isinstance(prefix, types.UniTuple) and \
            isinstance(prefix.dtype, types.UnicodeType):

        def startswith_tuple_impl(s, prefix, start=None, end=None):
            for item in prefix:
                if s.startswith(item, start, end):
                    return True
            return False

        return startswith_tuple_impl

    elif isinstance(prefix, types.UnicodeCharSeq):
        def startswith_char_seq_impl(s, prefix, start=None, end=None):
            return s.startswith(str(prefix), start, end)

        return startswith_char_seq_impl

    elif isinstance(prefix, types.UnicodeType):
        def startswith_unicode_impl(s, prefix, start=None, end=None):
            length, prefix_length = len(s), len(prefix)
            if start is None:
                start = 0
            if end is None:
                end = length

            start, end = _adjust_indices(length, start, end)
            if end - start < prefix_length:
                return False

            if prefix_length == 0:
                return True

            s_slice = s[start:end]

            return _cmp_region(s_slice, 0, prefix, 0, prefix_length) == 0

        return startswith_unicode_impl

    else:
        raise TypingError(
            "The arg 'prefix' should be a string or a tuple of strings")


@overload_method(types.UnicodeType, 'endswith')
def unicode_endswith(s, substr, start=None, end=None):
    if not (start is None or isinstance(start, (types.Omitted,
                                                types.Integer,
                                                types.NoneType))):
        raise TypingError('The arg must be a Integer or None')

    if not (end is None or isinstance(end, (types.Omitted,
                                            types.Integer,
                                            types.NoneType))):
        raise TypingError('The arg must be a Integer or None')

    if isinstance(substr, (types.Tuple, types.UniTuple)):
        def endswith_impl(s, substr, start=None, end=None):
            for item in substr:
                if s.endswith(item, start, end) is True:
                    return True

            return False
        return endswith_impl

    if isinstance(substr, types.UnicodeType):
        def endswith_impl(s, substr, start=None, end=None):
            length = len(s)
            sub_length = len(substr)
            if start is None:
                start = 0
            if end is None:
                end = length

            start, end = _adjust_indices(length, start, end)
            if end - start < sub_length:
                return False

            if sub_length == 0:
                return True

            s = s[start:end]
            offset = len(s) - sub_length

            return _cmp_region(s, offset, substr, 0, sub_length) == 0
        return endswith_impl

    if isinstance(substr, types.UnicodeCharSeq):
        def endswith_impl(s, substr, start=None, end=None):
            return s.endswith(str(substr), start, end)
        return endswith_impl


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11519-L11595    # noqa: E501
@overload_method(types.UnicodeType, 'expandtabs')
def unicode_expandtabs(data, tabsize=8):
    """Implements str.expandtabs()"""
    thety = tabsize
    # if the type is omitted, the concrete type is the value
    if isinstance(tabsize, types.Omitted):
        thety = tabsize.value
    # if the type is optional, the concrete type is the captured type
    elif isinstance(tabsize, types.Optional):
        thety = tabsize.type

    accepted = (types.Integer, int)
    if thety is not None and not isinstance(thety, accepted):
        raise TypingError(
            '"tabsize" must be {}, not {}'.format(accepted, tabsize))

    def expandtabs_impl(data, tabsize=8):
        length = len(data)
        j = line_pos = 0
        found = False
        for i in range(length):
            code_point = _get_code_point(data, i)
            if code_point == _Py_TAB:
                found = True
                if tabsize > 0:
                    # cannot overflow
                    incr = tabsize - (line_pos % tabsize)
                    if j > sys.maxsize - incr:
                        raise OverflowError('new string is too long')
                    line_pos += incr
                    j += incr
            else:
                if j > sys.maxsize - 1:
                    raise OverflowError('new string is too long')
                line_pos += 1
                j += 1
                if code_point in (_Py_LINEFEED, _Py_CARRIAGE_RETURN):
                    line_pos = 0

        if not found:
            return data

        res = _empty_string(data._kind, j, data._is_ascii)
        j = line_pos = 0
        for i in range(length):
            code_point = _get_code_point(data, i)
            if code_point == _Py_TAB:
                if tabsize > 0:
                    incr = tabsize - (line_pos % tabsize)
                    line_pos += incr
                    for idx in range(j, j + incr):
                        _set_code_point(res, idx, _Py_SPACE)
                    j += incr
            else:
                line_pos += 1
                _set_code_point(res, j, code_point)
                j += 1
                if code_point in (_Py_LINEFEED, _Py_CARRIAGE_RETURN):
                    line_pos = 0

        return res

    return expandtabs_impl


@overload_method(types.UnicodeType, 'split')
def unicode_split(a, sep=None, maxsplit=-1):
    if not (maxsplit == -1 or
            isinstance(maxsplit, (types.Omitted, types.Integer,
                                  types.IntegerLiteral))):
        return None  # fail typing if maxsplit is not an integer

    if isinstance(sep, types.UnicodeCharSeq):
        def split_impl(a, sep=None, maxsplit=-1):
            return a.split(str(sep), maxsplit=maxsplit)
        return split_impl

    if isinstance(sep, types.UnicodeType):
        def split_impl(a, sep=None, maxsplit=-1):
            a_len = len(a)
            sep_len = len(sep)

            if sep_len == 0:
                raise ValueError('empty separator')

            parts = []
            last = 0
            idx = 0

            if sep_len == 1 and maxsplit == -1:
                sep_code_point = _get_code_point(sep, 0)
                for idx in range(a_len):
                    if _get_code_point(a, idx) == sep_code_point:
                        parts.append(a[last:idx])
                        last = idx + 1
            else:
                split_count = 0

                while idx < a_len and (maxsplit == -1 or
                                       split_count < maxsplit):
                    if _cmp_region(a, idx, sep, 0, sep_len) == 0:
                        parts.append(a[last:idx])
                        idx += sep_len
                        last = idx
                        split_count += 1
                    else:
                        idx += 1

            if last <= a_len:
                parts.append(a[last:])

            return parts
        return split_impl
    elif sep is None or isinstance(sep, types.NoneType) or \
            getattr(sep, 'value', False) is None:
        def split_whitespace_impl(a, sep=None, maxsplit=-1):
            a_len = len(a)

            parts = []
            last = 0
            idx = 0
            split_count = 0
            in_whitespace_block = True

            for idx in range(a_len):
                code_point = _get_code_point(a, idx)
                is_whitespace = _PyUnicode_IsSpace(code_point)
                if in_whitespace_block:
                    if is_whitespace:
                        pass  # keep consuming space
                    else:
                        last = idx  # this is the start of the next string
                        in_whitespace_block = False
                else:
                    if not is_whitespace:
                        pass  # keep searching for whitespace transition
                    else:
                        parts.append(a[last:idx])
                        in_whitespace_block = True
                        split_count += 1
                        if maxsplit != -1 and split_count == maxsplit:
                            break

            if last <= a_len and not in_whitespace_block:
                parts.append(a[last:])

            return parts
        return split_whitespace_impl


def generate_rsplit_whitespace_impl(isspace_func):
    """Generate whitespace rsplit func based on either ascii or unicode"""

    def rsplit_whitespace_impl(data, sep=None, maxsplit=-1):
        # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/stringlib/split.h#L192-L240    # noqa: E501
        if maxsplit < 0:
            maxsplit = sys.maxsize

        result = []
        i = len(data) - 1
        while maxsplit > 0:
            while i >= 0:
                code_point = _get_code_point(data, i)
                if not isspace_func(code_point):
                    break
                i -= 1
            if i < 0:
                break
            j = i
            i -= 1
            while i >= 0:
                code_point = _get_code_point(data, i)
                if isspace_func(code_point):
                    break
                i -= 1
            result.append(data[i + 1:j + 1])
            maxsplit -= 1

        if i >= 0:
            # Only occurs when maxsplit was reached
            # Skip any remaining whitespace and copy to beginning of string
            while i >= 0:
                code_point = _get_code_point(data, i)
                if not isspace_func(code_point):
                    break
                i -= 1
            if i >= 0:
                result.append(data[0:i + 1])

        return result[::-1]

    return rsplit_whitespace_impl


unicode_rsplit_whitespace_impl = register_jitable(
    generate_rsplit_whitespace_impl(_PyUnicode_IsSpace))
ascii_rsplit_whitespace_impl = register_jitable(
    generate_rsplit_whitespace_impl(_Py_ISSPACE))


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L13095-L13108    # noqa: E501
@overload_method(types.UnicodeType, 'rsplit')
def unicode_rsplit(data, sep=None, maxsplit=-1):
    """Implements str.unicode_rsplit()"""

    def _unicode_rsplit_check_type(ty, name, accepted):
        """Check object belongs to one of specified types"""
        thety = ty
        # if the type is omitted, the concrete type is the value
        if isinstance(ty, types.Omitted):
            thety = ty.value
        # if the type is optional, the concrete type is the captured type
        elif isinstance(ty, types.Optional):
            thety = ty.type

        if thety is not None and not isinstance(thety, accepted):
            raise TypingError(
                '"{}" must be {}, not {}'.format(name, accepted, ty))

    _unicode_rsplit_check_type(sep, 'sep', (types.UnicodeType,
                                            types.UnicodeCharSeq,
                                            types.NoneType))
    _unicode_rsplit_check_type(maxsplit, 'maxsplit', (types.Integer, int))

    if sep is None or isinstance(sep, (types.NoneType, types.Omitted)):

        def rsplit_whitespace_impl(data, sep=None, maxsplit=-1):
            if data._is_ascii:
                return ascii_rsplit_whitespace_impl(data, sep, maxsplit)
            return unicode_rsplit_whitespace_impl(data, sep, maxsplit)

        return rsplit_whitespace_impl

    def rsplit_impl(data, sep=None, maxsplit=-1):
        sep = str(sep)
        # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/stringlib/split.h#L286-L333    # noqa: E501
        if data._kind < sep._kind or len(data) < len(sep):
            return [data]

        def _rsplit_char(data, ch, maxsplit):
            # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/stringlib/split.h#L242-L284    # noqa: E501
            result = []
            ch_code_point = _get_code_point(ch, 0)
            i = j = len(data) - 1
            while i >= 0 and maxsplit > 0:
                data_code_point = _get_code_point(data, i)
                if data_code_point == ch_code_point:
                    result.append(data[i + 1 : j + 1])
                    j = i = i - 1
                    maxsplit -= 1
                i -= 1
            if j >= -1:
                result.append(data[0 : j + 1])

            return result[::-1]

        if maxsplit < 0:
            maxsplit = sys.maxsize

        sep_length = len(sep)

        if sep_length == 0:
            raise ValueError('empty separator')
        if sep_length == 1:
            return _rsplit_char(data, sep, maxsplit)

        result = []
        j = len(data)
        while maxsplit > 0:
            pos = data.rfind(sep, start=0, end=j)
            if pos < 0:
                break
            result.append(data[pos + sep_length:j])
            j = pos
            maxsplit -= 1

        result.append(data[0:j])

        return result[::-1]

    return rsplit_impl


@overload_method(types.UnicodeType, 'center')
def unicode_center(string, width, fillchar=' '):
    if not isinstance(width, types.Integer):
        raise TypingError('The width must be an Integer')

    if isinstance(fillchar, types.UnicodeCharSeq):
        def center_impl(string, width, fillchar=' '):
            return string.center(width, str(fillchar))
        return center_impl

    if not (fillchar == ' ' or
            isinstance(fillchar, (types.Omitted, types.UnicodeType))):
        raise TypingError('The fillchar must be a UnicodeType')

    def center_impl(string, width, fillchar=' '):
        str_len = len(string)
        fillchar_len = len(fillchar)

        if fillchar_len != 1:
            raise ValueError('The fill character must be exactly one '
                             'character long')

        if width <= str_len:
            return string

        allmargin = width - str_len
        lmargin = (allmargin // 2) + (allmargin & width & 1)
        rmargin = allmargin - lmargin

        l_string = fillchar * lmargin
        if lmargin == rmargin:
            return l_string + string + l_string
        else:
            return l_string + string + (fillchar * rmargin)

    return center_impl


def gen_unicode_Xjust(STRING_FIRST):
    def unicode_Xjust(string, width, fillchar=' '):
        if not isinstance(width, types.Integer):
            raise TypingError('The width must be an Integer')

        if isinstance(fillchar, types.UnicodeCharSeq):
            if STRING_FIRST:
                def ljust_impl(string, width, fillchar=' '):
                    return string.ljust(width, str(fillchar))
                return ljust_impl
            else:
                def rjust_impl(string, width, fillchar=' '):
                    return string.rjust(width, str(fillchar))
                return rjust_impl

        if not (fillchar == ' ' or
                isinstance(fillchar, (types.Omitted, types.UnicodeType))):
            raise TypingError('The fillchar must be a UnicodeType')

        def impl(string, width, fillchar=' '):
            str_len = len(string)
            fillchar_len = len(fillchar)

            if fillchar_len != 1:
                raise ValueError('The fill character must be exactly one '
                                 'character long')

            if width <= str_len:
                return string

            newstr = (fillchar * (width - str_len))
            if STRING_FIRST:
                return string + newstr
            else:
                return newstr + string

        return impl

    return unicode_Xjust


overload_method(types.UnicodeType, 'rjust')(gen_unicode_Xjust(False))
overload_method(types.UnicodeType, 'ljust')(gen_unicode_Xjust(True))


def generate_splitlines_func(is_line_break_func):
    """Generate splitlines performer based on ascii or unicode line breaks."""
    def impl(data, keepends):
        # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/stringlib/split.h#L335-L389    # noqa: E501
        length = len(data)
        result = []
        i = j = 0
        while i < length:
            # find a line and append it
            while i < length:
                code_point = _get_code_point(data, i)
                if is_line_break_func(code_point):
                    break
                i += 1

            # skip the line break reading CRLF as one line break
            eol = i
            if i < length:
                if i + 1 < length:
                    cur_cp = _get_code_point(data, i)
                    next_cp = _get_code_point(data, i + 1)
                    if _Py_ISCARRIAGERETURN(cur_cp) and _Py_ISLINEFEED(next_cp):
                        i += 1
                i += 1
                if keepends:
                    eol = i

            result.append(data[j:eol])
            j = i

        return result

    return impl


_ascii_splitlines = register_jitable(generate_splitlines_func(_Py_ISLINEBREAK))
_unicode_splitlines = register_jitable(generate_splitlines_func(
    _PyUnicode_IsLineBreak))


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L10196-L10229    # noqa: E501
@overload_method(types.UnicodeType, 'splitlines')
def unicode_splitlines(data, keepends=False):
    """Implements str.splitlines()"""
    thety = keepends
    # if the type is omitted, the concrete type is the value
    if isinstance(keepends, types.Omitted):
        thety = keepends.value
    # if the type is optional, the concrete type is the captured type
    elif isinstance(keepends, types.Optional):
        thety = keepends.type

    accepted = (types.Integer, int, types.Boolean, bool)
    if thety is not None and not isinstance(thety, accepted):
        raise TypingError(
            '"{}" must be {}, not {}'.format('keepends', accepted, keepends))

    def splitlines_impl(data, keepends=False):
        if data._is_ascii:
            return _ascii_splitlines(data, keepends)

        return _unicode_splitlines(data, keepends)

    return splitlines_impl


@register_jitable
def join_list(sep, parts):
    parts_len = len(parts)
    if parts_len == 0:
        return ''

    # Precompute size and char_width of result
    sep_len = len(sep)
    length = (parts_len - 1) * sep_len
    kind = sep._kind
    is_ascii = sep._is_ascii
    for p in parts:
        length += len(p)
        kind = _pick_kind(kind, p._kind)
        is_ascii = _pick_ascii(is_ascii, p._is_ascii)

    result = _empty_string(kind, length, is_ascii)

    # populate string
    part = parts[0]
    _strncpy(result, 0, part, 0, len(part))
    dst_offset = len(part)
    for idx in range(1, parts_len):
        _strncpy(result, dst_offset, sep, 0, sep_len)
        dst_offset += sep_len
        part = parts[idx]
        _strncpy(result, dst_offset, part, 0, len(part))
        dst_offset += len(part)

    return result


@overload_method(types.UnicodeType, 'join')
def unicode_join(sep, parts):

    if isinstance(parts, types.List):
        if isinstance(parts.dtype, types.UnicodeType):
            def join_list_impl(sep, parts):
                return join_list(sep, parts)
            return join_list_impl
        elif isinstance(parts.dtype, types.UnicodeCharSeq):
            def join_list_impl(sep, parts):
                _parts = [str(p) for p in parts]
                return join_list(sep, _parts)
            return join_list_impl
        else:
            pass  # lists of any other type not supported
    elif isinstance(parts, types.IterableType):
        def join_iter_impl(sep, parts):
            parts_list = [p for p in parts]
            return sep.join(parts_list)
        return join_iter_impl
    elif isinstance(parts, types.UnicodeType):
        # Temporary workaround until UnicodeType is iterable
        def join_str_impl(sep, parts):
            parts_list = [parts[i] for i in range(len(parts))]
            return join_list(sep, parts_list)
        return join_str_impl


@overload_method(types.UnicodeType, 'zfill')
def unicode_zfill(string, width):
    if not isinstance(width, types.Integer):
        raise TypingError("<width> must be an Integer")

    def zfill_impl(string, width):

        str_len = len(string)

        if width <= str_len:
            return string

        first_char = string[0] if str_len else ''
        padding = '0' * (width - str_len)

        if first_char in ['+', '-']:
            newstr = first_char + padding + string[1:]
        else:
            newstr = padding + string

        return newstr

    return zfill_impl


# ------------------------------------------------------------------------------
# Strip functions
# ------------------------------------------------------------------------------
@register_jitable
def unicode_strip_left_bound(string, chars):
    str_len = len(string)

    i = 0
    if chars is not None:
        for i in range(str_len):
            if string[i] not in chars:
                return i
    else:
        for i in range(str_len):
            if not _PyUnicode_IsSpace(string[i]):
                return i

    return str_len


@register_jitable
def unicode_strip_right_bound(string, chars):
    str_len = len(string)
    i = 0
    if chars is not None:
        for i in range(str_len - 1, -1, -1):
            if string[i] not in chars:
                i += 1
                break
    else:
        for i in range(str_len - 1, -1, -1):
            if not _PyUnicode_IsSpace(string[i]):
                i += 1
                break

    return i


def unicode_strip_types_check(chars):
    if isinstance(chars, types.Optional):
        chars = chars.type  # catch optional type with invalid non-None type
    if not (chars is None or isinstance(chars, (types.Omitted,
                                                types.UnicodeType,
                                                types.NoneType))):
        raise TypingError('The arg must be a UnicodeType or None')


def _count_args_types_check(arg):
    if isinstance(arg, types.Optional):
        arg = arg.type
    if not (arg is None or isinstance(arg, (types.Omitted,
                                            types.Integer,
                                            types.NoneType))):
        raise TypingError("The slice indices must be an Integer or None")


@overload_method(types.UnicodeType, 'lstrip')
def unicode_lstrip(string, chars=None):

    if isinstance(chars, types.UnicodeCharSeq):
        def lstrip_impl(string, chars=None):
            return string.lstrip(str(chars))
        return lstrip_impl

    unicode_strip_types_check(chars)

    def lstrip_impl(string, chars=None):
        return string[unicode_strip_left_bound(string, chars):]
    return lstrip_impl


@overload_method(types.UnicodeType, 'rstrip')
def unicode_rstrip(string, chars=None):

    if isinstance(chars, types.UnicodeCharSeq):
        def rstrip_impl(string, chars=None):
            return string.rstrip(str(chars))
        return rstrip_impl

    unicode_strip_types_check(chars)

    def rstrip_impl(string, chars=None):
        return string[:unicode_strip_right_bound(string, chars)]
    return rstrip_impl


@overload_method(types.UnicodeType, 'strip')
def unicode_strip(string, chars=None):

    if isinstance(chars, types.UnicodeCharSeq):
        def strip_impl(string, chars=None):
            return string.strip(str(chars))
        return strip_impl

    unicode_strip_types_check(chars)

    def strip_impl(string, chars=None):
        lb = unicode_strip_left_bound(string, chars)
        rb = unicode_strip_right_bound(string, chars)
        return string[lb:rb]
    return strip_impl


# ------------------------------------------------------------------------------
# Slice functions
# ------------------------------------------------------------------------------

@register_jitable
def normalize_str_idx(idx, length, is_start=True):
    """
    Parameters
    ----------
    idx : int or None
        the index
    length : int
        the string length
    is_start : bool; optional with defaults to True
        Is it the *start* or the *stop* of the slice?

    Returns
    -------
    norm_idx : int
        normalized index
    """
    if idx is None:
        if is_start:
            return 0
        else:
            return length
    elif idx < 0:
        idx += length

    if idx < 0 or idx >= length:
        raise IndexError("string index out of range")

    return idx


@register_jitable
def _normalize_slice_idx_count(arg, slice_len, default):
    """
    Used for unicode_count

    If arg < -slice_len, returns 0 (prevents circle)

    If arg is within slice, e.g -slice_len <= arg < slice_len
    returns its real index via arg % slice_len

    If arg > slice_len, returns arg (in this case count must
    return 0 if it is start index)
    """

    if arg is None:
        return default
    if -slice_len <= arg < slice_len:
        return arg % slice_len
    return 0 if arg < 0 else arg


@intrinsic
def _normalize_slice(typingctx, sliceobj, length):
    """Fix slice object.
    """
    sig = sliceobj(sliceobj, length)

    def codegen(context, builder, sig, args):
        [slicetype, lengthtype] = sig.args
        [sliceobj, length] = args
        slice = context.make_helper(builder, slicetype, sliceobj)
        slicing.guard_invalid_slice(context, builder, slicetype, slice)
        slicing.fix_slice(builder, slice, length)
        return slice._getvalue()

    return sig, codegen


@intrinsic
def _slice_span(typingctx, sliceobj):
    """Compute the span from the given slice object.
    """
    sig = types.intp(sliceobj)

    def codegen(context, builder, sig, args):
        [slicetype] = sig.args
        [sliceobj] = args
        slice = context.make_helper(builder, slicetype, sliceobj)
        result_size = slicing.get_slice_length(builder, slice)
        return result_size

    return sig, codegen


@register_jitable(_nrt=False)
def _strncpy(dst, dst_offset, src, src_offset, n):
    if src._kind == dst._kind:
        byte_width = _kind_to_byte_width(src._kind)
        src_byte_offset = byte_width * src_offset
        dst_byte_offset = byte_width * dst_offset
        nbytes = n * byte_width
        memcpy_region(dst._data, dst_byte_offset, src._data,
                      src_byte_offset, nbytes, align=1)
    else:
        for i in range(n):
            _set_code_point(dst, dst_offset + i,
                            _get_code_point(src, src_offset + i))


@intrinsic
def _get_str_slice_view(typingctx, src_t, start_t, length_t):
    """Create a slice of a unicode string using a view of its data to avoid
    extra allocation.
    """
    assert src_t == types.unicode_type

    def codegen(context, builder, sig, args):
        src, start, length = args
        in_str = cgutils.create_struct_proxy(
            types.unicode_type)(context, builder, value=src)
        view_str = cgutils.create_struct_proxy(
            types.unicode_type)(context, builder)
        view_str.meminfo = in_str.meminfo
        view_str.kind = in_str.kind
        view_str.is_ascii = in_str.is_ascii
        view_str.length = length
        # hash value -1 to indicate "need to compute hash"
        view_str.hash = context.get_constant(_Py_hash_t, -1)
        # get a pointer to start of slice data
        bw_typ = context.typing_context.resolve_value_type(_kind_to_byte_width)
        bw_sig = bw_typ.get_call_type(
            context.typing_context, (types.int32,), {})
        bw_impl = context.get_function(bw_typ, bw_sig)
        byte_width = bw_impl(builder, (in_str.kind,))
        offset = builder.mul(start, byte_width)
        view_str.data = builder.gep(in_str.data, [offset])
        # Set parent pyobject to NULL
        view_str.parent = cgutils.get_null_value(view_str.parent.type)
        # incref original string
        if context.enable_nrt:
            context.nrt.incref(builder, sig.args[0], src)
        return view_str._getvalue()

    sig = types.unicode_type(types.unicode_type, types.intp, types.intp)
    return sig, codegen


@overload(operator.getitem)
def unicode_getitem(s, idx):
    if isinstance(s, types.UnicodeType):
        if isinstance(idx, types.Integer):
            def getitem_char(s, idx):
                idx = normalize_str_idx(idx, len(s))
                cp = _get_code_point(s, idx)
                kind = _codepoint_to_kind(cp)
                if kind == s._kind:
                    return _get_str_slice_view(s, idx, 1)
                else:
                    is_ascii = _codepoint_is_ascii(cp)
                    ret = _empty_string(kind, 1, is_ascii)
                    _set_code_point(ret, 0, cp)
                    return ret
            return getitem_char
        elif isinstance(idx, types.SliceType):
            def getitem_slice(s, idx):
                slice_idx = _normalize_slice(idx, len(s))
                span = _slice_span(slice_idx)

                cp = _get_code_point(s, slice_idx.start)
                kind = _codepoint_to_kind(cp)
                is_ascii = _codepoint_is_ascii(cp)

                # Check slice to see if it's homogeneous in kind
                for i in range(slice_idx.start + slice_idx.step,
                               slice_idx.stop, slice_idx.step):
                    cp = _get_code_point(s, i)
                    is_ascii &= _codepoint_is_ascii(cp)
                    new_kind = _codepoint_to_kind(cp)
                    if kind != new_kind:
                        kind = _pick_kind(kind, new_kind)
                    # TODO: it might be possible to break here if the kind
                    # is PY_UNICODE_4BYTE_KIND but there are potentially
                    # strings coming from other internal functions that are
                    # this wide and also actually ASCII (i.e. kind is larger
                    # than actually required for storing the code point), so
                    # it's necessary to continue.

                if slice_idx.step == 1 and kind == s._kind:
                    # Can return a view, the slice has the same kind as the
                    # string itself and it's a stride slice 1.
                    return _get_str_slice_view(s, slice_idx.start, span)
                else:
                    # It's heterogeneous in kind OR stride != 1
                    ret = _empty_string(kind, span, is_ascii)
                    cur = slice_idx.start
                    for i in range(span):
                        _set_code_point(ret, i, _get_code_point(s, cur))
                        cur += slice_idx.step
                    return ret

            return getitem_slice


# ------------------------------------------------------------------------------
# String operations
# ------------------------------------------------------------------------------


@overload(operator.add)
@overload(operator.iadd)
def unicode_concat(a, b):
    if isinstance(a, types.UnicodeType) and isinstance(b, types.UnicodeType):
        def concat_impl(a, b):
            new_length = a._length + b._length
            new_kind = _pick_kind(a._kind, b._kind)
            new_ascii = _pick_ascii(a._is_ascii, b._is_ascii)
            result = _empty_string(new_kind, new_length, new_ascii)
            for i in range(len(a)):
                _set_code_point(result, i, _get_code_point(a, i))
            for j in range(len(b)):
                _set_code_point(result, len(a) + j, _get_code_point(b, j))
            return result
        return concat_impl

    if isinstance(a, types.UnicodeType) and isinstance(b, types.UnicodeCharSeq):
        def concat_impl(a, b):
            return a + str(b)
        return concat_impl


@register_jitable
def _repeat_impl(str_arg, mult_arg):
    if str_arg == '' or mult_arg < 1:
        return ''
    elif mult_arg == 1:
        return str_arg
    else:
        new_length = str_arg._length * mult_arg
        new_kind = str_arg._kind
        result = _empty_string(new_kind, new_length, str_arg._is_ascii)
        # make initial copy into result
        len_a = len(str_arg)
        _strncpy(result, 0, str_arg, 0, len_a)
        # loop through powers of 2 for efficient copying
        copy_size = len_a
        while 2 * copy_size <= new_length:
            _strncpy(result, copy_size, result, 0, copy_size)
            copy_size *= 2

        if not 2 * copy_size == new_length:
            # if copy_size not an exact multiple it then needs
            # to complete the rest of the copies
            rest = new_length - copy_size
            _strncpy(result, copy_size, result, copy_size - rest, rest)
            return result


@overload(operator.mul)
def unicode_repeat(a, b):
    if isinstance(a, types.UnicodeType) and isinstance(b, types.Integer):
        def wrap(a, b):
            return _repeat_impl(a, b)
        return wrap
    elif isinstance(a, types.Integer) and isinstance(b, types.UnicodeType):
        def wrap(a, b):
            return _repeat_impl(b, a)
        return wrap


@overload(operator.not_)
def unicode_not(a):
    if isinstance(a, types.UnicodeType):
        def impl(a):
            return len(a) == 0
        return impl


@overload_method(types.UnicodeType, 'replace')
def unicode_replace(s, old_str, new_str, count=-1):
    thety = count
    if isinstance(count, types.Omitted):
        thety = count.value
    elif isinstance(count, types.Optional):
        thety = count.type

    if not isinstance(thety, (int, types.Integer)):
        raise TypingError('Unsupported parameters. The parameters '
                          'must be Integer. Given count: {}'.format(count))

    if not isinstance(old_str, (types.UnicodeType, types.NoneType)):
        raise TypingError('The object must be a UnicodeType.'
                          ' Given: {}'.format(old_str))

    if not isinstance(new_str, types.UnicodeType):
        raise TypingError('The object must be a UnicodeType.'
                          ' Given: {}'.format(new_str))

    def impl(s, old_str, new_str, count=-1):
        if count == 0:
            return s
        if old_str == '':
            schars = list(s)
            if count == -1:
                return new_str + new_str.join(schars) + new_str
            split_result = [new_str]
            min_count = min(len(schars), count)
            for i in range(min_count):
                split_result.append(schars[i])
                if i + 1 != min_count:
                    split_result.append(new_str)
                else:
                    split_result.append(''.join(schars[(i + 1):]))
            if count > len(schars):
                split_result.append(new_str)
            return ''.join(split_result)
        schars = s.split(old_str, count)
        result = new_str.join(schars)
        return result

    return impl

# ------------------------------------------------------------------------------
# String `is*()` methods
# ------------------------------------------------------------------------------


# generates isalpha/isalnum
def gen_isAlX(ascii_func, unicode_func):
    def unicode_isAlX(data):

        def impl(data):
            length = len(data)
            if length == 0:
                return False

            if length == 1:
                code_point = _get_code_point(data, 0)
                if data._is_ascii:
                    return ascii_func(code_point)
                else:
                    return unicode_func(code_point)

            if data._is_ascii:
                for i in range(length):
                    code_point = _get_code_point(data, i)
                    if not ascii_func(code_point):
                        return False

            for i in range(length):
                code_point = _get_code_point(data, i)
                if not unicode_func(code_point):
                    return False

            return True

        return impl
    return unicode_isAlX


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11928-L11964    # noqa: E501
overload_method(types.UnicodeType, 'isalpha')(gen_isAlX(_Py_ISALPHA,
                                                        _PyUnicode_IsAlpha))

_unicode_is_alnum = register_jitable(lambda x:
                                     (_PyUnicode_IsNumeric(x) or
                                      _PyUnicode_IsAlpha(x)))

# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11975-L12006    # noqa: E501
overload_method(types.UnicodeType, 'isalnum')(gen_isAlX(_Py_ISALNUM,
                                                        _unicode_is_alnum))


def _is_upper(is_lower, is_upper, is_title):
    # impl is an approximate translation of:
    # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11794-L11827    # noqa: E501
    # mixed with:
    # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/bytes_methods.c#L218-L242    # noqa: E501
    def impl(a):
        l = len(a)
        if l == 1:
            return is_upper(_get_code_point(a, 0))
        if l == 0:
            return False
        cased = False
        for idx in range(l):
            code_point = _get_code_point(a, idx)
            if is_lower(code_point) or is_title(code_point):
                return False
            elif (not cased and is_upper(code_point)):
                cased = True
        return cased
    return impl


_always_false = register_jitable(lambda x: False)
_ascii_is_upper = register_jitable(_is_upper(_Py_ISLOWER, _Py_ISUPPER,
                                             _always_false))
_unicode_is_upper = register_jitable(_is_upper(_PyUnicode_IsLowercase,
                                               _PyUnicode_IsUppercase,
                                               _PyUnicode_IsTitlecase))


@overload_method(types.UnicodeType, 'isupper')
def unicode_isupper(a):
    """
    Implements .isupper()
    """
    def impl(a):
        if a._is_ascii:
            return _ascii_is_upper(a)
        else:
            return _unicode_is_upper(a)
    return impl


@overload_method(types.UnicodeType, 'isascii')
def unicode_isascii(data):
    """Implements UnicodeType.isascii()"""

    def impl(data):
        return data._is_ascii
    return impl


@overload_method(types.UnicodeType, 'istitle')
def unicode_istitle(data):
    """
    Implements UnicodeType.istitle()
    The algorithm is an approximate translation from CPython:
    https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11829-L11885 # noqa: E501
    """

    def impl(data):
        length = len(data)
        if length == 1:
            char = _get_code_point(data, 0)
            return _PyUnicode_IsUppercase(char) or _PyUnicode_IsTitlecase(char)

        if length == 0:
            return False

        cased = False
        previous_is_cased = False
        for idx in range(length):
            char = _get_code_point(data, idx)
            if _PyUnicode_IsUppercase(char) or _PyUnicode_IsTitlecase(char):
                if previous_is_cased:
                    return False
                previous_is_cased = True
                cased = True
            elif _PyUnicode_IsLowercase(char):
                if not previous_is_cased:
                    return False
                previous_is_cased = True
                cased = True
            else:
                previous_is_cased = False

        return cased
    return impl


@overload_method(types.UnicodeType, 'islower')
def unicode_islower(data):
    """
    impl is an approximate translation of:
    https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L11900-L11933    # noqa: E501
    mixed with:
    https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/bytes_methods.c#L131-L156    # noqa: E501
    """

    def impl(data):
        length = len(data)
        if length == 1:
            return _PyUnicode_IsLowercase(_get_code_point(data, 0))
        if length == 0:
            return False

        cased = False
        for idx in range(length):
            cp = _get_code_point(data, idx)
            if _PyUnicode_IsUppercase(cp) or _PyUnicode_IsTitlecase(cp):
                return False
            elif not cased and _PyUnicode_IsLowercase(cp):
                cased = True
        return cased
    return impl


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12126-L12161    # noqa: E501
@overload_method(types.UnicodeType, 'isidentifier')
def unicode_isidentifier(data):
    """Implements UnicodeType.isidentifier()"""

    def impl(data):
        length = len(data)
        if length == 0:
            return False

        first_cp = _get_code_point(data, 0)
        if not _PyUnicode_IsXidStart(first_cp) and first_cp != 0x5F:
            return False

        for i in range(1, length):
            code_point = _get_code_point(data, i)
            if not _PyUnicode_IsXidContinue(code_point):
                return False

        return True

    return impl


# generator for simple unicode "isX" methods
def gen_isX(_PyUnicode_IS_func, empty_is_false=True):
    def unicode_isX(data):
        def impl(data):
            length = len(data)
            if length == 1:
                return _PyUnicode_IS_func(_get_code_point(data, 0))

            if empty_is_false and length == 0:
                return False

            for i in range(length):
                code_point = _get_code_point(data, i)
                if not _PyUnicode_IS_func(code_point):
                    return False

            return True

        return impl
    return unicode_isX


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L11896-L11925    # noqa: E501
overload_method(types.UnicodeType, 'isspace')(gen_isX(_PyUnicode_IsSpace))

# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12096-L12124    # noqa: E501
overload_method(types.UnicodeType, 'isnumeric')(gen_isX(_PyUnicode_IsNumeric))

# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12056-L12085    # noqa: E501
overload_method(types.UnicodeType, 'isdigit')(gen_isX(_PyUnicode_IsDigit))

# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12017-L12045    # noqa: E501
overload_method(types.UnicodeType, 'isdecimal')(
    gen_isX(_PyUnicode_IsDecimalDigit))

# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L12188-L12213    # noqa: E501
overload_method(types.UnicodeType, 'isprintable')(
    gen_isX(_PyUnicode_IsPrintable, False))

# ------------------------------------------------------------------------------
# String methods that apply a transformation to the characters themselves
# ------------------------------------------------------------------------------


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L9863-L9908    # noqa: E501
def case_operation(ascii_func, unicode_func):
    """Generate common case operation performer."""
    def impl(data):
        length = len(data)
        if length == 0:
            return _empty_string(data._kind, length, data._is_ascii)

        if data._is_ascii:
            res = _empty_string(data._kind, length, 1)
            ascii_func(data, res)
            return res

        # https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L9863-L9908    # noqa: E501
        tmp = _empty_string(PY_UNICODE_4BYTE_KIND, 3 * length, data._is_ascii)
        # maxchar should be inside of a list to be pass as argument by reference
        maxchars = [0]
        newlength = unicode_func(data, length, tmp, maxchars)
        maxchar = maxchars[0]
        newkind = _codepoint_to_kind(maxchar)
        res = _empty_string(newkind, newlength, _codepoint_is_ascii(maxchar))
        for i in range(newlength):
            _set_code_point(res, i, _get_code_point(tmp, i))

        return res

    return impl


# https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L9856-L9883    # noqa: E501
@register_jitable
def _handle_capital_sigma(data, length, idx):
    """This is a translation of the function that handles the capital sigma."""
    c = 0
    j = idx - 1
    while j >= 0:
        c = _get_code_point(data, j)
        if not _PyUnicode_IsCaseIgnorable(c):
            break
        j -= 1
    final_sigma = (j >= 0 and _PyUnicode_IsCased(c))
    if final_sigma:
        j = idx + 1
        while j < length:
            c = _get_code_point(data, j)
            if not _PyUnicode_IsCaseIgnorable(c):
                break
            j += 1
        final_sigma = (j == length or (not _PyUnicode_IsCased(c)))

    return 0x3c2 if final_sigma else 0x3c3


# https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L9885-L9895    # noqa: E501
@register_jitable
def _lower_ucs4(code_point, data, length, idx, mapped):
    """This is a translation of the function that lowers a character."""
    if code_point == 0x3A3:
        mapped[0] = _handle_capital_sigma(data, length, idx)
        return 1
    return _PyUnicode_ToLowerFull(code_point, mapped)


# https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L9946-L9965    # noqa: E501
def _gen_unicode_upper_or_lower(lower):
    def _do_upper_or_lower(data, length, res, maxchars):
        k = 0
        for idx in range(length):
            mapped = np.zeros(3, dtype=_Py_UCS4)
            code_point = _get_code_point(data, idx)
            if lower:
                n_res = _lower_ucs4(code_point, data, length, idx, mapped)
            else:
                # might be needed if call _do_upper_or_lower in unicode_upper
                n_res = _PyUnicode_ToUpperFull(code_point, mapped)
            for m in mapped[:n_res]:
                maxchars[0] = max(maxchars[0], m)
                _set_code_point(res, k, m)
                k += 1
        return k
    return _do_upper_or_lower


_unicode_upper = register_jitable(_gen_unicode_upper_or_lower(False))
_unicode_lower = register_jitable(_gen_unicode_upper_or_lower(True))


def _gen_ascii_upper_or_lower(func):
    def _ascii_upper_or_lower(data, res):
        for idx in range(len(data)):
            code_point = _get_code_point(data, idx)
            _set_code_point(res, idx, func(code_point))
    return _ascii_upper_or_lower


_ascii_upper = register_jitable(_gen_ascii_upper_or_lower(_Py_TOUPPER))
_ascii_lower = register_jitable(_gen_ascii_upper_or_lower(_Py_TOLOWER))


@overload_method(types.UnicodeType, 'lower')
def unicode_lower(data):
    """Implements .lower()"""
    return case_operation(_ascii_lower, _unicode_lower)


@overload_method(types.UnicodeType, 'upper')
def unicode_upper(data):
    """Implements .upper()"""
    return case_operation(_ascii_upper, _unicode_upper)


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L9819-L9834    # noqa: E501
@register_jitable
def _unicode_casefold(data, length, res, maxchars):
    k = 0
    mapped = np.zeros(3, dtype=_Py_UCS4)
    for idx in range(length):
        mapped.fill(0)
        code_point = _get_code_point(data, idx)
        n_res = _PyUnicode_ToFoldedFull(code_point, mapped)
        for m in mapped[:n_res]:
            maxchar = maxchars[0]
            maxchars[0] = max(maxchar, m)
            _set_code_point(res, k, m)
            k += 1

    return k


@register_jitable
def _ascii_casefold(data, res):
    for idx in range(len(data)):
        code_point = _get_code_point(data, idx)
        _set_code_point(res, idx, _Py_TOLOWER(code_point))


@overload_method(types.UnicodeType, 'casefold')
def unicode_casefold(data):
    """Implements str.casefold()"""
    return case_operation(_ascii_casefold, _unicode_casefold)


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L9737-L9759    # noqa: E501
@register_jitable
def _unicode_capitalize(data, length, res, maxchars):
    k = 0
    maxchar = 0
    mapped = np.zeros(3, dtype=_Py_UCS4)
    code_point = _get_code_point(data, 0)

    n_res = _PyUnicode_ToTitleFull(code_point, mapped)

    for m in mapped[:n_res]:
        maxchar = max(maxchar, m)
        _set_code_point(res, k, m)
        k += 1
    for idx in range(1, length):
        mapped.fill(0)
        code_point = _get_code_point(data, idx)
        n_res = _lower_ucs4(code_point, data, length, idx, mapped)
        for m in mapped[:n_res]:
            maxchar = max(maxchar, m)
            _set_code_point(res, k, m)
            k += 1
    maxchars[0] = maxchar
    return k


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/bytes_methods.c#L361-L382    # noqa: E501
@register_jitable
def _ascii_capitalize(data, res):
    code_point = _get_code_point(data, 0)
    _set_code_point(res, 0, _Py_TOUPPER(code_point))
    for idx in range(1, len(data)):
        code_point = _get_code_point(data, idx)
        _set_code_point(res, idx, _Py_TOLOWER(code_point))


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L10765-L10774    # noqa: E501
@overload_method(types.UnicodeType, 'capitalize')
def unicode_capitalize(data):
    return case_operation(_ascii_capitalize, _unicode_capitalize)


# https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L9996-L10021    # noqa: E501
@register_jitable
def _unicode_title(data, length, res, maxchars):
    """This is a translation of the function that titles a unicode string."""
    k = 0
    previous_cased = False
    mapped = np.empty(3, dtype=_Py_UCS4)
    for idx in range(length):
        mapped.fill(0)
        code_point = _get_code_point(data, idx)
        if previous_cased:
            n_res = _lower_ucs4(code_point, data, length, idx, mapped)
        else:
            n_res = _PyUnicode_ToTitleFull(_Py_UCS4(code_point), mapped)
        for m in mapped[:n_res]:
            maxchar, = maxchars
            maxchars[0] = max(maxchar, m)
            _set_code_point(res, k, m)
            k += 1
        previous_cased = _PyUnicode_IsCased(_Py_UCS4(code_point))
    return k


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/bytes_methods.c#L332-L352    # noqa: E501
@register_jitable
def _ascii_title(data, res):
    """ Does .title() on an ASCII string """
    previous_is_cased = False
    for idx in range(len(data)):
        code_point = _get_code_point(data, idx)
        if _Py_ISLOWER(code_point):
            if not previous_is_cased:
                code_point = _Py_TOUPPER(code_point)
            previous_is_cased = True
        elif _Py_ISUPPER(code_point):
            if previous_is_cased:
                code_point = _Py_TOLOWER(code_point)
            previous_is_cased = True
        else:
            previous_is_cased = False
        _set_code_point(res, idx, code_point)


# https://github.com/python/cpython/blob/201c8f79450628241574fba940e08107178dc3a5/Objects/unicodeobject.c#L10023-L10069    # noqa: E501
@overload_method(types.UnicodeType, 'title')
def unicode_title(data):
    """Implements str.title()"""
    # https://docs.python.org/3/library/stdtypes.html#str.title
    return case_operation(_ascii_title, _unicode_title)


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/bytes_methods.c#L391-L408    # noqa: E501
@register_jitable
def _ascii_swapcase(data, res):
    for idx in range(len(data)):
        code_point = _get_code_point(data, idx)
        if _Py_ISUPPER(code_point):
            code_point = _Py_TOLOWER(code_point)
        elif _Py_ISLOWER(code_point):
            code_point = _Py_TOUPPER(code_point)
        _set_code_point(res, idx, code_point)


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L9761-L9784    # noqa: E501
@register_jitable
def _unicode_swapcase(data, length, res, maxchars):
    k = 0
    maxchar = 0
    mapped = np.empty(3, dtype=_Py_UCS4)
    for idx in range(length):
        mapped.fill(0)
        code_point = _get_code_point(data, idx)
        if _PyUnicode_IsUppercase(code_point):
            n_res = _lower_ucs4(code_point, data, length, idx, mapped)
        elif _PyUnicode_IsLowercase(code_point):
            n_res = _PyUnicode_ToUpperFull(code_point, mapped)
        else:
            n_res = 1
            mapped[0] = code_point
        for m in mapped[:n_res]:
            maxchar = max(maxchar, m)
            _set_code_point(res, k, m)
            k += 1
    maxchars[0] = maxchar
    return k


@overload_method(types.UnicodeType, 'swapcase')
def unicode_swapcase(data):
    return case_operation(_ascii_swapcase, _unicode_swapcase)


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Python/bltinmodule.c#L1781-L1824    # noqa: E501
@overload(ord)
def ol_ord(c):
    if isinstance(c, types.UnicodeType):
        def impl(c):
            lc = len(c)
            if lc != 1:
                # CPython does TypeError
                raise TypeError("ord() expected a character")
            return _get_code_point(c, 0)
        return impl


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L2005-L2028    # noqa: E501
# This looks a bit different to the cpython implementation but, with the
# exception of a latin1 fast path is logically the same. It finds the "kind" of
# the codepoint `ch`, creates a length 1 string of that kind and then injects
# the code point into the zero position of that string. Cpython does similar but
# branches for each kind (this is encapsulated in Numba's _set_code_point).
@register_jitable
def _unicode_char(ch):
    assert ch <= _MAX_UNICODE
    kind = _codepoint_to_kind(ch)
    ret = _empty_string(kind, 1, kind == PY_UNICODE_1BYTE_KIND)
    _set_code_point(ret, 0, ch)
    return ret


_out_of_range_msg = "chr() arg not in range(0x%hx)" % _MAX_UNICODE


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Objects/unicodeobject.c#L3045-L3055    # noqa: E501
@register_jitable
def _PyUnicode_FromOrdinal(ordinal):
    if (ordinal < 0 or ordinal > _MAX_UNICODE):
        raise ValueError(_out_of_range_msg)

    return _unicode_char(_Py_UCS4(ordinal))


# https://github.com/python/cpython/blob/1d4b6ba19466aba0eb91c4ba01ba509acf18c723/Python/bltinmodule.c#L715-L720    # noqa: E501
@overload(chr)
def ol_chr(i):
    if isinstance(i, types.Integer):
        def impl(i):
            return _PyUnicode_FromOrdinal(i)
        return impl


@overload_method(types.UnicodeType, "__str__")
def unicode_str(s):
    return lambda s: s


@overload_method(types.UnicodeType, "__repr__")
def unicode_repr(s):
    # Can't use f-string as the impl ends up calling str and then repr, which
    # then recurses somewhere in imports.
    return lambda s: "'" + s + "'"


@overload_method(types.Integer, "__str__")
def integer_str(n):

    ten = n(10)

    def impl(n):
        flag = False
        if n < 0:
            n = -n
            flag = True
        if n == 0:
            return '0'
        length = flag + 1 + int(np.floor(np.log10(n)))
        kind = PY_UNICODE_1BYTE_KIND
        char_width = _kind_to_byte_width(kind)
        s = _malloc_string(kind, char_width, length, True)
        if flag:
            _set_code_point(s, 0, ord('-'))
        idx = length - 1
        while n > 0:
            n, digit = divmod(n, ten)
            c = ord('0') + digit
            _set_code_point(s, idx, c)
            idx -= 1
        return s
    return impl


@overload_method(types.Integer, "__repr__")
def integer_repr(n):
    return lambda n: n.__str__()


@overload_method(types.Boolean, "__repr__")
@overload_method(types.Boolean, "__str__")
def boolean_str(b):
    return lambda b: 'True' if b else 'False'


# ------------------------------------------------------------------------------
# iteration
# ------------------------------------------------------------------------------


@lower_builtin('getiter', types.UnicodeType)
def getiter_unicode(context, builder, sig, args):
    [ty] = sig.args
    [data] = args

    iterobj = context.make_helper(builder, sig.return_type)

    # set the index to zero
    zero = context.get_constant(types.uintp, 0)
    indexptr = cgutils.alloca_once_value(builder, zero)

    iterobj.index = indexptr

    # wire in the unicode type data
    iterobj.data = data

    # incref as needed
    if context.enable_nrt:
        context.nrt.incref(builder, ty, data)

    res = iterobj._getvalue()
    return impl_ret_new_ref(context, builder, sig.return_type, res)


@lower_builtin('iternext', types.UnicodeIteratorType)
# a new ref counted object is put into result._yield so set the new_ref to True!
@iternext_impl(RefType.NEW)
def iternext_unicode(context, builder, sig, args, result):
    [iterty] = sig.args
    [iter] = args

    tyctx = context.typing_context

    # get ref to unicode.__getitem__
    fnty = tyctx.resolve_value_type(operator.getitem)
    getitem_sig = fnty.get_call_type(tyctx, (types.unicode_type, types.uintp),
                                     {})
    getitem_impl = context.get_function(fnty, getitem_sig)

    # get ref to unicode.__len__
    fnty = tyctx.resolve_value_type(len)
    len_sig = fnty.get_call_type(tyctx, (types.unicode_type,), {})
    len_impl = context.get_function(fnty, len_sig)

    # grab unicode iterator struct
    iterobj = context.make_helper(builder, iterty, value=iter)

    # find the length of the string
    strlen = len_impl(builder, (iterobj.data,))

    # find the current index
    index = builder.load(iterobj.index)

    # see if the index is in range
    is_valid = builder.icmp_unsigned('<', index, strlen)
    result.set_valid(is_valid)

    with builder.if_then(is_valid):
        # return value at index
        gotitem = getitem_impl(builder, (iterobj.data, index,))
        result.yield_(gotitem)

        # bump index for next cycle
        nindex = cgutils.increment_index(builder, index)
        builder.store(nindex, iterobj.index)
