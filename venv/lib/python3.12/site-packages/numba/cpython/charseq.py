"""Implements operations on bytes and str (unicode) array items."""
import operator
import numpy as np
from llvmlite import ir

from numba.core import types, cgutils, config
from numba.core.extending import (overload, intrinsic, overload_method,
                                  lower_cast, register_jitable)
from numba.core.cgutils import is_nonelike
from numba.cpython import unicode

# bytes and str arrays items are of type CharSeq and UnicodeCharSeq,
# respectively.  See numpy/types/npytypes.py for CharSeq,
# UnicodeCharSeq definitions.  The corresponding data models are
# defined in numpy/datamodel/models.py. Boxing/unboxing of item types
# are defined in numpy/targets/boxing.py, see box_unicodecharseq,
# unbox_unicodecharseq, box_charseq, unbox_charseq.

s1_dtype = np.dtype('S1')
assert s1_dtype.itemsize == 1
if config.USE_LEGACY_TYPE_SYSTEM:
    bytes_type = types.Bytes(types.uint8, 1, "C", readonly=True)
else:
    bytes_type = types.Bytes(types.c_uint8, 1, "C", readonly=True)

# Currently, NumPy supports only UTF-32 arrays but this may change in
# future and the approach used here for supporting str arrays may need
# a revision depending on how NumPy will support UTF-8 and UTF-16
# arrays.
u1_dtype = np.dtype('U1')
unicode_byte_width = u1_dtype.itemsize
unicode_uint = {1: np.uint8, 2: np.uint16, 4: np.uint32}[unicode_byte_width]
unicode_kind = {1: unicode.PY_UNICODE_1BYTE_KIND,
                2: unicode.PY_UNICODE_2BYTE_KIND,
                4: unicode.PY_UNICODE_4BYTE_KIND}[unicode_byte_width]


# this is modified version of numba.unicode.make_deref_codegen
def make_deref_codegen(bitsize):
    def codegen(context, builder, signature, args):
        data, idx = args
        rawptr = cgutils.alloca_once_value(builder, value=data)
        ptr = builder.bitcast(rawptr, ir.IntType(bitsize).as_pointer())
        ch = builder.load(builder.gep(ptr, [idx]))
        return builder.zext(ch, ir.IntType(32))
    return codegen


@intrinsic
def deref_uint8(typingctx, data, offset):
    sig = types.uint32(data, types.intp)
    return sig, make_deref_codegen(8)


@intrinsic
def deref_uint16(typingctx, data, offset):
    sig = types.uint32(data, types.intp)
    return sig, make_deref_codegen(16)


@intrinsic
def deref_uint32(typingctx, data, offset):
    sig = types.uint32(data, types.intp)
    return sig, make_deref_codegen(32)


@register_jitable(_nrt=False)
def charseq_get_code(a, i):
    """Access i-th item of CharSeq object via code value
    """
    return deref_uint8(a, i)


@register_jitable
def charseq_get_value(a, i):
    """Access i-th item of CharSeq object via code value.

    null code is interpreted as IndexError
    """
    code = charseq_get_code(a, i)
    if code == 0:
        raise IndexError('index out of range')
    return code


@register_jitable(_nrt=False)
def unicode_charseq_get_code(a, i):
    """Access i-th item of UnicodeCharSeq object via code value
    """
    if unicode_byte_width == 4:
        return deref_uint32(a, i)
    elif unicode_byte_width == 2:
        return deref_uint16(a, i)
    elif unicode_byte_width == 1:
        return deref_uint8(a, i)
    else:
        raise NotImplementedError(
            'unicode_charseq_get_code: unicode_byte_width not in [1, 2, 4]')


@register_jitable
def unicode_get_code(a, i):
    """Access i-th item of UnicodeType object.
    """
    return unicode._get_code_point(a, i)


@register_jitable
def bytes_get_code(a, i):
    """Access i-th item of Bytes object.
        """
    return a[i]


def _get_code_impl(a):
    if isinstance(a, types.CharSeq):
        return charseq_get_code
    elif isinstance(a, types.Bytes):
        return bytes_get_code
    elif isinstance(a, types.UnicodeCharSeq):
        return unicode_charseq_get_code
    elif isinstance(a, types.UnicodeType):
        return unicode_get_code


def _same_kind(a, b):
    for t in [(types.CharSeq, types.Bytes),
              (types.UnicodeCharSeq, types.UnicodeType)]:
        if isinstance(a, t) and isinstance(b, t):
            return True
    return False


def _is_bytes(a):
    return isinstance(a, (types.CharSeq, types.Bytes))


def is_default(x, default):
    return x == default or isinstance(x, types.Omitted)


@register_jitable
def unicode_charseq_get_value(a, i):
    """Access i-th item of UnicodeCharSeq object via unicode value

    null code is interpreted as IndexError
    """
    code = unicode_charseq_get_code(a, i)
    if code == 0:
        raise IndexError('index out of range')
    # Return numpy equivalent of `chr(code)`
    return np.array(code, unicode_uint).view(u1_dtype)[()]


#
# CAST
#
# Currently, the following casting operations are supported:
#   Bytes -> CharSeq                 (ex: a=np.array(b'abc'); a[()] = b'123')
#   UnicodeType -> UnicodeCharSeq    (ex: a=np.array('abc'); a[()] = '123')
#   CharSeq -> Bytes                 (ex: a=np.array(b'abc'); b = bytes(a[()]))
#   UnicodeType -> Bytes             (ex: str('123')._to_bytes())
#
# The following casting operations can be implemented when required:
#   Bytes -> UnicodeCharSeq   (ex: a=np.array('abc'); a[()] = b'123')
#   UnicodeType -> CharSeq    (ex: a=np.array(b'abc'); a[()] = '123')
#   UnicodeType -> Bytes      (ex: bytes('123', 'utf8'))
#


@lower_cast(types.Bytes, types.CharSeq)
def bytes_to_charseq(context, builder, fromty, toty, val):
    barr = cgutils.create_struct_proxy(fromty)(context, builder, value=val)
    src = builder.bitcast(barr.data, ir.IntType(8).as_pointer())
    src_length = barr.nitems

    lty = context.get_value_type(toty)
    dstint_t = ir.IntType(8)
    dst_ptr = cgutils.alloca_once(builder, lty)
    dst = builder.bitcast(dst_ptr, dstint_t.as_pointer())

    dst_length = ir.Constant(src_length.type, toty.count)
    is_shorter_value = builder.icmp_unsigned('<', src_length, dst_length)
    count = builder.select(is_shorter_value, src_length, dst_length)
    with builder.if_then(is_shorter_value):
        cgutils.memset(builder,
                       dst,
                       ir.Constant(src_length.type,
                                   toty.count), 0)
    with cgutils.for_range(builder, count) as loop:
        in_ptr = builder.gep(src, [loop.index])
        in_val = builder.zext(builder.load(in_ptr), dstint_t)
        builder.store(in_val, builder.gep(dst, [loop.index]))

    return builder.load(dst_ptr)


def _make_constant_bytes(context, builder, nbytes):
    bstr_ctor = cgutils.create_struct_proxy(bytes_type)
    bstr = bstr_ctor(context, builder)

    if isinstance(nbytes, int):
        nbytes = ir.Constant(bstr.nitems.type, nbytes)

    bstr.meminfo = context.nrt.meminfo_alloc(builder, nbytes)
    bstr.nitems = nbytes
    bstr.itemsize = ir.Constant(bstr.itemsize.type, 1)
    bstr.data = context.nrt.meminfo_data(builder, bstr.meminfo)
    bstr.parent = cgutils.get_null_value(bstr.parent.type)
    # bstr.shape and bstr.strides are not used
    bstr.shape = cgutils.get_null_value(bstr.shape.type)
    bstr.strides = cgutils.get_null_value(bstr.strides.type)
    return bstr


@lower_cast(types.CharSeq, types.Bytes)
def charseq_to_bytes(context, builder, fromty, toty, val):
    bstr = _make_constant_bytes(context, builder, val.type.count)
    rawptr = cgutils.alloca_once_value(builder, value=val)
    ptr = builder.bitcast(rawptr, bstr.data.type)
    cgutils.memcpy(builder, bstr.data, ptr, bstr.nitems)
    return bstr


@lower_cast(types.UnicodeType, types.Bytes)
def unicode_to_bytes_cast(context, builder, fromty, toty, val):
    uni_str = cgutils.create_struct_proxy(fromty)(context, builder, value=val)
    src1 = builder.bitcast(uni_str.data, ir.IntType(8).as_pointer())
    notkind1 = builder.icmp_unsigned('!=', uni_str.kind,
                                     ir.Constant(uni_str.kind.type, 1))
    src_length = uni_str.length

    with builder.if_then(notkind1):
        context.call_conv.return_user_exc(
            builder, ValueError,
            ("cannot cast higher than 8-bit unicode_type to bytes",))

    bstr = _make_constant_bytes(context, builder, src_length)
    cgutils.memcpy(builder, bstr.data, src1, bstr.nitems)
    return bstr


@intrinsic
def _unicode_to_bytes(typingctx, s):
    # used in _to_bytes method
    assert s == types.unicode_type
    sig = bytes_type(s)

    def codegen(context, builder, signature, args):
        return unicode_to_bytes_cast(
            context, builder, s, bytes_type, args[0])._getvalue()
    return sig, codegen


@lower_cast(types.UnicodeType, types.UnicodeCharSeq)
def unicode_to_unicode_charseq(context, builder, fromty, toty, val):
    uni_str = cgutils.create_struct_proxy(fromty)(context, builder, value=val)
    src1 = builder.bitcast(uni_str.data, ir.IntType(8).as_pointer())
    src2 = builder.bitcast(uni_str.data, ir.IntType(16).as_pointer())
    src4 = builder.bitcast(uni_str.data, ir.IntType(32).as_pointer())
    kind1 = builder.icmp_unsigned('==', uni_str.kind,
                                  ir.Constant(uni_str.kind.type, 1))
    kind2 = builder.icmp_unsigned('==', uni_str.kind,
                                  ir.Constant(uni_str.kind.type, 2))
    kind4 = builder.icmp_unsigned('==', uni_str.kind,
                                  ir.Constant(uni_str.kind.type, 4))
    src_length = uni_str.length

    lty = context.get_value_type(toty)
    dstint_t = ir.IntType(8 * unicode_byte_width)
    dst_ptr = cgutils.alloca_once(builder, lty)
    dst = builder.bitcast(dst_ptr, dstint_t.as_pointer())

    dst_length = ir.Constant(src_length.type, toty.count)
    is_shorter_value = builder.icmp_unsigned('<', src_length, dst_length)
    count = builder.select(is_shorter_value, src_length, dst_length)
    with builder.if_then(is_shorter_value):
        cgutils.memset(builder,
                       dst,
                       ir.Constant(src_length.type,
                                   toty.count * unicode_byte_width), 0)

    with builder.if_then(kind1):
        with cgutils.for_range(builder, count) as loop:
            in_ptr = builder.gep(src1, [loop.index])
            in_val = builder.zext(builder.load(in_ptr), dstint_t)
            builder.store(in_val, builder.gep(dst, [loop.index]))

    with builder.if_then(kind2):
        if unicode_byte_width >= 2:
            with cgutils.for_range(builder, count) as loop:
                in_ptr = builder.gep(src2, [loop.index])
                in_val = builder.zext(builder.load(in_ptr), dstint_t)
                builder.store(in_val, builder.gep(dst, [loop.index]))
        else:
            context.call_conv.return_user_exc(
                builder, ValueError,
                ("cannot cast 16-bit unicode_type to %s-bit %s"
                 % (unicode_byte_width * 8, toty)))

    with builder.if_then(kind4):
        if unicode_byte_width >= 4:
            with cgutils.for_range(builder, count) as loop:
                in_ptr = builder.gep(src4, [loop.index])
                in_val = builder.zext(builder.load(in_ptr), dstint_t)
                builder.store(in_val, builder.gep(dst, [loop.index]))
        else:
            context.call_conv.return_user_exc(
                builder, ValueError,
                ("cannot cast 32-bit unicode_type to %s-bit %s"
                 % (unicode_byte_width * 8, toty)))

    return builder.load(dst_ptr)

#
#   Operations on bytes/str array items
#
# Implementation note: while some operations need
# CharSeq/UnicodeCharSeq specific implementations (getitem, len, str,
# etc), many operations can be supported by casting
# CharSeq/UnicodeCharSeq objects to Bytes/UnicodeType objects and
# re-use existing operations.
#
# However, in numba more operations are implemented for UnicodeType
# than for Bytes objects, hence the support for operations with bytes
# array items will be less complete than for str arrays. Although, in
# some cases (hash, contains, etc) the UnicodeType implementations can
# be reused for Bytes objects via using `_to_str` method.
#


@overload(operator.getitem)
def charseq_getitem(s, i):
    get_value = None
    if isinstance(i, types.Integer):
        if isinstance(s, types.CharSeq):
            get_value = charseq_get_value
        if isinstance(s, types.UnicodeCharSeq):
            get_value = unicode_charseq_get_value
    if get_value is not None:
        max_i = s.count
        msg = 'index out of range [0, %s]' % (max_i - 1)

        def getitem_impl(s, i):
            if i < max_i and i >= 0:
                return get_value(s, i)
            raise IndexError(msg)
        return getitem_impl


@overload(len)
def charseq_len(s):
    if isinstance(s, (types.CharSeq, types.UnicodeCharSeq)):
        get_code = _get_code_impl(s)
        n = s.count
        if n == 0:
            def len_impl(s):
                return 0
            return len_impl
        else:
            def len_impl(s):
                # return the index of the last non-null value (numpy
                # behavior)
                i = n
                code = 0
                while code == 0:
                    i = i - 1
                    if i < 0:
                        break
                    code = get_code(s, i)
                return i + 1
            return len_impl


@overload(operator.add)
@overload(operator.iadd)
def charseq_concat(a, b):
    if not _same_kind(a, b):
        return
    if (isinstance(a, types.UnicodeCharSeq) and
            isinstance(b, types.UnicodeType)):
        def impl(a, b):
            return str(a) + b
        return impl
    if (isinstance(b, types.UnicodeCharSeq) and
            isinstance(a, types.UnicodeType)):
        def impl(a, b):
            return a + str(b)
        return impl
    if (isinstance(a, types.UnicodeCharSeq) and
            isinstance(b, types.UnicodeCharSeq)):
        def impl(a, b):
            return str(a) + str(b)
        return impl
    if (isinstance(a, (types.CharSeq, types.Bytes)) and
            isinstance(b, (types.CharSeq, types.Bytes))):
        def impl(a, b):
            return (a._to_str() + b._to_str())._to_bytes()
        return impl


@overload(operator.mul)
def charseq_repeat(a, b):
    if isinstance(a, types.UnicodeCharSeq):
        def wrap(a, b):
            return str(a) * b
        return wrap
    if isinstance(b, types.UnicodeCharSeq):
        def wrap(a, b):
            return a * str(b)
        return wrap
    if isinstance(a, (types.CharSeq, types.Bytes)):
        def wrap(a, b):
            return (a._to_str() * b)._to_bytes()
        return wrap
    if isinstance(b, (types.CharSeq, types.Bytes)):
        def wrap(a, b):
            return (a * b._to_str())._to_bytes()
        return wrap


@overload(operator.not_)
def charseq_not(a):
    if isinstance(a, (types.UnicodeCharSeq, types.CharSeq, types.Bytes)):
        def impl(a):
            return len(a) == 0
        return impl


@overload(operator.eq)
def charseq_eq(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        def eq_impl(a, b):
            n = len(a)
            if n != len(b):
                return False
            for i in range(n):
                if left_code(a, i) != right_code(b, i):
                    return False
            return True
        return eq_impl


@overload(operator.ne)
def charseq_ne(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        def ne_impl(a, b):
            return not (a == b)
        return ne_impl


@overload(operator.lt)
def charseq_lt(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        def lt_impl(a, b):
            na = len(a)
            nb = len(b)
            n = min(na, nb)
            for i in range(n):
                ca, cb = left_code(a, i), right_code(b, i)
                if ca != cb:
                    return ca < cb
            return na < nb
        return lt_impl


@overload(operator.gt)
def charseq_gt(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        def gt_impl(a, b):
            return b < a
        return gt_impl


@overload(operator.le)
def charseq_le(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        def le_impl(a, b):
            return not (a > b)
        return le_impl


@overload(operator.ge)
def charseq_ge(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        def ge_impl(a, b):
            return not (a < b)
        return ge_impl


@overload(operator.contains)
def charseq_contains(a, b):
    if not _same_kind(a, b):
        return
    left_code = _get_code_impl(a)
    right_code = _get_code_impl(b)
    if left_code is not None and right_code is not None:
        if _is_bytes(a):
            def contains_impl(a, b):
                # Ideally, `return bytes(b) in bytes(a)` would be used
                # here, but numba Bytes does not implement
                # contains. So, using `unicode_type` implementation
                # here:
                return b._to_str() in a._to_str()
        else:
            def contains_impl(a, b):
                return str(b) in str(a)
        return contains_impl


@overload_method(types.UnicodeCharSeq, 'isascii')
@overload_method(types.CharSeq, 'isascii')
@overload_method(types.Bytes, 'isascii')
def charseq_isascii(s):
    get_code = _get_code_impl(s)

    def impl(s):
        for i in range(len(s)):
            if get_code(s, i) > 127:
                return False
        return True
    return impl


@overload_method(types.UnicodeCharSeq, '_get_kind')
@overload_method(types.CharSeq, '_get_kind')
def charseq_get_kind(s):
    get_code = _get_code_impl(s)

    def impl(s):
        max_code = 0
        for i in range(len(s)):
            code = get_code(s, i)
            if code > max_code:
                max_code = code
        if max_code > 0xffff:
            return unicode.PY_UNICODE_4BYTE_KIND
        if max_code > 0xff:
            return unicode.PY_UNICODE_2BYTE_KIND
        return unicode.PY_UNICODE_1BYTE_KIND
    return impl


@overload_method(types.UnicodeType, '_to_bytes')
def unicode_to_bytes_mth(s):
    """Convert unicode_type object to Bytes object.

    Note: The usage of _to_bytes method can be eliminated once all
    Python bytes operations are implemented for numba Bytes objects.

    """
    def impl(s):
        return _unicode_to_bytes(s)
    return impl


@overload_method(types.CharSeq, '_to_str')
@overload_method(types.Bytes, '_to_str')
def charseq_to_str_mth(s):
    """Convert bytes array item or bytes instance to UTF-8 str.

    Note: The usage of _to_str method can be eliminated once all
    Python bytes operations are implemented for numba Bytes objects.
    """
    get_code = _get_code_impl(s)

    def tostr_impl(s):
        n = len(s)
        is_ascii = s.isascii()
        result = unicode._empty_string(
            unicode.PY_UNICODE_1BYTE_KIND, n, is_ascii)
        for i in range(n):
            code = get_code(s, i)
            unicode._set_code_point(result, i, code)
        return result
    return tostr_impl


@overload_method(types.UnicodeCharSeq, "__str__")
def charseq_str(s):
    get_code = _get_code_impl(s)

    def str_impl(s):
        n = len(s)
        kind = s._get_kind()
        is_ascii = kind == 1 and s.isascii()
        result = unicode._empty_string(kind, n, is_ascii)
        for i in range(n):
            code = get_code(s, i)
            unicode._set_code_point(result, i, code)
        return result

    return str_impl


@overload(bytes)
def charseq_bytes(s):
    if isinstance(s, types.CharSeq):
        return lambda s: s


@overload_method(types.UnicodeCharSeq, '__hash__')
def unicode_charseq_hash(s):
    def impl(s):
        return hash(str(s))
    return impl


@overload_method(types.CharSeq, '__hash__')
def charseq_hash(s):
    def impl(s):
        # Ideally, `return hash(bytes(s))` would be used here but
        # numba Bytes does not implement hash (yet). However, for a
        # UTF-8 string `s`, we have hash(bytes(s)) == hash(s), hence,
        # we can convert CharSeq object to unicode_type and reuse its
        # hash implementation:
        return hash(s._to_str())
    return impl


@overload_method(types.UnicodeCharSeq, 'isupper')
def unicode_charseq_isupper(s):
    def impl(s):
        # workaround unicode_type.isupper bug: it returns int value
        return not not str(s).isupper()
    return impl


@overload_method(types.CharSeq, 'isupper')
def charseq_isupper(s):
    def impl(s):
        # return bytes(s).isupper()  # TODO: implement isupper for Bytes
        return not not s._to_str().isupper()
    return impl


@overload_method(types.UnicodeCharSeq, 'upper')
def unicode_charseq_upper(s):
    def impl(s):
        return str(s).upper()
    return impl


@overload_method(types.CharSeq, 'upper')
def charseq_upper(s):
    def impl(s):
        # return bytes(s).upper()  # TODO: implement upper for Bytes
        return s._to_str().upper()._to_bytes()
    return impl


@overload_method(types.UnicodeCharSeq, 'find')
@overload_method(types.CharSeq, 'find')
@overload_method(types.Bytes, 'find')
def unicode_charseq_find(a, b):
    if isinstance(a, types.UnicodeCharSeq):
        if isinstance(b, types.UnicodeCharSeq):
            def impl(a, b):
                return str(a).find(str(b))
            return impl
        if isinstance(b, types.UnicodeType):
            def impl(a, b):
                return str(a).find(b)
            return impl
    if isinstance(a, types.CharSeq):
        if isinstance(b, (types.CharSeq, types.Bytes)):
            def impl(a, b):
                return a._to_str().find(b._to_str())
            return impl
    if isinstance(a, types.UnicodeType):
        if isinstance(b, types.UnicodeCharSeq):
            def impl(a, b):
                return a.find(str(b))
            return impl
    if isinstance(a, types.Bytes):
        if isinstance(b, types.CharSeq):
            def impl(a, b):
                return a._to_str().find(b._to_str())
            return impl


@overload_method(types.UnicodeCharSeq, 'rfind')
@overload_method(types.CharSeq, 'rfind')
@overload_method(types.Bytes, 'rfind')
def unicode_charseq_rfind(a, b):
    if isinstance(a, types.UnicodeCharSeq):
        if isinstance(b, types.UnicodeCharSeq):
            def impl(a, b):
                return str(a).rfind(str(b))
            return impl
        if isinstance(b, types.UnicodeType):
            def impl(a, b):
                return str(a).rfind(b)
            return impl
    if isinstance(a, types.CharSeq):
        if isinstance(b, (types.CharSeq, types.Bytes)):
            def impl(a, b):
                return a._to_str().rfind(b._to_str())
            return impl
    if isinstance(a, types.UnicodeType):
        if isinstance(b, types.UnicodeCharSeq):
            def impl(a, b):
                return a.rfind(str(b))
            return impl
    if isinstance(a, types.Bytes):
        if isinstance(b, types.CharSeq):
            def impl(a, b):
                return a._to_str().rfind(b._to_str())
            return impl


@overload_method(types.UnicodeCharSeq, 'startswith')
@overload_method(types.CharSeq, 'startswith')
@overload_method(types.Bytes, 'startswith')
def unicode_charseq_startswith(a, b):
    if isinstance(a, types.UnicodeCharSeq):
        if isinstance(b, types.UnicodeCharSeq):
            def impl(a, b):
                return str(a).startswith(str(b))
            return impl
        if isinstance(b, types.UnicodeType):
            def impl(a, b):
                return str(a).startswith(b)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if isinstance(b, (types.CharSeq, types.Bytes)):
            def impl(a, b):
                return a._to_str().startswith(b._to_str())
            return impl


@overload_method(types.UnicodeCharSeq, 'endswith')
@overload_method(types.CharSeq, 'endswith')
@overload_method(types.Bytes, 'endswith')
def unicode_charseq_endswith(a, b):
    if isinstance(a, types.UnicodeCharSeq):
        if isinstance(b, types.UnicodeCharSeq):
            def impl(a, b):
                return str(a).endswith(str(b))
            return impl
        if isinstance(b, types.UnicodeType):
            def impl(a, b):
                return str(a).endswith(b)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if isinstance(b, (types.CharSeq, types.Bytes)):
            def impl(a, b):
                return a._to_str().endswith(b._to_str())
            return impl


@register_jitable
def _map_bytes(seq):
    return [s._to_bytes() for s in seq]


@overload_method(types.UnicodeCharSeq, 'split')
@overload_method(types.CharSeq, 'split')
@overload_method(types.Bytes, 'split')
def unicode_charseq_split(a, sep=None, maxsplit=-1):
    if not (maxsplit == -1 or
            isinstance(maxsplit, (types.Omitted, types.Integer,
                                  types.IntegerLiteral))):
        return None
    if isinstance(a, types.UnicodeCharSeq):
        if isinstance(sep, types.UnicodeCharSeq):
            def impl(a, sep=None, maxsplit=-1):
                return str(a).split(sep=str(sep), maxsplit=maxsplit)
            return impl
        if isinstance(sep, types.UnicodeType):
            def impl(a, sep=None, maxsplit=-1):
                return str(a).split(sep=sep, maxsplit=maxsplit)
            return impl
        if is_nonelike(sep):
            if is_default(maxsplit, -1):
                def impl(a, sep=None, maxsplit=-1):
                    return str(a).split()
            else:
                def impl(a, sep=None, maxsplit=-1):
                    return str(a).split(maxsplit=maxsplit)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if isinstance(sep, (types.CharSeq, types.Bytes)):
            def impl(a, sep=None, maxsplit=-1):
                return _map_bytes(a._to_str().split(sep._to_str(),
                                                    maxsplit=maxsplit))
            return impl
        if is_nonelike(sep):
            if is_default(maxsplit, -1):
                def impl(a, sep=None, maxsplit=-1):
                    return _map_bytes(a._to_str().split())
            else:
                def impl(a, sep=None, maxsplit=-1):
                    return _map_bytes(a._to_str().split(maxsplit=maxsplit))
            return impl

# NOT IMPLEMENTED: rsplit


@overload_method(types.UnicodeCharSeq, 'ljust')
@overload_method(types.CharSeq, 'ljust')
@overload_method(types.Bytes, 'ljust')
def unicode_charseq_ljust(a, width, fillchar=' '):
    if isinstance(a, types.UnicodeCharSeq):
        if is_default(fillchar, ' '):
            def impl(a, width, fillchar=' '):
                return str(a).ljust(width)
            return impl
        elif isinstance(fillchar, types.UnicodeCharSeq):
            def impl(a, width, fillchar=' '):
                return str(a).ljust(width, str(fillchar))
            return impl
        elif isinstance(fillchar, types.UnicodeType):
            def impl(a, width, fillchar=' '):
                return str(a).ljust(width, fillchar)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if is_default(fillchar, ' ') or is_default(fillchar, b' '):
            def impl(a, width, fillchar=' '):
                return a._to_str().ljust(width)._to_bytes()
            return impl
        elif isinstance(fillchar, (types.CharSeq, types.Bytes)):
            def impl(a, width, fillchar=' '):
                return a._to_str().ljust(width, fillchar._to_str())._to_bytes()
            return impl


@overload_method(types.UnicodeCharSeq, 'rjust')
@overload_method(types.CharSeq, 'rjust')
@overload_method(types.Bytes, 'rjust')
def unicode_charseq_rjust(a, width, fillchar=' '):
    if isinstance(a, types.UnicodeCharSeq):
        if is_default(fillchar, ' '):
            def impl(a, width, fillchar=' '):
                return str(a).rjust(width)
            return impl
        elif isinstance(fillchar, types.UnicodeCharSeq):
            def impl(a, width, fillchar=' '):
                return str(a).rjust(width, str(fillchar))
            return impl
        elif isinstance(fillchar, types.UnicodeType):
            def impl(a, width, fillchar=' '):
                return str(a).rjust(width, fillchar)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if is_default(fillchar, ' ') or is_default(fillchar, b' '):
            def impl(a, width, fillchar=' '):
                return a._to_str().rjust(width)._to_bytes()
            return impl
        elif isinstance(fillchar, (types.CharSeq, types.Bytes)):
            def impl(a, width, fillchar=' '):
                return a._to_str().rjust(width, fillchar._to_str())._to_bytes()
            return impl


@overload_method(types.UnicodeCharSeq, 'center')
@overload_method(types.CharSeq, 'center')
@overload_method(types.Bytes, 'center')
def unicode_charseq_center(a, width, fillchar=' '):
    if isinstance(a, types.UnicodeCharSeq):
        if is_default(fillchar, ' '):
            def impl(a, width, fillchar=' '):
                return str(a).center(width)
            return impl
        elif isinstance(fillchar, types.UnicodeCharSeq):
            def impl(a, width, fillchar=' '):
                return str(a).center(width, str(fillchar))
            return impl
        elif isinstance(fillchar, types.UnicodeType):
            def impl(a, width, fillchar=' '):
                return str(a).center(width, fillchar)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if is_default(fillchar, ' ') or is_default(fillchar, b' '):
            def impl(a, width, fillchar=' '):
                return a._to_str().center(width)._to_bytes()
            return impl
        elif isinstance(fillchar, (types.CharSeq, types.Bytes)):
            def impl(a, width, fillchar=' '):
                return a._to_str().center(width, fillchar._to_str())._to_bytes()
            return impl


@overload_method(types.UnicodeCharSeq, 'zfill')
@overload_method(types.CharSeq, 'zfill')
@overload_method(types.Bytes, 'zfill')
def unicode_charseq_zfill(a, width):
    if isinstance(a, types.UnicodeCharSeq):
        def impl(a, width):
            return str(a).zfill(width)
        return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        def impl(a, width):
            return a._to_str().zfill(width)._to_bytes()
        return impl


@overload_method(types.UnicodeCharSeq, 'lstrip')
@overload_method(types.CharSeq, 'lstrip')
@overload_method(types.Bytes, 'lstrip')
def unicode_charseq_lstrip(a, chars=None):
    if isinstance(a, types.UnicodeCharSeq):
        if is_nonelike(chars):
            def impl(a, chars=None):
                return str(a).lstrip()
            return impl
        elif isinstance(chars, types.UnicodeCharSeq):
            def impl(a, chars=None):
                return str(a).lstrip(str(chars))
            return impl
        elif isinstance(chars, types.UnicodeType):
            def impl(a, chars=None):
                return str(a).lstrip(chars)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if is_nonelike(chars):
            def impl(a, chars=None):
                return a._to_str().lstrip()._to_bytes()
            return impl
        elif isinstance(chars, (types.CharSeq, types.Bytes)):
            def impl(a, chars=None):
                return a._to_str().lstrip(chars._to_str())._to_bytes()
            return impl


@overload_method(types.UnicodeCharSeq, 'rstrip')
@overload_method(types.CharSeq, 'rstrip')
@overload_method(types.Bytes, 'rstrip')
def unicode_charseq_rstrip(a, chars=None):
    if isinstance(a, types.UnicodeCharSeq):
        if is_nonelike(chars):
            def impl(a, chars=None):
                return str(a).rstrip()
            return impl
        elif isinstance(chars, types.UnicodeCharSeq):
            def impl(a, chars=None):
                return str(a).rstrip(str(chars))
            return impl
        elif isinstance(chars, types.UnicodeType):
            def impl(a, chars=None):
                return str(a).rstrip(chars)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if is_nonelike(chars):
            def impl(a, chars=None):
                return a._to_str().rstrip()._to_bytes()
            return impl
        elif isinstance(chars, (types.CharSeq, types.Bytes)):
            def impl(a, chars=None):
                return a._to_str().rstrip(chars._to_str())._to_bytes()
            return impl


@overload_method(types.UnicodeCharSeq, 'strip')
@overload_method(types.CharSeq, 'strip')
@overload_method(types.Bytes, 'strip')
def unicode_charseq_strip(a, chars=None):
    if isinstance(a, types.UnicodeCharSeq):
        if is_nonelike(chars):
            def impl(a, chars=None):
                return str(a).strip()
            return impl
        elif isinstance(chars, types.UnicodeCharSeq):
            def impl(a, chars=None):
                return str(a).strip(str(chars))
            return impl
        elif isinstance(chars, types.UnicodeType):
            def impl(a, chars=None):
                return str(a).strip(chars)
            return impl
    if isinstance(a, (types.CharSeq, types.Bytes)):
        if is_nonelike(chars):
            def impl(a, chars=None):
                return a._to_str().strip()._to_bytes()
            return impl
        elif isinstance(chars, (types.CharSeq, types.Bytes)):
            def impl(a, chars=None):
                return a._to_str().strip(chars._to_str())._to_bytes()
            return impl


@overload_method(types.UnicodeCharSeq, 'join')
@overload_method(types.CharSeq, 'join')
@overload_method(types.Bytes, 'join')
def unicode_charseq_join(a, parts):

    if isinstance(a, types.UnicodeCharSeq):
        # assuming parts contains UnicodeCharSeq or UnicodeType objects
        def impl(a, parts):
            _parts = [str(p) for p in parts]
            return str(a).join(_parts)
        return impl

    if isinstance(a, (types.CharSeq, types.Bytes)):
        # assuming parts contains CharSeq or Bytes objects
        def impl(a, parts):
            _parts = [p._to_str() for p in parts]
            return a._to_str().join(_parts)._to_bytes()
        return impl
