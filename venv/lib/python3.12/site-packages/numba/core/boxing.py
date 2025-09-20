"""
Boxing and unboxing of native Numba values to / from CPython objects.
"""

from llvmlite import ir

from numba.core import types, cgutils
from numba.core.pythonapi import box, unbox, reflect, NativeValue
from numba.core.errors import NumbaNotImplementedError, TypingError
from numba.core.typing.typeof import typeof, Purpose

from numba.cpython import setobj, listobj
from numba.np import numpy_support
from contextlib import contextmanager, ExitStack


#
# Scalar types
#

@box(types.Boolean)
def box_bool(typ, val, c):
    return c.pyapi.bool_from_bool(val)

@unbox(types.Boolean)
def unbox_boolean(typ, obj, c):
    istrue = c.pyapi.object_istrue(obj)
    zero = ir.Constant(istrue.type, 0)
    val = c.builder.icmp_signed('!=', istrue, zero)
    return NativeValue(val, is_error=c.pyapi.c_api_error())


@box(types.IntegerLiteral)
@box(types.BooleanLiteral)
def box_literal_integer(typ, val, c):
    val = c.context.cast(c.builder, val, typ, typ.literal_type)
    return c.box(typ.literal_type, val)


@box(types.Integer)
def box_integer(typ, val, c):
    if typ.signed:
        ival = c.builder.sext(val, c.pyapi.longlong)
        return c.pyapi.long_from_longlong(ival)
    else:
        ullval = c.builder.zext(val, c.pyapi.ulonglong)
        return c.pyapi.long_from_ulonglong(ullval)

@unbox(types.Integer)
def unbox_integer(typ, obj, c):
    ll_type = c.context.get_argument_type(typ)
    val = cgutils.alloca_once(c.builder, ll_type)
    longobj = c.pyapi.number_long(obj)
    with c.pyapi.if_object_ok(longobj):
        if typ.signed:
            llval = c.pyapi.long_as_longlong(longobj)
        else:
            llval = c.pyapi.long_as_ulonglong(longobj)
        c.pyapi.decref(longobj)
        c.builder.store(c.builder.trunc(llval, ll_type), val)
    return NativeValue(c.builder.load(val),
                       is_error=c.pyapi.c_api_error())


@box(types.Float)
def box_float(typ, val, c):
    if typ == types.float32:
        dbval = c.builder.fpext(val, c.pyapi.double)
    else:
        assert typ == types.float64
        dbval = val
    return c.pyapi.float_from_double(dbval)

@unbox(types.Float)
def unbox_float(typ, obj, c):
    fobj = c.pyapi.number_float(obj)
    dbval = c.pyapi.float_as_double(fobj)
    c.pyapi.decref(fobj)
    if typ == types.float32:
        val = c.builder.fptrunc(dbval,
                                c.context.get_argument_type(typ))
    else:
        assert typ == types.float64
        val = dbval
    return NativeValue(val, is_error=c.pyapi.c_api_error())


@box(types.Complex)
def box_complex(typ, val, c):
    cval = c.context.make_complex(c.builder, typ, value=val)

    if typ == types.complex64:
        freal = c.builder.fpext(cval.real, c.pyapi.double)
        fimag = c.builder.fpext(cval.imag, c.pyapi.double)
    else:
        assert typ == types.complex128
        freal, fimag = cval.real, cval.imag
    return c.pyapi.complex_from_doubles(freal, fimag)

@unbox(types.Complex)
def unbox_complex(typ, obj, c):
    # First unbox to complex128, since that's what CPython gives us
    c128 = c.context.make_complex(c.builder, types.complex128)
    ok = c.pyapi.complex_adaptor(obj, c128._getpointer())
    failed = cgutils.is_false(c.builder, ok)

    with cgutils.if_unlikely(c.builder, failed):
        c.pyapi.err_set_string("PyExc_TypeError",
                               "conversion to %s failed" % (typ,))

    if typ == types.complex64:
        # Downcast to complex64 if necessary
        cplx = c.context.make_complex(c.builder, typ)
        cplx.real = c.context.cast(c.builder, c128.real,
                                   types.float64, types.float32)
        cplx.imag = c.context.cast(c.builder, c128.imag,
                                   types.float64, types.float32)
    else:
        assert typ == types.complex128
        cplx = c128
    return NativeValue(cplx._getvalue(), is_error=failed)


@box(types.NoneType)
def box_none(typ, val, c):
    return c.pyapi.make_none()

@unbox(types.NoneType)
@unbox(types.EllipsisType)
def unbox_none(typ, val, c):
    return NativeValue(c.context.get_dummy_value())


@box(types.NPDatetime)
def box_npdatetime(typ, val, c):
    return c.pyapi.create_np_datetime(val, typ.unit_code)

@unbox(types.NPDatetime)
def unbox_npdatetime(typ, obj, c):
    val = c.pyapi.extract_np_datetime(obj)
    return NativeValue(val, is_error=c.pyapi.c_api_error())


@box(types.NPTimedelta)
def box_nptimedelta(typ, val, c):
    return c.pyapi.create_np_timedelta(val, typ.unit_code)

@unbox(types.NPTimedelta)
def unbox_nptimedelta(typ, obj, c):
    val = c.pyapi.extract_np_timedelta(obj)
    return NativeValue(val, is_error=c.pyapi.c_api_error())


@box(types.RawPointer)
def box_raw_pointer(typ, val, c):
    """
    Convert a raw pointer to a Python int.
    """
    ll_intp = c.context.get_value_type(types.uintp)
    addr = c.builder.ptrtoint(val, ll_intp)
    return c.box(types.uintp, addr)


@box(types.EnumMember)
def box_enum(typ, val, c):
    """
    Fetch an enum member given its native value.
    """
    valobj = c.box(typ.dtype, val)
    # Call the enum class with the value object
    cls_obj = c.pyapi.unserialize(c.pyapi.serialize_object(typ.instance_class))
    return c.pyapi.call_function_objargs(cls_obj, (valobj,))


@unbox(types.EnumMember)
def unbox_enum(typ, obj, c):
    """
    Convert an enum member's value to its native value.
    """
    valobj = c.pyapi.object_getattr_string(obj, "value")
    return c.unbox(typ.dtype, valobj)


@box(types.UndefVar)
def box_undefvar(typ, val, c):
    """This type cannot be boxed, there's no Python equivalent"""
    msg = ("UndefVar type cannot be boxed, there is no Python equivalent of "
           "this type.")
    raise TypingError(msg)

#
# Composite types
#

@box(types.Record)
def box_record(typ, val, c):
    # Note we will create a copy of the record
    # This is the only safe way.
    size = ir.Constant(ir.IntType(32), val.type.pointee.count)
    ptr = c.builder.bitcast(val, ir.PointerType(ir.IntType(8)))
    return c.pyapi.recreate_record(ptr, size, typ.dtype, c.env_manager)


@unbox(types.Record)
def unbox_record(typ, obj, c):
    buf = c.pyapi.alloca_buffer()
    ptr = c.pyapi.extract_record_data(obj, buf)
    is_error = cgutils.is_null(c.builder, ptr)

    ltyp = c.context.get_value_type(typ)
    val = c.builder.bitcast(ptr, ltyp)

    def cleanup():
        c.pyapi.release_buffer(buf)
    return NativeValue(val, cleanup=cleanup, is_error=is_error)


@box(types.UnicodeCharSeq)
def box_unicodecharseq(typ, val, c):
    # XXX could kind be determined from strptr?
    unicode_kind = {
        1: c.pyapi.py_unicode_1byte_kind,
        2: c.pyapi.py_unicode_2byte_kind,
        4: c.pyapi.py_unicode_4byte_kind}[numpy_support.sizeof_unicode_char]
    kind = c.context.get_constant(types.int32, unicode_kind)
    rawptr = cgutils.alloca_once_value(c.builder, value=val)
    strptr = c.builder.bitcast(rawptr, c.pyapi.cstring)

    fullsize = c.context.get_constant(types.intp, typ.count)
    zero = fullsize.type(0)
    one = fullsize.type(1)
    step = fullsize.type(numpy_support.sizeof_unicode_char)
    count = cgutils.alloca_once_value(c.builder, zero)
    with cgutils.loop_nest(c.builder, [fullsize], fullsize.type) as [idx]:
        # Get char at idx
        ch = c.builder.load(c.builder.gep(strptr, [c.builder.mul(idx, step)]))
        # If the char is a non-null-byte, store the next index as count
        with c.builder.if_then(cgutils.is_not_null(c.builder, ch)):
            c.builder.store(c.builder.add(idx, one), count)
    strlen = c.builder.load(count)
    return c.pyapi.string_from_kind_and_data(kind, strptr, strlen)


@unbox(types.UnicodeCharSeq)
def unbox_unicodecharseq(typ, obj, c):
    lty = c.context.get_value_type(typ)

    ok, buffer, size, kind, is_ascii, hashv = \
        c.pyapi.string_as_string_size_and_kind(obj)

    # If conversion is ok, copy the buffer to the output storage.
    with cgutils.if_likely(c.builder, ok):
        # Check if the returned string size fits in the charseq
        storage_size = ir.Constant(size.type, typ.count)
        size_fits = c.builder.icmp_unsigned("<=", size, storage_size)

        # Allow truncation of string
        size = c.builder.select(size_fits, size, storage_size)

        # Initialize output to zero bytes
        null_string = ir.Constant(lty, None)
        outspace  = cgutils.alloca_once_value(c.builder, null_string)

        # We don't need to set the NULL-terminator because the storage
        # is already zero-filled.
        cgutils.memcpy(c.builder,
                        c.builder.bitcast(outspace, buffer.type),
                        buffer, size)

    ret = c.builder.load(outspace)
    return NativeValue(ret, is_error=c.builder.not_(ok))


@box(types.Bytes)
def box_bytes(typ, val, c):
    obj = c.context.make_helper(c.builder, typ, val)
    ret = c.pyapi.bytes_from_string_and_size(obj.data, obj.nitems)
    c.context.nrt.decref(c.builder, typ, val)
    return ret


@box(types.CharSeq)
def box_charseq(typ, val, c):
    rawptr = cgutils.alloca_once_value(c.builder, value=val)
    strptr = c.builder.bitcast(rawptr, c.pyapi.cstring)
    fullsize = c.context.get_constant(types.intp, typ.count)
    zero = fullsize.type(0)
    one = fullsize.type(1)
    count = cgutils.alloca_once_value(c.builder, zero)

    # Find the length of the string, mimicking Numpy's behaviour:
    # search for the last non-null byte in the underlying storage
    # (e.g. b'A\0\0B\0\0\0' will return the logical string b'A\0\0B')
    with cgutils.loop_nest(c.builder, [fullsize], fullsize.type) as [idx]:
        # Get char at idx
        ch = c.builder.load(c.builder.gep(strptr, [idx]))
        # If the char is a non-null-byte, store the next index as count
        with c.builder.if_then(cgutils.is_not_null(c.builder, ch)):
            c.builder.store(c.builder.add(idx, one), count)

    strlen = c.builder.load(count)
    return c.pyapi.bytes_from_string_and_size(strptr, strlen)


@unbox(types.CharSeq)
def unbox_charseq(typ, obj, c):
    lty = c.context.get_value_type(typ)
    ok, buffer, size = c.pyapi.string_as_string_and_size(obj)

    # If conversion is ok, copy the buffer to the output storage.
    with cgutils.if_likely(c.builder, ok):
        # Check if the returned string size fits in the charseq
        storage_size = ir.Constant(size.type, typ.count)
        size_fits = c.builder.icmp_unsigned("<=", size, storage_size)

        # Allow truncation of string
        size = c.builder.select(size_fits, size, storage_size)

        # Initialize output to zero bytes
        null_string = ir.Constant(lty, None)
        outspace  = cgutils.alloca_once_value(c.builder, null_string)

        # We don't need to set the NULL-terminator because the storage
        # is already zero-filled.
        cgutils.memcpy(c.builder,
                       c.builder.bitcast(outspace, buffer.type),
                       buffer, size)

    ret = c.builder.load(outspace)
    return NativeValue(ret, is_error=c.builder.not_(ok))


@box(types.Optional)
def box_optional(typ, val, c):
    optval = c.context.make_helper(c.builder, typ, val)
    ret = cgutils.alloca_once_value(c.builder, c.pyapi.borrow_none())
    with c.builder.if_else(optval.valid) as (then, otherwise):
        with then:
            validres = c.box(typ.type, optval.data)
            c.builder.store(validres, ret)
        with otherwise:
            c.builder.store(c.pyapi.make_none(), ret)
    return c.builder.load(ret)


@unbox(types.Optional)
def unbox_optional(typ, obj, c):
    """
    Convert object *obj* to a native optional structure.
    """
    noneval = c.context.make_optional_none(c.builder, typ.type)
    is_not_none = c.builder.icmp_signed('!=', obj, c.pyapi.borrow_none())

    retptr = cgutils.alloca_once(c.builder, noneval.type)
    errptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)

    with c.builder.if_else(is_not_none) as (then, orelse):
        with then:
            native = c.unbox(typ.type, obj)
            just = c.context.make_optional_value(c.builder,
                                                 typ.type, native.value)
            c.builder.store(just, retptr)
            c.builder.store(native.is_error, errptr)

        with orelse:
            c.builder.store(noneval, retptr)

    if native.cleanup is not None:
        def cleanup():
            with c.builder.if_then(is_not_none):
                native.cleanup()
    else:
        cleanup = None

    ret = c.builder.load(retptr)
    return NativeValue(ret, is_error=c.builder.load(errptr),
                       cleanup=cleanup)


@unbox(types.SliceType)
def unbox_slice(typ, obj, c):
    """
    Convert object *obj* to a native slice structure.
    """
    from numba.cpython import slicing
    ok, start, stop, step = c.pyapi.slice_as_ints(obj)
    sli = c.context.make_helper(c.builder, typ)
    sli.start = start
    sli.stop = stop
    sli.step = step
    return NativeValue(sli._getvalue(), is_error=c.builder.not_(ok))

@box(types.SliceLiteral)
def box_slice_literal(typ, val, c):
    # Check for integer overflows at compile time.
    slice_lit = typ.literal_value
    for field_name in ("start", "stop", "step"):
        field_obj = getattr(slice_lit, field_name)
        if isinstance(field_obj, int):
            try:
                typeof(field_obj, Purpose)
            except ValueError as e:
                raise ValueError((
                    f"Unable to create literal slice. "
                    f"Error encountered with {field_name} "
                    f"attribute. {str(e)}")
                )

    py_ctor, py_args = typ.literal_value.__reduce__()
    serialized_ctor = c.pyapi.serialize_object(py_ctor)
    serialized_args = c.pyapi.serialize_object(py_args)
    ctor = c.pyapi.unserialize(serialized_ctor)
    args = c.pyapi.unserialize(serialized_args)
    obj = c.pyapi.call(ctor, args)
    c.pyapi.decref(ctor)
    c.pyapi.decref(args)
    return obj

@unbox(types.StringLiteral)
def unbox_string_literal(typ, obj, c):
    # A string literal is a dummy value
    return NativeValue(c.context.get_dummy_value())

#
# Collections
#

# NOTE: boxing functions are supposed to steal any NRT references in
# the given native value.

@box(types.Array)
def box_array(typ, val, c):
    nativearycls = c.context.make_array(typ)
    nativeary = nativearycls(c.context, c.builder, value=val)
    if c.context.enable_nrt:
        np_dtype = numpy_support.as_dtype(typ.dtype)
        dtypeptr = c.env_manager.read_const(c.env_manager.add_const(np_dtype))
        newary = c.pyapi.nrt_adapt_ndarray_to_python(typ, val, dtypeptr)
        # Steals NRT ref
        c.context.nrt.decref(c.builder, typ, val)
        return newary
    else:
        parent = nativeary.parent
        c.pyapi.incref(parent)
        return parent


@unbox(types.Buffer)
def unbox_buffer(typ, obj, c):
    """
    Convert a Py_buffer-providing object to a native array structure.
    """
    buf = c.pyapi.alloca_buffer()
    res = c.pyapi.get_buffer(obj, buf)
    is_error = cgutils.is_not_null(c.builder, res)

    nativearycls = c.context.make_array(typ)
    nativeary = nativearycls(c.context, c.builder)
    aryptr = nativeary._getpointer()

    with cgutils.if_likely(c.builder, c.builder.not_(is_error)):
        ptr = c.builder.bitcast(aryptr, c.pyapi.voidptr)
        if c.context.enable_nrt:
            c.pyapi.nrt_adapt_buffer_from_python(buf, ptr)
        else:
            c.pyapi.numba_buffer_adaptor(buf, ptr)

    def cleanup():
        c.pyapi.release_buffer(buf)

    return NativeValue(c.builder.load(aryptr), is_error=is_error,
                       cleanup=cleanup)

@unbox(types.Array)
def unbox_array(typ, obj, c):
    """
    Convert a Numpy array object to a native array structure.
    """
    # This is necessary because unbox_buffer() does not work on some
    # dtypes, e.g. datetime64 and timedelta64.
    # TODO check matching dtype.
    #      currently, mismatching dtype will still work and causes
    #      potential memory corruption
    nativearycls = c.context.make_array(typ)
    nativeary = nativearycls(c.context, c.builder)
    aryptr = nativeary._getpointer()

    ptr = c.builder.bitcast(aryptr, c.pyapi.voidptr)
    if c.context.enable_nrt:
        errcode = c.pyapi.nrt_adapt_ndarray_from_python(obj, ptr)
    else:
        errcode = c.pyapi.numba_array_adaptor(obj, ptr)

    # TODO: here we have minimal typechecking by the itemsize.
    #       need to do better
    try:
        expected_itemsize = numpy_support.as_dtype(typ.dtype).itemsize
    except NumbaNotImplementedError:
        # Don't check types that can't be `as_dtype()`-ed
        itemsize_mismatch = cgutils.false_bit
    else:
        expected_itemsize = nativeary.itemsize.type(expected_itemsize)
        itemsize_mismatch = c.builder.icmp_unsigned(
            '!=',
            nativeary.itemsize,
            expected_itemsize,
            )

    failed = c.builder.or_(
        cgutils.is_not_null(c.builder, errcode),
        itemsize_mismatch,
    )
    # Handle error
    with c.builder.if_then(failed, likely=False):
        c.pyapi.err_set_string("PyExc_TypeError",
                               "can't unbox array from PyObject into "
                               "native value.  The object maybe of a "
                               "different type")
    return NativeValue(c.builder.load(aryptr), is_error=failed)


@box(types.Tuple)
@box(types.UniTuple)
def box_tuple(typ, val, c):
    """
    Convert native array or structure *val* to a tuple object.
    """
    tuple_val = c.pyapi.tuple_new(typ.count)

    for i, dtype in enumerate(typ):
        item = c.builder.extract_value(val, i)
        obj = c.box(dtype, item)
        c.pyapi.tuple_setitem(tuple_val, i, obj)

    return tuple_val

@box(types.NamedTuple)
@box(types.NamedUniTuple)
def box_namedtuple(typ, val, c):
    """
    Convert native array or structure *val* to a namedtuple object.
    """
    cls_obj = c.pyapi.unserialize(c.pyapi.serialize_object(typ.instance_class))
    tuple_obj = box_tuple(typ, val, c)
    obj = c.pyapi.call(cls_obj, tuple_obj)
    c.pyapi.decref(cls_obj)
    c.pyapi.decref(tuple_obj)
    return obj


@unbox(types.BaseTuple)
def unbox_tuple(typ, obj, c):
    """
    Convert tuple *obj* to a native array (if homogeneous) or structure.
    """
    n = len(typ)
    values = []
    cleanups = []
    lty = c.context.get_value_type(typ)

    is_error_ptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)
    value_ptr = cgutils.alloca_once(c.builder, lty)

    # Issue #1638: need to check the tuple size
    actual_size = c.pyapi.tuple_size(obj)
    size_matches = c.builder.icmp_unsigned('==', actual_size,
                                            ir.Constant(actual_size.type, n))
    with c.builder.if_then(c.builder.not_(size_matches), likely=False):
        c.pyapi.err_format(
            "PyExc_ValueError",
            "size mismatch for tuple, expected %d element(s) but got %%zd" % (n,),
            actual_size)
        c.builder.store(cgutils.true_bit, is_error_ptr)

    # We unbox the items even if not `size_matches`, to avoid issues with
    # the generated IR (instruction doesn't dominate all uses)
    for i, eltype in enumerate(typ):
        elem = c.pyapi.tuple_getitem(obj, i)
        native = c.unbox(eltype, elem)
        values.append(native.value)
        with c.builder.if_then(native.is_error, likely=False):
            c.builder.store(cgutils.true_bit, is_error_ptr)
        if native.cleanup is not None:
            cleanups.append(native.cleanup)

    value = c.context.make_tuple(c.builder, typ, values)
    c.builder.store(value, value_ptr)

    if cleanups:
        with c.builder.if_then(size_matches, likely=True):
            def cleanup():
                for func in reversed(cleanups):
                    func()
    else:
        cleanup = None

    return NativeValue(c.builder.load(value_ptr), cleanup=cleanup,
                       is_error=c.builder.load(is_error_ptr))


@box(types.List)
def box_list(typ, val, c):
    """
    Convert native list *val* to a list object.
    """
    list = listobj.ListInstance(c.context, c.builder, typ, val)
    obj = list.parent
    res = cgutils.alloca_once_value(c.builder, obj)
    with c.builder.if_else(cgutils.is_not_null(c.builder, obj)) as (has_parent, otherwise):
        with has_parent:
            # List is actually reflected => return the original object
            # (note not all list instances whose *type* is reflected are
            #  actually reflected; see numba.tests.test_lists for an example)
            c.pyapi.incref(obj)

        with otherwise:
            # Build a new Python list
            nitems = list.size
            obj = c.pyapi.list_new(nitems)
            with c.builder.if_then(cgutils.is_not_null(c.builder, obj),
                                   likely=True):
                with cgutils.for_range(c.builder, nitems) as loop:
                    item = list.getitem(loop.index)
                    list.incref_value(item)
                    itemobj = c.box(typ.dtype, item)
                    c.pyapi.list_setitem(obj, loop.index, itemobj)

            c.builder.store(obj, res)

    # Steal NRT ref
    c.context.nrt.decref(c.builder, typ, val)
    return c.builder.load(res)


class _NumbaTypeHelper(object):
    """A helper for acquiring `numba.typeof` for type checking.

    Usage
    -----

        # `c` is the boxing context.
        with _NumbaTypeHelper(c) as nth:
            # This contextmanager maintains the lifetime of the `numba.typeof`
            # function.
            the_numba_type = nth.typeof(some_object)
            # Do work on the type object
            do_checks(the_numba_type)
            # Cleanup
            c.pyapi.decref(the_numba_type)
        # At this point *nth* should not be used.
    """
    def __init__(self, c):
        self.c = c

    def __enter__(self):
        c = self.c
        numba_name = c.context.insert_const_string(c.builder.module, 'numba')
        numba_mod = c.pyapi.import_module(numba_name)
        typeof_fn = c.pyapi.object_getattr_string(numba_mod, 'typeof')
        self.typeof_fn = typeof_fn
        c.pyapi.decref(numba_mod)
        return self

    def __exit__(self, *args, **kwargs):
        c = self.c
        c.pyapi.decref(self.typeof_fn)

    def typeof(self, obj):
        res = self.c.pyapi.call_function_objargs(self.typeof_fn, [obj])
        return res


def _python_list_to_native(typ, obj, c, size, listptr, errorptr):
    """
    Construct a new native list from a Python list.
    """
    def check_element_type(nth, itemobj, expected_typobj):
        typobj = nth.typeof(itemobj)
        # Check if *typobj* is NULL
        with c.builder.if_then(
                cgutils.is_null(c.builder, typobj),
                likely=False,
                ):
            c.builder.store(cgutils.true_bit, errorptr)
            loop.do_break()
        # Mandate that objects all have the same exact type
        type_mismatch = c.builder.icmp_signed('!=', typobj, expected_typobj)

        with c.builder.if_then(type_mismatch, likely=False):
            c.builder.store(cgutils.true_bit, errorptr)
            c.pyapi.err_format(
                "PyExc_TypeError",
                "can't unbox heterogeneous list: %S != %S",
                expected_typobj, typobj,
                )
            c.pyapi.decref(typobj)
            loop.do_break()
        c.pyapi.decref(typobj)

    # Allocate a new native list
    ok, list = listobj.ListInstance.allocate_ex(c.context, c.builder, typ, size)
    with c.builder.if_else(ok, likely=True) as (if_ok, if_not_ok):
        with if_ok:
            list.size = size
            zero = ir.Constant(size.type, 0)
            with c.builder.if_then(c.builder.icmp_signed('>', size, zero),
                                   likely=True):
                # Traverse Python list and unbox objects into native list
                with _NumbaTypeHelper(c) as nth:
                    # Note: *expected_typobj* can't be NULL
                    expected_typobj = nth.typeof(c.pyapi.list_getitem(obj, zero))
                    with cgutils.for_range(c.builder, size) as loop:
                        itemobj = c.pyapi.list_getitem(obj, loop.index)
                        check_element_type(nth, itemobj, expected_typobj)
                        # XXX we don't call native cleanup for each
                        # list element, since that would require keeping
                        # of which unboxings have been successful.
                        native = c.unbox(typ.dtype, itemobj)
                        with c.builder.if_then(native.is_error, likely=False):
                            c.builder.store(cgutils.true_bit, errorptr)
                            loop.do_break()
                        # The reference is borrowed so incref=False
                        list.setitem(loop.index, native.value, incref=False)
                    c.pyapi.decref(expected_typobj)
            if typ.reflected:
                list.parent = obj
            # Stuff meminfo pointer into the Python object for
            # later reuse.
            with c.builder.if_then(c.builder.not_(c.builder.load(errorptr)),
                                                  likely=False):
                c.pyapi.object_set_private_data(obj, list.meminfo)
            list.set_dirty(False)
            c.builder.store(list.value, listptr)

        with if_not_ok:
            c.builder.store(cgutils.true_bit, errorptr)

    # If an error occurred, drop the whole native list
    with c.builder.if_then(c.builder.load(errorptr)):
        c.context.nrt.decref(c.builder, typ, list.value)


@unbox(types.List)
def unbox_list(typ, obj, c):
    """
    Convert list *obj* to a native list.

    If list was previously unboxed, we reuse the existing native list
    to ensure consistency.
    """
    size = c.pyapi.list_size(obj)

    errorptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)
    listptr = cgutils.alloca_once(c.builder, c.context.get_value_type(typ))

    # See if the list was previously unboxed, if so, re-use the meminfo.
    ptr = c.pyapi.object_get_private_data(obj)

    with c.builder.if_else(cgutils.is_not_null(c.builder, ptr)) \
        as (has_meminfo, otherwise):

        with has_meminfo:
            # List was previously unboxed => reuse meminfo
            list = listobj.ListInstance.from_meminfo(c.context, c.builder, typ, ptr)
            list.size = size
            if typ.reflected:
                list.parent = obj
            c.builder.store(list.value, listptr)

        with otherwise:
            _python_list_to_native(typ, obj, c, size, listptr, errorptr)

    def cleanup():
        # Clean up the associated pointer, as the meminfo is now invalid.
        c.pyapi.object_reset_private_data(obj)

    return NativeValue(c.builder.load(listptr),
                       is_error=c.builder.load(errorptr),
                       cleanup=cleanup)


@reflect(types.List)
def reflect_list(typ, val, c):
    """
    Reflect the native list's contents into the Python object.
    """
    if not typ.reflected:
        return
    if typ.dtype.reflected:
        msg = "cannot reflect element of reflected container: {}\n".format(typ)
        raise TypeError(msg)

    list = listobj.ListInstance(c.context, c.builder, typ, val)
    with c.builder.if_then(list.dirty, likely=False):
        obj = list.parent
        size = c.pyapi.list_size(obj)
        new_size = list.size
        diff = c.builder.sub(new_size, size)
        diff_gt_0 = c.builder.icmp_signed('>=', diff,
                                          ir.Constant(diff.type, 0))
        with c.builder.if_else(diff_gt_0) as (if_grow, if_shrink):
            # XXX no error checking below
            with if_grow:
                # First overwrite existing items
                with cgutils.for_range(c.builder, size) as loop:
                    item = list.getitem(loop.index)
                    list.incref_value(item)
                    itemobj = c.box(typ.dtype, item)
                    c.pyapi.list_setitem(obj, loop.index, itemobj)
                # Then add missing items
                with cgutils.for_range(c.builder, diff) as loop:
                    idx = c.builder.add(size, loop.index)
                    item = list.getitem(idx)
                    list.incref_value(item)
                    itemobj = c.box(typ.dtype, item)
                    c.pyapi.list_append(obj, itemobj)
                    c.pyapi.decref(itemobj)

            with if_shrink:
                # First delete list tail
                c.pyapi.list_setslice(obj, new_size, size, None)
                # Then overwrite remaining items
                with cgutils.for_range(c.builder, new_size) as loop:
                    item = list.getitem(loop.index)
                    list.incref_value(item)
                    itemobj = c.box(typ.dtype, item)
                    c.pyapi.list_setitem(obj, loop.index, itemobj)

        # Mark the list clean, in case it is reflected twice
        list.set_dirty(False)


def _python_set_to_native(typ, obj, c, size, setptr, errorptr):
    """
    Construct a new native set from a Python set.
    """
    # Allocate a new native set
    ok, inst = setobj.SetInstance.allocate_ex(c.context, c.builder, typ, size)
    with c.builder.if_else(ok, likely=True) as (if_ok, if_not_ok):
        with if_ok:
            # Traverse Python set and unbox objects into native set
            typobjptr = cgutils.alloca_once_value(c.builder,
                                                  ir.Constant(c.pyapi.pyobj, None))

            with c.pyapi.set_iterate(obj) as loop:
                itemobj = loop.value
                # Mandate that objects all have the same exact type
                typobj = c.pyapi.get_type(itemobj)
                expected_typobj = c.builder.load(typobjptr)

                with c.builder.if_else(
                    cgutils.is_null(c.builder, expected_typobj),
                    likely=False) as (if_first, if_not_first):
                    with if_first:
                        # First iteration => store item type
                        c.builder.store(typobj, typobjptr)
                    with if_not_first:
                        # Otherwise, check item type
                        type_mismatch = c.builder.icmp_signed('!=', typobj,
                                                              expected_typobj)
                        with c.builder.if_then(type_mismatch, likely=False):
                            c.builder.store(cgutils.true_bit, errorptr)
                            c.pyapi.err_set_string("PyExc_TypeError",
                                                   "can't unbox heterogeneous set")
                            loop.do_break()

                # XXX we don't call native cleanup for each set element,
                # since that would require keeping track
                # of which unboxings have been successful.
                native = c.unbox(typ.dtype, itemobj)
                with c.builder.if_then(native.is_error, likely=False):
                    c.builder.store(cgutils.true_bit, errorptr)
                inst.add_pyapi(c.pyapi, native.value, do_resize=False)

            if typ.reflected:
                inst.parent = obj
            # Associate meminfo pointer with the Python object for later reuse.
            with c.builder.if_then(c.builder.not_(c.builder.load(errorptr)),
                                   likely=False):
                c.pyapi.object_set_private_data(obj, inst.meminfo)
            inst.set_dirty(False)
            c.builder.store(inst.value, setptr)

        with if_not_ok:
            c.builder.store(cgutils.true_bit, errorptr)

    # If an error occurred, drop the whole native set
    with c.builder.if_then(c.builder.load(errorptr)):
        c.context.nrt.decref(c.builder, typ, inst.value)


@unbox(types.Set)
def unbox_set(typ, obj, c):
    """
    Convert set *obj* to a native set.

    If set was previously unboxed, we reuse the existing native set
    to ensure consistency.
    """
    size = c.pyapi.set_size(obj)

    errorptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)
    setptr = cgutils.alloca_once(c.builder, c.context.get_value_type(typ))

    # See if the set was previously unboxed, if so, re-use the meminfo.
    ptr = c.pyapi.object_get_private_data(obj)

    with c.builder.if_else(cgutils.is_not_null(c.builder, ptr)) \
        as (has_meminfo, otherwise):

        with has_meminfo:
            # Set was previously unboxed => reuse meminfo
            inst = setobj.SetInstance.from_meminfo(c.context, c.builder, typ, ptr)
            if typ.reflected:
                inst.parent = obj
            c.builder.store(inst.value, setptr)

        with otherwise:
            _python_set_to_native(typ, obj, c, size, setptr, errorptr)

    def cleanup():
        # Clean up the associated pointer, as the meminfo is now invalid.
        c.pyapi.object_reset_private_data(obj)

    return NativeValue(c.builder.load(setptr),
                       is_error=c.builder.load(errorptr),
                       cleanup=cleanup)


def _native_set_to_python_list(typ, payload, c):
    """
    Create a Python list from a native set's items.
    """
    nitems = payload.used
    listobj = c.pyapi.list_new(nitems)
    ok = cgutils.is_not_null(c.builder, listobj)
    with c.builder.if_then(ok, likely=True):
        index = cgutils.alloca_once_value(c.builder,
                                          ir.Constant(nitems.type, 0))
        with payload._iterate() as loop:
            i = c.builder.load(index)
            item = loop.entry.key
            c.context.nrt.incref(c.builder, typ.dtype, item)
            itemobj = c.box(typ.dtype, item)
            c.pyapi.list_setitem(listobj, i, itemobj)
            i = c.builder.add(i, ir.Constant(i.type, 1))
            c.builder.store(i, index)

    return ok, listobj


@box(types.Set)
def box_set(typ, val, c):
    """
    Convert native set *val* to a set object.
    """
    inst = setobj.SetInstance(c.context, c.builder, typ, val)
    obj = inst.parent
    res = cgutils.alloca_once_value(c.builder, obj)

    with c.builder.if_else(cgutils.is_not_null(c.builder, obj)) as (has_parent, otherwise):
        with has_parent:
            # Set is actually reflected => return the original object
            # (note not all set instances whose *type* is reflected are
            #  actually reflected; see numba.tests.test_sets for an example)
            c.pyapi.incref(obj)

        with otherwise:
            # Build a new Python list and then create a set from that
            payload = inst.payload
            ok, listobj = _native_set_to_python_list(typ, payload, c)
            with c.builder.if_then(ok, likely=True):
                obj = c.pyapi.set_new(listobj)
                c.pyapi.decref(listobj)
                c.builder.store(obj, res)

    # Steal NRT ref
    c.context.nrt.decref(c.builder, typ, val)
    return c.builder.load(res)

@reflect(types.Set)
def reflect_set(typ, val, c):
    """
    Reflect the native set's contents into the Python object.
    """
    if not typ.reflected:
        return
    inst = setobj.SetInstance(c.context, c.builder, typ, val)
    payload = inst.payload

    with c.builder.if_then(payload.dirty, likely=False):
        obj = inst.parent
        # XXX errors are not dealt with below
        c.pyapi.set_clear(obj)

        # Build a new Python list and then update the set with that
        ok, listobj = _native_set_to_python_list(typ, payload, c)
        with c.builder.if_then(ok, likely=True):
            c.pyapi.set_update(obj, listobj)
            c.pyapi.decref(listobj)

        # Mark the set clean, in case it is reflected twice
        inst.set_dirty(False)


#
# Other types
#

@box(types.Generator)
def box_generator(typ, val, c):
    return c.pyapi.from_native_generator(val, typ, c.env_manager.env_ptr)

@unbox(types.Generator)
def unbox_generator(typ, obj, c):
    return c.pyapi.to_native_generator(obj, typ)


@box(types.DType)
def box_dtype(typ, val, c):
    np_dtype = numpy_support.as_dtype(typ.dtype)
    return c.pyapi.unserialize(c.pyapi.serialize_object(np_dtype))

@unbox(types.DType)
def unbox_dtype(typ, val, c):
    return NativeValue(c.context.get_dummy_value())


@box(types.NumberClass)
def box_number_class(typ, val, c):
    np_dtype = numpy_support.as_dtype(typ.dtype)
    return c.pyapi.unserialize(c.pyapi.serialize_object(np_dtype))

@unbox(types.NumberClass)
def unbox_number_class(typ, val, c):
    return NativeValue(c.context.get_dummy_value())


@box(types.PyObject)
@box(types.Object)
def box_pyobject(typ, val, c):
    return val

@unbox(types.PyObject)
@unbox(types.Object)
def unbox_pyobject(typ, obj, c):
    return NativeValue(obj)


@unbox(types.ExternalFunctionPointer)
def unbox_funcptr(typ, obj, c):
    if typ.get_pointer is None:
        raise NotImplementedError(typ)

    # Call get_pointer() on the object to get the raw pointer value
    ptrty = c.context.get_function_pointer_type(typ)
    ret = cgutils.alloca_once_value(c.builder,
                                    ir.Constant(ptrty, None),
                                    name='fnptr')
    ser = c.pyapi.serialize_object(typ.get_pointer)
    get_pointer = c.pyapi.unserialize(ser)
    with cgutils.if_likely(c.builder,
                           cgutils.is_not_null(c.builder, get_pointer)):
        intobj = c.pyapi.call_function_objargs(get_pointer, (obj,))
        c.pyapi.decref(get_pointer)
        with cgutils.if_likely(c.builder,
                               cgutils.is_not_null(c.builder, intobj)):
            ptr = c.pyapi.long_as_voidptr(intobj)
            c.pyapi.decref(intobj)
            c.builder.store(c.builder.bitcast(ptr, ptrty), ret)
    return NativeValue(c.builder.load(ret), is_error=c.pyapi.c_api_error())

@box(types.DeferredType)
def box_deferred(typ, val, c):
    out = c.pyapi.from_native_value(typ.get(),
                                    c.builder.extract_value(val, [0]),
                                    env_manager=c.env_manager)
    return out


@unbox(types.DeferredType)
def unbox_deferred(typ, obj, c):
    native_value = c.pyapi.to_native_value(typ.get(), obj)
    model = c.context.data_model_manager[typ]
    res = model.set(c.builder, model.make_uninitialized(), native_value.value)
    return NativeValue(res, is_error=native_value.is_error,
                       cleanup=native_value.cleanup)


@unbox(types.Dispatcher)
def unbox_dispatcher(typ, obj, c):
    # In native code, Dispatcher types can be casted to FunctionType.
    return NativeValue(obj)


@box(types.Dispatcher)
def box_pyobject(typ, val, c):
    c.pyapi.incref(val)
    return val


def unbox_unsupported(typ, obj, c):
    c.pyapi.err_set_string("PyExc_TypeError",
                           "can't unbox {!r} type".format(typ))
    res = c.context.get_constant_null(typ)
    return NativeValue(res, is_error=cgutils.true_bit)


def box_unsupported(typ, val, c):
    msg = "cannot convert native %s to Python object" % (typ,)
    c.pyapi.err_set_string("PyExc_TypeError", msg)
    res = c.pyapi.get_null_object()
    return res


@box(types.Literal)
def box_literal(typ, val, c):
    # Const type contains the python object of the constant value,
    # which we can directly return.
    retval = typ.literal_value
    # Serialize the value into the IR
    return c.pyapi.unserialize(c.pyapi.serialize_object(retval))


@box(types.MemInfoPointer)
def box_meminfo_pointer(typ, val, c):
    return c.pyapi.nrt_meminfo_as_pyobject(val)


@unbox(types.MemInfoPointer)
def unbox_meminfo_pointer(typ, obj, c):
    res = c.pyapi.nrt_meminfo_from_pyobject(obj)
    errored = cgutils.is_null(c.builder, res)
    return NativeValue(res, is_error=errored)

@unbox(types.TypeRef)
def unbox_typeref(typ, val, c):
    return NativeValue(c.context.get_dummy_value(), is_error=cgutils.false_bit)


@box(types.LiteralStrKeyDict)
def box_LiteralStrKeyDict(typ, val, c):
    return box_unsupported(typ, val, c)


# Original implementation at: https://github.com/numba/numba/issues/4499#issuecomment-1063138477
@unbox(types.NumPyRandomBitGeneratorType)
def unbox_numpy_random_bitgenerator(typ, obj, c):
    """
    The bit_generator instance has a `.ctypes` attr which is a namedtuple
    with the following members (types):
    * state_address (Python int)
    * state (ctypes.c_void_p)
    * next_uint64 (ctypes.CFunctionType instance)
    * next_uint32 (ctypes.CFunctionType instance)
    * next_double (ctypes.CFunctionType instance)
    * bit_generator (ctypes.c_void_p)
    """

    is_error_ptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)
    extra_refs = []

    def clear_extra_refs():
        for _ref in extra_refs:
            c.pyapi.decref(_ref)

    def handle_failure():
        c.builder.store(cgutils.true_bit, is_error_ptr)
        clear_extra_refs()

    with ExitStack() as stack:

        def object_getattr_safely(obj, attr):
            attr_obj = c.pyapi.object_getattr_string(obj, attr)
            extra_refs.append(attr_obj)
            return attr_obj

        struct_ptr = cgutils.create_struct_proxy(typ)(c.context, c.builder)
        struct_ptr.parent = obj

        # Get the .ctypes attr
        ctypes_binding = object_getattr_safely(obj, 'ctypes')
        with cgutils.early_exit_if_null(c.builder, stack, ctypes_binding):
            handle_failure()

        # Look up the "state_address" member and wire it into the struct
        interface_state_address = object_getattr_safely(
            ctypes_binding, 'state_address')
        with cgutils.early_exit_if_null(c.builder, stack, interface_state_address):
            handle_failure()

        setattr(struct_ptr, 'state_address',
                c.unbox(types.uintp, interface_state_address).value)

        # Look up the "state" member and wire it into the struct
        interface_state = object_getattr_safely(ctypes_binding, 'state')
        with cgutils.early_exit_if_null(c.builder, stack, interface_state):
            handle_failure()

        interface_state_value = object_getattr_safely(
            interface_state, 'value')
        with cgutils.early_exit_if_null(c.builder, stack, interface_state_value):
            handle_failure()
        setattr(
            struct_ptr,
            'state',
            c.unbox(
                types.uintp,
                interface_state_value).value)

        # Want to store callable function pointers to these CFunctionTypes, so
        # import ctypes and use it to cast the CFunctionTypes to c_void_p and
        # store the results.
        # First find ctypes.cast, and ctypes.c_void_p
        ctypes_name = c.context.insert_const_string(c.builder.module, 'ctypes')
        ctypes_module = c.pyapi.import_module(ctypes_name)
        extra_refs.append(ctypes_module)
        with cgutils.early_exit_if_null(c.builder, stack, ctypes_module):
            handle_failure()

        ct_cast = object_getattr_safely(ctypes_module, 'cast')
        with cgutils.early_exit_if_null(c.builder, stack, ct_cast):
            handle_failure()

        ct_voidptr_ty = object_getattr_safely(ctypes_module, 'c_void_p')
        with cgutils.early_exit_if_null(c.builder, stack, ct_voidptr_ty):
            handle_failure()

        # This wires in the fnptrs referred to by name
        def wire_in_fnptrs(name):
            # Find the CFunctionType function
            interface_next_fn = c.pyapi.object_getattr_string(
                ctypes_binding, name)

            extra_refs.append(interface_next_fn)
            with cgutils.early_exit_if_null(c.builder, stack, interface_next_fn):
                handle_failure()

            # Want to do ctypes.cast(CFunctionType, ctypes.c_void_p), create an
            # args tuple for that.
            args = c.pyapi.tuple_pack([interface_next_fn, ct_voidptr_ty])
            with cgutils.early_exit_if_null(c.builder, stack, args):
                handle_failure()
            extra_refs.append(args)

            # Call ctypes.cast()
            interface_next_fn_casted = c.pyapi.call(ct_cast, args)
            extra_refs.append(interface_next_fn_casted)

            # Fetch the .value attr on the resulting ctypes.c_void_p for storage
            # in the function pointer slot.
            interface_next_fn_casted_value = object_getattr_safely(
                interface_next_fn_casted, 'value')
            with cgutils.early_exit_if_null(c.builder, stack, interface_next_fn_casted_value):
                handle_failure()

            # Wire up
            setattr(struct_ptr, f'fnptr_{name}',
                    c.unbox(types.uintp, interface_next_fn_casted_value).value)


        wire_in_fnptrs('next_double')
        wire_in_fnptrs('next_uint64')
        wire_in_fnptrs('next_uint32')

        clear_extra_refs()

    return NativeValue(struct_ptr._getvalue(), is_error=c.builder.load(is_error_ptr))

_bit_gen_type = types.NumPyRandomBitGeneratorType('bit_generator')

@unbox(types.NumPyRandomGeneratorType)
def unbox_numpy_random_generator(typ, obj, c):
    """
    Here we're creating a NumPyRandomGeneratorType StructModel with following fields:
    * ('bit_generator', _bit_gen_type): The unboxed BitGenerator associated with
                                        this Generator object instance.
    * ('parent', types.pyobject): Pointer to the original Generator PyObject.
    * ('meminfo', types.MemInfoPointer(types.voidptr)): The information about the memory
        stored at the pointer (to the original Generator PyObject). This is useful for
        keeping track of reference counts within the Python runtime. Helps prevent cases
        where deletion happens in Python runtime without NRT being awareness of it.
    """
    is_error_ptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)

    with ExitStack() as stack:
        struct_ptr = cgutils.create_struct_proxy(typ)(c.context, c.builder)
        bit_gen_inst = c.pyapi.object_getattr_string(obj, 'bit_generator')
        with cgutils.early_exit_if_null(c.builder, stack, bit_gen_inst):
            c.builder.store(cgutils.true_bit, is_error_ptr)
        unboxed = c.unbox(_bit_gen_type, bit_gen_inst).value
        struct_ptr.bit_generator = unboxed
        struct_ptr.parent = obj
        NULL = cgutils.voidptr_t(None)
        struct_ptr.meminfo = c.pyapi.nrt_meminfo_new_from_pyobject(
            NULL,  # there's no data
            obj,   # the python object, the call to nrt_meminfo_new_from_pyobject
                # will py_incref
        )
        c.pyapi.decref(bit_gen_inst)

    return NativeValue(struct_ptr._getvalue(), is_error=c.builder.load(is_error_ptr))


@box(types.NumPyRandomGeneratorType)
def box_numpy_random_generator(typ, val, c):
    inst = c.context.make_helper(c.builder, typ, val)
    obj = inst.parent
    res = cgutils.alloca_once_value(c.builder, obj)
    c.pyapi.incref(obj)
    # Steal NRT ref
    c.context.nrt.decref(c.builder, typ, val)
    return c.builder.load(res)
