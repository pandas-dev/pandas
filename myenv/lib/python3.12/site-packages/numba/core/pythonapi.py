from collections import namedtuple
import contextlib
import pickle
import hashlib
import sys

from llvmlite import ir
from llvmlite.ir import Constant

import ctypes
from numba import _helperlib
from numba.core import (
    types, utils, config, lowering, cgutils, imputils, serialize,
)
from numba.core.utils import PYVERSION

PY_UNICODE_1BYTE_KIND = _helperlib.py_unicode_1byte_kind
PY_UNICODE_2BYTE_KIND = _helperlib.py_unicode_2byte_kind
PY_UNICODE_4BYTE_KIND = _helperlib.py_unicode_4byte_kind
if PYVERSION in ((3, 9), (3, 10), (3, 11)):
    PY_UNICODE_WCHAR_KIND = _helperlib.py_unicode_wchar_kind


class _Registry(object):

    def __init__(self):
        self.functions = {}

    def register(self, typeclass):
        assert issubclass(typeclass, types.Type)
        def decorator(func):
            if typeclass in self.functions:
                raise KeyError("duplicate registration for %s" % (typeclass,))
            self.functions[typeclass] = func
            return func
        return decorator

    def lookup(self, typeclass, default=None):
        assert issubclass(typeclass, types.Type)
        for cls in typeclass.__mro__:
            func = self.functions.get(cls)
            if func is not None:
                return func
        return default

# Registries of boxing / unboxing implementations
_boxers = _Registry()
_unboxers = _Registry()
_reflectors = _Registry()

box = _boxers.register
unbox = _unboxers.register
reflect = _reflectors.register

class _BoxContext(namedtuple("_BoxContext",
                  ("context", "builder", "pyapi", "env_manager"))):
    """
    The facilities required by boxing implementations.
    """
    __slots__ = ()

    def box(self, typ, val):
        return self.pyapi.from_native_value(typ, val, self.env_manager)


class _UnboxContext(namedtuple("_UnboxContext",
                    ("context", "builder", "pyapi"))):
    """
    The facilities required by unboxing implementations.
    """
    __slots__ = ()

    def unbox(self, typ, obj):
        return self.pyapi.to_native_value(typ, obj)


class _ReflectContext(namedtuple("_ReflectContext",
                      ("context", "builder", "pyapi", "env_manager",
                       "is_error"))):
    """
    The facilities required by reflection implementations.
    """
    __slots__ = ()

    # XXX the error bit is currently unused by consumers (e.g. PyCallWrapper)
    def set_error(self):
        self.builder.store(self.is_error, cgutils.true_bit)

    def box(self, typ, val):
        return self.pyapi.from_native_value(typ, val, self.env_manager)

    def reflect(self, typ, val):
        return self.pyapi.reflect_native_value(typ, val, self.env_manager)


class NativeValue(object):
    """
    Encapsulate the result of converting a Python object to a native value,
    recording whether the conversion was successful and how to cleanup.
    """

    def __init__(self, value, is_error=None, cleanup=None):
        self.value = value
        self.is_error = is_error if is_error is not None else cgutils.false_bit
        self.cleanup = cleanup


class EnvironmentManager(object):

    def __init__(self, pyapi, env, env_body, env_ptr):
        assert isinstance(env, lowering.Environment)
        self.pyapi = pyapi
        self.env = env
        self.env_body = env_body
        self.env_ptr = env_ptr

    def add_const(self, const):
        """
        Add a constant to the environment, return its index.
        """
        # All constants are frozen inside the environment
        if isinstance(const, str):
            const = sys.intern(const)
        for index, val in enumerate(self.env.consts):
            if val is const:
                break
        else:
            index = len(self.env.consts)
            self.env.consts.append(const)
        return index

    def read_const(self, index):
        """
        Look up constant number *index* inside the environment body.
        A borrowed reference is returned.

        The returned LLVM value may have NULL value at runtime which indicates
        an error at runtime.
        """
        assert index < len(self.env.consts)

        builder = self.pyapi.builder
        consts = self.env_body.consts
        ret = cgutils.alloca_once(builder, self.pyapi.pyobj, zfill=True)
        with builder.if_else(cgutils.is_not_null(builder, consts)) as \
                (br_not_null, br_null):
            with br_not_null:
                getitem = self.pyapi.list_getitem(consts, index)
                builder.store(getitem, ret)
            with br_null:
                # This can happen when the Environment is accidentally released
                # and has subsequently been garbage collected.
                self.pyapi.err_set_string(
                    "PyExc_RuntimeError",
                    "`env.consts` is NULL in `read_const`",
                )
        return builder.load(ret)


_IteratorLoop = namedtuple('_IteratorLoop', ('value', 'do_break'))


class PythonAPI(object):
    """
    Code generation facilities to call into the CPython C API (and related
    helpers).
    """

    def __init__(self, context, builder):
        """
        Note: Maybe called multiple times when lowering a function
        """
        self.context = context
        self.builder = builder

        self.module = builder.basic_block.function.module
        # A unique mapping of serialized objects in this module
        try:
            self.module.__serialized
        except AttributeError:
            self.module.__serialized = {}

        # Initialize types
        self.pyobj = self.context.get_argument_type(types.pyobject)
        self.pyobjptr = self.pyobj.as_pointer()
        self.voidptr = ir.PointerType(ir.IntType(8))
        self.long = ir.IntType(ctypes.sizeof(ctypes.c_long) * 8)
        self.ulong = self.long
        self.longlong = ir.IntType(ctypes.sizeof(ctypes.c_ulonglong) * 8)
        self.ulonglong = self.longlong
        self.double = ir.DoubleType()
        self.py_ssize_t = self.context.get_value_type(types.intp)
        self.cstring = ir.PointerType(ir.IntType(8))
        self.gil_state = ir.IntType(_helperlib.py_gil_state_size * 8)
        self.py_buffer_t = ir.ArrayType(ir.IntType(8), _helperlib.py_buffer_size)
        self.py_hash_t = self.py_ssize_t
        self.py_unicode_1byte_kind = _helperlib.py_unicode_1byte_kind
        self.py_unicode_2byte_kind = _helperlib.py_unicode_2byte_kind
        self.py_unicode_4byte_kind = _helperlib.py_unicode_4byte_kind

    def get_env_manager(self, env, env_body, env_ptr):
        return EnvironmentManager(self, env, env_body, env_ptr)

    def emit_environment_sentry(self, envptr, return_pyobject=False,
                                debug_msg=''):
        """Emits LLVM code to ensure the `envptr` is not NULL
        """
        is_null = cgutils.is_null(self.builder, envptr)
        with cgutils.if_unlikely(self.builder, is_null):
            if return_pyobject:
                fnty = self.builder.function.type.pointee
                assert fnty.return_type == self.pyobj
                self.err_set_string(
                    "PyExc_RuntimeError", f"missing Environment: {debug_msg}",
                )
                self.builder.ret(self.get_null_object())
            else:
                self.context.call_conv.return_user_exc(
                    self.builder, RuntimeError,
                    (f"missing Environment: {debug_msg}",),
                )

    # ------ Python API -----

    #
    # Basic object API
    #

    def incref(self, obj):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="Py_IncRef")
        self.builder.call(fn, [obj])

    def decref(self, obj):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="Py_DecRef")
        self.builder.call(fn, [obj])

    def get_type(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="numba_py_type")
        return self.builder.call(fn, [obj])

    #
    # Argument unpacking
    #

    def parse_tuple_and_keywords(self, args, kws, fmt, keywords, *objs):
        charptr = ir.PointerType(ir.IntType(8))
        charptrary = ir.PointerType(charptr)
        argtypes = [self.pyobj, self.pyobj, charptr, charptrary]
        fnty = ir.FunctionType(ir.IntType(32), argtypes, var_arg=True)
        fn = self._get_function(fnty, name="PyArg_ParseTupleAndKeywords")
        return self.builder.call(fn, [args, kws, fmt, keywords] + list(objs))

    def parse_tuple(self, args, fmt, *objs):
        charptr = ir.PointerType(ir.IntType(8))
        argtypes = [self.pyobj, charptr]
        fnty = ir.FunctionType(ir.IntType(32), argtypes, var_arg=True)
        fn = self._get_function(fnty, name="PyArg_ParseTuple")
        return self.builder.call(fn, [args, fmt] + list(objs))

    def unpack_tuple(self, args, name, n_min, n_max, *objs):
        charptr = ir.PointerType(ir.IntType(8))
        argtypes = [self.pyobj, charptr, self.py_ssize_t, self.py_ssize_t]
        fnty = ir.FunctionType(ir.IntType(32), argtypes, var_arg=True)
        fn = self._get_function(fnty, name="PyArg_UnpackTuple")
        n_min = Constant(self.py_ssize_t, int(n_min))
        n_max = Constant(self.py_ssize_t, int(n_max))
        if isinstance(name, str):
            name = self.context.insert_const_string(self.builder.module, name)
        return self.builder.call(fn, [args, name, n_min, n_max] + list(objs))

    #
    # Exception and errors
    #

    def err_occurred(self):
        fnty = ir.FunctionType(self.pyobj, ())
        fn = self._get_function(fnty, name="PyErr_Occurred")
        return self.builder.call(fn, ())

    def err_clear(self):
        fnty = ir.FunctionType(ir.VoidType(), ())
        fn = self._get_function(fnty, name="PyErr_Clear")
        return self.builder.call(fn, ())

    def err_set_string(self, exctype, msg):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj, self.cstring])
        fn = self._get_function(fnty, name="PyErr_SetString")
        if isinstance(exctype, str):
            exctype = self.get_c_object(exctype)
        if isinstance(msg, str):
            msg = self.context.insert_const_string(self.module, msg)
        return self.builder.call(fn, (exctype, msg))

    def err_format(self, exctype, msg, *format_args):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj, self.cstring], var_arg=True)
        fn = self._get_function(fnty, name="PyErr_Format")
        if isinstance(exctype, str):
            exctype = self.get_c_object(exctype)
        if isinstance(msg, str):
            msg = self.context.insert_const_string(self.module, msg)
        return self.builder.call(fn, (exctype, msg) + tuple(format_args))

    def raise_object(self, exc=None):
        """
        Raise an arbitrary exception (type or value or (type, args)
        or None - if reraising).  A reference to the argument is consumed.
        """
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="numba_do_raise")
        if exc is None:
            exc = self.make_none()
        return self.builder.call(fn, (exc,))

    def err_set_object(self, exctype, excval):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyErr_SetObject")
        if isinstance(exctype, str):
            exctype = self.get_c_object(exctype)
        return self.builder.call(fn, (exctype, excval))

    def err_set_none(self, exctype):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="PyErr_SetNone")
        if isinstance(exctype, str):
            exctype = self.get_c_object(exctype)
        return self.builder.call(fn, (exctype,))

    def err_write_unraisable(self, obj):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="PyErr_WriteUnraisable")
        return self.builder.call(fn, (obj,))

    def err_fetch(self, pty, pval, ptb):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobjptr] * 3)
        fn = self._get_function(fnty, name="PyErr_Fetch")
        return self.builder.call(fn, (pty, pval, ptb))

    def err_restore(self, ty, val, tb):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj] * 3)
        fn = self._get_function(fnty, name="PyErr_Restore")
        return self.builder.call(fn, (ty, val, tb))

    @contextlib.contextmanager
    def err_push(self, keep_new=False):
        """
        Temporarily push the current error indicator while the code
        block is executed.  If *keep_new* is True and the code block
        raises a new error, the new error is kept, otherwise the old
        error indicator is restored at the end of the block.
        """
        pty, pval, ptb = [cgutils.alloca_once(self.builder, self.pyobj)
                          for i in range(3)]
        self.err_fetch(pty, pval, ptb)
        yield
        ty = self.builder.load(pty)
        val = self.builder.load(pval)
        tb = self.builder.load(ptb)
        if keep_new:
            new_error = cgutils.is_not_null(self.builder, self.err_occurred())
            with self.builder.if_else(new_error, likely=False) as (if_error, if_ok):
                with if_error:
                    # Code block raised an error, keep it
                    self.decref(ty)
                    self.decref(val)
                    self.decref(tb)
                with if_ok:
                    # Restore previous error
                    self.err_restore(ty, val, tb)
        else:
            self.err_restore(ty, val, tb)

    def get_c_object(self, name):
        """
        Get a Python object through its C-accessible *name*
        (e.g. "PyExc_ValueError").  The underlying variable must be
        a `PyObject *`, and the value of that pointer is returned.
        """
        # A LLVM global variable is implicitly a pointer to the declared
        # type, so fix up by using pyobj.pointee.
        return self.context.get_c_value(self.builder, self.pyobj.pointee, name,
                                        dllimport=True)

    def raise_missing_global_error(self, name):
        msg = "global name '%s' is not defined" % name
        cstr = self.context.insert_const_string(self.module, msg)
        self.err_set_string("PyExc_NameError", cstr)

    def raise_missing_name_error(self, name):
        msg = "name '%s' is not defined" % name
        cstr = self.context.insert_const_string(self.module, msg)
        self.err_set_string("PyExc_NameError", cstr)

    def fatal_error(self, msg):
        fnty = ir.FunctionType(ir.VoidType(), [self.cstring])
        fn = self._get_function(fnty, name="Py_FatalError")
        fn.attributes.add("noreturn")
        cstr = self.context.insert_const_string(self.module, msg)
        self.builder.call(fn, (cstr,))

    #
    # Concrete dict API
    #

    def dict_getitem_string(self, dic, name):
        """Lookup name inside dict

        Returns a borrowed reference
        """
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.cstring])
        fn = self._get_function(fnty, name="PyDict_GetItemString")
        cstr = self.context.insert_const_string(self.module, name)
        return self.builder.call(fn, [dic, cstr])

    def dict_getitem(self, dic, name):
        """Lookup name inside dict

        Returns a borrowed reference
        """
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyDict_GetItem")
        return self.builder.call(fn, [dic, name])

    def dict_new(self, presize=0):
        if presize == 0:
            fnty = ir.FunctionType(self.pyobj, ())
            fn = self._get_function(fnty, name="PyDict_New")
            return self.builder.call(fn, ())
        else:
            fnty = ir.FunctionType(self.pyobj, [self.py_ssize_t])
            fn = self._get_function(fnty, name="_PyDict_NewPresized")
            return self.builder.call(fn,
                                     [Constant(self.py_ssize_t, int(presize))])

    def dict_setitem(self, dictobj, nameobj, valobj):
        fnty = ir.FunctionType(ir.IntType(32), (self.pyobj, self.pyobj,
                                                self.pyobj))
        fn = self._get_function(fnty, name="PyDict_SetItem")
        return self.builder.call(fn, (dictobj, nameobj, valobj))

    def dict_setitem_string(self, dictobj, name, valobj):
        fnty = ir.FunctionType(ir.IntType(32), (self.pyobj, self.cstring,
                                                self.pyobj))
        fn = self._get_function(fnty, name="PyDict_SetItemString")
        cstr = self.context.insert_const_string(self.module, name)
        return self.builder.call(fn, (dictobj, cstr, valobj))

    def dict_pack(self, keyvalues):
        """
        Args
        -----
        keyvalues: iterable of (str, llvm.Value of PyObject*)
        """
        dictobj = self.dict_new()
        with self.if_object_ok(dictobj):
            for k, v in keyvalues:
                self.dict_setitem_string(dictobj, k, v)
        return dictobj

    #
    # Concrete number APIs
    #

    def float_from_double(self, fval):
        fnty = ir.FunctionType(self.pyobj, [self.double])
        fn = self._get_function(fnty, name="PyFloat_FromDouble")
        return self.builder.call(fn, [fval])

    def number_as_ssize_t(self, numobj):
        fnty = ir.FunctionType(self.py_ssize_t, [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_AsSsize_t")
        # We don't want any clipping, so pass OverflowError as the 2nd arg
        exc_class = self.get_c_object("PyExc_OverflowError")
        return self.builder.call(fn, [numobj, exc_class])

    def number_long(self, numobj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_Long")
        return self.builder.call(fn, [numobj])

    def long_as_ulonglong(self, numobj):
        fnty = ir.FunctionType(self.ulonglong, [self.pyobj])
        fn = self._get_function(fnty, name="PyLong_AsUnsignedLongLong")
        return self.builder.call(fn, [numobj])

    def long_as_longlong(self, numobj):
        fnty = ir.FunctionType(self.ulonglong, [self.pyobj])
        fn = self._get_function(fnty, name="PyLong_AsLongLong")
        return self.builder.call(fn, [numobj])

    def long_as_voidptr(self, numobj):
        """
        Convert the given Python integer to a void*.  This is recommended
        over number_as_ssize_t as it isn't affected by signedness.
        """
        fnty = ir.FunctionType(self.voidptr, [self.pyobj])
        fn = self._get_function(fnty, name="PyLong_AsVoidPtr")
        return self.builder.call(fn, [numobj])

    def _long_from_native_int(self, ival, func_name, native_int_type,
                              signed):
        fnty = ir.FunctionType(self.pyobj, [native_int_type])
        fn = self._get_function(fnty, name=func_name)
        resptr = cgutils.alloca_once(self.builder, self.pyobj)
        fn = self._get_function(fnty, name=func_name)
        self.builder.store(self.builder.call(fn, [ival]), resptr)

        return self.builder.load(resptr)

    def long_from_long(self, ival):
        func_name = "PyLong_FromLong"
        fnty = ir.FunctionType(self.pyobj, [self.long])
        fn = self._get_function(fnty, name=func_name)
        return self.builder.call(fn, [ival])

    def long_from_ulong(self, ival):
        return self._long_from_native_int(ival, "PyLong_FromUnsignedLong",
                                          self.long, signed=False)

    def long_from_ssize_t(self, ival):
        return self._long_from_native_int(ival, "PyLong_FromSsize_t",
                                          self.py_ssize_t, signed=True)

    def long_from_longlong(self, ival):
        return self._long_from_native_int(ival, "PyLong_FromLongLong",
                                          self.longlong, signed=True)

    def long_from_ulonglong(self, ival):
        return self._long_from_native_int(ival, "PyLong_FromUnsignedLongLong",
                                          self.ulonglong, signed=False)

    def long_from_signed_int(self, ival):
        """
        Return a Python integer from any native integer value.
        """
        bits = ival.type.width
        if bits <= self.long.width:
            return self.long_from_long(self.builder.sext(ival, self.long))
        elif bits <= self.longlong.width:
            return self.long_from_longlong(self.builder.sext(ival, self.longlong))
        else:
            raise OverflowError("integer too big (%d bits)" % (bits))

    def long_from_unsigned_int(self, ival):
        """
        Same as long_from_signed_int, but for unsigned values.
        """
        bits = ival.type.width
        if bits <= self.ulong.width:
            return self.long_from_ulong(self.builder.zext(ival, self.ulong))
        elif bits <= self.ulonglong.width:
            return self.long_from_ulonglong(self.builder.zext(ival, self.ulonglong))
        else:
            raise OverflowError("integer too big (%d bits)" % (bits))

    def _get_number_operator(self, name):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_%s" % name)
        return fn

    def _call_number_operator(self, name, lhs, rhs, inplace=False):
        if inplace:
            name = "InPlace" + name
        fn = self._get_number_operator(name)
        return self.builder.call(fn, [lhs, rhs])

    def number_add(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Add", lhs, rhs, inplace=inplace)

    def number_subtract(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Subtract", lhs, rhs, inplace=inplace)

    def number_multiply(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Multiply", lhs, rhs, inplace=inplace)

    def number_truedivide(self, lhs, rhs, inplace=False):
        return self._call_number_operator("TrueDivide", lhs, rhs, inplace=inplace)

    def number_floordivide(self, lhs, rhs, inplace=False):
        return self._call_number_operator("FloorDivide", lhs, rhs, inplace=inplace)

    def number_remainder(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Remainder", lhs, rhs, inplace=inplace)

    def number_matrix_multiply(self, lhs, rhs, inplace=False):
        return self._call_number_operator("MatrixMultiply", lhs, rhs, inplace=inplace)

    def number_lshift(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Lshift", lhs, rhs, inplace=inplace)

    def number_rshift(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Rshift", lhs, rhs, inplace=inplace)

    def number_and(self, lhs, rhs, inplace=False):
        return self._call_number_operator("And", lhs, rhs, inplace=inplace)

    def number_or(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Or", lhs, rhs, inplace=inplace)

    def number_xor(self, lhs, rhs, inplace=False):
        return self._call_number_operator("Xor", lhs, rhs, inplace=inplace)

    def number_power(self, lhs, rhs, inplace=False):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj] * 3)
        fname = "PyNumber_InPlacePower" if inplace else "PyNumber_Power"
        fn = self._get_function(fnty, fname)
        return self.builder.call(fn, [lhs, rhs, self.borrow_none()])

    def number_negative(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_Negative")
        return self.builder.call(fn, (obj,))

    def number_positive(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_Positive")
        return self.builder.call(fn, (obj,))

    def number_float(self, val):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_Float")
        return self.builder.call(fn, [val])

    def number_invert(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyNumber_Invert")
        return self.builder.call(fn, (obj,))

    def float_as_double(self, fobj):
        fnty = ir.FunctionType(self.double, [self.pyobj])
        fn = self._get_function(fnty, name="PyFloat_AsDouble")
        return self.builder.call(fn, [fobj])

    def bool_from_bool(self, bval):
        """
        Get a Python bool from a LLVM boolean.
        """
        longval = self.builder.zext(bval, self.long)
        return self.bool_from_long(longval)

    def bool_from_long(self, ival):
        fnty = ir.FunctionType(self.pyobj, [self.long])
        fn = self._get_function(fnty, name="PyBool_FromLong")
        return self.builder.call(fn, [ival])

    def complex_from_doubles(self, realval, imagval):
        fnty = ir.FunctionType(self.pyobj, [ir.DoubleType(), ir.DoubleType()])
        fn = self._get_function(fnty, name="PyComplex_FromDoubles")
        return self.builder.call(fn, [realval, imagval])

    def complex_real_as_double(self, cobj):
        fnty = ir.FunctionType(ir.DoubleType(), [self.pyobj])
        fn = self._get_function(fnty, name="PyComplex_RealAsDouble")
        return self.builder.call(fn, [cobj])

    def complex_imag_as_double(self, cobj):
        fnty = ir.FunctionType(ir.DoubleType(), [self.pyobj])
        fn = self._get_function(fnty, name="PyComplex_ImagAsDouble")
        return self.builder.call(fn, [cobj])

    #
    # Concrete slice API
    #
    def slice_as_ints(self, obj):
        """
        Read the members of a slice of integers.

        Returns a (ok, start, stop, step) tuple where ok is a boolean and
        the following members are pointer-sized ints.
        """
        pstart = cgutils.alloca_once(self.builder, self.py_ssize_t)
        pstop = cgutils.alloca_once(self.builder, self.py_ssize_t)
        pstep = cgutils.alloca_once(self.builder, self.py_ssize_t)
        fnty = ir.FunctionType(ir.IntType(32),
                               [self.pyobj] + [self.py_ssize_t.as_pointer()] * 3)
        fn = self._get_function(fnty, name="numba_unpack_slice")
        res = self.builder.call(fn, (obj, pstart, pstop, pstep))
        start = self.builder.load(pstart)
        stop = self.builder.load(pstop)
        step = self.builder.load(pstep)
        return cgutils.is_null(self.builder, res), start, stop, step

    #
    # List and sequence APIs
    #

    def sequence_getslice(self, obj, start, stop):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.py_ssize_t,
                                            self.py_ssize_t])
        fn = self._get_function(fnty, name="PySequence_GetSlice")
        return self.builder.call(fn, (obj, start, stop))

    def sequence_tuple(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PySequence_Tuple")
        return self.builder.call(fn, [obj])

    def sequence_concat(self, obj1, obj2):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PySequence_Concat")
        return self.builder.call(fn, [obj1, obj2])

    def list_new(self, szval):
        fnty = ir.FunctionType(self.pyobj, [self.py_ssize_t])
        fn = self._get_function(fnty, name="PyList_New")
        return self.builder.call(fn, [szval])

    def list_size(self, lst):
        fnty = ir.FunctionType(self.py_ssize_t, [self.pyobj])
        fn = self._get_function(fnty, name="PyList_Size")
        return self.builder.call(fn, [lst])

    def list_append(self, lst, val):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyList_Append")
        return self.builder.call(fn, [lst, val])

    def list_setitem(self, lst, idx, val):
        """
        Warning: Steals reference to ``val``
        """
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.py_ssize_t,
                                                self.pyobj])
        fn = self._get_function(fnty, name="PyList_SetItem")
        return self.builder.call(fn, [lst, idx, val])

    def list_getitem(self, lst, idx):
        """
        Returns a borrowed reference.
        """
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.py_ssize_t])
        fn = self._get_function(fnty, name="PyList_GetItem")
        if isinstance(idx, int):
            idx = self.context.get_constant(types.intp, idx)
        return self.builder.call(fn, [lst, idx])

    def list_setslice(self, lst, start, stop, obj):
        if obj is None:
            obj = self.get_null_object()
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.py_ssize_t,
                                                self.py_ssize_t, self.pyobj])
        fn = self._get_function(fnty, name="PyList_SetSlice")
        return self.builder.call(fn, (lst, start, stop, obj))


    #
    # Concrete tuple API
    #

    def tuple_getitem(self, tup, idx):
        """
        Borrow reference
        """
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.py_ssize_t])
        fn = self._get_function(fnty, name="PyTuple_GetItem")
        idx = self.context.get_constant(types.intp, idx)
        return self.builder.call(fn, [tup, idx])

    def tuple_pack(self, items):
        fnty = ir.FunctionType(self.pyobj, [self.py_ssize_t], var_arg=True)
        fn = self._get_function(fnty, name="PyTuple_Pack")
        n = self.context.get_constant(types.intp, len(items))
        args = [n]
        args.extend(items)
        return self.builder.call(fn, args)

    def tuple_size(self, tup):
        fnty = ir.FunctionType(self.py_ssize_t, [self.pyobj])
        fn = self._get_function(fnty, name="PyTuple_Size")
        return self.builder.call(fn, [tup])

    def tuple_new(self, count):
        fnty = ir.FunctionType(self.pyobj, [ir.IntType(32)])
        fn = self._get_function(fnty, name='PyTuple_New')
        return self.builder.call(fn, [self.context.get_constant(types.int32,
                                                                count)])

    def tuple_setitem(self, tuple_val, index, item):
        """
        Steals a reference to `item`.
        """
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, ir.IntType(32), self.pyobj])
        setitem_fn = self._get_function(fnty, name='PyTuple_SetItem')
        index = self.context.get_constant(types.int32, index)
        self.builder.call(setitem_fn, [tuple_val, index, item])

    #
    # Concrete set API
    #

    def set_new(self, iterable=None):
        if iterable is None:
            iterable = self.get_null_object()
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PySet_New")
        return self.builder.call(fn, [iterable])

    def set_add(self, set, value):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PySet_Add")
        return self.builder.call(fn, [set, value])

    def set_clear(self, set):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj])
        fn = self._get_function(fnty, name="PySet_Clear")
        return self.builder.call(fn, [set])

    def set_size(self, set):
        fnty = ir.FunctionType(self.py_ssize_t, [self.pyobj])
        fn = self._get_function(fnty, name="PySet_Size")
        return self.builder.call(fn, [set])

    def set_update(self, set, iterable):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="_PySet_Update")
        return self.builder.call(fn, [set, iterable])

    def set_next_entry(self, set, posptr, keyptr, hashptr):
        fnty = ir.FunctionType(ir.IntType(32),
                               [self.pyobj, self.py_ssize_t.as_pointer(),
                                self.pyobj.as_pointer(), self.py_hash_t.as_pointer()])
        fn = self._get_function(fnty, name="_PySet_NextEntry")
        return self.builder.call(fn, (set, posptr, keyptr, hashptr))

    @contextlib.contextmanager
    def set_iterate(self, set):
        builder = self.builder

        hashptr = cgutils.alloca_once(builder, self.py_hash_t, name="hashptr")
        keyptr = cgutils.alloca_once(builder, self.pyobj, name="keyptr")
        posptr = cgutils.alloca_once_value(builder,
                                           Constant(self.py_ssize_t, 0),
                                           name="posptr")

        bb_body = builder.append_basic_block("bb_body")
        bb_end = builder.append_basic_block("bb_end")

        builder.branch(bb_body)
        def do_break():
            builder.branch(bb_end)

        with builder.goto_block(bb_body):
            r = self.set_next_entry(set, posptr, keyptr, hashptr)
            finished = cgutils.is_null(builder, r)
            with builder.if_then(finished, likely=False):
                builder.branch(bb_end)
            yield _IteratorLoop(builder.load(keyptr), do_break)
            builder.branch(bb_body)

        builder.position_at_end(bb_end)

    #
    # GIL APIs
    #

    def gil_ensure(self):
        """
        Ensure the GIL is acquired.
        The returned value must be consumed by gil_release().
        """
        gilptrty = ir.PointerType(self.gil_state)
        fnty = ir.FunctionType(ir.VoidType(), [gilptrty])
        fn = self._get_function(fnty, "numba_gil_ensure")
        gilptr = cgutils.alloca_once(self.builder, self.gil_state)
        self.builder.call(fn, [gilptr])
        return gilptr

    def gil_release(self, gil):
        """
        Release the acquired GIL by gil_ensure().
        Must be paired with a gil_ensure().
        """
        gilptrty = ir.PointerType(self.gil_state)
        fnty = ir.FunctionType(ir.VoidType(), [gilptrty])
        fn = self._get_function(fnty, "numba_gil_release")
        return self.builder.call(fn, [gil])

    def save_thread(self):
        """
        Release the GIL and return the former thread state
        (an opaque non-NULL pointer).
        """
        fnty = ir.FunctionType(self.voidptr, [])
        fn = self._get_function(fnty, name="PyEval_SaveThread")
        return self.builder.call(fn, [])

    def restore_thread(self, thread_state):
        """
        Restore the given thread state by reacquiring the GIL.
        """
        fnty = ir.FunctionType(ir.VoidType(), [self.voidptr])
        fn = self._get_function(fnty, name="PyEval_RestoreThread")
        self.builder.call(fn, [thread_state])

    #
    # Generic object private data (a way of associating an arbitrary void *
    # pointer to an arbitrary Python object).
    #

    def object_get_private_data(self, obj):
        fnty = ir.FunctionType(self.voidptr, [self.pyobj])
        fn = self._get_function(fnty, name="numba_get_pyobject_private_data")
        return self.builder.call(fn, (obj,))

    def object_set_private_data(self, obj, ptr):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj, self.voidptr])
        fn = self._get_function(fnty, name="numba_set_pyobject_private_data")
        return self.builder.call(fn, (obj, ptr))

    def object_reset_private_data(self, obj):
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="numba_reset_pyobject_private_data")
        return self.builder.call(fn, (obj,))


    #
    # Other APIs (organize them better!)
    #

    def import_module_noblock(self, modname):
        fnty = ir.FunctionType(self.pyobj, [self.cstring])
        fn = self._get_function(fnty, name="PyImport_ImportModuleNoBlock")
        return self.builder.call(fn, [modname])

    def call_function_objargs(self, callee, objargs):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj], var_arg=True)
        fn = self._get_function(fnty, name="PyObject_CallFunctionObjArgs")
        args = [callee] + list(objargs)
        args.append(self.context.get_constant_null(types.pyobject))
        return self.builder.call(fn, args)

    def call_method(self, callee, method, objargs=()):
        cname = self.context.insert_const_string(self.module, method)
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.cstring, self.cstring],
                             var_arg=True)
        fn = self._get_function(fnty, name="PyObject_CallMethod")
        fmt = 'O' * len(objargs)
        cfmt = self.context.insert_const_string(self.module, fmt)
        args = [callee, cname, cfmt]
        if objargs:
            args.extend(objargs)
        args.append(self.context.get_constant_null(types.pyobject))
        return self.builder.call(fn, args)

    def call(self, callee, args=None, kws=None):
        if args_was_none := args is None:
            args = self.tuple_new(0)
        if kws is None:
            kws = self.get_null_object()
        fnty = ir.FunctionType(self.pyobj, [self.pyobj] * 3)
        fn = self._get_function(fnty, name="PyObject_Call")
        result = self.builder.call(fn, (callee, args, kws))
        if args_was_none:
            self.decref(args)
        return result

    def object_type(self, obj):
        """Emit a call to ``PyObject_Type(obj)`` to get the type of ``obj``.
        """
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyObject_Type")
        return self.builder.call(fn, (obj,))

    def object_istrue(self, obj):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj])
        fn = self._get_function(fnty, name="PyObject_IsTrue")
        return self.builder.call(fn, [obj])

    def object_not(self, obj):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj])
        fn = self._get_function(fnty, name="PyObject_Not")
        return self.builder.call(fn, [obj])

    def object_richcompare(self, lhs, rhs, opstr):
        """
        Refer to Python source Include/object.h for macros definition
        of the opid.
        """
        ops = ['<', '<=', '==', '!=', '>', '>=']
        if opstr in ops:
            opid = ops.index(opstr)
            fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj, ir.IntType(32)])
            fn = self._get_function(fnty, name="PyObject_RichCompare")
            lopid = self.context.get_constant(types.int32, opid)
            return self.builder.call(fn, (lhs, rhs, lopid))
        elif opstr == 'is':
            bitflag = self.builder.icmp_unsigned('==', lhs, rhs)
            return self.bool_from_bool(bitflag)
        elif opstr == 'is not':
            bitflag = self.builder.icmp_unsigned('!=', lhs, rhs)
            return self.bool_from_bool(bitflag)
        elif opstr in ('in', 'not in'):
            fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj])
            fn = self._get_function(fnty, name="PySequence_Contains")
            status = self.builder.call(fn, (rhs, lhs))
            negone = self.context.get_constant(types.int32, -1)
            is_good = self.builder.icmp_unsigned('!=', status, negone)
            # Stack allocate output and initialize to Null
            outptr = cgutils.alloca_once_value(self.builder,
                                               Constant(self.pyobj, None))
            # If PySequence_Contains returns non-error value
            with cgutils.if_likely(self.builder, is_good):
                if opstr == 'not in':
                    status = self.builder.not_(status)
                # Store the status as a boolean object
                truncated = self.builder.trunc(status, ir.IntType(1))
                self.builder.store(self.bool_from_bool(truncated),
                                   outptr)

            return self.builder.load(outptr)
        else:
            raise NotImplementedError("Unknown operator {op!r}".format(
                op=opstr))

    def iter_next(self, iterobj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyIter_Next")
        return self.builder.call(fn, [iterobj])

    def object_getiter(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyObject_GetIter")
        return self.builder.call(fn, [obj])

    def object_getattr_string(self, obj, attr):
        cstr = self.context.insert_const_string(self.module, attr)
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.cstring])
        fn = self._get_function(fnty, name="PyObject_GetAttrString")
        return self.builder.call(fn, [obj, cstr])

    def object_getattr(self, obj, attr):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyObject_GetAttr")
        return self.builder.call(fn, [obj, attr])

    def object_setattr_string(self, obj, attr, val):
        cstr = self.context.insert_const_string(self.module, attr)
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.cstring, self.pyobj])
        fn = self._get_function(fnty, name="PyObject_SetAttrString")
        return self.builder.call(fn, [obj, cstr, val])

    def object_setattr(self, obj, attr, val):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyObject_SetAttr")
        return self.builder.call(fn, [obj, attr, val])

    def object_delattr_string(self, obj, attr):
        # PyObject_DelAttrString() is actually a C macro calling
        # PyObject_SetAttrString() with value == NULL.
        return self.object_setattr_string(obj, attr, self.get_null_object())

    def object_delattr(self, obj, attr):
        # PyObject_DelAttr() is actually a C macro calling
        # PyObject_SetAttr() with value == NULL.
        return self.object_setattr(obj, attr, self.get_null_object())

    def object_getitem(self, obj, key):
        """
        Return obj[key]
        """
        fnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyObject_GetItem")
        return self.builder.call(fn, (obj, key))

    def object_setitem(self, obj, key, val):
        """
        obj[key] = val
        """
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyObject_SetItem")
        return self.builder.call(fn, (obj, key, val))

    def object_delitem(self, obj, key):
        """
        del obj[key]
        """
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.pyobj])
        fn = self._get_function(fnty, name="PyObject_DelItem")
        return self.builder.call(fn, (obj, key))

    def string_as_string(self, strobj):
        fnty = ir.FunctionType(self.cstring, [self.pyobj])
        fname = "PyUnicode_AsUTF8"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [strobj])

    def string_as_string_and_size(self, strobj):
        """
        Returns a tuple of ``(ok, buffer, length)``.
        The ``ok`` is i1 value that is set if ok.
        The ``buffer`` is a i8* of the output buffer.
        The ``length`` is a i32/i64 (py_ssize_t) of the length of the buffer.
        """

        p_length = cgutils.alloca_once(self.builder, self.py_ssize_t)
        fnty = ir.FunctionType(self.cstring, [self.pyobj,
                                              self.py_ssize_t.as_pointer()])
        fname = "PyUnicode_AsUTF8AndSize"
        fn = self._get_function(fnty, name=fname)

        buffer = self.builder.call(fn, [strobj, p_length])
        ok = self.builder.icmp_unsigned('!=',
                                        Constant(buffer.type, None),
                                        buffer)
        return (ok, buffer, self.builder.load(p_length))

    def string_as_string_size_and_kind(self, strobj):
        """
        Returns a tuple of ``(ok, buffer, length, kind)``.
        The ``ok`` is i1 value that is set if ok.
        The ``buffer`` is a i8* of the output buffer.
        The ``length`` is a i32/i64 (py_ssize_t) of the length of the buffer.
        The ``kind`` is a i32 (int32) of the Unicode kind constant
        The ``hash`` is a long/uint64_t (py_hash_t) of the Unicode constant hash
        """
        p_length = cgutils.alloca_once(self.builder, self.py_ssize_t)
        p_kind = cgutils.alloca_once(self.builder, ir.IntType(32))
        p_ascii = cgutils.alloca_once(self.builder, ir.IntType(32))
        p_hash = cgutils.alloca_once(self.builder, self.py_hash_t)
        fnty = ir.FunctionType(self.cstring, [self.pyobj,
                                              self.py_ssize_t.as_pointer(),
                                              ir.IntType(32).as_pointer(),
                                              ir.IntType(32).as_pointer(),
                                              self.py_hash_t.as_pointer()])
        fname = "numba_extract_unicode"
        fn = self._get_function(fnty, name=fname)

        buffer = self.builder.call(
            fn, [strobj, p_length, p_kind, p_ascii, p_hash])
        ok = self.builder.icmp_unsigned('!=',
                                        Constant(buffer.type, None),
                                        buffer)
        return (ok, buffer, self.builder.load(p_length),
                self.builder.load(p_kind), self.builder.load(p_ascii),
                self.builder.load(p_hash))

    def string_from_string_and_size(self, string, size):
        fnty = ir.FunctionType(self.pyobj, [self.cstring, self.py_ssize_t])
        fname = "PyString_FromStringAndSize"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [string, size])

    def string_from_string(self, string):
        fnty = ir.FunctionType(self.pyobj, [self.cstring])
        fname = "PyUnicode_FromString"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [string])

    def string_from_kind_and_data(self, kind, string, size):
        fnty = ir.FunctionType(self.pyobj, [ir.IntType(32), self.cstring, self.py_ssize_t])
        fname = "PyUnicode_FromKindAndData"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [kind, string, size])

    def bytes_as_string(self, obj):
        fnty = ir.FunctionType(self.cstring, [self.pyobj])
        fname = "PyBytes_AsString"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [obj])

    def bytes_as_string_and_size(self, obj, p_buffer, p_length):
        fnty = ir.FunctionType(
            ir.IntType(32),
            [self.pyobj, self.cstring.as_pointer(), self.py_ssize_t.as_pointer()],
        )
        fname = "PyBytes_AsStringAndSize"
        fn = self._get_function(fnty, name=fname)
        result = self.builder.call(fn, [obj, p_buffer, p_length])
        ok = self.builder.icmp_signed("!=", Constant(result.type, -1), result)
        return ok

    def bytes_from_string_and_size(self, string, size):
        fnty = ir.FunctionType(self.pyobj, [self.cstring, self.py_ssize_t])
        fname = "PyBytes_FromStringAndSize"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [string, size])

    def object_hash(self, obj):
        fnty = ir.FunctionType(self.py_hash_t, [self.pyobj, ])
        fname = "PyObject_Hash"
        fn = self._get_function(fnty, name=fname)
        return self.builder.call(fn, [obj,])

    def object_str(self, obj):
        fnty = ir.FunctionType(self.pyobj, [self.pyobj])
        fn = self._get_function(fnty, name="PyObject_Str")
        return self.builder.call(fn, [obj])

    def make_none(self):
        obj = self.borrow_none()
        self.incref(obj)
        return obj

    def borrow_none(self):
        return self.get_c_object("_Py_NoneStruct")

    def sys_write_stdout(self, fmt, *args):
        fnty = ir.FunctionType(ir.VoidType(), [self.cstring], var_arg=True)
        fn = self._get_function(fnty, name="PySys_FormatStdout")
        return self.builder.call(fn, (fmt,) + args)

    def object_dump(self, obj):
        """
        Dump a Python object on C stderr.  For debugging purposes.
        """
        fnty = ir.FunctionType(ir.VoidType(), [self.pyobj])
        fn = self._get_function(fnty, name="_PyObject_Dump")
        return self.builder.call(fn, (obj,))

    #
    # NRT (Numba runtime) APIs
    #

    def nrt_adapt_ndarray_to_python(self, aryty, ary, dtypeptr):
        assert self.context.enable_nrt, "NRT required"

        intty = ir.IntType(32)
        # Embed the Python type of the array (maybe subclass) in the LLVM IR.
        serial_aryty_pytype = self.unserialize(self.serialize_object(aryty.box_type))

        fnty = ir.FunctionType(self.pyobj,
                               [self.voidptr, self.pyobj, intty, intty, self.pyobj])
        fn = self._get_function(fnty, name="NRT_adapt_ndarray_to_python_acqref")
        fn.args[0].add_attribute('nocapture')

        ndim = self.context.get_constant(types.int32, aryty.ndim)
        writable = self.context.get_constant(types.int32, int(aryty.mutable))

        aryptr = cgutils.alloca_once_value(self.builder, ary)
        return self.builder.call(fn, [self.builder.bitcast(aryptr,
                                                           self.voidptr),
                                      serial_aryty_pytype,
                                      ndim, writable, dtypeptr])

    def nrt_meminfo_new_from_pyobject(self, data, pyobj):
        """
        Allocate a new MemInfo with data payload borrowed from a python
        object.
        """
        mod = self.builder.module
        fnty = ir.FunctionType(
            cgutils.voidptr_t,
            [cgutils.voidptr_t, cgutils.voidptr_t],
            )
        fn = cgutils.get_or_insert_function(
            mod,
            fnty,
            "NRT_meminfo_new_from_pyobject",
            )
        fn.args[0].add_attribute('nocapture')
        fn.args[1].add_attribute('nocapture')
        fn.return_value.add_attribute("noalias")
        return self.builder.call(fn, [data, pyobj])

    def nrt_meminfo_as_pyobject(self, miptr):
        mod = self.builder.module
        fnty = ir.FunctionType(
            self.pyobj,
            [cgutils.voidptr_t]
        )
        fn = cgutils.get_or_insert_function(
            mod,
            fnty,
            'NRT_meminfo_as_pyobject',
        )
        fn.return_value.add_attribute("noalias")
        return self.builder.call(fn, [miptr])

    def nrt_meminfo_from_pyobject(self, miobj):
        mod = self.builder.module
        fnty = ir.FunctionType(
            cgutils.voidptr_t,
            [self.pyobj]
        )
        fn = cgutils.get_or_insert_function(
            mod,
            fnty,
            'NRT_meminfo_from_pyobject',
        )
        fn.return_value.add_attribute("noalias")
        return self.builder.call(fn, [miobj])

    def nrt_adapt_ndarray_from_python(self, ary, ptr):
        assert self.context.enable_nrt
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.voidptr])
        fn = self._get_function(fnty, name="NRT_adapt_ndarray_from_python")
        fn.args[0].add_attribute('nocapture')
        fn.args[1].add_attribute('nocapture')
        return self.builder.call(fn, (ary, ptr))

    def nrt_adapt_buffer_from_python(self, buf, ptr):
        assert self.context.enable_nrt
        fnty = ir.FunctionType(ir.VoidType(), [ir.PointerType(self.py_buffer_t),
                                               self.voidptr])
        fn = self._get_function(fnty, name="NRT_adapt_buffer_from_python")
        fn.args[0].add_attribute('nocapture')
        fn.args[1].add_attribute('nocapture')
        return self.builder.call(fn, (buf, ptr))

    # ------ utils -----

    def _get_function(self, fnty, name):
        return cgutils.get_or_insert_function(self.module, fnty, name)

    def alloca_obj(self):
        return self.builder.alloca(self.pyobj)

    def alloca_buffer(self):
        """
        Return a pointer to a stack-allocated, zero-initialized Py_buffer.
        """
        # Treat the buffer as an opaque array of bytes
        ptr = cgutils.alloca_once_value(self.builder,
                                        Constant(self.py_buffer_t, None))
        return ptr

    @contextlib.contextmanager
    def if_object_ok(self, obj):
        with cgutils.if_likely(self.builder,
                               cgutils.is_not_null(self.builder, obj)):
            yield

    def print_object(self, obj):
        strobj = self.object_str(obj)
        cstr = self.string_as_string(strobj)
        fmt = self.context.insert_const_string(self.module, "%s")
        self.sys_write_stdout(fmt, cstr)
        self.decref(strobj)

    def print_string(self, text):
        fmt = self.context.insert_const_string(self.module, text)
        self.sys_write_stdout(fmt)

    def get_null_object(self):
        return Constant(self.pyobj, None)

    def return_none(self):
        none = self.make_none()
        self.builder.ret(none)

    def list_pack(self, items):
        n = len(items)
        seq = self.list_new(self.context.get_constant(types.intp, n))
        with self.if_object_ok(seq):
            for i in range(n):
                idx = self.context.get_constant(types.intp, i)
                self.incref(items[i])
                self.list_setitem(seq, idx, items[i])
        return seq

    def unserialize(self, structptr):
        """
        Unserialize some data.  *structptr* should be a pointer to
        a {i8* data, i32 length, i8* hashbuf, i8* func_ptr, i32 alloc_flag}
        structure.
        """
        fnty = ir.FunctionType(self.pyobj,
                             (self.voidptr, ir.IntType(32), self.voidptr))
        fn = self._get_function(fnty, name="numba_unpickle")
        ptr = self.builder.extract_value(self.builder.load(structptr), 0)
        n = self.builder.extract_value(self.builder.load(structptr), 1)
        hashed = self.builder.extract_value(self.builder.load(structptr), 2)
        return self.builder.call(fn, (ptr, n, hashed))

    def build_dynamic_excinfo_struct(self, struct_gv, exc_args):
        """
        Serialize some data at runtime. Returns a pointer to a python tuple
        (bytes_data, hash) where the first element is the serialized data as
        bytes and the second its hash.
        """
        fnty = ir.FunctionType(self.pyobj, (self.pyobj, self.pyobj))
        fn = self._get_function(fnty, name="numba_runtime_build_excinfo_struct")
        return self.builder.call(fn, (struct_gv, exc_args))

    def serialize_uncached(self, obj):
        """
        Same as serialize_object(), but don't create a global variable,
        simply return a literal for structure:
        {i8* data, i32 length, i8* hashbuf, i8* func_ptr, i32 alloc_flag}
        """
        # First make the array constant
        data = serialize.dumps(obj)
        assert len(data) < 2**31
        name = ".const.pickledata.%s" % (id(obj) if config.DIFF_IR == 0 else "DIFF_IR")
        bdata = cgutils.make_bytearray(data)
        # Make SHA1 hash on the pickled content
        # NOTE: update buffer size in numba_unpickle() when changing the
        #       hash algorithm.
        hashed = cgutils.make_bytearray(hashlib.sha1(data).digest())
        arr = self.context.insert_unique_const(self.module, name, bdata)
        hasharr = self.context.insert_unique_const(
            self.module, f"{name}.sha1", hashed,
        )
        # Then populate the structure constant
        struct = Constant.literal_struct([
            arr.bitcast(self.voidptr),
            Constant(ir.IntType(32), arr.type.pointee.count),
            hasharr.bitcast(self.voidptr),
            cgutils.get_null_value(self.voidptr),
            Constant(ir.IntType(32), 0),
            ])
        return struct

    def serialize_object(self, obj):
        """
        Serialize the given object in the bitcode, and return it
        as a pointer to a
        {i8* data, i32 length, i8* hashbuf, i8* fn_ptr, i32 alloc_flag},
        structure constant (suitable for passing to unserialize()).
        """
        try:
            gv = self.module.__serialized[obj]
        except KeyError:
            struct = self.serialize_uncached(obj)
            name = ".const.picklebuf.%s" % (id(obj) if config.DIFF_IR == 0 else "DIFF_IR")
            gv = self.context.insert_unique_const(self.module, name, struct)
            # Make the id() (and hence the name) unique while populating the module.
            self.module.__serialized[obj] = gv
        return gv

    def c_api_error(self):
        return cgutils.is_not_null(self.builder, self.err_occurred())

    def to_native_value(self, typ, obj):
        """
        Unbox the Python object as the given Numba type.
        A NativeValue instance is returned.
        """
        from numba.core.boxing import unbox_unsupported

        impl = _unboxers.lookup(typ.__class__, unbox_unsupported)
        c = _UnboxContext(self.context, self.builder, self)
        return impl(typ, obj, c)

    def from_native_return(self, typ, val, env_manager):
        assert not isinstance(typ, types.Optional), "callconv should have " \
                                                    "prevented the return of " \
                                                    "optional value"
        out = self.from_native_value(typ, val, env_manager)
        return out

    def from_native_value(self, typ, val, env_manager=None):
        """
        Box the native value of the given Numba type.  A Python object
        pointer is returned (NULL if an error occurred).
        This method steals any native (NRT) reference embedded in *val*.
        """
        from numba.core.boxing import box_unsupported

        impl = _boxers.lookup(typ.__class__, box_unsupported)

        c = _BoxContext(self.context, self.builder, self, env_manager)
        return impl(typ, val, c)

    def reflect_native_value(self, typ, val, env_manager=None):
        """
        Reflect the native value onto its Python original, if any.
        An error bit (as an LLVM value) is returned.
        """
        impl = _reflectors.lookup(typ.__class__)
        if impl is None:
            # Reflection isn't needed for most types
            return cgutils.false_bit

        is_error = cgutils.alloca_once_value(self.builder, cgutils.false_bit)
        c = _ReflectContext(self.context, self.builder, self, env_manager,
                            is_error)
        impl(typ, val, c)
        return self.builder.load(c.is_error)

    def to_native_generator(self, obj, typ):
        """
        Extract the generator structure pointer from a generator *obj*
        (a _dynfunc.Generator instance).
        """
        gen_ptr_ty = ir.PointerType(self.context.get_data_type(typ))
        value = self.context.get_generator_state(self.builder, obj, gen_ptr_ty)
        return NativeValue(value)

    def from_native_generator(self, val, typ, env=None):
        """
        Make a Numba generator (a _dynfunc.Generator instance) from a
        generator structure pointer *val*.
        *env* is an optional _dynfunc.Environment instance to be wrapped
        in the generator.
        """
        llty = self.context.get_data_type(typ)
        assert not llty.is_pointer
        gen_struct_size = self.context.get_abi_sizeof(llty)

        gendesc = self.context.get_generator_desc(typ)

        # This is the PyCFunctionWithKeywords generated by PyCallWrapper
        genfnty = ir.FunctionType(self.pyobj, [self.pyobj, self.pyobj, self.pyobj])
        genfn = self._get_function(genfnty, name=gendesc.llvm_cpython_wrapper_name)

        # This is the raw finalizer generated by _lower_generator_finalize_func()
        finalizerty = ir.FunctionType(ir.VoidType(), [self.voidptr])
        if typ.has_finalizer:
            finalizer = self._get_function(finalizerty, name=gendesc.llvm_finalizer_name)
        else:
            finalizer = Constant(ir.PointerType(finalizerty), None)

        # PyObject *numba_make_generator(state_size, initial_state, nextfunc, finalizer, env)
        fnty = ir.FunctionType(self.pyobj, [self.py_ssize_t,
                                            self.voidptr,
                                            ir.PointerType(genfnty),
                                            ir.PointerType(finalizerty),
                                            self.voidptr])
        fn = self._get_function(fnty, name="numba_make_generator")

        state_size = Constant(self.py_ssize_t, gen_struct_size)
        initial_state = self.builder.bitcast(val, self.voidptr)
        if env is None:
            env = self.get_null_object()
        env = self.builder.bitcast(env, self.voidptr)

        return self.builder.call(fn,
                                 (state_size, initial_state, genfn, finalizer, env))

    def numba_array_adaptor(self, ary, ptr):
        assert not self.context.enable_nrt
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, self.voidptr])
        fn = self._get_function(fnty, name="numba_adapt_ndarray")
        fn.args[0].add_attribute('nocapture')
        fn.args[1].add_attribute('nocapture')
        return self.builder.call(fn, (ary, ptr))

    def numba_buffer_adaptor(self, buf, ptr):
        fnty = ir.FunctionType(ir.VoidType(),
                             [ir.PointerType(self.py_buffer_t), self.voidptr])
        fn = self._get_function(fnty, name="numba_adapt_buffer")
        fn.args[0].add_attribute('nocapture')
        fn.args[1].add_attribute('nocapture')
        return self.builder.call(fn, (buf, ptr))

    def complex_adaptor(self, cobj, cmplx):
        fnty = ir.FunctionType(ir.IntType(32), [self.pyobj, cmplx.type])
        fn = self._get_function(fnty, name="numba_complex_adaptor")
        return self.builder.call(fn, [cobj, cmplx])

    def extract_record_data(self, obj, pbuf):
        fnty = ir.FunctionType(self.voidptr,
                               [self.pyobj, ir.PointerType(self.py_buffer_t)])
        fn = self._get_function(fnty, name="numba_extract_record_data")
        return self.builder.call(fn, [obj, pbuf])

    def get_buffer(self, obj, pbuf):
        fnty = ir.FunctionType(ir.IntType(32),
                               [self.pyobj, ir.PointerType(self.py_buffer_t)])
        fn = self._get_function(fnty, name="numba_get_buffer")
        return self.builder.call(fn, [obj, pbuf])

    def release_buffer(self, pbuf):
        fnty = ir.FunctionType(ir.VoidType(), [ir.PointerType(self.py_buffer_t)])
        fn = self._get_function(fnty, name="numba_release_buffer")
        return self.builder.call(fn, [pbuf])

    def extract_np_datetime(self, obj):
        fnty = ir.FunctionType(ir.IntType(64), [self.pyobj])
        fn = self._get_function(fnty, name="numba_extract_np_datetime")
        return self.builder.call(fn, [obj])

    def extract_np_timedelta(self, obj):
        fnty = ir.FunctionType(ir.IntType(64), [self.pyobj])
        fn = self._get_function(fnty, name="numba_extract_np_timedelta")
        return self.builder.call(fn, [obj])

    def create_np_datetime(self, val, unit_code):
        unit_code = Constant(ir.IntType(32), int(unit_code))
        fnty = ir.FunctionType(self.pyobj, [ir.IntType(64), ir.IntType(32)])
        fn = self._get_function(fnty, name="numba_create_np_datetime")
        return self.builder.call(fn, [val, unit_code])

    def create_np_timedelta(self, val, unit_code):
        unit_code = Constant(ir.IntType(32), int(unit_code))
        fnty = ir.FunctionType(self.pyobj, [ir.IntType(64), ir.IntType(32)])
        fn = self._get_function(fnty, name="numba_create_np_timedelta")
        return self.builder.call(fn, [val, unit_code])

    def recreate_record(self, pdata, size, dtype, env_manager):
        fnty = ir.FunctionType(self.pyobj, [ir.PointerType(ir.IntType(8)),
                                            ir.IntType(32), self.pyobj])
        fn = self._get_function(fnty, name="numba_recreate_record")
        dtypeaddr = env_manager.read_const(env_manager.add_const(dtype))
        return self.builder.call(fn, [pdata, size, dtypeaddr])

    def string_from_constant_string(self, string):
        cstr = self.context.insert_const_string(self.module, string)
        sz = self.context.get_constant(types.intp, len(string))
        return self.string_from_string_and_size(cstr, sz)

    def call_jit_code(self, func, sig, args):
        """Calls into Numba jitted code and propagate error using the Python
        calling convention.

        Parameters
        ----------
        func : function
            The Python function to be compiled. This function is compiled
            in nopython-mode.
        sig : numba.typing.Signature
            The function signature for *func*.
        args : Sequence[llvmlite.binding.Value]
            LLVM values to use as arguments.

        Returns
        -------
        (is_error, res) :  2-tuple of llvmlite.binding.Value.
            is_error : true iff *func* raised an exception.
            res : Returned value from *func* iff *is_error* is false.

        If *is_error* is true, this method will adapt the nopython exception
        into a Python exception. Caller should return NULL to Python to
        indicate an error.
        """
        # Compile *func*
        builder = self.builder
        cres = self.context.compile_subroutine(builder, func, sig)
        got_retty = cres.signature.return_type
        retty = sig.return_type
        if got_retty != retty:
            # This error indicates an error in *func* or the caller of this
            # method.
            raise errors.LoweringError(
                f'mismatching signature {got_retty} != {retty}.\n'
            )
        # Call into *func*
        status, res = self.context.call_internal_no_propagate(
            builder, cres.fndesc, sig, args,
        )
        # Post-call handling for *func*
        is_error_ptr = cgutils.alloca_once(builder, cgutils.bool_t, zfill=True)
        res_type = self.context.get_value_type(sig.return_type)
        res_ptr = cgutils.alloca_once(builder, res_type, zfill=True)

        # Handle error and adapt the nopython exception into cpython exception
        with builder.if_else(status.is_error) as (has_err, no_err):
            with has_err:
                builder.store(status.is_error, is_error_ptr)
                # Set error state in the Python interpreter
                self.context.call_conv.raise_error(builder, self, status)
            with no_err:
                # Handle returned value
                res = imputils.fix_returning_optional(
                    self.context, builder, sig, status, res,
                )
                builder.store(res, res_ptr)

        is_error = builder.load(is_error_ptr)
        res = builder.load(res_ptr)
        return is_error, res


class ObjModeUtils:
    """Internal utils for calling objmode dispatcher from within NPM code.
    """
    def __init__(self, pyapi):
        self.pyapi = pyapi

    def load_dispatcher(self, fnty, argtypes):
        builder = self.pyapi.builder
        tyctx = self.pyapi.context
        m = builder.module

        # Add a global variable to cache the objmode dispatcher
        gv = ir.GlobalVariable(
            m, self.pyapi.pyobj,
            name=m.get_unique_name("cached_objmode_dispatcher"),
        )
        gv.initializer = gv.type.pointee(None)
        gv.linkage = 'internal'

        # Make a basic-block to common exit
        bb_end = builder.append_basic_block("bb_end")

        if serialize.is_serialiable(fnty.dispatcher):
            serialized_dispatcher = self.pyapi.serialize_object(
                (fnty.dispatcher, tuple(argtypes)),
            )
            compile_args = self.pyapi.unserialize(serialized_dispatcher)
            # unserialize (unpickling) can fail
            failed_unser = cgutils.is_null(builder, compile_args)
            with builder.if_then(failed_unser):
                # early exit. `gv` is still null.
                builder.branch(bb_end)

        cached = builder.load(gv)
        with builder.if_then(cgutils.is_null(builder, cached)):
            if serialize.is_serialiable(fnty.dispatcher):
                cls = type(self)
                compiler = self.pyapi.unserialize(
                    self.pyapi.serialize_object(cls._call_objmode_dispatcher)
                )
                callee = self.pyapi.call_function_objargs(
                    compiler, [compile_args],
                )
                # Clean up
                self.pyapi.decref(compiler)
                self.pyapi.decref(compile_args)
            else:
                entry_pt = fnty.dispatcher.compile(tuple(argtypes))
                callee = tyctx.add_dynamic_addr(
                    builder, id(entry_pt), info="with_objectmode",
                )
            # Incref the dispatcher and cache it
            self.pyapi.incref(callee)
            builder.store(callee, gv)
        # Jump to the exit block
        builder.branch(bb_end)
        # Define the exit block
        builder.position_at_end(bb_end)
        callee = builder.load(gv)
        return callee

    @staticmethod
    def _call_objmode_dispatcher(compile_args):
        dispatcher, argtypes = compile_args
        entrypt = dispatcher.compile(argtypes)
        return entrypt
