"""
Calling conventions for Numba-compiled functions.
"""

from collections import namedtuple
from collections.abc import Iterable
import itertools
import hashlib

from llvmlite import ir

from numba.core import types, cgutils, errors
from numba.core.base import PYOBJECT, GENERIC_POINTER


TryStatus = namedtuple('TryStatus', ['in_try', 'excinfo'])


Status = namedtuple("Status",
                    ("code",
                     # If the function returned ok (a value or None)
                     "is_ok",
                     # If the function returned None
                     "is_none",
                     # If the function errored out (== not is_ok)
                     "is_error",
                     # If the generator exited with StopIteration
                     "is_stop_iteration",
                     # If the function errored with an already set exception
                     "is_python_exc",
                     # If the function errored with a user exception
                     "is_user_exc",
                     # The pointer to the exception info structure (for user
                     # exceptions)
                     "excinfoptr",
                     ))

int32_t = ir.IntType(32)
int64_t = ir.IntType(64)
errcode_t = int32_t


def _const_int(code):
    return ir.Constant(errcode_t, code)


RETCODE_OK = _const_int(0)
RETCODE_EXC = _const_int(-1)
RETCODE_NONE = _const_int(-2)
# StopIteration
RETCODE_STOPIT = _const_int(-3)

FIRST_USEREXC = 1

RETCODE_USEREXC = _const_int(FIRST_USEREXC)


class BaseCallConv(object):

    def __init__(self, context):
        self.context = context

    def return_optional_value(self, builder, retty, valty, value):
        if valty == types.none:
            # Value is none
            self.return_native_none(builder)

        elif retty == valty:
            # Value is an optional, need a runtime switch
            optval = self.context.make_helper(builder, retty, value=value)

            validbit = cgutils.as_bool_bit(builder, optval.valid)
            with builder.if_then(validbit):
                retval = self.context.get_return_value(builder, retty.type,
                                                       optval.data)
                self.return_value(builder, retval)

            self.return_native_none(builder)

        elif not isinstance(valty, types.Optional):
            # Value is not an optional, need a cast
            if valty != retty.type:
                value = self.context.cast(builder, value, fromty=valty,
                                          toty=retty.type)
            retval = self.context.get_return_value(builder, retty.type, value)
            self.return_value(builder, retval)

        else:
            raise NotImplementedError("returning {0} for {1}".format(valty,
                                                                     retty))

    def return_native_none(self, builder):
        self._return_errcode_raw(builder, RETCODE_NONE)

    def return_exc(self, builder):
        self._return_errcode_raw(builder, RETCODE_EXC)

    def return_stop_iteration(self, builder):
        self._return_errcode_raw(builder, RETCODE_STOPIT)

    def get_return_type(self, ty):
        """
        Get the actual type of the return argument for Numba type *ty*.
        """
        restype = self.context.data_model_manager[ty].get_return_type()
        return restype.as_pointer()

    def init_call_helper(self, builder):
        """
        Initialize and return a call helper object for the given builder.
        """
        ch = self._make_call_helper(builder)
        builder.__call_helper = ch
        return ch

    def _get_call_helper(self, builder):
        return builder.__call_helper

    def unpack_exception(self, builder, pyapi, status):
        return pyapi.unserialize(status.excinfoptr)

    def raise_error(self, builder, pyapi, status):
        """
        Given a non-ok *status*, raise the corresponding Python exception.
        """
        bbend = builder.function.append_basic_block()

        with builder.if_then(status.is_user_exc):
            # Unserialize user exception.
            # Make sure another error may not interfere.
            pyapi.err_clear()
            exc = self.unpack_exception(builder, pyapi, status)
            with cgutils.if_likely(builder,
                                   cgutils.is_not_null(builder, exc)):
                pyapi.raise_object(exc)  # steals ref
            builder.branch(bbend)

        with builder.if_then(status.is_stop_iteration):
            pyapi.err_set_none("PyExc_StopIteration")
            builder.branch(bbend)

        with builder.if_then(status.is_python_exc):
            # Error already raised => nothing to do
            builder.branch(bbend)

        pyapi.err_set_string("PyExc_SystemError",
                             "unknown error when calling native function")
        builder.branch(bbend)

        builder.position_at_end(bbend)

    def decode_arguments(self, builder, argtypes, func):
        """
        Get the decoded (unpacked) Python arguments with *argtypes*
        from LLVM function *func*.  A tuple of LLVM values is returned.
        """
        raw_args = self.get_arguments(func)
        arginfo = self._get_arg_packer(argtypes)
        return arginfo.from_arguments(builder, raw_args)

    def _get_arg_packer(self, argtypes):
        """
        Get an argument packer for the given argument types.
        """
        return self.context.get_arg_packer(argtypes)


class MinimalCallConv(BaseCallConv):
    """
    A minimal calling convention, suitable for e.g. GPU targets.
    The implemented function signature is:

        retcode_t (<Python return type>*, ... <Python arguments>)

    The return code will be one of the RETCODE_* constants or a
    function-specific user exception id (>= RETCODE_USEREXC).

    Caller is responsible for allocating a slot for the return value
    (passed as a pointer in the first argument).
    """

    def _make_call_helper(self, builder):
        return _MinimalCallHelper()

    def return_value(self, builder, retval):
        retptr = builder.function.args[0]
        assert retval.type == retptr.type.pointee, \
            (str(retval.type), str(retptr.type.pointee))
        builder.store(retval, retptr)
        self._return_errcode_raw(builder, RETCODE_OK)

    def return_user_exc(self, builder, exc, exc_args=None, loc=None,
                        func_name=None):
        if exc is not None and not issubclass(exc, BaseException):
            raise TypeError("exc should be None or exception class, got %r"
                            % (exc,))
        if exc_args is not None and not isinstance(exc_args, tuple):
            raise TypeError("exc_args should be None or tuple, got %r"
                            % (exc_args,))

        # Build excinfo struct
        if loc is not None:
            fname = loc._raw_function_name()
            if fname is None:
                # could be exec(<string>) or REPL, try func_name
                fname = func_name

            locinfo = (fname, loc.filename, loc.line)
            if None in locinfo:
                locinfo = None
        else:
            locinfo = None

        call_helper = self._get_call_helper(builder)
        exc_id = call_helper._add_exception(exc, exc_args, locinfo)
        self._return_errcode_raw(builder, _const_int(exc_id))

    def return_status_propagate(self, builder, status):
        self._return_errcode_raw(builder, status.code)

    def _return_errcode_raw(self, builder, code):
        if isinstance(code, int):
            code = _const_int(code)
        builder.ret(code)

    def _get_return_status(self, builder, code):
        """
        Given a return *code*, get a Status instance.
        """
        norm = builder.icmp_signed('==', code, RETCODE_OK)
        none = builder.icmp_signed('==', code, RETCODE_NONE)
        ok = builder.or_(norm, none)
        err = builder.not_(ok)
        exc = builder.icmp_signed('==', code, RETCODE_EXC)
        is_stop_iteration = builder.icmp_signed('==', code, RETCODE_STOPIT)
        is_user_exc = builder.icmp_signed('>=', code, RETCODE_USEREXC)

        status = Status(code=code,
                        is_ok=ok,
                        is_error=err,
                        is_python_exc=exc,
                        is_none=none,
                        is_user_exc=is_user_exc,
                        is_stop_iteration=is_stop_iteration,
                        excinfoptr=None)
        return status

    def get_function_type(self, restype, argtypes):
        """
        Get the implemented Function type for *restype* and *argtypes*.
        """
        arginfo = self._get_arg_packer(argtypes)
        argtypes = list(arginfo.argument_types)
        resptr = self.get_return_type(restype)
        fnty = ir.FunctionType(errcode_t, [resptr] + argtypes)
        return fnty

    def decorate_function(self, fn, args, fe_argtypes, noalias=False):
        """
        Set names and attributes of function arguments.
        """
        assert not noalias
        arginfo = self._get_arg_packer(fe_argtypes)
        arginfo.assign_names(self.get_arguments(fn),
                             ['arg.' + a for a in args])
        fn.args[0].name = ".ret"

    def get_arguments(self, func):
        """
        Get the Python-level arguments of LLVM *func*.
        """
        return func.args[1:]

    def call_function(self, builder, callee, resty, argtys, args):
        """
        Call the Numba-compiled *callee*.
        """
        retty = callee.args[0].type.pointee
        retvaltmp = cgutils.alloca_once(builder, retty)
        # initialize return value
        builder.store(cgutils.get_null_value(retty), retvaltmp)

        arginfo = self._get_arg_packer(argtys)
        args = arginfo.as_arguments(builder, args)
        realargs = [retvaltmp] + list(args)
        code = builder.call(callee, realargs)
        status = self._get_return_status(builder, code)
        retval = builder.load(retvaltmp)
        out = self.context.get_returned_value(builder, resty, retval)
        return status, out


class _MinimalCallHelper(object):
    """
    A call helper object for the "minimal" calling convention.
    User exceptions are represented as integer codes and stored in
    a mapping for retrieval from the caller.
    """

    def __init__(self):
        self.exceptions = {}

    def _add_exception(self, exc, exc_args, locinfo):
        """
        Add a new user exception to this helper. Returns an integer that can be
        used to refer to the added exception in future.

        Parameters
        ----------
        exc :
            exception type
        exc_args : None or tuple
            exception args
        locinfo : tuple
            location information
        """
        exc_id = len(self.exceptions) + FIRST_USEREXC
        self.exceptions[exc_id] = exc, exc_args, locinfo
        return exc_id

    def get_exception(self, exc_id):
        """
        Get information about a user exception. Returns a tuple of
        (exception type, exception args, location information).

        Parameters
        ----------
        id : integer
            The ID of the exception to look up
        """
        try:
            return self.exceptions[exc_id]
        except KeyError:
            msg = "unknown error %d in native function" % exc_id
            exc = SystemError
            exc_args = (msg,)
            locinfo = None
            return exc, exc_args, locinfo


# The structure type constructed by PythonAPI.serialize_uncached()
# i.e a {i8* pickle_buf, i32 pickle_bufsz, i8* hash_buf, i8* fn, i32 alloc_flag}
PICKLE_BUF_IDX = 0
PICKLE_BUFSZ_IDX = 1
HASH_BUF_IDX = 2
UNWRAP_FUNC_IDX = 3
ALLOC_FLAG_IDX = 4
excinfo_t = ir.LiteralStructType(
    [GENERIC_POINTER, int32_t, GENERIC_POINTER, GENERIC_POINTER, int32_t])
excinfo_ptr_t = ir.PointerType(excinfo_t)


class CPUCallConv(BaseCallConv):
    """
    The calling convention for CPU targets.
    The implemented function signature is:

        retcode_t (<Python return type>*, excinfo **, ... <Python arguments>)

    The return code will be one of the RETCODE_* constants.
    If RETCODE_USEREXC, the exception info pointer will be filled with
    a pointer to a constant struct describing the raised exception.

    Caller is responsible for allocating slots for the return value
    and the exception info pointer (passed as first and second arguments,
    respectively).
    """
    _status_ids = itertools.count(1)

    def _make_call_helper(self, builder):
        return None

    def return_value(self, builder, retval):
        retptr = self._get_return_argument(builder.function)
        assert retval.type == retptr.type.pointee, \
            (str(retval.type), str(retptr.type.pointee))
        builder.store(retval, retptr)
        self._return_errcode_raw(builder, RETCODE_OK)

    def build_excinfo_struct(self, exc, exc_args, loc, func_name):
        # Build excinfo struct
        if loc is not None:
            fname = loc._raw_function_name()
            if fname is None:
                # could be exec(<string>) or REPL, try func_name
                fname = func_name

            locinfo = (fname, loc.filename, loc.line)
            if None in locinfo:
                locinfo = None
        else:
            locinfo = None

        exc = (exc, exc_args, locinfo)
        return exc

    def set_static_user_exc(self, builder, exc, exc_args=None, loc=None,
                            func_name=None):
        if exc is not None and not issubclass(exc, BaseException):
            raise TypeError("exc should be None or exception class, got %r"
                            % (exc,))
        if exc_args is not None and not isinstance(exc_args, tuple):
            raise TypeError("exc_args should be None or tuple, got %r"
                            % (exc_args,))
        # None is indicative of no args, set the exc_args to an empty tuple
        # as PyObject_CallObject(exc, exc_args) requires the second argument to
        # be a tuple (or nullptr, but doing this makes it consistent)
        if exc_args is None:
            exc_args = tuple()

        # An exception in Numba is defined as the excinfo_t struct defined
        # above. Some arguments in this struct are not used, depending on
        # which kind of exception is being raised. A static exception uses
        # only the first three members whilst a dynamic exception uses all
        # members:
        #
        #             static exc - last 2 args are NULL and 0
        #             vvv  vvv  vvv
        # excinfo_t: {i8*, i32, i8*, i8*, i32}
        #                       ^^^  ^^^  ^^^
        #                       dynamic exc only - first 2 args are used for
        #                                          static info
        #
        # Comment below details how the struct is used in the case of a dynamic
        # exception. For dynamic exceptions, see
        # CPUCallConv::set_dynamic_user_exc
        #
        # {i8*, ___, ___, ___, ___}
        #   ^  serialized info about the exception (loc, kind, compile time
        #                                           args)
        #
        # {___, i32, ___, ___, ___}
        #        ^  len(serialized_exception)
        #
        # {___, ___, i8*, ___, ___}
        #             ^  Store a list of native values in a dynamic exception.
        #                Or a hash(serialized_exception) in a static exc.
        #
        # {___, ___, ___, i8*, ___}
        #                  ^  "NULL" as this member is not used in a static exc
        #
        # {___, ___, ___, ___, i32}
        #                       ^  Number of dynamic args in the exception. For
        #                          static exceptions, this value is "0"

        pyapi = self.context.get_python_api(builder)
        exc = self.build_excinfo_struct(exc, exc_args, loc, func_name)
        struct_gv = pyapi.serialize_object(exc)
        excptr = self._get_excinfo_argument(builder.function)
        store = builder.store(struct_gv, excptr)
        md = builder.module.add_metadata([ir.IntType(1)(1)])
        store.set_metadata("numba_exception_output", md)

    def return_user_exc(self, builder, exc, exc_args=None, loc=None,
                        func_name=None):
        try_info = getattr(builder, '_in_try_block', False)
        self.set_static_user_exc(builder, exc, exc_args=exc_args,
                                 loc=loc, func_name=func_name)
        self.check_try_status(builder)
        if try_info:
            # This is a hack for old-style impl.
            # We will branch directly to the exception handler.
            builder.branch(try_info['target'])
        else:
            # Return from the current function
            self._return_errcode_raw(builder, RETCODE_USEREXC)

    def unpack_dynamic_exception(self, builder, pyapi, status):
        excinfo_ptr = status.excinfoptr

        # load the serialized exception buffer from the module and create
        # a python bytes object
        picklebuf = builder.extract_value(
            builder.load(excinfo_ptr), PICKLE_BUF_IDX)
        picklebuf_sz = builder.extract_value(
            builder.load(excinfo_ptr), PICKLE_BUFSZ_IDX)
        static_exc_bytes = pyapi.bytes_from_string_and_size(
            picklebuf, builder.sext(picklebuf_sz, pyapi.py_ssize_t))

        # Load dynamic args (i8*) and the unwrap function
        dyn_args = builder.extract_value(
            builder.load(excinfo_ptr), HASH_BUF_IDX)
        func_ptr = builder.extract_value(
            builder.load(excinfo_ptr), UNWRAP_FUNC_IDX)

        # Convert the unwrap function to a function pointer and call it.
        # Function returns a python tuple with dynamic arguments converted to
        # CPython objects
        fnty = ir.FunctionType(PYOBJECT, [GENERIC_POINTER])
        fn = builder.bitcast(func_ptr, fnty.as_pointer())
        py_tuple = builder.call(fn, [dyn_args])

        # We check at this stage if creating the Python tuple was successful
        # or not. Note the exception is raised by calling PyErr_SetString
        # directly as the current function is the CPython wrapper.
        failed = cgutils.is_null(builder, py_tuple)
        with cgutils.if_unlikely(builder, failed):
            msg = ('Error creating Python tuple from runtime exception '
                   'arguments')
            pyapi.err_set_string("PyExc_RuntimeError", msg)
            # Return NULL to indicate an error was raised
            fnty = builder.function.function_type
            if not isinstance(fnty.return_type, ir.VoidType):
                # in some ufuncs, the return type is void
                builder.ret(cgutils.get_null_value(fnty.return_type))
            else:
                builder.ret_void()

        # merge static and dynamic variables
        excinfo = pyapi.build_dynamic_excinfo_struct(static_exc_bytes, py_tuple)

        # At this point, one can free the entire excinfo_ptr struct
        if self.context.enable_nrt:
            # One can safely emit a free instruction as it is only executed
            # if its in a dynamic exception branch
            self.context.nrt.free(
                builder, builder.bitcast(excinfo_ptr, pyapi.voidptr))
        return excinfo

    def unpack_exception(self, builder, pyapi, status):
        # Emit code that checks the alloc flag (last excinfo member)
        # if alloc_flag > 0:
        #     (dynamic) unpack the exception to retrieve runtime information
        #               and merge with static info
        # else:
        #     (static) unserialize the exception using pythonapi.unserialize

        excinfo_ptr = status.excinfoptr
        alloc_flag = builder.extract_value(builder.load(excinfo_ptr),
                                           ALLOC_FLAG_IDX)
        gt = builder.icmp_signed('>', alloc_flag, int32_t(0))
        with builder.if_else(gt) as (then, otherwise):
            with then:
                dyn_exc = self.unpack_dynamic_exception(builder, pyapi, status)
                bb_then = builder.block
            with otherwise:
                static_exc = pyapi.unserialize(excinfo_ptr)
                bb_else = builder.block
        phi = builder.phi(static_exc.type)
        phi.add_incoming(dyn_exc, bb_then)
        phi.add_incoming(static_exc, bb_else)
        return phi

    def emit_unwrap_dynamic_exception_fn(self, module, st_type, nb_types):
        # Create a function that converts a list of runtime arguments to a tuple
        # of PyObjects. i.e.:
        #
        #   @njit('void(float, int64)')
        #   def func(a, b):
        #       raise ValueError(a, 123, b)
        #
        # The last three arguments of the exception info structure will hold:
        #   {___, ___, i8*, i8*, i32}
        #               ^ A ptr to a {f32, i64} struct
        #                    ^ function ptr that converts i8* -> {f32, i64}* ->
        #                      python tuple
        #                          ^ Number of dynamic arguments = 2
        #

        _hash = hashlib.sha1(str(st_type).encode()).hexdigest()
        name = f'__excinfo_unwrap_args{_hash}'
        if name in module.globals:
            return module.globals.get(name)

        fnty = ir.FunctionType(GENERIC_POINTER, [GENERIC_POINTER])
        fn = ir.Function(module, fnty, name)
        # Linkage is changed to linkonce_odr for PYCC. External is the default.
        fn.linkage = "external"

        # prevent the function from being inlined
        fn.attributes.add('nounwind')
        fn.attributes.add('noinline')

        bb_entry = fn.append_basic_block('')
        builder = ir.IRBuilder(bb_entry)
        pyapi = self.context.get_python_api(builder)

        # i8* -> {native arg1 type, native arg2 type, ...}
        st_type_ptr = st_type.as_pointer()
        st_ptr = builder.bitcast(fn.args[0], st_type_ptr)
        # compile time values are stored as None
        nb_types = [typ for typ in nb_types if typ is not None]

        # convert native values into CPython objects
        objs = []
        env_manager = self.context.get_env_manager(builder,
                                                   return_pyobject=True)
        for i, typ in enumerate(nb_types):
            val = builder.extract_value(builder.load(st_ptr), i)
            obj = pyapi.from_native_value(typ, val, env_manager=env_manager)

            # If object cannot be boxed, raise an exception
            if obj == cgutils.get_null_value(obj.type):
                # When not supported, abort compilation
                msg = f'Cannot convert native {typ} to a Python object.'
                raise errors.TypingError(msg)

            objs.append(obj)

        # at this point, a pointer to the list of runtime values can be freed
        self.context.nrt.free(builder,
                              self._get_return_argument(builder.function))

        # Create a tuple of CPython objects
        tup = pyapi.tuple_pack(objs)
        builder.ret(tup)

        return fn

    def emit_wrap_args_insts(self, builder, pyapi, struct_type, exc_args):
        """
        Create an anonymous struct containing the given LLVM *values*.
        """
        st_size = pyapi.py_ssize_t(self.context.get_abi_sizeof(struct_type))

        st_ptr = builder.bitcast(
            self.context.nrt.allocate(builder, st_size),
            struct_type.as_pointer())

        # skip compile-time values
        exc_args = [arg for arg in exc_args if isinstance(arg, ir.Value)]

        zero = int32_t(0)
        for idx, arg in enumerate(exc_args):
            builder.store(arg, builder.gep(st_ptr, [zero, int32_t(idx)]))

        return st_ptr

    def set_dynamic_user_exc(self, builder, exc, exc_args, nb_types, loc=None,
                             func_name=None):
        """
        Compute the required bits to emit an exception with dynamic (runtime)
        values
        """
        if not issubclass(exc, BaseException):
            raise TypeError("exc should be an exception class, got %r"
                            % (exc,))
        if exc_args is not None and not isinstance(exc_args, tuple):
            raise TypeError("exc_args should be None or tuple, got %r"
                            % (exc_args,))

        # An exception in Numba is defined as the excinfo_t struct defined
        # above. Some arguments in this struct are not used, depending on
        # which kind of exception is being raised. A static exception uses
        # only the first three members whilst a dynamic exception uses all
        # members:
        #
        #             static exc - last 2 args are NULL and 0
        #             vvv  vvv  vvv
        # excinfo_t: {i8*, i32, i8*, i8*, i32}
        #                       ^^^  ^^^  ^^^
        #                       dynamic exc only - first 2 args are used for
        #                                          static info
        #
        # Comment below details how the struct is used in the case of a dynamic
        # exception. For static exception, see CPUCallConv::set_static_user_exc
        #
        # {i8*, ___, ___, ___, ___}
        #   ^  serialized info about the exception (loc, kind, compile time
        #                                           args)
        #
        # {___, i32, ___, ___, ___}
        #        ^  len(serialized_exception)
        #
        # {___, ___, i8*, ___, ___}
        #             ^  Store a list of native values in a dynamic exception.
        #                Or a hash(serialized_exception) in a static exc.
        #
        # {___, ___, ___, i8*, ___}
        #                  ^  Pointer to function that convert native values
        #                     into PyObject*. NULL in the case of a static
        #                     exception
        #
        # {___, ___, ___, ___, i32}
        #                       ^  Number of dynamic args in the exception.
        #                          Default is "0"
        #
        # The following code will:
        # 1) Serialize compile time information and store them in the first
        #    two args {i8*, i32, ___, ___, ___} of excinfo_t
        # 2) Emit the required code for converting native values to CPython
        #    objects. Those objects are stored in the last three args
        #    {___, ___, i8*, i8*, i32} of excinfo_t
        # 3) Allocate a new excinfo_t struct
        # 4) Fill excinfo_t struct and copy the pointer to the excinfo** arg

        # serialize comp. time args
        pyapi = self.context.get_python_api(builder)
        dummy = self.context.get_dummy_value()
        exc_args_static = tuple(
            [dummy if isinstance(arg, ir.Value) else arg for arg in exc_args])
        exc = self.build_excinfo_struct(exc, exc_args_static, loc, func_name)
        excinfo_pp = self._get_excinfo_argument(builder.function)
        struct_gv = builder.load(pyapi.serialize_object(exc))

        # Create the struct for runtime args and emit a function to convert it
        # into a Python tuple
        struct_type = ir.LiteralStructType([arg.type for arg in exc_args if
                                            isinstance(arg, ir.Value)])
        st_ptr = self.emit_wrap_args_insts(builder, pyapi, struct_type,
                                           exc_args)
        unwrap_fn = self.emit_unwrap_dynamic_exception_fn(
            builder.module, struct_type, nb_types)

        # allocate the excinfo struct
        exc_size = pyapi.py_ssize_t(self.context.get_abi_sizeof(excinfo_t))
        excinfo_p = builder.bitcast(
            self.context.nrt.allocate(builder, exc_size),
            excinfo_ptr_t)

        # fill the args
        zero = int32_t(0)
        exc_fields = (builder.extract_value(struct_gv, PICKLE_BUF_IDX),
                      builder.extract_value(struct_gv, PICKLE_BUFSZ_IDX),
                      builder.bitcast(st_ptr, GENERIC_POINTER),
                      builder.bitcast(unwrap_fn, GENERIC_POINTER),
                      int32_t(len(struct_type)))
        for idx, arg in enumerate(exc_fields):
            builder.store(arg, builder.gep(excinfo_p, [zero, int32_t(idx)]))
        builder.store(excinfo_p, excinfo_pp)

    def return_dynamic_user_exc(self, builder, exc, exc_args, nb_types,
                                loc=None, func_name=None):
        """
        Same as ::return_user_exc but for dynamic exceptions
        """
        self.set_dynamic_user_exc(builder, exc, exc_args, nb_types,
                                  loc=loc, func_name=func_name)
        self._return_errcode_raw(builder, RETCODE_USEREXC)

    def _get_try_state(self, builder):
        try:
            return builder.__eh_try_state
        except AttributeError:
            ptr = cgutils.alloca_once(
                builder, cgutils.intp_t, name='try_state', zfill=True,
            )
            builder.__eh_try_state = ptr
            return ptr

    def check_try_status(self, builder):
        try_state_ptr = self._get_try_state(builder)
        try_depth = builder.load(try_state_ptr)
        # try_depth > 0
        in_try = builder.icmp_unsigned('>', try_depth, try_depth.type(0))

        excinfoptr = self._get_excinfo_argument(builder.function)
        excinfo = builder.load(excinfoptr)

        return TryStatus(in_try=in_try, excinfo=excinfo)

    def set_try_status(self, builder):
        try_state_ptr = self._get_try_state(builder)
        # Increment try depth
        old = builder.load(try_state_ptr)
        new = builder.add(old, old.type(1))
        builder.store(new, try_state_ptr)

    def unset_try_status(self, builder):
        try_state_ptr = self._get_try_state(builder)
        # Decrement try depth
        old = builder.load(try_state_ptr)
        new = builder.sub(old, old.type(1))
        builder.store(new, try_state_ptr)

        # Needs to reset the exception state so that the exception handler
        # will run normally.
        excinfoptr = self._get_excinfo_argument(builder.function)
        null = cgutils.get_null_value(excinfoptr.type.pointee)
        builder.store(null, excinfoptr)

    def return_status_propagate(self, builder, status):
        trystatus = self.check_try_status(builder)
        excptr = self._get_excinfo_argument(builder.function)
        builder.store(status.excinfoptr, excptr)
        with builder.if_then(builder.not_(trystatus.in_try)):
            self._return_errcode_raw(builder, status.code)

    def _return_errcode_raw(self, builder, code):
        builder.ret(code)

    def _get_return_status(self, builder, code, excinfoptr):
        """
        Given a return *code* and *excinfoptr*, get a Status instance.
        """
        norm = builder.icmp_signed('==', code, RETCODE_OK)
        none = builder.icmp_signed('==', code, RETCODE_NONE)
        exc = builder.icmp_signed('==', code, RETCODE_EXC)
        is_stop_iteration = builder.icmp_signed('==', code, RETCODE_STOPIT)
        ok = builder.or_(norm, none)
        err = builder.not_(ok)
        is_user_exc = builder.icmp_signed('>=', code, RETCODE_USEREXC)
        excinfoptr = builder.select(is_user_exc, excinfoptr,
                                    ir.Constant(excinfo_ptr_t, ir.Undefined))

        status = Status(code=code,
                        is_ok=ok,
                        is_error=err,
                        is_python_exc=exc,
                        is_none=none,
                        is_user_exc=is_user_exc,
                        is_stop_iteration=is_stop_iteration,
                        excinfoptr=excinfoptr)
        return status

    def get_function_type(self, restype, argtypes):
        """
        Get the implemented Function type for *restype* and *argtypes*.
        """
        arginfo = self._get_arg_packer(argtypes)
        argtypes = list(arginfo.argument_types)
        resptr = self.get_return_type(restype)
        fnty = ir.FunctionType(errcode_t,
                               [resptr, ir.PointerType(excinfo_ptr_t)]
                               + argtypes)
        return fnty

    def decorate_function(self, fn, args, fe_argtypes, noalias=False):
        """
        Set names of function arguments, and add useful attributes to them.
        """
        arginfo = self._get_arg_packer(fe_argtypes)
        arginfo.assign_names(self.get_arguments(fn),
                             ['arg.' + a for a in args])
        retarg = self._get_return_argument(fn)
        retarg.name = "retptr"
        retarg.add_attribute("nocapture")
        retarg.add_attribute("noalias")
        excarg = self._get_excinfo_argument(fn)
        excarg.name = "excinfo"
        excarg.add_attribute("nocapture")
        excarg.add_attribute("noalias")

        if noalias:
            args = self.get_arguments(fn)
            for a in args:
                if isinstance(a.type, ir.PointerType):
                    a.add_attribute("nocapture")
                    a.add_attribute("noalias")

        # Add metadata to mark functions that may need NRT
        # thus disabling aggressive refct pruning in removerefctpass.py
        def type_may_always_need_nrt(ty):
            # Returns True if it's a non-Array type that is contains MemInfo
            if not isinstance(ty, types.Array):
                dmm = self.context.data_model_manager
                if dmm[ty].contains_nrt_meminfo():
                    return True
            return False

        args_may_always_need_nrt = any(
            map(type_may_always_need_nrt, fe_argtypes)
        )

        if args_may_always_need_nrt:
            nmd = fn.module.add_named_metadata(
                'numba_args_may_always_need_nrt',
            )
            nmd.add(fn.module.add_metadata([fn]))

    def get_arguments(self, func):
        """
        Get the Python-level arguments of LLVM *func*.
        """
        return func.args[2:]

    def _get_return_argument(self, func):
        return func.args[0]

    def _get_excinfo_argument(self, func):
        return func.args[1]

    def call_function(self, builder, callee, resty, argtys, args,
                      attrs=None):
        """
        Call the Numba-compiled *callee*.
        Parameters:
        -----------
        attrs: LLVM style string or iterable of individual attributes, default
               is None which specifies no attributes. Examples:
               LLVM style string: "noinline fast"
               Equivalent iterable: ("noinline", "fast")
        """
        retty = self.get_return_type(resty).pointee
        actual_retty = self._get_return_argument(callee.function_type).pointee
        if retty != actual_retty:
            m = f"Function type returns {actual_retty} but resty={retty}"
            raise ValueError(m)

        retvaltmp = cgutils.alloca_once(builder, retty)
        # initialize return value to zeros
        builder.store(cgutils.get_null_value(retty), retvaltmp)

        excinfoptr = cgutils.alloca_once(builder, ir.PointerType(excinfo_t),
                                         name="excinfo")

        arginfo = self._get_arg_packer(argtys)
        args = list(arginfo.as_arguments(builder, args))
        realargs = [retvaltmp, excinfoptr] + args
        # deal with attrs, it's fine to specify a load in a string like
        # "noinline fast" as per LLVM or equally as an iterable of individual
        # attributes.
        if attrs is None:
            _attrs = ()
        elif isinstance(attrs, Iterable) and not isinstance(attrs, str):
            _attrs = tuple(attrs)
        else:
            raise TypeError("attrs must be an iterable of strings or None")
        code = builder.call(callee, realargs, attrs=_attrs)
        status = self._get_return_status(builder, code,
                                         builder.load(excinfoptr))
        retval = builder.load(retvaltmp)
        out = self.context.get_returned_value(builder, resty, retval)
        return status, out


class ErrorModel(object):

    def __init__(self, call_conv):
        self.call_conv = call_conv

    def fp_zero_division(self, builder, exc_args=None, loc=None):
        if self.raise_on_fp_zero_division:
            self.call_conv.return_user_exc(builder, ZeroDivisionError, exc_args,
                                           loc)
            return True
        else:
            return False


class PythonErrorModel(ErrorModel):
    """
    The Python error model.  Any invalid FP input raises an exception.
    """
    raise_on_fp_zero_division = True


class NumpyErrorModel(ErrorModel):
    """
    In the Numpy error model, floating-point errors don't raise an
    exception.  The FPU exception state is inspected by Numpy at the
    end of a ufunc's execution and a warning is raised if appropriate.

    Note there's no easy way to set the FPU exception state from LLVM.
    Instructions known to set an FP exception can be optimized away:
        https://llvm.org/bugs/show_bug.cgi?id=6050
        http://lists.llvm.org/pipermail/llvm-dev/2014-September/076918.html
        http://lists.llvm.org/pipermail/llvm-commits/Week-of-Mon-20140929/237997.html
    """
    raise_on_fp_zero_division = False


error_models = {
    'python': PythonErrorModel,
    'numpy': NumpyErrorModel,
}


def create_error_model(model_name, context):
    """
    Create an error model instance for the given target context.
    """
    return error_models[model_name](context.call_conv)
