from llvmlite.ir import Constant, IRBuilder
import llvmlite.ir

from numba.core import types, config, cgutils


class _ArgManager(object):
    """
    A utility class to handle argument unboxing and cleanup
    """
    def __init__(self, context, builder, api, env_manager, endblk, nargs):
        self.context = context
        self.builder = builder
        self.api = api
        self.env_manager = env_manager
        self.arg_count = 0  # how many function arguments have been processed
        self.cleanups = []
        self.nextblk = endblk

    def add_arg(self, obj, ty):
        """
        Unbox argument and emit code that handles any error during unboxing.
        Args are cleaned up in reverse order of the parameter list, and
        cleanup begins as soon as unboxing of any argument fails. E.g. failure
        on arg2 will result in control flow going through:

            arg2.err -> arg1.err -> arg0.err -> arg.end (returns)
        """
        # Unbox argument
        native = self.api.to_native_value(ty, obj)

        # If an error occurred, go to the cleanup block for
        # the previous argument
        with cgutils.if_unlikely(self.builder, native.is_error):
            self.builder.branch(self.nextblk)

        # Define the cleanup function for the argument
        def cleanup_arg():
            # Native value reflection
            self.api.reflect_native_value(ty, native.value, self.env_manager)

            # Native value cleanup
            if native.cleanup is not None:
                native.cleanup()

            # NRT cleanup
            # (happens after the native value cleanup as the latter
            #  may need the native value)
            if self.context.enable_nrt:
                self.context.nrt.decref(self.builder, ty, native.value)

        self.cleanups.append(cleanup_arg)

        # Write the on-error cleanup block for this argument
        cleanupblk = self.builder.append_basic_block(
            "arg%d.err" % self.arg_count)
        with self.builder.goto_block(cleanupblk):
            cleanup_arg()
            # Go to next cleanup block
            self.builder.branch(self.nextblk)

        self.nextblk = cleanupblk
        self.arg_count += 1
        return native.value

    def emit_cleanup(self):
        """
        Emit the cleanup code after returning from the wrapped function.
        """
        for dtor in self.cleanups:
            dtor()


class _GilManager(object):
    """
    A utility class to handle releasing the GIL and then re-acquiring it
    again.
    """

    def __init__(self, builder, api, argman):
        self.builder = builder
        self.api = api
        self.argman = argman
        self.thread_state = api.save_thread()

    def emit_cleanup(self):
        self.api.restore_thread(self.thread_state)
        self.argman.emit_cleanup()


class PyCallWrapper(object):
    def __init__(self, context, module, func, fndesc, env, call_helper,
                 release_gil):
        self.context = context
        self.module = module
        self.func = func
        self.fndesc = fndesc
        self.env = env
        self.release_gil = release_gil

    def build(self):
        wrapname = self.fndesc.llvm_cpython_wrapper_name

        # This is the signature of PyCFunctionWithKeywords
        # (see CPython's methodobject.h)
        pyobj = self.context.get_argument_type(types.pyobject)
        wrapty = llvmlite.ir.FunctionType(pyobj, [pyobj, pyobj, pyobj])
        wrapper = llvmlite.ir.Function(self.module, wrapty, name=wrapname)

        builder = IRBuilder(wrapper.append_basic_block('entry'))

        # - `closure` will receive the `self` pointer stored in the
        #   PyCFunction object (see _dynfunc.c)
        # - `args` and `kws` will receive the tuple and dict objects
        #   of positional and keyword arguments, respectively.
        closure, args, kws = wrapper.args
        closure.name = 'py_closure'
        args.name = 'py_args'
        kws.name = 'py_kws'

        api = self.context.get_python_api(builder)
        self.build_wrapper(api, builder, closure, args, kws)

        return wrapper, api

    def build_wrapper(self, api, builder, closure, args, kws):
        nargs = len(self.fndesc.argtypes)

        objs = [api.alloca_obj() for _ in range(nargs)]
        parseok = api.unpack_tuple(args, self.fndesc.qualname,
                                   nargs, nargs, *objs)

        pred = builder.icmp_unsigned(
            '==',
            parseok,
            Constant(parseok.type, None))
        with cgutils.if_unlikely(builder, pred):
            builder.ret(api.get_null_object())

        # Block that returns after erroneous argument unboxing/cleanup
        endblk = builder.append_basic_block("arg.end")
        with builder.goto_block(endblk):
            builder.ret(api.get_null_object())

        # Get the Environment object
        env_manager = self.get_env(api, builder)

        cleanup_manager = _ArgManager(self.context, builder, api,
                                      env_manager, endblk, nargs)

        # Compute the arguments to the compiled Numba function.
        innerargs = []
        for obj, ty in zip(objs, self.fndesc.argtypes):
            if isinstance(ty, types.Omitted):
                # It's an omitted value => ignore dummy Python object
                innerargs.append(None)
            else:
                val = cleanup_manager.add_arg(builder.load(obj), ty)
                innerargs.append(val)

        if self.release_gil:
            cleanup_manager = _GilManager(builder, api, cleanup_manager)

        # We elect to not inline the top level user function into the call
        # wrapper, this incurs an overhead of a function call, however, it
        # increases optimisation stability in that the optimised user function
        # is what will actually be run and it is this function that all the
        # inspection tools "see". Further, this makes optimisation "stable" in
        # that calling the user function from e.g. C or from this wrapper will
        # result in the same code executing, were inlining permitted this may
        # not be the case as the inline could trigger additional optimisation
        # as the function goes into the wrapper, this resulting in the executing
        # instruction stream being different from that of the instruction stream
        # present in the user function.
        status, retval = self.context.call_conv.call_function(
            builder, self.func, self.fndesc.restype, self.fndesc.argtypes,
            innerargs, attrs=('noinline',))
        # Do clean up
        self.debug_print(builder, "# callwrapper: emit_cleanup")
        cleanup_manager.emit_cleanup()
        self.debug_print(builder, "# callwrapper: emit_cleanup end")

        # Determine return status
        with builder.if_then(status.is_ok, likely=True):
            # Ok => return boxed Python value
            with builder.if_then(status.is_none):
                api.return_none()

            retty = self._simplified_return_type()
            obj = api.from_native_return(retty, retval, env_manager)
            builder.ret(obj)

        # Error out
        self.context.call_conv.raise_error(builder, api, status)
        builder.ret(api.get_null_object())

    def get_env(self, api, builder):
        """Get the Environment object which is declared as a global
        in the module of the wrapped function.
        """
        envname = self.context.get_env_name(self.fndesc)
        gvptr = self.context.declare_env_global(builder.module, envname)
        envptr = builder.load(gvptr)

        env_body = self.context.get_env_body(builder, envptr)

        api.emit_environment_sentry(envptr, return_pyobject=True,
                                    debug_msg=self.fndesc.env_name)
        env_manager = api.get_env_manager(self.env, env_body, envptr)
        return env_manager

    def _simplified_return_type(self):
        """
        The NPM callconv has already converted simplified optional types.
        We can simply use the value type from it.
        """
        restype = self.fndesc.restype
        # Optional type
        if isinstance(restype, types.Optional):
            return restype.type
        else:
            return restype

    def debug_print(self, builder, msg):
        if config.DEBUG_JIT:
            self.context.debug_print(builder, "DEBUGJIT: {0}".format(msg))
