# -*- coding: utf-8 -*-


import collections
import functools
import sys
import types as pytypes
import uuid
import weakref
from contextlib import ExitStack
from abc import abstractmethod

from numba import _dispatcher
from numba.core import (
    utils, types, errors, typing, serialize, config, compiler, sigutils
)
from numba.core.compiler_lock import global_compiler_lock
from numba.core.typeconv.rules import default_type_manager
from numba.core.typing.templates import fold_arguments
from numba.core.typing.typeof import Purpose, typeof
from numba.core.bytecode import get_code_object
from numba.core.caching import NullCache, FunctionCache
from numba.core import entrypoints
import numba.core.event as ev


class OmittedArg(object):
    """
    A placeholder for omitted arguments with a default value.
    """

    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return "omitted arg(%r)" % (self.value,)

    @property
    def _numba_type_(self):
        return types.Omitted(self.value)


class _FunctionCompiler(object):
    def __init__(self, py_func, targetdescr, targetoptions, locals,
                 pipeline_class):
        self.py_func = py_func
        self.targetdescr = targetdescr
        self.targetoptions = targetoptions
        self.locals = locals
        self.pysig = utils.pysignature(self.py_func)
        self.pipeline_class = pipeline_class
        # Remember key=(args, return_type) combinations that will fail
        # compilation to avoid compilation attempt on them.  The values are
        # the exceptions.
        self._failed_cache = {}

    def fold_argument_types(self, args, kws):
        """
        Given positional and named argument types, fold keyword arguments
        and resolve defaults by inserting types.Omitted() instances.

        A (pysig, argument types) tuple is returned.
        """
        def normal_handler(index, param, value):
            return value

        def default_handler(index, param, default):
            return types.Omitted(default)

        def stararg_handler(index, param, values):
            return types.StarArgTuple(values)
        # For now, we take argument values from the @jit function
        args = fold_arguments(self.pysig, args, kws,
                              normal_handler,
                              default_handler,
                              stararg_handler)
        return self.pysig, args

    def compile(self, args, return_type):
        status, retval = self._compile_cached(args, return_type)
        if status:
            return retval
        else:
            raise retval

    def _compile_cached(self, args, return_type):
        key = tuple(args), return_type
        try:
            return False, self._failed_cache[key]
        except KeyError:
            pass

        try:
            retval = self._compile_core(args, return_type)
        except errors.TypingError as e:
            self._failed_cache[key] = e
            return False, e
        else:
            return True, retval

    def _compile_core(self, args, return_type):
        flags = compiler.Flags()
        self.targetdescr.options.parse_as_flags(flags, self.targetoptions)
        flags = self._customize_flags(flags)

        impl = self._get_implementation(args, {})
        cres = compiler.compile_extra(self.targetdescr.typing_context,
                                      self.targetdescr.target_context,
                                      impl,
                                      args=args, return_type=return_type,
                                      flags=flags, locals=self.locals,
                                      pipeline_class=self.pipeline_class)
        # Check typing error if object mode is used
        if cres.typing_error is not None and not flags.enable_pyobject:
            raise cres.typing_error
        return cres

    def get_globals_for_reduction(self):
        return serialize._get_function_globals_for_reduction(self.py_func)

    def _get_implementation(self, args, kws):
        return self.py_func

    def _customize_flags(self, flags):
        return flags


class _GeneratedFunctionCompiler(_FunctionCompiler):

    def __init__(self, py_func, targetdescr, targetoptions, locals,
                 pipeline_class):
        super(_GeneratedFunctionCompiler, self).__init__(
            py_func, targetdescr, targetoptions, locals, pipeline_class)
        self.impls = set()

    def get_globals_for_reduction(self):
        # This will recursively get the globals used by any nested
        # implementation function.
        return serialize._get_function_globals_for_reduction(self.py_func)

    def _get_implementation(self, args, kws):
        impl = self.py_func(*args, **kws)
        # Check the generating function and implementation signatures are
        # compatible, otherwise compiling would fail later.
        pysig = utils.pysignature(self.py_func)
        implsig = utils.pysignature(impl)
        ok = len(pysig.parameters) == len(implsig.parameters)
        if ok:
            for pyparam, implparam in zip(pysig.parameters.values(),
                                          implsig.parameters.values()):
                # We allow the implementation to omit default values, but
                # if it mentions them, they should have the same value...
                if (pyparam.name != implparam.name or
                    pyparam.kind != implparam.kind or
                    (implparam.default is not implparam.empty and
                     implparam.default != pyparam.default)):
                    ok = False
        if not ok:
            raise TypeError("generated implementation %s should be compatible "
                            "with signature '%s', but has signature '%s'"
                            % (impl, pysig, implsig))
        self.impls.add(impl)
        return impl


_CompileStats = collections.namedtuple(
    '_CompileStats', ('cache_path', 'cache_hits', 'cache_misses'))


class CompilingCounter(object):
    """
    A simple counter that increment in __enter__ and decrement in __exit__.
    """

    def __init__(self):
        self.counter = 0

    def __enter__(self):
        assert self.counter >= 0
        self.counter += 1

    def __exit__(self, *args, **kwargs):
        self.counter -= 1
        assert self.counter >= 0

    def __bool__(self):
        return self.counter > 0

    __nonzero__ = __bool__


class _DispatcherBase(_dispatcher.Dispatcher):
    """
    Common base class for dispatcher Implementations.
    """

    __numba__ = "py_func"

    def __init__(self, arg_count, py_func, pysig, can_fallback,
                 exact_match_required):
        self._tm = default_type_manager

        # A mapping of signatures to compile results
        self.overloads = collections.OrderedDict()

        self.py_func = py_func
        # other parts of Numba assume the old Python 2 name for code object
        self.func_code = get_code_object(py_func)
        # but newer python uses a different name
        self.__code__ = self.func_code
        # a place to keep an active reference to the types of the active call
        self._types_active_call = []
        # Default argument values match the py_func
        self.__defaults__ = py_func.__defaults__

        argnames = tuple(pysig.parameters)
        default_values = self.py_func.__defaults__ or ()
        defargs = tuple(OmittedArg(val) for val in default_values)
        try:
            lastarg = list(pysig.parameters.values())[-1]
        except IndexError:
            has_stararg = False
        else:
            has_stararg = lastarg.kind == lastarg.VAR_POSITIONAL
        _dispatcher.Dispatcher.__init__(self, self._tm.get_pointer(),
                                        arg_count, self._fold_args,
                                        argnames, defargs,
                                        can_fallback,
                                        has_stararg,
                                        exact_match_required)

        self.doc = py_func.__doc__
        self._compiling_counter = CompilingCounter()
        weakref.finalize(self, self._make_finalizer())

    def _compilation_chain_init_hook(self):
        """
        This will be called ahead of any part of compilation taking place (this
        even includes being ahead of working out the types of the arguments).
        This permits activities such as initialising extension entry points so
        that the compiler knows about additional externally defined types etc
        before it does anything.
        """
        entrypoints.init_all()

    def _reset_overloads(self):
        self._clear()
        self.overloads.clear()

    def _make_finalizer(self):
        """
        Return a finalizer function that will release references to
        related compiled functions.
        """
        overloads = self.overloads
        targetctx = self.targetctx

        # Early-bind utils.shutting_down() into the function's local namespace
        # (see issue #689)
        def finalizer(shutting_down=utils.shutting_down):
            # The finalizer may crash at shutdown, skip it (resources
            # will be cleared by the process exiting, anyway).
            if shutting_down():
                return
            # This function must *not* hold any reference to self:
            # we take care to bind the necessary objects in the closure.
            for cres in overloads.values():
                try:
                    targetctx.remove_user_function(cres.entry_point)
                except KeyError:
                    pass

        return finalizer

    @property
    def signatures(self):
        """
        Returns a list of compiled function signatures.
        """
        return list(self.overloads)

    @property
    def nopython_signatures(self):
        return [cres.signature for cres in self.overloads.values()
                if not cres.objectmode]

    def disable_compile(self, val=True):
        """Disable the compilation of new signatures at call time.
        """
        # If disabling compilation then there must be at least one signature
        assert (not val) or len(self.signatures) > 0
        self._can_compile = not val

    def add_overload(self, cres):
        args = tuple(cres.signature.args)
        sig = [a._code for a in args]
        self._insert(sig, cres.entry_point, cres.objectmode)
        self.overloads[args] = cres

    def fold_argument_types(self, args, kws):
        return self._compiler.fold_argument_types(args, kws)

    def get_call_template(self, args, kws):
        """
        Get a typing.ConcreteTemplate for this dispatcher and the given
        *args* and *kws* types.  This allows to resolve the return type.

        A (template, pysig, args, kws) tuple is returned.
        """
        # XXX how about a dispatcher template class automating the
        # following?

        # Fold keyword arguments and resolve default values
        pysig, args = self._compiler.fold_argument_types(args, kws)
        kws = {}
        # Ensure an overload is available
        if self._can_compile:
            self.compile(tuple(args))

        # Create function type for typing
        func_name = self.py_func.__name__
        name = "CallTemplate({0})".format(func_name)
        # The `key` isn't really used except for diagnosis here,
        # so avoid keeping a reference to `cfunc`.
        call_template = typing.make_concrete_template(
            name, key=func_name, signatures=self.nopython_signatures)
        return call_template, pysig, args, kws

    def get_overload(self, sig):
        """
        Return the compiled function for the given signature.
        """
        args, return_type = sigutils.normalize_signature(sig)
        return self.overloads[tuple(args)].entry_point

    @property
    def is_compiling(self):
        """
        Whether a specialization is currently being compiled.
        """
        return self._compiling_counter

    def _compile_for_args(self, *args, **kws):
        """
        For internal use.  Compile a specialized version of the function
        for the given *args* and *kws*, and return the resulting callable.
        """
        assert not kws
        # call any initialisation required for the compilation chain (e.g.
        # extension point registration).
        self._compilation_chain_init_hook()

        def error_rewrite(e, issue_type):
            """
            Rewrite and raise Exception `e` with help supplied based on the
            specified issue_type.
            """
            if config.SHOW_HELP:
                help_msg = errors.error_extras[issue_type]
                e.patch_message('\n'.join((str(e).rstrip(), help_msg)))
            if config.FULL_TRACEBACKS:
                raise e
            else:
                raise e.with_traceback(None)

        argtypes = []
        for a in args:
            if isinstance(a, OmittedArg):
                argtypes.append(types.Omitted(a.value))
            else:
                argtypes.append(self.typeof_pyval(a))

        return_val = None
        try:
            return_val = self.compile(tuple(argtypes))
        except errors.ForceLiteralArg as e:
            # Received request for compiler re-entry with the list of arguments
            # indicated by e.requested_args.
            # First, check if any of these args are already Literal-ized
            already_lit_pos = [i for i in e.requested_args
                               if isinstance(args[i], types.Literal)]
            if already_lit_pos:
                # Abort compilation if any argument is already a Literal.
                # Letting this continue will cause infinite compilation loop.
                m = ("Repeated literal typing request.\n"
                     "{}.\n"
                     "This is likely caused by an error in typing. "
                     "Please see nested and suppressed exceptions.")
                info = ', '.join('Arg #{} is {}'.format(i, args[i])
                                 for i in sorted(already_lit_pos))
                raise errors.CompilerError(m.format(info))
            # Convert requested arguments into a Literal.
            args = [(types.literal
                     if i in e.requested_args
                     else lambda x: x)(args[i])
                    for i, v in enumerate(args)]
            # Re-enter compilation with the Literal-ized arguments
            return_val = self._compile_for_args(*args)

        except errors.TypingError as e:
            # Intercept typing error that may be due to an argument
            # that failed inferencing as a Numba type
            failed_args = []
            for i, arg in enumerate(args):
                val = arg.value if isinstance(arg, OmittedArg) else arg
                try:
                    tp = typeof(val, Purpose.argument)
                except ValueError as typeof_exc:
                    failed_args.append((i, str(typeof_exc)))
                else:
                    if tp is None:
                        failed_args.append(
                            (i, f"cannot determine Numba type of value {val}"))
            if failed_args:
                # Patch error message to ease debugging
                args_str = "\n".join(
                    f"- argument {i}: {err}" for i, err in failed_args
                )
                msg = (f"{str(e).rstrip()} \n\nThis error may have been caused "
                       f"by the following argument(s):\n{args_str}\n")
                e.patch_message(msg)

            error_rewrite(e, 'typing')
        except errors.UnsupportedError as e:
            # Something unsupported is present in the user code, add help info
            error_rewrite(e, 'unsupported_error')
        except (errors.NotDefinedError, errors.RedefinedError,
                errors.VerificationError) as e:
            # These errors are probably from an issue with either the code
            # supplied being syntactically or otherwise invalid
            error_rewrite(e, 'interpreter')
        except errors.ConstantInferenceError as e:
            # this is from trying to infer something as constant when it isn't
            # or isn't supported as a constant
            error_rewrite(e, 'constant_inference')
        except Exception as e:
            if config.SHOW_HELP:
                if hasattr(e, 'patch_message'):
                    help_msg = errors.error_extras['reportable']
                    e.patch_message('\n'.join((str(e).rstrip(), help_msg)))
            # ignore the FULL_TRACEBACKS config, this needs reporting!
            raise e
        finally:
            self._types_active_call = []
        return return_val

    def inspect_llvm(self, signature=None):
        """Get the LLVM intermediate representation generated by compilation.

        Parameters
        ----------
        signature : tuple of numba types, optional
            Specify a signature for which to obtain the LLVM IR. If None, the
            IR is returned for all available signatures.

        Returns
        -------
        llvm : dict[signature, str] or str
            Either the LLVM IR string for the specified signature, or, if no
            signature was given, a dictionary mapping signatures to LLVM IR
            strings.
        """
        if signature is not None:
            lib = self.overloads[signature].library
            return lib.get_llvm_str()

        return dict((sig, self.inspect_llvm(sig)) for sig in self.signatures)

    def inspect_asm(self, signature=None):
        """Get the generated assembly code.

        Parameters
        ----------
        signature : tuple of numba types, optional
            Specify a signature for which to obtain the assembly code. If
            None, the assembly code is returned for all available signatures.

        Returns
        -------
        asm : dict[signature, str] or str
            Either the assembly code for the specified signature, or, if no
            signature was given, a dictionary mapping signatures to assembly
            code.
        """
        if signature is not None:
            lib = self.overloads[signature].library
            return lib.get_asm_str()

        return dict((sig, self.inspect_asm(sig)) for sig in self.signatures)

    def inspect_types(self, file=None, signature=None,
                      pretty=False, style='default', **kwargs):
        """Print/return Numba intermediate representation (IR)-annotated code.

        Parameters
        ----------
        file : file-like object, optional
            File to which to print. Defaults to sys.stdout if None. Must be
            None if ``pretty=True``.
        signature : tuple of numba types, optional
            Print/return the intermediate representation for only the given
            signature. If None, the IR is printed for all available signatures.
        pretty : bool, optional
            If True, an Annotate object will be returned that can render the
            IR with color highlighting in Jupyter and IPython. ``file`` must
            be None if ``pretty`` is True. Additionally, the ``pygments``
            library must be installed for ``pretty=True``.
        style : str, optional
            Choose a style for rendering. Ignored if ``pretty`` is ``False``.
            This is directly consumed by ``pygments`` formatters. To see a
            list of available styles, import ``pygments`` and run
            ``list(pygments.styles.get_all_styles())``.

        Returns
        -------
        annotated : Annotate object, optional
            Only returned if ``pretty=True``, otherwise this function is only
            used for its printing side effect. If ``pretty=True``, an Annotate
            object is returned that can render itself in Jupyter and IPython.
        """
        overloads = self.overloads
        if signature is not None:
            overloads = {signature: self.overloads[signature]}

        if not pretty:
            if file is None:
                file = sys.stdout

            for ver, res in overloads.items():
                print("%s %s" % (self.py_func.__name__, ver), file=file)
                print('-' * 80, file=file)
                print(res.type_annotation, file=file)
                print('=' * 80, file=file)
        else:
            if file is not None:
                raise ValueError("`file` must be None if `pretty=True`")
            from numba.core.annotations.pretty_annotate import Annotate
            return Annotate(self, signature=signature, style=style)

    def inspect_cfg(self, signature=None, show_wrapper=None, **kwargs):
        """
        For inspecting the CFG of the function.

        By default the CFG of the user function is shown.  The *show_wrapper*
        option can be set to "python" or "cfunc" to show the python wrapper
        function or the *cfunc* wrapper function, respectively.

        Parameters accepted in kwargs
        -----------------------------
        filename : string, optional
            the name of the output file, if given this will write the output to
            filename
        view : bool, optional
            whether to immediately view the optional output file
        highlight : bool, set, dict, optional
            what, if anything, to highlight, options are:
            { incref : bool, # highlight NRT_incref calls
              decref : bool, # highlight NRT_decref calls
              returns : bool, # highlight exits which are normal returns
              raises : bool, # highlight exits which are from raise
              meminfo : bool, # highlight calls to NRT*meminfo
              branches : bool, # highlight true/false branches
             }
            Default is True which sets all of the above to True. Supplying a set
            of strings is also accepted, these are interpreted as key:True with
            respect to the above dictionary. e.g. {'incref', 'decref'} would
            switch on highlighting on increfs and decrefs.
        interleave: bool, set, dict, optional
            what, if anything, to interleave in the LLVM IR, options are:
            { python: bool # interleave python source code with the LLVM IR
              lineinfo: bool # interleave line information markers with the LLVM
                             # IR
            }
            Default is True which sets all of the above to True. Supplying a set
            of strings is also accepted, these are interpreted as key:True with
            respect to the above dictionary. e.g. {'python',} would
            switch on interleaving of python source code in the LLVM IR.
        strip_ir : bool, optional
            Default is False. If set to True all LLVM IR that is superfluous to
            that requested in kwarg `highlight` will be removed.
        show_key : bool, optional
            Default is True. Create a "key" for the highlighting in the rendered
            CFG.
        fontsize : int, optional
            Default is 8. Set the fontsize in the output to this value.
        """
        if signature is not None:
            cres = self.overloads[signature]
            lib = cres.library
            if show_wrapper == 'python':
                fname = cres.fndesc.llvm_cpython_wrapper_name
            elif show_wrapper == 'cfunc':
                fname = cres.fndesc.llvm_cfunc_wrapper_name
            else:
                fname = cres.fndesc.mangled_name
            return lib.get_function_cfg(fname, py_func=self.py_func, **kwargs)

        return dict((sig, self.inspect_cfg(sig, show_wrapper=show_wrapper))
                    for sig in self.signatures)

    def inspect_disasm_cfg(self, signature=None):
        """
        For inspecting the CFG of the disassembly of the function.

        Requires python package: r2pipe
        Requires radare2 binary on $PATH.
        Notebook rendering requires python package: graphviz

        signature : tuple of Numba types, optional
            Print/return the disassembly CFG for only the given signatures.
            If None, the IR is printed for all available signatures.
        """
        if signature is not None:
            cres = self.overloads[signature]
            lib = cres.library
            return lib.get_disasm_cfg(cres.fndesc.mangled_name)

        return dict((sig, self.inspect_disasm_cfg(sig))
                    for sig in self.signatures)

    def get_annotation_info(self, signature=None):
        """
        Gets the annotation information for the function specified by
        signature. If no signature is supplied a dictionary of signature to
        annotation information is returned.
        """
        signatures = self.signatures if signature is None else [signature]
        out = collections.OrderedDict()
        for sig in signatures:
            cres = self.overloads[sig]
            ta = cres.type_annotation
            key = (ta.func_id.filename + ':' + str(ta.func_id.firstlineno + 1),
                   ta.signature)
            out[key] = ta.annotate_raw()[key]
        return out

    def _explain_ambiguous(self, *args, **kws):
        """
        Callback for the C _Dispatcher object.
        """
        assert not kws, "kwargs not handled"
        args = tuple([self.typeof_pyval(a) for a in args])
        # The order here must be deterministic for testing purposes, which
        # is ensured by the OrderedDict.
        sigs = self.nopython_signatures
        # This will raise
        self.typingctx.resolve_overload(self.py_func, sigs, args, kws,
                                        allow_ambiguous=False)

    def _explain_matching_error(self, *args, **kws):
        """
        Callback for the C _Dispatcher object.
        """
        assert not kws, "kwargs not handled"
        args = [self.typeof_pyval(a) for a in args]
        msg = ("No matching definition for argument type(s) %s"
               % ', '.join(map(str, args)))
        raise TypeError(msg)

    def _search_new_conversions(self, *args, **kws):
        """
        Callback for the C _Dispatcher object.
        Search for approximately matching signatures for the given arguments,
        and ensure the corresponding conversions are registered in the C++
        type manager.
        """
        assert not kws, "kwargs not handled"
        args = [self.typeof_pyval(a) for a in args]
        found = False
        for sig in self.nopython_signatures:
            conv = self.typingctx.install_possible_conversions(args, sig.args)
            if conv:
                found = True
        return found

    def __repr__(self):
        return "%s(%s)" % (type(self).__name__, self.py_func)

    def typeof_pyval(self, val):
        """
        Resolve the Numba type of Python value *val*.
        This is called from numba._dispatcher as a fallback if the native code
        cannot decide the type.
        """
        # Not going through the resolve_argument_type() indirection
        # can save a couple µs.
        try:
            tp = typeof(val, Purpose.argument)
        except ValueError:
            tp = types.pyobject
        else:
            if tp is None:
                tp = types.pyobject
        self._types_active_call.append(tp)
        return tp

    def _callback_add_timer(self, duration, cres, lock_name):
        md = cres.metadata
        # md can be None when code is loaded from cache
        if md is not None:
            timers = md.setdefault("timers", {})
            if lock_name not in timers:
                # Only write if the metadata does not exist
                timers[lock_name] = duration
            else:
                msg = f"'{lock_name} metadata is already defined."
                raise AssertionError(msg)

    def _callback_add_compiler_timer(self, duration, cres):
        return self._callback_add_timer(duration, cres,
                                        lock_name="compiler_lock")

    def _callback_add_llvm_timer(self, duration, cres):
        return self._callback_add_timer(duration, cres,
                                        lock_name="llvm_lock")


class _MemoMixin:
    __uuid = None
    # A {uuid -> instance} mapping, for deserialization
    _memo = weakref.WeakValueDictionary()
    # hold refs to last N functions deserialized, retaining them in _memo
    # regardless of whether there is another reference
    _recent = collections.deque(maxlen=config.FUNCTION_CACHE_SIZE)

    @property
    def _uuid(self):
        """
        An instance-specific UUID, to avoid multiple deserializations of
        a given instance.

        Note: this is lazily-generated, for performance reasons.
        """
        u = self.__uuid
        if u is None:
            u = str(uuid.uuid4())
            self._set_uuid(u)
        return u

    def _set_uuid(self, u):
        assert self.__uuid is None
        self.__uuid = u
        self._memo[u] = self
        self._recent.append(self)


class Dispatcher(serialize.ReduceMixin, _MemoMixin, _DispatcherBase):
    """
    Implementation of user-facing dispatcher objects (i.e. created using
    the @jit decorator).
    This is an abstract base class. Subclasses should define the targetdescr
    class attribute.
    """
    _fold_args = True

    __numba__ = 'py_func'

    def __init__(self, py_func, locals={}, targetoptions={},
                 pipeline_class=compiler.Compiler):
        """
        Parameters
        ----------
        py_func: function object to be compiled
        locals: dict, optional
            Mapping of local variable names to Numba types.  Used to override
            the types deduced by the type inference engine.
        targetoptions: dict, optional
            Target-specific config options.
        pipeline_class: type numba.compiler.CompilerBase
            The compiler pipeline type.
        """
        self.typingctx = self.targetdescr.typing_context
        self.targetctx = self.targetdescr.target_context

        pysig = utils.pysignature(py_func)
        arg_count = len(pysig.parameters)
        can_fallback = not targetoptions.get('nopython', False)

        _DispatcherBase.__init__(self, arg_count, py_func, pysig, can_fallback,
                                 exact_match_required=False)

        functools.update_wrapper(self, py_func)

        self.targetoptions = targetoptions
        self.locals = locals
        self._cache = NullCache()
        compiler_class = _FunctionCompiler
        self._compiler = compiler_class(py_func, self.targetdescr,
                                        targetoptions, locals, pipeline_class)
        self._cache_hits = collections.Counter()
        self._cache_misses = collections.Counter()

        self._type = types.Dispatcher(self)
        self.typingctx.insert_global(self, self._type)

    def dump(self, tab=''):
        print(f'{tab}DUMP {type(self).__name__}[{self.py_func.__name__}'
              f', type code={self._type._code}]')
        for cres in self.overloads.values():
            cres.dump(tab=tab + '  ')
        print(f'{tab}END DUMP {type(self).__name__}[{self.py_func.__name__}]')

    @property
    def _numba_type_(self):
        return types.Dispatcher(self)

    def enable_caching(self):
        self._cache = FunctionCache(self.py_func)

    def __get__(self, obj, objtype=None):
        '''Allow a JIT function to be bound as a method to an object'''
        if obj is None:  # Unbound method
            return self
        else:  # Bound method
            return pytypes.MethodType(self, obj)

    def _reduce_states(self):
        """
        Reduce the instance for pickling.  This will serialize
        the original function as well the compilation options and
        compiled signatures, but not the compiled code itself.

        NOTE: part of ReduceMixin protocol
        """
        if self._can_compile:
            sigs = []
        else:
            sigs = [cr.signature for cr in self.overloads.values()]

        return dict(
            uuid=str(self._uuid),
            py_func=self.py_func,
            locals=self.locals,
            targetoptions=self.targetoptions,
            can_compile=self._can_compile,
            sigs=sigs,
        )

    @classmethod
    def _rebuild(cls, uuid, py_func, locals, targetoptions,
                 can_compile, sigs):
        """
        Rebuild an Dispatcher instance after it was __reduce__'d.

        NOTE: part of ReduceMixin protocol
        """
        try:
            return cls._memo[uuid]
        except KeyError:
            pass
        self = cls(py_func, locals, targetoptions)
        # Make sure this deserialization will be merged with subsequent ones
        self._set_uuid(uuid)
        for sig in sigs:
            self.compile(sig)
        self._can_compile = can_compile
        return self

    def compile(self, sig):
        with ExitStack() as scope:
            cres = None

            def cb_compiler(dur):
                if cres is not None:
                    self._callback_add_compiler_timer(dur, cres)

            def cb_llvm(dur):
                if cres is not None:
                    self._callback_add_llvm_timer(dur, cres)

            scope.enter_context(ev.install_timer("numba:compiler_lock",
                                                 cb_compiler))
            scope.enter_context(ev.install_timer("numba:llvm_lock", cb_llvm))
            scope.enter_context(global_compiler_lock)

            if not self._can_compile:
                raise RuntimeError("compilation disabled")
            # Use counter to track recursion compilation depth
            with self._compiling_counter:
                args, return_type = sigutils.normalize_signature(sig)
                # Don't recompile if signature already exists
                existing = self.overloads.get(tuple(args))
                if existing is not None:
                    return existing.entry_point
                # Try to load from disk cache
                cres = self._cache.load_overload(sig, self.targetctx)
                if cres is not None:
                    self._cache_hits[sig] += 1
                    # XXX fold this in add_overload()? (also see compiler.py)
                    if not cres.objectmode:
                        self.targetctx.insert_user_function(cres.entry_point,
                                                            cres.fndesc,
                                                            [cres.library])
                    self.add_overload(cres)
                    return cres.entry_point

                self._cache_misses[sig] += 1
                ev_details = dict(
                    dispatcher=self,
                    args=args,
                    return_type=return_type,
                )
                with ev.trigger_event("numba:compile", data=ev_details):
                    try:
                        cres = self._compiler.compile(args, return_type)
                    except errors.ForceLiteralArg as e:
                        def folded(args, kws):
                            return self._compiler.fold_argument_types(args,
                                                                      kws)[1]
                        raise e.bind_fold_arguments(folded)
                    self.add_overload(cres)
                self._cache.save_overload(sig, cres)
                return cres.entry_point

    def get_compile_result(self, sig):
        """Compile (if needed) and return the compilation result with the
        given signature.

        Returns ``CompileResult``.
        Raises ``NumbaError`` if the signature is incompatible.
        """
        atypes = tuple(sig.args)
        if atypes not in self.overloads:
            if self._can_compile:
                # Compiling may raise any NumbaError
                self.compile(atypes)
            else:
                msg = f"{sig} not available and compilation disabled"
                raise errors.TypingError(msg)
        return self.overloads[atypes]

    def recompile(self):
        """
        Recompile all signatures afresh.
        """
        sigs = list(self.overloads)
        old_can_compile = self._can_compile
        # Ensure the old overloads are disposed of,
        # including compiled functions.
        self._make_finalizer()()
        self._reset_overloads()
        self._cache.flush()
        self._can_compile = True
        try:
            for sig in sigs:
                self.compile(sig)
        finally:
            self._can_compile = old_can_compile

    @property
    def stats(self):
        return _CompileStats(
            cache_path=self._cache.cache_path,
            cache_hits=self._cache_hits,
            cache_misses=self._cache_misses,
        )

    def parallel_diagnostics(self, signature=None, level=1):
        """
        Print parallel diagnostic information for the given signature. If no
        signature is present it is printed for all known signatures. level is
        used to adjust the verbosity, level=1 (default) is minimal verbosity,
        and 2, 3, and 4 provide increasing levels of verbosity.
        """
        def dump(sig):
            ol = self.overloads[sig]
            pfdiag = ol.metadata.get('parfor_diagnostics', None)
            if pfdiag is None:
                msg = "No parfors diagnostic available, is 'parallel=True' set?"
                raise ValueError(msg)
            pfdiag.dump(level)
        if signature is not None:
            dump(signature)
        else:
            [dump(sig) for sig in self.signatures]

    def get_metadata(self, signature=None):
        """
        Obtain the compilation metadata for a given signature.
        """
        if signature is not None:
            return self.overloads[signature].metadata
        else:
            return dict(
                (sig,self.overloads[sig].metadata) for sig in self.signatures
            )

    def get_function_type(self):
        """Return unique function type of dispatcher when possible, otherwise
        return None.

        A Dispatcher instance has unique function type when it
        contains exactly one compilation result and its compilation
        has been disabled (via its disable_compile method).
        """
        if not self._can_compile and len(self.overloads) == 1:
            cres = tuple(self.overloads.values())[0]
            return types.FunctionType(cres.signature)


class LiftedCode(serialize.ReduceMixin, _MemoMixin, _DispatcherBase):
    """
    Implementation of the hidden dispatcher objects used for lifted code
    (a lifted loop is really compiled as a separate function).
    """
    _fold_args = False
    can_cache = False

    def __init__(self, func_ir, typingctx, targetctx, flags, locals):
        self.func_ir = func_ir
        self.lifted_from = None

        self.typingctx = typingctx
        self.targetctx = targetctx
        self.flags = flags
        self.locals = locals

        _DispatcherBase.__init__(self, self.func_ir.arg_count,
                                 self.func_ir.func_id.func,
                                 self.func_ir.func_id.pysig,
                                 can_fallback=True,
                                 exact_match_required=False)

    def _reduce_states(self):
        """
        Reduce the instance for pickling.  This will serialize
        the original function as well the compilation options and
        compiled signatures, but not the compiled code itself.

        NOTE: part of ReduceMixin protocol
        """
        return dict(
            uuid=self._uuid, func_ir=self.func_ir, flags=self.flags,
            locals=self.locals, extras=self._reduce_extras(),
        )

    def _reduce_extras(self):
        """
        NOTE: sub-class can override to add extra states
        """
        return {}

    @classmethod
    def _rebuild(cls, uuid, func_ir, flags, locals, extras):
        """
        Rebuild an Dispatcher instance after it was __reduce__'d.

        NOTE: part of ReduceMixin protocol
        """
        try:
            return cls._memo[uuid]
        except KeyError:
            pass

        # NOTE: We are assuming that this is must be cpu_target, which is true
        #       for now.
        # TODO: refactor this to not assume on `cpu_target`

        from numba.core import registry
        typingctx = registry.cpu_target.typing_context
        targetctx = registry.cpu_target.target_context

        self = cls(func_ir, typingctx, targetctx, flags, locals, **extras)
        self._set_uuid(uuid)
        return self

    def get_source_location(self):
        """Return the starting line number of the loop.
        """
        return self.func_ir.loc.line

    def _pre_compile(self, args, return_type, flags):
        """Pre-compile actions
        """
        pass

    @abstractmethod
    def compile(self, sig):
        """Lifted code should implement a compilation method that will return
        a CompileResult.entry_point for the given signature."""
        pass

    def _get_dispatcher_for_current_target(self):
        # Lifted code does not honor the target switch currently.
        # No work has been done to check if this can be allowed.
        return self


class LiftedLoop(LiftedCode):
    def _pre_compile(self, args, return_type, flags):
        assert not flags.enable_looplift, "Enable looplift flags is on"

    def compile(self, sig):
        with ExitStack() as scope:
            cres = None

            def cb_compiler(dur):
                if cres is not None:
                    self._callback_add_compiler_timer(dur, cres)

            def cb_llvm(dur):
                if cres is not None:
                    self._callback_add_llvm_timer(dur, cres)

            scope.enter_context(ev.install_timer("numba:compiler_lock",
                                                 cb_compiler))
            scope.enter_context(ev.install_timer("numba:llvm_lock", cb_llvm))
            scope.enter_context(global_compiler_lock)

            # Use counter to track recursion compilation depth
            with self._compiling_counter:
                # XXX this is mostly duplicated from Dispatcher.
                flags = self.flags
                args, return_type = sigutils.normalize_signature(sig)

                # Don't recompile if signature already exists
                # (e.g. if another thread compiled it before we got the lock)
                existing = self.overloads.get(tuple(args))
                if existing is not None:
                    return existing.entry_point

                self._pre_compile(args, return_type, flags)

                # copy the flags, use nopython first
                npm_loop_flags = flags.copy()
                npm_loop_flags.force_pyobject = False

                pyobject_loop_flags = flags.copy()
                pyobject_loop_flags.force_pyobject = True

                # Clone IR to avoid (some of the) mutation in the rewrite pass
                cloned_func_ir = self.func_ir.copy()

                ev_details = dict(
                    dispatcher=self,
                    args=args,
                    return_type=return_type,
                )
                with ev.trigger_event("numba:compile", data=ev_details):
                    # this emulates "object mode fall-back", try nopython, if it
                    # fails, then try again in object mode.
                    try:
                        cres = compiler.compile_ir(typingctx=self.typingctx,
                                                   targetctx=self.targetctx,
                                                   func_ir=cloned_func_ir,
                                                   args=args,
                                                   return_type=return_type,
                                                   flags=npm_loop_flags,
                                                   locals=self.locals,
                                                   lifted=(),
                                                   lifted_from=self.lifted_from,
                                                   is_lifted_loop=True,)
                    except errors.TypingError:
                        cres = compiler.compile_ir(typingctx=self.typingctx,
                                                   targetctx=self.targetctx,
                                                   func_ir=cloned_func_ir,
                                                   args=args,
                                                   return_type=return_type,
                                                   flags=pyobject_loop_flags,
                                                   locals=self.locals,
                                                   lifted=(),
                                                   lifted_from=self.lifted_from,
                                                   is_lifted_loop=True,)
                    # Check typing error if object mode is used
                    if (cres.typing_error is not None):
                        raise cres.typing_error
                    self.add_overload(cres)
                return cres.entry_point


class LiftedWith(LiftedCode):

    can_cache = True

    def _reduce_extras(self):
        return dict(output_types=self.output_types)

    @property
    def _numba_type_(self):
        return types.Dispatcher(self)

    def get_call_template(self, args, kws):
        """
        Get a typing.ConcreteTemplate for this dispatcher and the given
        *args* and *kws* types.  This enables the resolving of the return type.

        A (template, pysig, args, kws) tuple is returned.
        """
        # Ensure an overload is available
        if self._can_compile:
            self.compile(tuple(args))

        pysig = None
        # Create function type for typing
        func_name = self.py_func.__name__
        name = "CallTemplate({0})".format(func_name)
        # The `key` isn't really used except for diagnosis here,
        # so avoid keeping a reference to `cfunc`.
        call_template = typing.make_concrete_template(
            name, key=func_name, signatures=self.nopython_signatures)
        return call_template, pysig, args, kws

    def compile(self, sig):
        # this is similar to LiftedLoop's compile but does not have the
        # "fallback" to object mode part.
        with ExitStack() as scope:
            cres = None

            def cb_compiler(dur):
                if cres is not None:
                    self._callback_add_compiler_timer(dur, cres)

            def cb_llvm(dur):
                if cres is not None:
                    self._callback_add_llvm_timer(dur, cres)

            scope.enter_context(ev.install_timer("numba:compiler_lock",
                                                 cb_compiler))
            scope.enter_context(ev.install_timer("numba:llvm_lock", cb_llvm))
            scope.enter_context(global_compiler_lock)

            # Use counter to track recursion compilation depth
            with self._compiling_counter:
                # XXX this is mostly duplicated from Dispatcher.
                flags = self.flags
                args, return_type = sigutils.normalize_signature(sig)

                # Don't recompile if signature already exists
                # (e.g. if another thread compiled it before we got the lock)
                existing = self.overloads.get(tuple(args))
                if existing is not None:
                    return existing.entry_point

                self._pre_compile(args, return_type, flags)

                # Clone IR to avoid (some of the) mutation in the rewrite pass
                cloned_func_ir = self.func_ir.copy()

                ev_details = dict(
                    dispatcher=self,
                    args=args,
                    return_type=return_type,
                )
                with ev.trigger_event("numba:compile", data=ev_details):
                    cres = compiler.compile_ir(typingctx=self.typingctx,
                                               targetctx=self.targetctx,
                                               func_ir=cloned_func_ir,
                                               args=args,
                                               return_type=return_type,
                                               flags=flags, locals=self.locals,
                                               lifted=(),
                                               lifted_from=self.lifted_from,
                                               is_lifted_loop=True,)

                    # Check typing error if object mode is used
                    if (cres.typing_error is not None and
                            not flags.enable_pyobject):
                        raise cres.typing_error
                    self.add_overload(cres)
                return cres.entry_point


class ObjModeLiftedWith(LiftedWith):
    def __init__(self, *args, **kwargs):
        self.output_types = kwargs.pop('output_types', None)
        super(LiftedWith, self).__init__(*args, **kwargs)
        if not self.flags.force_pyobject:
            raise ValueError("expecting `flags.force_pyobject`")
        if self.output_types is None:
            raise TypeError('`output_types` must be provided')
        # switch off rewrites, they have no effect
        self.flags.no_rewrites = True

    @property
    def _numba_type_(self):
        return types.ObjModeDispatcher(self)

    def get_call_template(self, args, kws):
        """
        Get a typing.ConcreteTemplate for this dispatcher and the given
        *args* and *kws* types.  This enables the resolving of the return type.

        A (template, pysig, args, kws) tuple is returned.
        """
        assert not kws
        self._legalize_arg_types(args)
        # Coerce to object mode
        args = [types.ffi_forced_object] * len(args)

        if self._can_compile:
            self.compile(tuple(args))

        signatures = [typing.signature(self.output_types, *args)]
        pysig = None
        func_name = self.py_func.__name__
        name = "CallTemplate({0})".format(func_name)
        call_template = typing.make_concrete_template(
            name, key=func_name, signatures=signatures)

        return call_template, pysig, args, kws

    def _legalize_arg_types(self, args):
        for i, a in enumerate(args, start=1):
            if isinstance(a, types.List):
                msg = (
                    'Does not support list type inputs into '
                    'with-context for arg {}'
                )
                raise errors.TypingError(msg.format(i))
            elif isinstance(a, types.Dispatcher):
                msg = (
                    'Does not support function type inputs into '
                    'with-context for arg {}'
                )
                raise errors.TypingError(msg.format(i))

    @global_compiler_lock
    def compile(self, sig):
        args, _ = sigutils.normalize_signature(sig)
        sig = (types.ffi_forced_object,) * len(args)
        return super().compile(sig)


# Initialize typeof machinery
_dispatcher.typeof_init(
    OmittedArg,
    dict((str(t), t._code) for t in types.number_domain))
