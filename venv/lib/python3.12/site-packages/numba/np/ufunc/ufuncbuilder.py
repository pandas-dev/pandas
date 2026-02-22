# -*- coding: utf-8 -*-

import inspect
import warnings
from contextlib import contextmanager

from numba.core import config, targetconfig
from numba.core.decorators import jit
from numba.core.descriptors import TargetDescriptor
from numba.core.extending import is_jitted
from numba.core.errors import NumbaDeprecationWarning
from numba.core.options import TargetOptions, include_default_options
from numba.core.registry import cpu_target
from numba.core.target_extension import dispatcher_registry, target_registry
from numba.core import utils, types, serialize, compiler, sigutils
from numba.np.numpy_support import as_dtype
from numba.np.ufunc import _internal
from numba.np.ufunc.sigparse import parse_signature
from numba.np.ufunc.wrappers import build_ufunc_wrapper, build_gufunc_wrapper
from numba.core.caching import FunctionCache, NullCache
from numba.core.compiler_lock import global_compiler_lock


_options_mixin = include_default_options(
    "nopython",
    "forceobj",
    "boundscheck",
    "fastmath",
    "writable_args"
)


class UFuncTargetOptions(_options_mixin, TargetOptions):

    def finalize(self, flags, options):
        if not flags.is_set("enable_pyobject"):
            flags.enable_pyobject = True

        if not flags.is_set("enable_looplift"):
            flags.enable_looplift = True

        flags.inherit_if_not_set("nrt", default=True)

        if not flags.is_set("debuginfo"):
            flags.debuginfo = config.DEBUGINFO_DEFAULT

        if not flags.is_set("boundscheck"):
            flags.boundscheck = flags.debuginfo

        flags.enable_pyobject_looplift = True

        flags.inherit_if_not_set("fastmath")


class UFuncTarget(TargetDescriptor):
    options = UFuncTargetOptions

    def __init__(self):
        super().__init__('ufunc')

    @property
    def typing_context(self):
        return cpu_target.typing_context

    @property
    def target_context(self):
        return cpu_target.target_context


ufunc_target = UFuncTarget()


class UFuncDispatcher(serialize.ReduceMixin):
    """
    An object handling compilation of various signatures for a ufunc.
    """
    targetdescr = ufunc_target

    def __init__(self, py_func, locals=None, targetoptions=None):
        if locals is None:
            locals = {}
        if targetoptions is None:
            targetoptions = {}
        self.py_func = py_func
        self.overloads = utils.UniqueDict()
        self.targetoptions = targetoptions
        self.locals = locals
        self.cache = NullCache()

    def _reduce_states(self):
        """
        NOTE: part of ReduceMixin protocol
        """
        return dict(
            pyfunc=self.py_func,
            locals=self.locals,
            targetoptions=self.targetoptions,
        )

    @classmethod
    def _rebuild(cls, pyfunc, locals, targetoptions):
        """
        NOTE: part of ReduceMixin protocol
        """
        return cls(py_func=pyfunc, locals=locals, targetoptions=targetoptions)

    def enable_caching(self):
        self.cache = FunctionCache(self.py_func)

    def compile(self, sig, locals=None, **targetoptions):
        if locals is None:
            locals = {}
        locs = self.locals.copy()
        locs.update(locals)

        topt = self.targetoptions.copy()
        topt.update(targetoptions)

        flags = compiler.Flags()
        self.targetdescr.options.parse_as_flags(flags, topt)

        flags.no_cpython_wrapper = True
        flags.error_model = "numpy"
        # Disable loop lifting
        # The feature requires a real
        #  python function
        flags.enable_looplift = False

        return self._compile_core(sig, flags, locals)

    def _compile_core(self, sig, flags, locals):
        """
        Trigger the compiler on the core function or load a previously
        compiled version from the cache.  Returns the CompileResult.
        """
        typingctx = self.targetdescr.typing_context
        targetctx = self.targetdescr.target_context

        @contextmanager
        def store_overloads_on_success():
            # use to ensure overloads are stored on success
            try:
                yield
            except Exception:
                raise
            else:
                exists = self.overloads.get(cres.signature)
                if exists is None:
                    self.overloads[cres.signature] = cres

        # Use cache and compiler in a critical section
        with global_compiler_lock:
            with targetconfig.ConfigStack().enter(flags.copy()):
                with store_overloads_on_success():
                    # attempt look up of existing
                    cres = self.cache.load_overload(sig, targetctx)
                    if cres is not None:
                        return cres

                    # Compile
                    args, return_type = sigutils.normalize_signature(sig)
                    cres = compiler.compile_extra(typingctx, targetctx,
                                                  self.py_func, args=args,
                                                  return_type=return_type,
                                                  flags=flags, locals=locals)

                    # cache lookup failed before so safe to save
                    self.cache.save_overload(sig, cres)

                    return cres


dispatcher_registry[target_registry['npyufunc']] = UFuncDispatcher


# Utility functions

def _compile_element_wise_function(nb_func, targetoptions, sig):
    # Do compilation
    # Return CompileResult to test
    cres = nb_func.compile(sig, **targetoptions)
    args, return_type = sigutils.normalize_signature(sig)
    return cres, args, return_type


def _finalize_ufunc_signature(cres, args, return_type):
    '''Given a compilation result, argument types, and a return type,
    build a valid Numba signature after validating that it doesn't
    violate the constraints for the compilation mode.
    '''
    if return_type is None:
        if cres.objectmode:
            # Object mode is used and return type is not specified
            raise TypeError("return type must be specified for object mode")
        else:
            return_type = cres.signature.return_type

    assert return_type != types.pyobject
    return return_type(*args)


def _build_element_wise_ufunc_wrapper(cres, signature):
    '''Build a wrapper for the ufunc loop entry point given by the
    compilation result object, using the element-wise signature.
    '''
    ctx = cres.target_context
    library = cres.library
    fname = cres.fndesc.llvm_func_name

    with global_compiler_lock:
        info = build_ufunc_wrapper(library, ctx, fname, signature,
                                   cres.objectmode, cres)
        ptr = info.library.get_pointer_to_function(info.name)
    # Get dtypes
    dtypenums = [as_dtype(a).num for a in signature.args]
    dtypenums.append(as_dtype(signature.return_type).num)
    return dtypenums, ptr, cres.environment


_identities = {
    0: _internal.PyUFunc_Zero,
    1: _internal.PyUFunc_One,
    None: _internal.PyUFunc_None,
    "reorderable": _internal.PyUFunc_ReorderableNone,
}


def parse_identity(identity):
    """
    Parse an identity value and return the corresponding low-level value
    for Numpy.
    """
    try:
        identity = _identities[identity]
    except KeyError:
        raise ValueError("Invalid identity value %r" % (identity,))
    return identity


@contextmanager
def _suppress_deprecation_warning_nopython_not_supplied():
    """This suppresses the NumbaDeprecationWarning that occurs through the use
    of `jit` without the `nopython` kwarg. This use of `jit` occurs in a few
    places in the `{g,}ufunc` mechanism in Numba, predominantly to wrap the
    "kernel" function."""
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                category=NumbaDeprecationWarning,
                                message=(".*The 'nopython' keyword argument "
                                         "was not supplied*"),)
        yield


# Class definitions

class _BaseUFuncBuilder(object):

    def add(self, sig=None):
        if hasattr(self, 'targetoptions'):
            targetoptions = self.targetoptions
        else:
            targetoptions = self.nb_func.targetoptions
        cres, args, return_type = _compile_element_wise_function(
            self.nb_func, targetoptions, sig)
        sig = self._finalize_signature(cres, args, return_type)
        self._sigs.append(sig)
        self._cres[sig] = cres
        return cres

    def disable_compile(self):
        """
        Disable the compilation of new signatures at call time.
        """
        # Override this for implementations that support lazy compilation


class UFuncBuilder(_BaseUFuncBuilder):

    def __init__(self, py_func, identity=None, cache=False, targetoptions=None):
        if targetoptions is None:
            targetoptions = {}
        if is_jitted(py_func):
            py_func = py_func.py_func
        self.py_func = py_func
        self.identity = parse_identity(identity)
        with _suppress_deprecation_warning_nopython_not_supplied():
            self.nb_func = jit(_target='npyufunc',
                               cache=cache,
                               **targetoptions)(py_func)
        self._sigs = []
        self._cres = {}

    def _finalize_signature(self, cres, args, return_type):
        '''Slated for deprecation, use ufuncbuilder._finalize_ufunc_signature()
        instead.
        '''
        return _finalize_ufunc_signature(cres, args, return_type)

    def build_ufunc(self):
        with global_compiler_lock:
            dtypelist = []
            ptrlist = []
            if not self.nb_func:
                raise TypeError("No definition")

            # Get signature in the order they are added
            keepalive = []
            cres = None
            for sig in self._sigs:
                cres = self._cres[sig]
                dtypenums, ptr, env = self.build(cres, sig)
                dtypelist.append(dtypenums)
                ptrlist.append(int(ptr))
                keepalive.append((cres.library, env))

            datlist = [None] * len(ptrlist)

            if cres is None:
                argspec = inspect.getfullargspec(self.py_func)
                inct = len(argspec.args)
            else:
                inct = len(cres.signature.args)
            outct = 1

            # Becareful that fromfunc does not provide full error checking yet.
            # If typenum is out-of-bound, we have nasty memory corruptions.
            # For instance, -1 for typenum will cause segfault.
            # If elements of type-list (2nd arg) is tuple instead,
            # there will also memory corruption. (Seems like code rewrite.)
            ufunc = _internal.fromfunc(
                self.py_func.__name__, self.py_func.__doc__,
                ptrlist, dtypelist, inct, outct, datlist,
                keepalive, self.identity,
            )

            return ufunc

    def build(self, cres, signature):
        '''Slated for deprecation, use
        ufuncbuilder._build_element_wise_ufunc_wrapper().
        '''
        return _build_element_wise_ufunc_wrapper(cres, signature)


class GUFuncBuilder(_BaseUFuncBuilder):

    # TODO handle scalar
    def __init__(self, py_func, signature, identity=None, cache=False,
                 targetoptions=None, writable_args=()):
        if targetoptions is None:
            targetoptions = {}
        self.py_func = py_func
        self.identity = parse_identity(identity)
        with _suppress_deprecation_warning_nopython_not_supplied():
            self.nb_func = jit(_target='npyufunc', cache=cache)(py_func)
        self.signature = signature
        self.sin, self.sout = parse_signature(signature)
        self.targetoptions = targetoptions
        self.cache = cache
        self._sigs = []
        self._cres = {}

        transform_arg = _get_transform_arg(py_func)
        self.writable_args = tuple([transform_arg(a) for a in writable_args])

    def _finalize_signature(self, cres, args, return_type):
        if not cres.objectmode and cres.signature.return_type != types.void:
            raise TypeError("gufunc kernel must have void return type")

        if return_type is None:
            return_type = types.void

        return return_type(*args)

    @global_compiler_lock
    def build_ufunc(self):
        type_list = []
        func_list = []
        if not self.nb_func:
            raise TypeError("No definition")

        # Get signature in the order they are added
        keepalive = []
        for sig in self._sigs:
            cres = self._cres[sig]
            dtypenums, ptr, env = self.build(cres)
            type_list.append(dtypenums)
            func_list.append(int(ptr))
            keepalive.append((cres.library, env))

        datalist = [None] * len(func_list)

        nin = len(self.sin)
        nout = len(self.sout)

        # Pass envs to fromfuncsig to bind to the lifetime of the ufunc object
        ufunc = _internal.fromfunc(
            self.py_func.__name__, self.py_func.__doc__,
            func_list, type_list, nin, nout, datalist,
            keepalive, self.identity, self.signature, self.writable_args
        )
        return ufunc

    def build(self, cres):
        """
        Returns (dtype numbers, function ptr, EnvironmentObject)
        """
        # Builder wrapper for ufunc entry point
        signature = cres.signature
        info = build_gufunc_wrapper(
            self.py_func, cres, self.sin, self.sout,
            cache=self.cache, is_parfors=False,
        )

        env = info.env
        ptr = info.library.get_pointer_to_function(info.name)
        # Get dtypes
        dtypenums = []
        for a in signature.args:
            if isinstance(a, types.Array):
                ty = a.dtype
            else:
                ty = a
            dtypenums.append(as_dtype(ty).num)
        return dtypenums, ptr, env


def _get_transform_arg(py_func):
    """Return function that transform arg into index"""
    args = inspect.getfullargspec(py_func).args
    pos_by_arg = {arg: i for i, arg in enumerate(args)}

    def transform_arg(arg):
        if isinstance(arg, int):
            return arg

        try:
            return pos_by_arg[arg]
        except KeyError:
            msg = (f"Specified writable arg {arg} not found in arg list "
                   f"{args} for function {py_func.__qualname__}")
            raise RuntimeError(msg)

    return transform_arg
