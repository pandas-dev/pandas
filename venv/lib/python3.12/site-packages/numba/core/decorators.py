"""
Define @jit and related decorators.
"""


import sys
import warnings
import inspect
import logging
from types import MappingProxyType

from numba.core.errors import DeprecationError, NumbaDeprecationWarning
from numba.stencils.stencil import stencil
from numba.core import config, extending, sigutils, registry

_logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Decorators

_msg_deprecated_signature_arg = ("Deprecated keyword argument `{0}`. "
                                 "Signatures should be passed as the first "
                                 "positional argument.")


def jit(signature_or_function=None, locals=MappingProxyType({}), cache=False,
        pipeline_class=None, boundscheck=None, **options):
    """
    This decorator is used to compile a Python function into native code.

    Args
    -----
    signature_or_function:
        The (optional) signature or list of signatures to be compiled.
        If not passed, required signatures will be compiled when the
        decorated function is called, depending on the argument values.
        As a convenience, you can directly pass the function to be compiled
        instead.

    locals: dict
        Mapping of local variable names to Numba types. Used to override the
        types deduced by Numba's type inference engine.

    pipeline_class: type numba.compiler.CompilerBase
            The compiler pipeline type for customizing the compilation stages.

    options:
        For a cpu target, valid options are:
            nopython: bool
                Set to True to disable the use of PyObjects and Python API
                calls. The default behavior is to allow the use of PyObjects
                and Python API. Default value is True.

            forceobj: bool
                Set to True to force the use of PyObjects for every value.
                Default value is False.

            looplift: bool
                Set to True to enable jitting loops in nopython mode while
                leaving surrounding code in object mode. This allows functions
                to allocate NumPy arrays and use Python objects, while the
                tight loops in the function can still be compiled in nopython
                mode. Any arrays that the tight loop uses should be created
                before the loop is entered. Default value is True.

            error_model: str
                The error-model affects divide-by-zero behavior.
                Valid values are 'python' and 'numpy'. The 'python' model
                raises exception.  The 'numpy' model sets the result to
                *+/-inf* or *nan*. Default value is 'python'.

            inline: str or callable
                The inline option will determine whether a function is inlined
                at into its caller if called. String options are 'never'
                (default) which will never inline, and 'always', which will
                always inline. If a callable is provided it will be called with
                the call expression node that is requesting inlining, the
                caller's IR and callee's IR as arguments, it is expected to
                return Truthy as to whether to inline.
                NOTE: This inlining is performed at the Numba IR level and is in
                no way related to LLVM inlining.

            boundscheck: bool or None
                Set to True to enable bounds checking for array indices. Out
                of bounds accesses will raise IndexError. The default is to
                not do bounds checking. If False, bounds checking is disabled,
                out of bounds accesses can produce garbage results or segfaults.
                However, enabling bounds checking will slow down typical
                functions, so it is recommended to only use this flag for
                debugging. You can also set the NUMBA_BOUNDSCHECK environment
                variable to 0 or 1 to globally override this flag. The default
                value is None, which under normal execution equates to False,
                but if debug is set to True then bounds checking will be
                enabled.

    Returns
    --------
    A callable usable as a compiled function.  Actual compiling will be
    done lazily if no explicit signatures are passed.

    Examples
    --------
    The function can be used in the following ways:

    1) jit(signatures, **targetoptions) -> jit(function)

        Equivalent to:

            d = dispatcher(function, targetoptions)
            for signature in signatures:
                d.compile(signature)

        Create a dispatcher object for a python function.  Then, compile
        the function with the given signature(s).

        Example:

            @jit("int32(int32, int32)")
            def foo(x, y):
                return x + y

            @jit(["int32(int32, int32)", "float32(float32, float32)"])
            def bar(x, y):
                return x + y

    2) jit(function, **targetoptions) -> dispatcher

        Create a dispatcher function object that specializes at call site.

        Examples:

            @jit
            def foo(x, y):
                return x + y

            @jit(nopython=True)
            def bar(x, y):
                return x + y

    """
    locals = dict(locals)
    forceobj = options.get('forceobj', False)
    if 'argtypes' in options:
        raise DeprecationError(_msg_deprecated_signature_arg.format('argtypes'))
    if 'restype' in options:
        raise DeprecationError(_msg_deprecated_signature_arg.format('restype'))
    nopython = options.get('nopython', None)
    if nopython is not None:
        assert type(nopython) is bool, "nopython option must be a bool"
    if nopython is True and forceobj:
        raise ValueError("Only one of 'nopython' or 'forceobj' can be True.")
    target = options.pop('_target', 'cpu')

    if nopython is False:
        msg = ("The keyword argument 'nopython=False' was supplied. From "
               "Numba 0.59.0 the default is True and supplying this argument "
               "has no effect.")
        warnings.warn(msg, NumbaDeprecationWarning)
    # nopython is True by default since 0.59.0, but if `forceobj` is set
    # `nopython` needs to set to False so that things like typing of args in the
    # dispatcher layer continues to work.
    if forceobj:
        options['nopython'] = False
    else:
        options['nopython'] = True

    options['boundscheck'] = boundscheck

    # Handle signature
    if signature_or_function is None:
        # No signature, no function
        pyfunc = None
        sigs = None
    elif isinstance(signature_or_function, list):
        # A list of signatures is passed
        pyfunc = None
        sigs = signature_or_function
    elif sigutils.is_signature(signature_or_function):
        # A single signature is passed
        pyfunc = None
        sigs = [signature_or_function]
    else:
        # A function is passed
        pyfunc = signature_or_function
        sigs = None

    dispatcher_args = {}
    if pipeline_class is not None:
        dispatcher_args['pipeline_class'] = pipeline_class
    wrapper = _jit(sigs, locals=locals, target=target, cache=cache,
                   targetoptions=options, **dispatcher_args)
    if pyfunc is not None:
        return wrapper(pyfunc)
    else:
        return wrapper


def _jit(sigs, locals, target, cache, targetoptions, **dispatcher_args):

    from numba.core.target_extension import resolve_dispatcher_from_str
    dispatcher = resolve_dispatcher_from_str(target)

    def wrapper(func):
        if extending.is_jitted(func):
            raise TypeError(
                "A jit decorator was called on an already jitted function "
                f"{func}.  If trying to access the original python "
                f"function, use the {func}.py_func attribute."
            )

        if not inspect.isfunction(func):
            raise TypeError(
                "The decorated object is not a function (got type "
                f"{type(func)})."
            )

        if config.ENABLE_CUDASIM and target == 'cuda':
            from numba import cuda
            return cuda.jit(func)
        if config.DISABLE_JIT and not target == 'npyufunc':
            return func
        disp = dispatcher(py_func=func, locals=locals,
                          targetoptions=targetoptions,
                          **dispatcher_args)
        if cache:
            disp.enable_caching()
        if sigs is not None:
            # Register the Dispatcher to the type inference mechanism,
            # even though the decorator hasn't returned yet.
            from numba.core import typeinfer
            with typeinfer.register_dispatcher(disp):
                for sig in sigs:
                    disp.compile(sig)
                disp.disable_compile()
        return disp

    return wrapper


def njit(*args, **kws):
    """
    Equivalent to jit(nopython=True)

    See documentation for jit function/decorator for full description.
    """
    if 'nopython' in kws:
        warnings.warn('nopython is set for njit and is ignored', RuntimeWarning)
    if 'forceobj' in kws:
        warnings.warn('forceobj is set for njit and is ignored', RuntimeWarning)
        del kws['forceobj']
    kws.update({'nopython': True})
    return jit(*args, **kws)


def cfunc(sig, locals=MappingProxyType({}), cache=False, pipeline_class=None, **options):
    """
    This decorator is used to compile a Python function into a C callback
    usable with foreign C libraries.

    Usage::
        @cfunc("float64(float64, float64)", nopython=True, cache=True)
        def add(a, b):
            return a + b

    """
    locals = dict(locals)
    sig = sigutils.normalize_signature(sig)

    def wrapper(func):
        from numba.core.ccallback import CFunc
        additional_args = {}
        if pipeline_class is not None:
            additional_args['pipeline_class'] = pipeline_class
        res = CFunc(func, sig, locals=locals, options=options, **additional_args)
        if cache:
            res.enable_caching()
        res.compile()
        return res

    return wrapper


def jit_module(**kwargs):
    """ Automatically ``jit``-wraps functions defined in a Python module

    Note that ``jit_module`` should only be called at the end of the module to
    be jitted. In addition, only functions which are defined in the module
    ``jit_module`` is called from are considered for automatic jit-wrapping.
    See the Numba documentation for more information about what can/cannot be
    jitted.

    :param kwargs: Keyword arguments to pass to ``jit`` such as ``nopython``
                   or ``error_model``.

    """
    # Get the module jit_module is being called from
    frame = inspect.stack()[1]
    module = inspect.getmodule(frame[0])
    # Replace functions in module with jit-wrapped versions
    for name, obj in module.__dict__.items():
        if inspect.isfunction(obj) and inspect.getmodule(obj) == module:
            _logger.debug("Auto decorating function {} from module {} with jit "
                          "and options: {}".format(obj, module.__name__, kwargs))
            module.__dict__[name] = jit(obj, **kwargs)
