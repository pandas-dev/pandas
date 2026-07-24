import os
import uuid
import weakref
import collections
import functools
from types import MappingProxyType

import numba
from numba.core import types, errors, utils, config

# Exported symbols
from numba.core.typing.typeof import typeof_impl  # noqa: F401
from numba.core.typing.asnumbatype import as_numba_type  # noqa: F401
from numba.core.typing.templates import infer, infer_getattr  # noqa: F401
from numba.core.imputils import (  # noqa: F401
    lower_builtin, lower_getattr, lower_getattr_generic,  # noqa: F401
    lower_setattr, lower_setattr_generic, lower_cast)  # noqa: F401
from numba.core.datamodel import models   # noqa: F401
from numba.core.datamodel import register_default as register_model  # noqa: F401, E501
from numba.core.pythonapi import box, unbox, reflect, NativeValue  # noqa: F401
from numba._helperlib import _import_cython_function  # noqa: F401
from numba.core.serialize import ReduceMixin


def type_callable(func):
    """
    Decorate a function as implementing typing for the callable *func*.
    *func* can be a callable object (probably a global) or a string
    denoting a built-in operation (such 'getitem' or '__array_wrap__')
    """
    from numba.core.typing.templates import (CallableTemplate, infer,
                                             infer_global)
    if not callable(func) and not isinstance(func, str):
        raise TypeError("`func` should be a function or string")
    try:
        func_name = func.__name__
    except AttributeError:
        func_name = str(func)

    def decorate(typing_func):
        def generic(self):
            return typing_func(self.context)

        name = "%s_CallableTemplate" % (func_name,)
        bases = (CallableTemplate,)
        class_dict = dict(key=func, generic=generic)
        template = type(name, bases, class_dict)
        infer(template)
        if callable(func):
            infer_global(func, types.Function(template))
        return typing_func

    return decorate


# By default, an *overload* does not have a cpython wrapper because it is not
# callable from python. It also has `nopython=True`, this has been default since
# its inception!
_overload_default_jit_options = {'no_cpython_wrapper': True,
                                 'nopython':True}


def overload(
    func,
    jit_options=MappingProxyType({}),
    strict=True,
    inline='never',
    prefer_literal=False,
    **kwargs,
):
    """
    A decorator marking the decorated function as typing and implementing
    *func* in nopython mode.

    The decorated function will have the same formal parameters as *func*
    and be passed the Numba types of those parameters.  It should return
    a function implementing *func* for the given types.

    Here is an example implementing len() for tuple types::

        @overload(len)
        def tuple_len(seq):
            if isinstance(seq, types.BaseTuple):
                n = len(seq)
                def len_impl(seq):
                    return n
                return len_impl

    Compiler options can be passed as an dictionary using the **jit_options**
    argument.

    Overloading strictness (that the typing and implementing signatures match)
    is enforced by the **strict** keyword argument, it is recommended that this
    is set to True (default).

    To handle a function that accepts imprecise types, an overload
    definition can return 2-tuple of ``(signature, impl_function)``, where
    the ``signature`` is a ``typing.Signature`` specifying the precise
    signature to be used; and ``impl_function`` is the same implementation
    function as in the simple case.

    If the kwarg inline determines whether the overload is inlined in the
    calling function and can be one of three values:
    * 'never' (default) - the overload is never inlined.
    * 'always' - the overload is always inlined.
    * a function that takes two arguments, both of which are instances of a
      namedtuple with fields:
        * func_ir
        * typemap
        * calltypes
        * signature
      The first argument holds the information from the caller, the second
      holds the information from the callee. The function should return Truthy
      to determine whether to inline, this essentially permitting custom
      inlining rules (typical use might be cost models).

    The *prefer_literal* option allows users to control if literal types should
    be tried first or last. The default (`False`) is to use non-literal types.
    Implementations that can specialize based on literal values should set the
    option to `True`. Note, this option maybe expanded in the near future to
    allow for more control (e.g. disabling non-literal types).

    **kwargs prescribes additional arguments passed through to the overload
    template. The only accepted key at present is 'target' which is a string
    corresponding to the target that this overload should be bound against.
    """
    from numba.core.typing.templates import make_overload_template, infer_global

    # set default options
    jit_options = dict(jit_options)
    opts = _overload_default_jit_options.copy()
    opts.update(jit_options)  # let user options override

    # TODO: abort now if the kwarg 'target' relates to an unregistered target,
    # this requires sorting out the circular imports first.

    def decorate(overload_func):
        template = make_overload_template(func, overload_func, opts, strict,
                                          inline, prefer_literal, **kwargs)
        infer(template)
        if callable(func):
            infer_global(func, types.Function(template))
        return overload_func

    return decorate


def register_jitable(*args, **kwargs):
    """
    Register a regular python function that can be executed by the python
    interpreter and can be compiled into a nopython function when referenced
    by other jit'ed functions.  Can be used as::

        @register_jitable
        def foo(x, y):
            return x + y

    Or, with compiler options::

        @register_jitable(_nrt=False) # disable runtime allocation
        def foo(x, y):
            return x + y

    """
    def wrap(fn):
        # It is just a wrapper for @overload
        inline = kwargs.pop('inline', 'never')

        @overload(fn, jit_options=kwargs, inline=inline, strict=False)
        def ov_wrap(*args, **kwargs):
            return fn
        return fn

    if kwargs:
        return wrap
    else:
        return wrap(*args)


def overload_attribute(typ, attr, **kwargs):
    """
    A decorator marking the decorated function as typing and implementing
    attribute *attr* for the given Numba type in nopython mode.

    *kwargs* are passed to the underlying `@overload` call.

    Here is an example implementing .nbytes for array types::

        @overload_attribute(types.Array, 'nbytes')
        def array_nbytes(arr):
            def get(arr):
                return arr.size * arr.itemsize
            return get
    """
    # TODO implement setters
    from numba.core.typing.templates import make_overload_attribute_template

    def decorate(overload_func):
        template = make_overload_attribute_template(
            typ, attr, overload_func,
            **kwargs
        )
        infer_getattr(template)
        overload(overload_func, **kwargs)(overload_func)
        return overload_func

    return decorate


def _overload_method_common(typ, attr, **kwargs):
    """Common code for overload_method and overload_classmethod
    """
    from numba.core.typing.templates import make_overload_method_template

    def decorate(overload_func):
        copied_kwargs = kwargs.copy() # avoid mutating parent dict
        template = make_overload_method_template(
            typ, attr, overload_func,
            inline=copied_kwargs.pop('inline', 'never'),
            prefer_literal=copied_kwargs.pop('prefer_literal', False),
            **copied_kwargs,
        )
        infer_getattr(template)
        overload(overload_func, **kwargs)(overload_func)
        return overload_func

    return decorate


def overload_method(typ, attr, **kwargs):
    """
    A decorator marking the decorated function as typing and implementing
    method *attr* for the given Numba type in nopython mode.

    *kwargs* are passed to the underlying `@overload` call.

    Here is an example implementing .take() for array types::

        @overload_method(types.Array, 'take')
        def array_take(arr, indices):
            if isinstance(indices, types.Array):
                def take_impl(arr, indices):
                    n = indices.shape[0]
                    res = np.empty(n, arr.dtype)
                    for i in range(n):
                        res[i] = arr[indices[i]]
                    return res
                return take_impl
    """
    return _overload_method_common(typ, attr, **kwargs)


def overload_classmethod(typ, attr, **kwargs):
    """
    A decorator marking the decorated function as typing and implementing
    classmethod *attr* for the given Numba type in nopython mode.


    Similar to ``overload_method``.


    Here is an example implementing a classmethod on the Array type to call
    ``np.arange()``::

        @overload_classmethod(types.Array, "make")
        def ov_make(cls, nitems):
            def impl(cls, nitems):
                return np.arange(nitems)
            return impl

    The above code will allow the following to work in jit-compiled code::

        @njit
        def foo(n):
            return types.Array.make(n)
    """
    return _overload_method_common(types.TypeRef(typ), attr, **kwargs)


def make_attribute_wrapper(typeclass, struct_attr, python_attr):
    """
    Make an automatic attribute wrapper exposing member named *struct_attr*
    as a read-only attribute named *python_attr*.
    The given *typeclass*'s model must be a StructModel subclass.
    """
    from numba.core.typing.templates import AttributeTemplate
    from numba.core.datamodel import default_manager
    from numba.core.datamodel.models import StructModel
    from numba.core.imputils import impl_ret_borrowed
    from numba.core import cgutils

    if not isinstance(typeclass, type) or not issubclass(typeclass, types.Type):
        raise TypeError("typeclass should be a Type subclass, got %s"
                        % (typeclass,))

    def get_attr_fe_type(typ):
        """
        Get the Numba type of member *struct_attr* in *typ*.
        """
        model = default_manager.lookup(typ)
        if not isinstance(model, StructModel):
            raise TypeError("make_struct_attribute_wrapper() needs a type "
                            "with a StructModel, but got %s" % (model,))
        return model.get_member_fe_type(struct_attr)

    @infer_getattr
    class StructAttribute(AttributeTemplate):
        key = typeclass

        def generic_resolve(self, typ, attr):
            if attr == python_attr:
                return get_attr_fe_type(typ)

    @lower_getattr(typeclass, python_attr)
    def struct_getattr_impl(context, builder, typ, val):
        val = cgutils.create_struct_proxy(typ)(context, builder, value=val)
        attrty = get_attr_fe_type(typ)
        attrval = getattr(val, struct_attr)
        return impl_ret_borrowed(context, builder, attrty, attrval)


class _Intrinsic(ReduceMixin):
    """
    Dummy callable for intrinsic
    """
    _memo: weakref.WeakValueDictionary = weakref.WeakValueDictionary()
    __cache_size = config.FUNCTION_CACHE_SIZE # type: ignore
    # hold refs to last N functions deserialized, retaining them in _memo
    # regardless of whether there is another reference
    _recent: collections.deque = collections.deque(maxlen=__cache_size)

    __uuid = None

    def __init__(self, name, defn, prefer_literal=False, **kwargs):
        self._ctor_kwargs = kwargs
        self._name = name
        self._defn = defn
        self._prefer_literal = prefer_literal
        functools.update_wrapper(self, defn)

    @property
    def _uuid(self):
        """
        An instance-specific UUID, to avoid multiple deserializations of
        a given instance.

        Note this is lazily-generated, for performance reasons.
        """
        u = self.__uuid
        if u is None:
            u = str(uuid.uuid1())
            self._set_uuid(u)
        return u

    def _set_uuid(self, u):
        assert self.__uuid is None
        self.__uuid = u
        self._memo[u] = self
        self._recent.append(self)

    def _register(self):
        # _ctor_kwargs
        from numba.core.typing.templates import (make_intrinsic_template,
                                                 infer_global)

        template = make_intrinsic_template(self, self._defn, self._name,
                                           prefer_literal=self._prefer_literal,
                                           kwargs=self._ctor_kwargs)
        infer(template)
        infer_global(self, types.Function(template))

    def __call__(self, *args, **kwargs):
        """
        This is only defined to pretend to be a callable from CPython.
        """
        msg = '{0} is not usable in pure-python'.format(self)
        raise NotImplementedError(msg)

    def __repr__(self):
        return "<intrinsic {0}>".format(self._name)

    def __deepcopy__(self, memo):
        # NOTE: Intrinsic are immutable and we don't need to copy.
        #       This is triggered from deepcopy of statements.
        return self

    def _reduce_states(self):
        """
        NOTE: part of ReduceMixin protocol
        """
        return dict(uuid=self._uuid, name=self._name, defn=self._defn)

    @classmethod
    def _rebuild(cls, uuid, name, defn):
        """
        NOTE: part of ReduceMixin protocol
        """
        try:
            return cls._memo[uuid]
        except KeyError:
            llc = cls(name=name, defn=defn)
            llc._register()
            llc._set_uuid(uuid)
            return llc


def intrinsic(*args, **kwargs):
    """
    A decorator marking the decorated function as typing and implementing
    *func* in nopython mode using the llvmlite IRBuilder API.  This is an escape
    hatch for expert users to build custom LLVM IR that will be inlined to
    the caller.

    The first argument to *func* is the typing context.  The rest of the
    arguments corresponds to the type of arguments of the decorated function.
    These arguments are also used as the formal argument of the decorated
    function.  If *func* has the signature ``foo(typing_context, arg0, arg1)``,
    the decorated function will have the signature ``foo(arg0, arg1)``.

    The return values of *func* should be a 2-tuple of expected type signature,
    and a code-generation function that will passed to ``lower_builtin``.
    For unsupported operation, return None.

    Here is an example implementing a ``cast_int_to_byte_ptr`` that cast
    any integer to a byte pointer::

        @intrinsic
        def cast_int_to_byte_ptr(typingctx, src):
            # check for accepted types
            if isinstance(src, types.Integer):
                # create the expected type signature
                result_type = types.CPointer(types.uint8)
                sig = result_type(types.uintp)
                # defines the custom code generation
                def codegen(context, builder, signature, args):
                    # llvm IRBuilder code here
                    [src] = args
                    rtype = signature.return_type
                    llrtype = context.get_value_type(rtype)
                    return builder.inttoptr(src, llrtype)
                return sig, codegen
    """
    # Make inner function for the actual work
    def _intrinsic(func):
        name = getattr(func, '__name__', str(func))
        llc = _Intrinsic(name, func, **kwargs)
        llc._register()
        return llc

    if not kwargs:
        # No option is given
        return _intrinsic(*args)
    else:
        # options are given, create a new callable to recv the
        # definition function
        def wrapper(func):
            return _intrinsic(func)
        return wrapper


def get_cython_function_address(module_name, function_name):
    """
    Get the address of a Cython function.

    Args
    ----
    module_name:
        Name of the Cython module
    function_name:
        Name of the Cython function

    Returns
    -------
    A Python int containing the address of the function

    """
    return _import_cython_function(module_name, function_name)


def include_path():
    """Returns the C include directory path.
    """
    include_dir = os.path.dirname(os.path.dirname(numba.__file__))
    path = os.path.abspath(include_dir)
    return path


def sentry_literal_args(pysig, literal_args, args, kwargs):
    """Ensures that the given argument types (in *args* and *kwargs*) are
    literally typed for a function with the python signature *pysig* and the
    list of literal argument names in *literal_args*.

    Alternatively, this is the same as::

        SentryLiteralArgs(literal_args).for_pysig(pysig).bind(*args, **kwargs)
    """
    boundargs = pysig.bind(*args, **kwargs)

    # Find literal argument positions and whether it is satisfied.
    request_pos = set()
    missing = False
    for i, (k, v) in enumerate(boundargs.arguments.items()):
        if k in literal_args:
            request_pos.add(i)
            if not isinstance(v, types.Literal):
                missing = True
    if missing:
        # Yes, there are missing required literal arguments
        e = errors.ForceLiteralArg(request_pos)

        # A helper function to fold arguments
        def folded(args, kwargs):
            out = pysig.bind(*args, **kwargs).arguments.values()
            return tuple(out)

        raise e.bind_fold_arguments(folded)


class SentryLiteralArgs(collections.namedtuple(
        '_SentryLiteralArgs', ['literal_args'])):
    """
    Parameters
    ----------
    literal_args : Sequence[str]
        A sequence of names for literal arguments

    Examples
    --------

    The following line:

    >>> SentryLiteralArgs(literal_args).for_pysig(pysig).bind(*args, **kwargs)

    is equivalent to:

    >>> sentry_literal_args(pysig, literal_args, args, kwargs)
    """
    def for_function(self, func):
        """Bind the sentry to the signature of *func*.

        Parameters
        ----------
        func : Function
            A python function.

        Returns
        -------
        obj : BoundLiteralArgs
        """
        return self.for_pysig(utils.pysignature(func))

    def for_pysig(self, pysig):
        """Bind the sentry to the given signature *pysig*.

        Parameters
        ----------
        pysig : inspect.Signature


        Returns
        -------
        obj : BoundLiteralArgs
        """
        return BoundLiteralArgs(
            pysig=pysig,
            literal_args=self.literal_args,
        )


class BoundLiteralArgs(collections.namedtuple(
        'BoundLiteralArgs', ['pysig', 'literal_args'])):
    """
    This class is usually created by SentryLiteralArgs.
    """
    def bind(self, *args, **kwargs):
        """Bind to argument types.
        """
        return sentry_literal_args(
            self.pysig,
            self.literal_args,
            args,
            kwargs,
        )


def is_jitted(function):
    """Returns True if a function is wrapped by one of the Numba @jit
    decorators, for example: numba.jit, numba.njit

    The purpose of this function is to provide a means to check if a function is
    already JIT decorated.
    """

    # don't want to export this so import locally
    from numba.core.dispatcher import Dispatcher
    return isinstance(function, Dispatcher)
