"""Core object proxy and function wrapper implementations."""

import inspect
import operator
import sys
import types


class WrapperNotInitializedError(ValueError):
    """
    Exception raised when a wrapper is in an inconsistent state: __init__ was
    called but __wrapped__ is not set. Inherits from ValueError only, so it is
    not silently swallowed by hasattr/getattr/except AttributeError patterns.
    """

    pass


class _ObjectProxyMethods:

    # We use properties to override the values of __module__ and
    # __doc__. If we add these in ObjectProxy, the derived class
    # __dict__ will still be setup to have string variants of these
    # attributes and the rules of descriptors means that they appear to
    # take precedence over the properties in the base class. To avoid
    # that, we copy the properties into the derived class type itself
    # via a meta class. In that way the properties will always take
    # precedence.
    #
    # Note that because these properties end up in the class __dict__,
    # type-level access (e.g. ObjectProxy.__module__) would return the
    # property object rather than a string, since CPython's
    # type.__module__ getter does a raw dict lookup without invoking the
    # descriptor protocol. The metaclass has its own __module__ and
    # __doc__ properties to handle type-level access correctly.

    @property
    def __module__(self):
        return self.__wrapped__.__module__

    @__module__.setter
    def __module__(self, value):
        self.__wrapped__.__module__ = value

    @property
    def __doc__(self):
        return self.__wrapped__.__doc__

    @__doc__.setter
    def __doc__(self, value):
        self.__wrapped__.__doc__ = value

    # We similar use a property for __dict__. We need __dict__ to be
    # explicit to ensure that vars() works as expected.

    @property
    def __dict__(self):
        return self.__wrapped__.__dict__

    # Need to also propagate the special __weakref__ attribute for case
    # where decorating classes which will define this. If do not define
    # it and use a function like inspect.getmembers() on a decorator
    # class it will fail. This can't be in the derived classes.

    @property
    def __weakref__(self):
        return self.__wrapped__.__weakref__


class _ObjectProxyDictBase:
    """Base class whose sole purpose is to provide a ``getset_descriptor``
    for ``__dict__`` that is valid for all ``ObjectProxy`` subclasses.
    The metaclass installs this descriptor as ``__self_dict__`` so that
    the real instance dictionary of the proxy can always be accessed,
    even though ``ObjectProxy`` replaces ``__dict__`` with a property
    that delegates to the wrapped object."""

    pass


_REAL_DICT_DESCRIPTOR = type.__dict__["__dict__"].__get__(_ObjectProxyDictBase)[
    "__dict__"
]


def _get_self_dict(self):
    return _REAL_DICT_DESCRIPTOR.__get__(self)


# Wrapping the descriptor in a read-only property ensures that
# ``proxy.__self_dict__ = value`` raises AttributeError rather than
# replacing the real instance dictionary (which would break the proxy).

_SELF_DICT_PROPERTY = property(_get_self_dict)


class _ObjectProxyMetaType(type):
    # Properties on the metaclass control type-level access to __module__
    # and __doc__ (e.g. ObjectProxy.__module__). Without these, the
    # instance-level properties copied from _ObjectProxyMethods into each
    # class dict would shadow the string values that type.__new__ sets,
    # causing type.__module__ to return a property object instead of a
    # string. The metaclass properties read from internal keys where the
    # real values are saved.

    @property
    def __module__(cls):
        return cls.__dict__.get("_cls_real_module", "builtins")

    @__module__.setter
    def __module__(cls, value):
        type.__setattr__(cls, "_cls_real_module", value)

    @property
    def __doc__(cls):
        return cls.__dict__.get("_cls_real_doc")

    @__doc__.setter
    def __doc__(cls, value):
        type.__setattr__(cls, "_cls_real_doc", value)

    def __new__(cls, name, bases, dictionary):
        # Copy our special properties into the class so that they
        # always take precedence over attributes of the same name added
        # during construction of a derived class. This is to save
        # duplicating the implementation for them in all derived classes.
        #
        # Because this overwrites the __module__ and __doc__ strings
        # that would normally be in the class dict with property objects,
        # we save the original values first and store them under internal
        # keys. The metaclass properties above read from these keys to
        # ensure type-level access (e.g. MyProxy.__module__) returns a
        # string rather than a property object.

        real_module = dictionary.get("__module__")
        real_doc = dictionary.get("__doc__")

        # If the subclass defines its own __dict__ property, preserve it
        # rather than overwriting it with the default delegating property
        # from _ObjectProxyMethods. When __dict__ is not explicitly
        # defined in the class body, it will not be present in the
        # dictionary at this point.

        custom_dict = dictionary.get("__dict__")

        dictionary.update(vars(_ObjectProxyMethods))

        if custom_dict is not None:
            dictionary["__dict__"] = custom_dict

        dictionary.setdefault("__self_dict__", _SELF_DICT_PROPERTY)

        klass = type.__new__(cls, name, bases, dictionary)

        if real_module is not None:
            type.__setattr__(klass, "_cls_real_module", real_module)
        if real_doc is not None:
            type.__setattr__(klass, "_cls_real_doc", real_doc)

        return klass


class ObjectProxy(_ObjectProxyDictBase, metaclass=_ObjectProxyMetaType):
    """A transparent object proxy that delegates attribute access to a
    wrapped object."""

    __class_getitem__ = classmethod(types.GenericAlias)

    def __init__(self, wrapped):
        """Create an object proxy around the given object."""

        if wrapped is None:
            try:
                callback = object.__getattribute__(self, "__wrapped_factory__")
            except AttributeError:
                callback = None

            if callback is not None:
                # If wrapped is none and class has a __wrapped_factory__
                # method, then we don't set __wrapped__ yet and instead will
                # defer creation of the wrapped object until it is first
                # needed.

                pass

            else:
                object.__setattr__(self, "__wrapped__", wrapped)
        else:
            object.__setattr__(self, "__wrapped__", wrapped)

        object.__setattr__(self, "__init_called__", True)

        # Python 3.2+ has the __qualname__ attribute, but it does not
        # allow it to be overridden using a property and it must instead
        # be an actual string object instead.

        try:
            object.__setattr__(self, "__qualname__", wrapped.__qualname__)
        except AttributeError:
            pass

        # Python 3.10 onwards also does not allow itself to be overridden
        # using a property and it must instead be set explicitly. Python
        # 3.14 onwards uses deferred evaluation of annotations via the
        # __annotate__ attribute, so we copy that instead to avoid
        # triggering eager evaluation which can fail if names referenced
        # in annotations have been shadowed.

        if sys.version_info >= (3, 14):
            try:
                object.__setattr__(self, "__annotate__", wrapped.__annotate__)
            except AttributeError:
                pass
        else:
            try:
                object.__setattr__(self, "__annotations__", wrapped.__annotations__)
            except AttributeError:
                pass

    @property
    def __object_proxy__(self):
        return ObjectProxy

    def __self_setattr__(self, name, value):
        object.__setattr__(self, name, value)

    @property
    def __name__(self):
        return self.__wrapped__.__name__

    @__name__.setter
    def __name__(self, value):
        self.__wrapped__.__name__ = value

    @property
    def __class__(self):
        return self.__wrapped__.__class__

    @__class__.setter
    def __class__(self, value):
        self.__wrapped__.__class__ = value

    def __dir__(self):
        return dir(self.__wrapped__)

    def __str__(self):
        return str(self.__wrapped__)

    def __bytes__(self):
        return bytes(self.__wrapped__)

    def __repr__(self):
        return f"<{type(self).__name__} at 0x{id(self):x} for {type(self.__wrapped__).__name__} at 0x{id(self.__wrapped__):x}>"

    def __format__(self, format_spec):
        return format(self.__wrapped__, format_spec)

    def __reversed__(self):
        return reversed(self.__wrapped__)

    def __round__(self, ndigits=None):
        return round(self.__wrapped__, ndigits)

    def __mro_entries__(self, bases):
        if not isinstance(self.__wrapped__, type) and hasattr(
            self.__wrapped__, "__mro_entries__"
        ):
            return self.__wrapped__.__mro_entries__(bases)
        return (self.__wrapped__,)

    def __lt__(self, other):
        return self.__wrapped__ < other

    def __le__(self, other):
        return self.__wrapped__ <= other

    def __eq__(self, other):
        return self.__wrapped__ == other

    def __ne__(self, other):
        return self.__wrapped__ != other

    def __gt__(self, other):
        return self.__wrapped__ > other

    def __ge__(self, other):
        return self.__wrapped__ >= other

    def __hash__(self):
        return hash(self.__wrapped__)

    def __bool__(self):
        return bool(self.__wrapped__)

    def __setattr__(self, name, value):
        if name.startswith("_self_"):
            object.__setattr__(self, name, value)

        elif name == "__wrapped__":
            object.__setattr__(self, name, value)

            try:
                object.__delattr__(self, "__qualname__")
            except AttributeError:
                pass
            try:
                object.__setattr__(self, "__qualname__", value.__qualname__)
            except AttributeError:
                pass
            if sys.version_info >= (3, 14):
                try:
                    object.__delattr__(self, "__annotate__")
                except AttributeError:
                    pass
                try:
                    object.__setattr__(self, "__annotate__", value.__annotate__)
                except AttributeError:
                    pass
            else:
                try:
                    object.__delattr__(self, "__annotations__")
                except AttributeError:
                    pass
                try:
                    object.__setattr__(self, "__annotations__", value.__annotations__)
                except AttributeError:
                    pass

            __wrapped_setattr_fixups__ = getattr(
                self, "__wrapped_setattr_fixups__", None
            )

            if __wrapped_setattr_fixups__ is not None:
                __wrapped_setattr_fixups__()

        elif name == "__qualname__":
            setattr(self.__wrapped__, name, value)
            object.__setattr__(self, name, value)

        elif name == "__annotations__":
            setattr(self.__wrapped__, name, value)
            object.__setattr__(self, name, value)

        elif name == "__annotate__":
            setattr(self.__wrapped__, name, value)
            object.__setattr__(self, name, value)

        elif hasattr(type(self), name):
            object.__setattr__(self, name, value)

        else:
            setattr(self.__wrapped__, name, value)

    def __getattr__(self, name):
        # If we need to lookup `__wrapped__` then the `__init__()` method
        # cannot have been called, or this is a lazy object proxy which is
        # deferring creation of the wrapped object until it is first needed.

        if name == "__wrapped__":
            # Note that we use existance of `__wrapped_factory__` to gate whether
            # we can attempt to initialize the wrapped object lazily, but it is
            # `__wrapped_get__` that we actually call to do the initialization.
            # This is so that we can handle multithreading correctly by having
            # `__wrapped_get__` use a lock to protect against multiple threads
            # trying to initialize the wrapped object at the same time.

            try:
                object.__getattribute__(self, "__wrapped_factory__")
            except AttributeError:
                pass
            else:
                return object.__getattribute__(self, "__wrapped_get__")()

            # If __init__ was called but __wrapped__ is not set, the wrapper
            # is in an inconsistent state. Raise WrapperNotInitializedError
            # (a ValueError, not AttributeError) so it is not silently
            # swallowed by hasattr/getattr patterns.

            try:
                object.__getattribute__(self, "__init_called__")
            except AttributeError:
                raise AttributeError(
                    f"'{type(self).__name__}' object has no attribute " f"'__wrapped__'"
                )

            raise WrapperNotInitializedError(
                "wrapper is in an inconsistent state: __wrapped__ is not set"
            )

        return getattr(self.__wrapped__, name)

    def __delattr__(self, name):
        if name.startswith("_self_"):
            object.__delattr__(self, name)

        elif name == "__wrapped__":
            raise TypeError("can't delete __wrapped__ attribute")

        elif name == "__qualname__":
            object.__delattr__(self, name)
            delattr(self.__wrapped__, name)

        elif name == "__annotations__":
            try:
                object.__delattr__(self, name)
            except AttributeError:
                pass
            delattr(self.__wrapped__, name)

        elif name == "__annotate__":
            try:
                object.__delattr__(self, name)
            except AttributeError:
                pass
            delattr(self.__wrapped__, name)

        elif hasattr(type(self), name):
            object.__delattr__(self, name)

        else:
            delattr(self.__wrapped__, name)

    def __add__(self, other):
        return self.__wrapped__ + other

    def __sub__(self, other):
        return self.__wrapped__ - other

    def __mul__(self, other):
        return self.__wrapped__ * other

    def __truediv__(self, other):
        return operator.truediv(self.__wrapped__, other)

    def __floordiv__(self, other):
        return self.__wrapped__ // other

    def __mod__(self, other):
        return self.__wrapped__ % other

    def __divmod__(self, other):
        return divmod(self.__wrapped__, other)

    def __pow__(self, other, *args):
        return pow(self.__wrapped__, other, *args)

    def __lshift__(self, other):
        return self.__wrapped__ << other

    def __rshift__(self, other):
        return self.__wrapped__ >> other

    def __and__(self, other):
        return self.__wrapped__ & other

    def __xor__(self, other):
        return self.__wrapped__ ^ other

    def __or__(self, other):
        return self.__wrapped__ | other

    def __radd__(self, other):
        return other + self.__wrapped__

    def __rsub__(self, other):
        return other - self.__wrapped__

    def __rmul__(self, other):
        return other * self.__wrapped__

    def __rtruediv__(self, other):
        return operator.truediv(other, self.__wrapped__)

    def __rfloordiv__(self, other):
        return other // self.__wrapped__

    def __rmod__(self, other):
        return other % self.__wrapped__

    def __rdivmod__(self, other):
        return divmod(other, self.__wrapped__)

    def __rpow__(self, other, *args):
        return pow(other, self.__wrapped__, *args)

    def __rlshift__(self, other):
        return other << self.__wrapped__

    def __rrshift__(self, other):
        return other >> self.__wrapped__

    def __rand__(self, other):
        return other & self.__wrapped__

    def __rxor__(self, other):
        return other ^ self.__wrapped__

    def __ror__(self, other):
        return other | self.__wrapped__

    def __iadd__(self, other):
        if hasattr(self.__wrapped__, "__iadd__"):
            self.__wrapped__ += other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ + other)

    def __isub__(self, other):
        if hasattr(self.__wrapped__, "__isub__"):
            self.__wrapped__ -= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ - other)

    def __imul__(self, other):
        if hasattr(self.__wrapped__, "__imul__"):
            self.__wrapped__ *= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ * other)

    def __itruediv__(self, other):
        if hasattr(self.__wrapped__, "__itruediv__"):
            self.__wrapped__ /= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ / other)

    def __ifloordiv__(self, other):
        if hasattr(self.__wrapped__, "__ifloordiv__"):
            self.__wrapped__ //= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ // other)

    def __imod__(self, other):
        if hasattr(self.__wrapped__, "__imod__"):
            self.__wrapped__ %= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ % other)

    def __ipow__(self, other):  # type: ignore[misc]
        if hasattr(self.__wrapped__, "__ipow__"):
            self.__wrapped__ **= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__**other)

    def __ilshift__(self, other):
        if hasattr(self.__wrapped__, "__ilshift__"):
            self.__wrapped__ <<= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ << other)

    def __irshift__(self, other):
        if hasattr(self.__wrapped__, "__irshift__"):
            self.__wrapped__ >>= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ >> other)

    def __iand__(self, other):
        if hasattr(self.__wrapped__, "__iand__"):
            self.__wrapped__ &= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ & other)

    def __ixor__(self, other):
        if hasattr(self.__wrapped__, "__ixor__"):
            self.__wrapped__ ^= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ ^ other)

    def __ior__(self, other):
        if hasattr(self.__wrapped__, "__ior__"):
            self.__wrapped__ |= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ | other)

    def __neg__(self):
        return -self.__wrapped__

    def __pos__(self):
        return +self.__wrapped__

    def __abs__(self):
        return abs(self.__wrapped__)

    def __invert__(self):
        return ~self.__wrapped__

    def __int__(self):
        return int(self.__wrapped__)

    def __float__(self):
        return float(self.__wrapped__)

    def __complex__(self):
        return complex(self.__wrapped__)

    def __index__(self):
        return operator.index(self.__wrapped__)

    def __matmul__(self, other):
        return self.__wrapped__ @ other

    def __rmatmul__(self, other):
        return other @ self.__wrapped__

    def __imatmul__(self, other):
        if hasattr(self.__wrapped__, "__imatmul__"):
            self.__wrapped__ @= other
            return self
        else:
            return self.__object_proxy__(self.__wrapped__ @ other)

    def __len__(self):
        return len(self.__wrapped__)

    def __contains__(self, value):
        return value in self.__wrapped__

    def __getitem__(self, key):
        return self.__wrapped__[key]

    def __setitem__(self, key, value):
        self.__wrapped__[key] = value

    def __delitem__(self, key):
        del self.__wrapped__[key]

    def __enter__(self):
        return self.__wrapped__.__enter__()

    def __exit__(self, *args, **kwargs):
        return self.__wrapped__.__exit__(*args, **kwargs)

    def __aenter__(self):
        return self.__wrapped__.__aenter__()

    def __aexit__(self, *args, **kwargs):
        return self.__wrapped__.__aexit__(*args, **kwargs)

    def __copy__(self):
        raise NotImplementedError("object proxy must define __copy__()")

    def __deepcopy__(self, memo):
        raise NotImplementedError("object proxy must define __deepcopy__()")

    def __reduce__(self):
        raise NotImplementedError("object proxy must define __reduce__()")

    def __instancecheck__(self, instance):
        return isinstance(instance, self.__wrapped__)

    def __subclasscheck__(self, subclass):
        if hasattr(subclass, "__wrapped__"):
            return issubclass(subclass.__wrapped__, self.__wrapped__)
        else:
            return issubclass(subclass, self.__wrapped__)


class CallableObjectProxy(ObjectProxy):
    """An object proxy for callable objects that also forwards calls."""

    def __call__(*args, **kwargs):
        def _unpack_self(self, *args):
            return self, args

        self, args = _unpack_self(*args)

        return self.__wrapped__(*args, **kwargs)


class PartialCallableObjectProxy(ObjectProxy):
    """A callable object proxy that supports partial application of arguments
    and keywords.
    """

    def __init__(*args, **kwargs):
        """Create a callable object proxy with partial application of the given
        arguments and keywords. This behaves the same as `functools.partial`, but
        implemented using the `ObjectProxy` class to provide better support for
        introspection.
        """

        def _unpack_self(self, *args):
            return self, args

        self, args = _unpack_self(*args)

        if len(args) < 1:
            raise TypeError("partial type takes at least one argument")

        wrapped, args = args[0], args[1:]

        if not callable(wrapped):
            raise TypeError("the first argument must be callable")

        # Explicit class in super() is used because the proxy overrides
        # __class__ and MRO-related methods to delegate to the wrapped
        # object, which can interfere with bare super().
        super(PartialCallableObjectProxy, self).__init__(wrapped)

        self._self_args = args
        self._self_kwargs = kwargs

    def __call__(*args, **kwargs):
        def _unpack_self(self, *args):
            return self, args

        self, args = _unpack_self(*args)

        _args = self._self_args + args

        _kwargs = dict(self._self_kwargs)
        _kwargs.update(kwargs)

        return self.__wrapped__(*_args, **_kwargs)


class _FunctionWrapperBase(ObjectProxy):

    def __init__(
        self,
        wrapped,
        instance,
        wrapper,
        enabled=None,
        binding="callable",
        parent=None,
        owner=None,
    ):

        # Explicit class in super() is used because the proxy overrides
        # __class__ and MRO-related methods to delegate to the wrapped
        # object, which can interfere with bare super().
        super(_FunctionWrapperBase, self).__init__(wrapped)

        object.__setattr__(self, "_self_instance", instance)
        object.__setattr__(self, "_self_wrapper", wrapper)
        object.__setattr__(self, "_self_enabled", enabled)
        object.__setattr__(self, "_self_binding", binding)
        object.__setattr__(self, "_self_parent", parent)
        object.__setattr__(self, "_self_owner", owner)

    def __get__(self, instance, owner):
        # This method handles both unbound and bound derived wrapper classes.
        # It is kept in the base class as the amount of common code makes it
        # impractical to split into the derived classes.
        #
        # The distinguishing attribute which determines whether we are being
        # called in an unbound or bound wrapper is the parent attribute. If
        # binding has never occurred, then the parent will be None.
        #
        # First therefore, is if we are called in an unbound wrapper. In this
        # case we perform the binding.
        #
        # We have two special cases to worry about here. These are where we are
        # decorating a class or builtin function as neither provide a __get__()
        # method to call. In this case we simply return self.
        #
        # Note that we otherwise still do binding even if instance is None and
        # accessing an unbound instance method from a class. This is because we
        # need to be able to later detect that specific case as we will need to
        # extract the instance from the first argument of those passed in.

        if self._self_parent is None:
            # Technically can probably just check for existence of __get__ on
            # the wrapped object, but this is more explicit.

            if self._self_binding == "builtin":
                return self

            if self._self_binding == "class":
                return self

            binder = getattr(self.__wrapped__, "__get__", None)

            if binder is None:
                return self

            descriptor = binder(instance, owner)

            return self.__bound_function_wrapper__(
                descriptor,
                instance,
                self._self_wrapper,
                self._self_enabled,
                self._self_binding,
                self,
                owner,
            )

        # Now we have the case of binding occurring a second time on what was
        # already a bound function. In this case we would usually return
        # ourselves again. This mirrors what Python does.
        #
        # The special case this time is where we were originally bound with an
        # instance of None and we were likely an instance method. In that case
        # we rebind against the original wrapped function from the parent again.

        if self._self_instance is None and self._self_binding in (
            "function",
            "instancemethod",
            "callable",
        ):
            descriptor = self._self_parent.__wrapped__.__get__(instance, owner)

            return self._self_parent.__bound_function_wrapper__(
                descriptor,
                instance,
                self._self_wrapper,
                self._self_enabled,
                self._self_binding,
                self._self_parent,
                owner,
            )

        return self

    def __call__(*args, **kwargs):
        def _unpack_self(self, *args):
            return self, args

        self, args = _unpack_self(*args)

        # If enabled has been specified, then evaluate it at this point
        # and if the wrapper is not to be executed, then simply return
        # the bound function rather than a bound wrapper for the bound
        # function. When evaluating enabled, if it is callable we call
        # it, otherwise we evaluate it as a boolean.

        if self._self_enabled is not None:
            if callable(self._self_enabled):
                if not self._self_enabled():
                    return self.__wrapped__(*args, **kwargs)
            elif not self._self_enabled:
                return self.__wrapped__(*args, **kwargs)

        # This can occur where initial function wrapper was applied to
        # a function that was already bound to an instance. In that case
        # we want to extract the instance from the function and use it.

        if self._self_binding in (
            "function",
            "instancemethod",
            "classmethod",
            "callable",
        ):
            if self._self_instance is None:
                instance = getattr(self.__wrapped__, "__self__", None)
                if instance is not None:
                    return self._self_wrapper(self.__wrapped__, instance, args, kwargs)

        # This is generally invoked when the wrapped function is being
        # called as a normal function and is not bound to a class as an
        # instance method. This is also invoked in the case where the
        # wrapped function was a method, but this wrapper was in turn
        # wrapped using the staticmethod decorator.

        return self._self_wrapper(self.__wrapped__, self._self_instance, args, kwargs)

    def __set_name__(self, owner, name):
        # This is a special method use to supply information to
        # descriptors about what the name of variable in a class
        # definition is. Not wanting to add this to ObjectProxy as not
        # sure of broader implications of doing that. Thus restrict to
        # FunctionWrapper used by decorators.

        if hasattr(self.__wrapped__, "__set_name__"):
            self.__wrapped__.__set_name__(owner, name)


_FUNCTION_WRAPPER_SLOTS = frozenset(
    (
        "_self_instance",
        "_self_wrapper",
        "_self_enabled",
        "_self_binding",
        "_self_parent",
        "_self_owner",
    )
)


class BoundFunctionWrapper(_FunctionWrapperBase):
    """A wrapper for bound methods, classmethods, and staticmethods."""

    def __setattr__(self, name, value):
        if name.startswith("_self_") and name not in _FUNCTION_WRAPPER_SLOTS:
            if self._self_parent is not None:
                object.__setattr__(self._self_parent, name, value)
                return
        super().__setattr__(name, value)

    def __getattr__(self, name):
        if self._self_parent is not None:
            try:
                return getattr(self._self_parent, name)
            except AttributeError:
                pass
        return super().__getattr__(name)

    def __call__(*args, **kwargs):
        def _unpack_self(self, *args):
            return self, args

        self, args = _unpack_self(*args)

        # If enabled has been specified, then evaluate it at this point and if
        # the wrapper is not to be executed, then simply return the bound
        # function rather than a bound wrapper for the bound function. When
        # evaluating enabled, if it is callable we call it, otherwise we
        # evaluate it as a boolean.

        if self._self_enabled is not None:
            if callable(self._self_enabled):
                if not self._self_enabled():
                    return self.__wrapped__(*args, **kwargs)
            elif not self._self_enabled:
                return self.__wrapped__(*args, **kwargs)

        # We need to do things different depending on whether we are likely
        # wrapping an instance method vs a static method or class method.

        if self._self_binding == "function":
            if self._self_instance is None and args:
                instance, newargs = args[0], args[1:]
                if isinstance(instance, self._self_owner):
                    wrapped = PartialCallableObjectProxy(self.__wrapped__, instance)
                    return self._self_wrapper(wrapped, instance, newargs, kwargs)

            return self._self_wrapper(
                self.__wrapped__, self._self_instance, args, kwargs
            )

        elif self._self_binding == "callable":
            if self._self_instance is None and args:
                # This situation can occur where someone is calling the
                # instancemethod via the class type and passing the instance as
                # the first argument. We need to shift the args before making
                # the call to the wrapper and effectively bind the instance to
                # the wrapped function using a partial so the wrapper doesn't
                # see anything as being different.

                instance, newargs = args[0], args[1:]
                if isinstance(instance, self._self_owner):
                    wrapped = PartialCallableObjectProxy(self.__wrapped__, instance)
                    return self._self_wrapper(wrapped, instance, newargs, kwargs)

            return self._self_wrapper(
                self.__wrapped__, self._self_instance, args, kwargs
            )

        else:
            # As in this case we would be dealing with a classmethod or
            # staticmethod, then _self_instance will only tell us whether
            # when calling the classmethod or staticmethod they did it via an
            # instance of the class it is bound to and not the case where
            # done by the class type itself. We thus ignore _self_instance
            # and use the __self__ attribute of the bound function instead.
            # For a classmethod, this means instance will be the class type
            # and for a staticmethod it will be None. This is probably the
            # more useful thing we can pass through even though we loose
            # knowledge of whether they were called on the instance vs the
            # class type, as it reflects what they have available in the
            # decoratored function.

            instance = getattr(self.__wrapped__, "__self__", None)

            return self._self_wrapper(self.__wrapped__, instance, args, kwargs)


class FunctionWrapper(_FunctionWrapperBase):
    """
    A wrapper for callable objects that can be used to apply decorators to
    functions, methods, classmethods, and staticmethods, or any other callable.
    It handles binding and unbinding of methods, and allows for the wrapper to
    be enabled or disabled.
    """

    __bound_function_wrapper__ = BoundFunctionWrapper

    def __init__(self, wrapped, wrapper, enabled=None):
        """
        Initialize the `FunctionWrapper` with the `wrapped` callable, the
        `wrapper` function, and an optional `enabled` argument. The `enabled`
        argument can be a boolean or a callable that returns a boolean. When a
        callable is provided, it will be called each time the wrapper is
        invoked to determine if the wrapper function should be executed or
        whether the wrapped function should be called directly. If `enabled`
        is not provided, the wrapper is enabled by default.
        """

        # What it is we are wrapping here could be anything. We need to
        # try and detect specific cases though. In particular, we need
        # to detect when we are given something that is a method of a
        # class. Further, we need to know when it is likely an instance
        # method, as opposed to a class or static method. This can
        # become problematic though as there isn't strictly a fool proof
        # method of knowing.
        #
        # The situations we could encounter when wrapping a method are:
        #
        # 1. The wrapper is being applied as part of a decorator which
        # is a part of the class definition. In this case what we are
        # given is the raw unbound function, classmethod or staticmethod
        # wrapper objects.
        #
        # The problem here is that we will not know we are being applied
        # in the context of the class being set up. This becomes
        # important later for the case of an instance method, because in
        # that case we just see it as a raw function and can't
        # distinguish it from wrapping a normal function outside of
        # a class context.
        #
        # 2. The wrapper is being applied when performing monkey
        # patching of the class type afterwards and the method to be
        # wrapped was retrieved direct from the __dict__ of the class
        # type. This is effectively the same as (1) above.
        #
        # 3. The wrapper is being applied when performing monkey
        # patching of the class type afterwards and the method to be
        # wrapped was retrieved from the class type. In this case
        # binding will have been performed where the instance against
        # which the method is bound will be None at that point.
        #
        # This case is a problem because we can no longer tell if the
        # method was a static method, plus if using Python3, we cannot
        # tell if it was an instance method as the concept of an
        # unnbound method no longer exists.
        #
        # 4. The wrapper is being applied when performing monkey
        # patching of an instance of a class. In this case binding will
        # have been performed where the instance was not None.
        #
        # This case is a problem because we can no longer tell if the
        # method was a static method.
        #
        # Overall, the best we can do is look at the original type of the
        # object which was wrapped prior to any binding being done and
        # see if it is an instance of classmethod or staticmethod. In
        # the case where other decorators are between us and them, if
        # they do not propagate the __class__  attribute so that the
        # isinstance() checks works, then likely this will do the wrong
        # thing where classmethod and staticmethod are used.
        #
        # Since it is likely to be very rare that anyone even puts
        # decorators around classmethod and staticmethod, likelihood of
        # that being an issue is very small, so we accept it and suggest
        # that those other decorators be fixed. It is also only an issue
        # if a decorator wants to actually do things with the arguments.
        #
        # As to not being able to identify static methods properly, we
        # just hope that that isn't something people are going to want
        # to wrap, or if they do suggest they do it the correct way by
        # ensuring that it is decorated in the class definition itself,
        # or patch it in the __dict__ of the class type.
        #
        # So to get the best outcome we can, whenever we aren't sure what
        # it is, we label it as a 'callable'. If it was already bound and
        # that is rebound later, we assume that it will be an instance
        # method and try and cope with the possibility that the 'self'
        # argument it being passed as an explicit argument and shuffle
        # the arguments around to extract 'self' for use as the instance.

        binding = None

        if isinstance(wrapped, _FunctionWrapperBase):
            binding = wrapped._self_binding

        if not binding:
            if inspect.isbuiltin(wrapped):
                binding = "builtin"

            elif inspect.isfunction(wrapped):
                binding = "function"

            elif inspect.isclass(wrapped):
                binding = "class"

            elif isinstance(wrapped, classmethod):
                binding = "classmethod"

            elif isinstance(wrapped, staticmethod):
                binding = "staticmethod"

            elif hasattr(wrapped, "__self__"):
                if inspect.isclass(wrapped.__self__):
                    binding = "classmethod"
                elif inspect.ismethod(wrapped):
                    binding = "instancemethod"
                else:
                    binding = "callable"

            else:
                binding = "callable"

        # Explicit class in super() is used because the proxy overrides
        # __class__ and MRO-related methods to delegate to the wrapped
        # object, which can interfere with bare super().
        super(FunctionWrapper, self).__init__(wrapped, None, wrapper, enabled, binding)
