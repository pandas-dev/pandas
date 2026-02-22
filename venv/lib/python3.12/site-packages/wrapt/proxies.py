"""Variants of ObjectProxy for different use cases."""

from collections.abc import Callable
from types import ModuleType

from .__wrapt__ import BaseObjectProxy
from .decorators import synchronized

# Define ObjectProxy which for compatibility adds `__iter__()` support which
# has been removed from `BaseObjectProxy`.


class ObjectProxy(BaseObjectProxy):
    """A generic object proxy which forwards special methods as needed.
    For backwards compatibility this class adds support for `__iter__()`. If
    you don't need backward compatibility for `__iter__()` support then it is
    preferable to use `BaseObjectProxy` directly. If you want automatic
    support for special dunder methods for callables, iterators, and async,
    then use `AutoObjectProxy`."""

    @property
    def __object_proxy__(self):
        return ObjectProxy

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)

    def __iter__(self):
        return iter(self.__wrapped__)


# Define variant of ObjectProxy which can automatically adjust to the wrapped
# object and add special dunder methods.


def __wrapper_call__(*args, **kwargs):
    def _unpack_self(self, *args):
        return self, args

    self, args = _unpack_self(*args)

    return self.__wrapped__(*args, **kwargs)


def __wrapper_iter__(self):
    return iter(self.__wrapped__)


def __wrapper_next__(self):
    return self.__wrapped__.__next__()


def __wrapper_aiter__(self):
    return self.__wrapped__.__aiter__()


async def __wrapper_anext__(self):
    return await self.__wrapped__.__anext__()


def __wrapper_length_hint__(self):
    return self.__wrapped__.__length_hint__()


def __wrapper_await__(self):
    return (yield from self.__wrapped__.__await__())


def __wrapper_get__(self, instance, owner):
    return self.__wrapped__.__get__(instance, owner)


def __wrapper_set__(self, instance, value):
    return self.__wrapped__.__set__(instance, value)


def __wrapper_delete__(self, instance):
    return self.__wrapped__.__delete__(instance)


def __wrapper_set_name__(self, owner, name):
    return self.__wrapped__.__set_name__(owner, name)


class AutoObjectProxy(BaseObjectProxy):
    """An object proxy which can automatically adjust to the wrapped object
    and add special dunder methods as needed. Note that this creates a new
    class for each instance, so it has much higher memory overhead than using
    `BaseObjectProxy` directly. If you know what special dunder methods you need
    then it is preferable to use `BaseObjectProxy` directly and add them to a
    subclass as needed. If you only need `__iter__()` support for backwards
    compatibility then use `ObjectProxy` instead.
    """

    def __new__(cls, wrapped):
        """Injects special dunder methods into a dynamically created subclass
        as needed based on the wrapped object.
        """

        namespace = {}

        wrapped_attrs = dir(wrapped)
        class_attrs = set(dir(cls))

        if callable(wrapped) and "__call__" not in class_attrs:
            namespace["__call__"] = __wrapper_call__

        if "__iter__" in wrapped_attrs and "__iter__" not in class_attrs:
            namespace["__iter__"] = __wrapper_iter__

        if "__next__" in wrapped_attrs and "__next__" not in class_attrs:
            namespace["__next__"] = __wrapper_next__

        if "__aiter__" in wrapped_attrs and "__aiter__" not in class_attrs:
            namespace["__aiter__"] = __wrapper_aiter__

        if "__anext__" in wrapped_attrs and "__anext__" not in class_attrs:
            namespace["__anext__"] = __wrapper_anext__

        if "__length_hint__" in wrapped_attrs and "__length_hint__" not in class_attrs:
            namespace["__length_hint__"] = __wrapper_length_hint__

        # Note that not providing compatibility with generator-based coroutines
        # (PEP 342) here as they are removed in Python 3.11+ and were deprecated
        # in 3.8.

        if "__await__" in wrapped_attrs and "__await__" not in class_attrs:
            namespace["__await__"] = __wrapper_await__

        if "__get__" in wrapped_attrs and "__get__" not in class_attrs:
            namespace["__get__"] = __wrapper_get__

        if "__set__" in wrapped_attrs and "__set__" not in class_attrs:
            namespace["__set__"] = __wrapper_set__

        if "__delete__" in wrapped_attrs and "__delete__" not in class_attrs:
            namespace["__delete__"] = __wrapper_delete__

        if "__set_name__" in wrapped_attrs and "__set_name__" not in class_attrs:
            namespace["__set_name__"] = __wrapper_set_name__

        name = cls.__name__

        if cls is AutoObjectProxy:
            name = BaseObjectProxy.__name__

        return super(AutoObjectProxy, cls).__new__(type(name, (cls,), namespace))

    def __wrapped_setattr_fixups__(self):
        """Adjusts special dunder methods on the class as needed based on the
        wrapped object, when `__wrapped__` is changed.
        """

        cls = type(self)
        class_attrs = set(dir(cls))

        if callable(self.__wrapped__):
            if "__call__" not in class_attrs:
                cls.__call__ = __wrapper_call__
        elif getattr(cls, "__call__", None) is __wrapper_call__:
            delattr(cls, "__call__")

        if hasattr(self.__wrapped__, "__iter__"):
            if "__iter__" not in class_attrs:
                cls.__iter__ = __wrapper_iter__
        elif getattr(cls, "__iter__", None) is __wrapper_iter__:
            delattr(cls, "__iter__")

        if hasattr(self.__wrapped__, "__next__"):
            if "__next__" not in class_attrs:
                cls.__next__ = __wrapper_next__
        elif getattr(cls, "__next__", None) is __wrapper_next__:
            delattr(cls, "__next__")

        if hasattr(self.__wrapped__, "__aiter__"):
            if "__aiter__" not in class_attrs:
                cls.__aiter__ = __wrapper_aiter__
        elif getattr(cls, "__aiter__", None) is __wrapper_aiter__:
            delattr(cls, "__aiter__")

        if hasattr(self.__wrapped__, "__anext__"):
            if "__anext__" not in class_attrs:
                cls.__anext__ = __wrapper_anext__
        elif getattr(cls, "__anext__", None) is __wrapper_anext__:
            delattr(cls, "__anext__")

        if hasattr(self.__wrapped__, "__length_hint__"):
            if "__length_hint__" not in class_attrs:
                cls.__length_hint__ = __wrapper_length_hint__
        elif getattr(cls, "__length_hint__", None) is __wrapper_length_hint__:
            delattr(cls, "__length_hint__")

        if hasattr(self.__wrapped__, "__await__"):
            if "__await__" not in class_attrs:
                cls.__await__ = __wrapper_await__
        elif getattr(cls, "__await__", None) is __wrapper_await__:
            delattr(cls, "__await__")

        if hasattr(self.__wrapped__, "__get__"):
            if "__get__" not in class_attrs:
                cls.__get__ = __wrapper_get__
        elif getattr(cls, "__get__", None) is __wrapper_get__:
            delattr(cls, "__get__")

        if hasattr(self.__wrapped__, "__set__"):
            if "__set__" not in class_attrs:
                cls.__set__ = __wrapper_set__
        elif getattr(cls, "__set__", None) is __wrapper_set__:
            delattr(cls, "__set__")

        if hasattr(self.__wrapped__, "__delete__"):
            if "__delete__" not in class_attrs:
                cls.__delete__ = __wrapper_delete__
        elif getattr(cls, "__delete__", None) is __wrapper_delete__:
            delattr(cls, "__delete__")

        if hasattr(self.__wrapped__, "__set_name__"):
            if "__set_name__" not in class_attrs:
                cls.__set_name__ = __wrapper_set_name__
        elif getattr(cls, "__set_name__", None) is __wrapper_set_name__:
            delattr(cls, "__set_name__")


class LazyObjectProxy(AutoObjectProxy):
    """An object proxy which can generate/create the wrapped object on demand
    when it is first needed.
    """

    def __new__(cls, callback=None, *, interface=...):
        """Injects special dunder methods into a dynamically created subclass
        as needed based on the wrapped object.
        """

        if interface is ...:
            interface = type(None)

        namespace = {}

        interface_attrs = dir(interface)
        class_attrs = set(dir(cls))

        if "__call__" in interface_attrs and "__call__" not in class_attrs:
            namespace["__call__"] = __wrapper_call__

        if "__iter__" in interface_attrs and "__iter__" not in class_attrs:
            namespace["__iter__"] = __wrapper_iter__

        if "__next__" in interface_attrs and "__next__" not in class_attrs:
            namespace["__next__"] = __wrapper_next__

        if "__aiter__" in interface_attrs and "__aiter__" not in class_attrs:
            namespace["__aiter__"] = __wrapper_aiter__

        if "__anext__" in interface_attrs and "__anext__" not in class_attrs:
            namespace["__anext__"] = __wrapper_anext__

        if (
            "__length_hint__" in interface_attrs
            and "__length_hint__" not in class_attrs
        ):
            namespace["__length_hint__"] = __wrapper_length_hint__

        # Note that not providing compatibility with generator-based coroutines
        # (PEP 342) here as they are removed in Python 3.11+ and were deprecated
        # in 3.8.

        if "__await__" in interface_attrs and "__await__" not in class_attrs:
            namespace["__await__"] = __wrapper_await__

        if "__get__" in interface_attrs and "__get__" not in class_attrs:
            namespace["__get__"] = __wrapper_get__

        if "__set__" in interface_attrs and "__set__" not in class_attrs:
            namespace["__set__"] = __wrapper_set__

        if "__delete__" in interface_attrs and "__delete__" not in class_attrs:
            namespace["__delete__"] = __wrapper_delete__

        if "__set_name__" in interface_attrs and "__set_name__" not in class_attrs:
            namespace["__set_name__"] = __wrapper_set_name__

        name = cls.__name__

        return super(AutoObjectProxy, cls).__new__(type(name, (cls,), namespace))

    def __init__(self, callback=None, *, interface=...):
        """Initialize the object proxy with wrapped object as `None` but due
        to presence of special `__wrapped_factory__` attribute addded first,
        this will actually trigger the deferred creation of the wrapped object
        when first needed.
        """

        if callback is not None:
            self.__wrapped_factory__ = callback

        super().__init__(None)

    __wrapped_initialized__ = False

    def __wrapped_factory__(self):
        return None

    def __wrapped_get__(self):
        """Gets the wrapped object, creating it if necessary."""

        # We synchronize on the class type, which will be unique to this instance
        # since we inherit from `AutoObjectProxy` which creates a new class
        # for each instance. If we synchronize on `self` or the method then
        # we can end up in infinite recursion via `__getattr__()`.

        with synchronized(type(self)):
            # We were called because `__wrapped__` was not set, but because of
            # multiple threads we may find that it has been set by the time
            # we get the lock. So check again now whether `__wrapped__` is set.
            # If it is then just return it, otherwise call the factory to
            # create it.

            if self.__wrapped_initialized__:
                return self.__wrapped__

            self.__wrapped__ = self.__wrapped_factory__()

            self.__wrapped_initialized__ = True

            return self.__wrapped__


def lazy_import(name, attribute=None, *, interface=...):
    """Lazily imports the module `name`, returning a `LazyObjectProxy` which
    will import the module when it is first needed. When `name is a dotted name,
    then the full dotted name is imported and the last module is taken as the
    target. If `attribute` is provided then it is used to retrieve an attribute
    from the module.
    """

    if attribute is not None:
        if interface is ...:
            interface = Callable
    else:
        if interface is ...:
            interface = ModuleType

    def _import():
        module = __import__(name, fromlist=[""])

        if attribute is not None:
            return getattr(module, attribute)

        return module

    return LazyObjectProxy(_import, interface=interface)
