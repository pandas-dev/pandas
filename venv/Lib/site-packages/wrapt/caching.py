"""Caching decorators. Currently provides ``lru_cache``, a drop-in
replacement for ``functools.lru_cache`` with correct handling of instance
methods: a separate cache is maintained per instance, stored as an
attribute on the instance itself so it is cleaned up with the instance
by the garbage collector. For plain functions, class methods, and static
methods a single shared cache is used, matching ``functools.lru_cache``.
"""

from functools import lru_cache as _functools_lru_cache
from functools import partial

from .__wrapt__ import BaseObjectProxy, BoundFunctionWrapper, FunctionWrapper
from .decorators import decorator
from .synchronization import synchronized

# Decorator that applies functools.lru_cache to the wrapped function.
# Unlike using functools.lru_cache directly, this works correctly with
# instance methods and class methods by maintaining a separate cache
# per instance or class. Each per-instance cache is stored as an
# attribute on the instance itself, so it is automatically cleaned up
# when the instance is garbage collected. For plain functions and
# static methods, a single cache is stored on the wrapper itself.
# The underlying caches are created lazily on first call, with
# double-checked locking via synchronized to avoid races in
# multi-threaded code.
#
# Custom FunctionWrapper and BoundFunctionWrapper subclasses are used
# so that cache_info() and cache_clear() are available directly on
# the decorated function. For bound methods, the BoundFunctionWrapper
# uses _self_instance to locate the per-instance cache.


class _BoundLRUCacheFunctionWrapper(BoundFunctionWrapper):

    def _is_instance_method(self):
        return self._self_binding == "function"

    def __call__(self, *args, **kwargs):
        parent = self._self_parent

        # For class methods and static methods, use a single shared
        # cache stored on the parent wrapper. The cache is created
        # from the bound __wrapped__ since the parent's __wrapped__
        # is the raw classmethod/staticmethod descriptor.

        if not self._is_instance_method():
            if parent._self_cache is None:
                with synchronized(parent):
                    if parent._self_cache is None:
                        parent._self_cache = _functools_lru_cache(
                            **parent._self_lru_kwargs
                        )(self.__wrapped__)

            return parent._self_cache(*args, **kwargs)

        # Instance method — per-instance cache stored as an attribute
        # on the instance so it is cleaned up with the instance by the
        # garbage collector.

        instance = self._self_instance
        cache_attr = parent._self_cache_attr

        cache = getattr(instance, cache_attr, None)

        if cache is None:
            with synchronized(parent):
                cache = getattr(instance, cache_attr, None)

                if cache is None:
                    cache = _functools_lru_cache(**parent._self_lru_kwargs)(
                        self.__wrapped__
                    )

                    # If the instance the method is bound to is a wrapt
                    # object proxy, a plain setattr() would fall through and
                    # store the cache on the wrapped object rather than the
                    # proxy. Use type() rather than isinstance() so the check
                    # sees the real proxy type and is not fooled by the proxy
                    # overriding __class__ to report the wrapped object's
                    # type. Proxies expose __self_setattr__() which stores the
                    # attribute on the proxy itself.

                    if issubclass(type(instance), BaseObjectProxy):
                        instance.__self_setattr__(cache_attr, cache)
                    else:
                        setattr(instance, cache_attr, cache)

        return cache(*args, **kwargs)

    def cache_info(self):
        """Return the cache statistics for this binding's cache, or
        ``None`` if the cache has not yet been created.
        """

        if not self._is_instance_method():
            return self._self_parent.cache_info()

        cache = getattr(self._self_instance, self._self_parent._self_cache_attr, None)

        if cache is not None:
            return cache.cache_info()

        return None

    def cache_clear(self):
        """Clear this binding's cache."""

        if not self._is_instance_method():
            self._self_parent.cache_clear()
            return

        cache = getattr(self._self_instance, self._self_parent._self_cache_attr, None)

        if cache is not None:
            cache.cache_clear()

    def cache_parameters(self):
        """Return the parameters used to create the cache."""

        if not self._is_instance_method():
            return self._self_parent.cache_parameters()

        cache = getattr(self._self_instance, self._self_parent._self_cache_attr, None)

        if cache is not None:
            return cache.cache_parameters()

        return None


class _LRUCacheFunctionWrapper(FunctionWrapper):

    __bound_function_wrapper__ = _BoundLRUCacheFunctionWrapper

    def __init__(self, wrapped, wrapper, **kwargs):
        super().__init__(wrapped, wrapper, **kwargs)

        # Extract the LRU cache configuration that was attached to the
        # wrapper function before it was passed to the decorator.

        self._self_lru_kwargs = wrapper._self_lru_kwargs
        self._self_cache = None

        # Use __func__ to get the name for classmethod/staticmethod
        # descriptors which lack __name__ on Python < 3.10.

        name = getattr(wrapped, "__name__", None)

        if name is None:
            name = wrapped.__func__.__name__

        # The cache attribute name must be unique per decorated method so
        # that a method overridden in a subclass does not share the same
        # per-instance cache slot as the method it overrides. If they shared
        # a slot, a subclass method calling super() would find the subclass
        # cache and re-enter its own body, recursing forever. The owning
        # class is not known here, since the decorator runs on the raw
        # function before the class exists, but each decorated method has its
        # own wrapper instance, so id(self) is a stable discriminator unique
        # to this method definition.

        self._self_cache_attr = "_lru_cache_" + name + "_" + str(id(self))

    def __call__(self, *args, **kwargs):
        # Plain function or static method — single cache stored
        # on the wrapper itself.

        if self._self_cache is None:
            with synchronized(self):
                if self._self_cache is None:
                    self._self_cache = _functools_lru_cache(**self._self_lru_kwargs)(
                        self.__wrapped__
                    )

        return self._self_cache(*args, **kwargs)

    def cache_info(self):
        """Return the cache statistics, or ``None`` if the cache has
        not yet been created.
        """

        if self._self_cache is not None:
            return self._self_cache.cache_info()

        return None

    def cache_clear(self):
        """Clear the cache and reset the statistics."""

        if self._self_cache is not None:
            self._self_cache.cache_clear()

    def cache_parameters(self):
        """Return the parameters used to create the cache, or ``None``
        if the cache has not yet been created.
        """

        if self._self_cache is not None:
            return self._self_cache.cache_parameters()

        return None


def lru_cache(func=None, /, **kwargs):
    """A decorator that applies ``functools.lru_cache`` to the wrapped
    function, with correct handling for instance methods, class methods,
    and static methods.

    For instance methods, a separate ``functools.lru_cache`` is
    maintained per instance. The cache is stored as an attribute on the
    instance itself, so it is automatically cleaned up when the instance
    is garbage collected. This means instances do not need to be
    hashable, each instance gets its own full ``maxsize`` budget, and
    no external mapping prevents garbage collection.

    For plain functions, class methods, and static methods, a single
    shared cache is used.

    All keyword arguments are passed through to ``functools.lru_cache``.

    Cache management methods ``cache_info()`` and ``cache_clear()`` are
    available directly on the decorated function. For bound methods,
    these operate on the per-instance cache for the bound instance.
    """

    if func is None:
        return partial(lru_cache, **kwargs)

    def _wrapper(wrapped, instance, args, _kwargs):
        return wrapped(*args, **_kwargs)

    _wrapper._self_lru_kwargs = kwargs

    return decorator(_wrapper, proxy=_LRUCacheFunctionWrapper)(func)
