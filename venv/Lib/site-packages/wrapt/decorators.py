"""This module implements decorators for implementing other decorators
as well as some commonly used decorators.

"""

from functools import partial
from inspect import isclass, signature

from .__wrapt__ import (
    BaseObjectProxy,
    BoundFunctionWrapper,
    CallableObjectProxy,
    FunctionWrapper,
)
from .arguments import formatargspec

# Adapter wrapper for the wrapped function which will overlay certain
# properties from the adapter function onto the wrapped function so that
# functions such as inspect.getfullargspec(), inspect.signature() and
# inspect.getsource() return the correct results one would expect.


class _AdapterFunctionCode(CallableObjectProxy):

    def __init__(self, wrapped_code, adapter_code):
        # Explicit class in super() is used because the proxy overrides
        # __class__ and MRO-related methods to delegate to the wrapped
        # object, which can interfere with bare super().
        super(_AdapterFunctionCode, self).__init__(wrapped_code)
        self._self_adapter_code = adapter_code

    @property
    def co_argcount(self):
        return self._self_adapter_code.co_argcount

    @property
    def co_code(self):
        return self._self_adapter_code.co_code

    @property
    def co_flags(self):
        return self._self_adapter_code.co_flags

    @property
    def co_kwonlyargcount(self):
        return self._self_adapter_code.co_kwonlyargcount

    @property
    def co_varnames(self):
        return self._self_adapter_code.co_varnames


class _AdapterFunctionSurrogate(CallableObjectProxy):

    def __init__(self, wrapped, adapter):
        # Explicit class in super() is used because the proxy overrides
        # __class__ and MRO-related methods to delegate to the wrapped
        # object, which can interfere with bare super().
        super(_AdapterFunctionSurrogate, self).__init__(wrapped)
        self._self_adapter = adapter

    @property
    def __code__(self):
        return _AdapterFunctionCode(
            self.__wrapped__.__code__, self._self_adapter.__code__
        )

    @property
    def __defaults__(self):
        return self._self_adapter.__defaults__

    @property
    def __kwdefaults__(self):
        return self._self_adapter.__kwdefaults__

    @property
    def __signature__(self):
        return signature(self._self_adapter)


class _BoundAdapterFunctionWrapper(BoundFunctionWrapper):

    @property
    def __func__(self):
        return _AdapterFunctionSurrogate(
            self.__wrapped__.__func__, self._self_parent._self_adapter
        )

    @property
    def __signature__(self):
        return signature(self._self_parent._self_adapter)


class _AdapterFunctionWrapper(FunctionWrapper):

    __bound_function_wrapper__ = _BoundAdapterFunctionWrapper

    def __init__(self, *args, **kwargs):
        adapter = kwargs.pop("adapter")
        # Explicit class in super() is used because the proxy overrides
        # __class__ and MRO-related methods to delegate to the wrapped
        # object, which can interfere with bare super().
        super(_AdapterFunctionWrapper, self).__init__(*args, **kwargs)
        self._self_surrogate = _AdapterFunctionSurrogate(self.__wrapped__, adapter)
        self._self_adapter = adapter

    @property
    def __code__(self):
        return self._self_surrogate.__code__

    @property
    def __defaults__(self):
        return self._self_surrogate.__defaults__

    @property
    def __kwdefaults__(self):
        return self._self_surrogate.__kwdefaults__

    @property
    def __signature__(self):
        return self._self_surrogate.__signature__


class AdapterFactory:
    def __call__(self, wrapped):
        raise NotImplementedError()


class _DelegatedAdapterFactory(AdapterFactory):
    def __init__(self, factory):
        # Explicit class in super() for consistency with proxy subclasses.
        super(_DelegatedAdapterFactory, self).__init__()
        self.factory = factory

    def __call__(self, wrapped):
        return self.factory(wrapped)


adapter_factory = _DelegatedAdapterFactory

# Decorator for creating other decorators. This decorator and the
# wrappers which they use are designed to properly preserve any name
# attributes, function signatures etc, in addition to the wrappers
# themselves acting like a transparent proxy for the original wrapped
# function so the wrapper is effectively indistinguishable from the
# original wrapped function.


def decorator(wrapper=None, /, *, enabled=None, adapter=None, proxy=FunctionWrapper):
    """
    The decorator should be supplied with a single positional argument
    which is the `wrapper` function to be used to implement the
    decorator. This may be preceded by a step whereby the keyword
    arguments are supplied to customise the behaviour of the
    decorator. The `adapter` argument is used to optionally denote a
    separate function which is notionally used by an adapter
    decorator. In that case parts of the function `__code__` and
    `__defaults__` attributes are used from the adapter function
    rather than those of the wrapped function. This allows for the
    argument specification from `inspect.getfullargspec()` and similar
    functions to be overridden with a prototype for a different
    function than what was wrapped. The `enabled` argument provides a
    way to enable/disable the use of the decorator. If the type of
    `enabled` is a boolean, then it is evaluated immediately and the
    wrapper not even applied if it is `False`. If not a boolean, it will
    be evaluated when the wrapper is called for an unbound wrapper,
    and when binding occurs for a bound wrapper. When being evaluated,
    if `enabled` is callable it will be called to obtain the value to
    be checked. If `False`, the wrapper will not be called and instead
    the original wrapped function will be called directly instead.
    The `proxy` argument provides a way of passing a custom version of
    the `FunctionWrapper` class used in decorating the function.
    """

    if wrapper is not None:
        # Helper function for creating wrapper of the appropriate
        # time when we need it down below.

        def _build(wrapped, wrapper, enabled=None, adapter=None):
            if adapter:
                if isinstance(adapter, AdapterFactory):
                    adapter = adapter(wrapped)

                if not callable(adapter):
                    ns = {}

                    # Check if the signature argument specification has
                    # annotations. If it does then we need to remember
                    # it but also drop it when attempting to manufacture
                    # a standin adapter function. This is necessary else
                    # it will try and look up any types referenced in
                    # the annotations in the empty namespace we use,
                    # which will fail.

                    annotations = {}

                    if not isinstance(adapter, str):
                        if len(adapter) == 7:
                            annotations = adapter[-1]
                            adapter = adapter[:-1]
                        adapter = formatargspec(*adapter)

                    exec(f"def adapter{adapter}: pass", ns, ns)
                    adapter = ns["adapter"]

                    # Override the annotations for the manufactured
                    # adapter function so they match the original
                    # adapter signature argument specification.

                    if annotations:
                        adapter.__annotations__ = annotations

                return _AdapterFunctionWrapper(
                    wrapped=wrapped, wrapper=wrapper, enabled=enabled, adapter=adapter
                )

            return proxy(wrapped=wrapped, wrapper=wrapper, enabled=enabled)

        # The wrapper has been provided so return the final decorator.
        # The decorator is itself one of our function wrappers so we
        # can determine when it is applied to functions, instance methods
        # or class methods. This allows us to bind the instance or class
        # method so the appropriate self or cls attribute is supplied
        # when it is finally called.

        def _wrapper(wrapped, instance, args, kwargs):
            # We first check for the case where the decorator was applied
            # to a class type.
            #
            #     @decorator
            #     class mydecoratorclass:
            #         def __init__(self, arg=None):
            #             self.arg = arg
            #         def __call__(self, wrapped, instance, args, kwargs):
            #             return wrapped(*args, **kwargs)
            #
            #     @mydecoratorclass(arg=1)
            #     def function():
            #         pass
            #
            # In this case an instance of the class is to be used as the
            # decorator wrapper function. If args was empty at this point,
            # then it means that there were optional keyword arguments
            # supplied to be used when creating an instance of the class
            # to be used as the wrapper function.

            if instance is None and isclass(wrapped) and not args:
                # We still need to be passed the target function to be
                # wrapped as yet, so we need to return a further function
                # to be able to capture it.

                def _capture(target_wrapped):
                    # Now have the target function to be wrapped and need
                    # to create an instance of the class which is to act
                    # as the decorator wrapper function. Before we do that,
                    # we need to first check that use of the decorator
                    # hadn't been disabled by a simple boolean. If it was,
                    # the target function to be wrapped is returned instead.

                    _enabled = enabled
                    if type(_enabled) is bool:
                        if not _enabled:
                            return target_wrapped
                        _enabled = None

                    # Now create an instance of the class which is to act
                    # as the decorator wrapper function. Any arguments had
                    # to be supplied as keyword only arguments so that is
                    # all we pass when creating it.

                    target_wrapper = wrapped(**kwargs)

                    # Finally build the wrapper itself and return it.

                    return _build(target_wrapped, target_wrapper, _enabled, adapter)

                return _capture

            # We should always have the target function to be wrapped at
            # this point as the first (and only) value in args.

            target_wrapped = args[0]

            # Need to now check that use of the decorator hadn't been
            # disabled by a simple boolean. If it was, then target
            # function to be wrapped is returned instead.

            _enabled = enabled
            if type(_enabled) is bool:
                if not _enabled:
                    return target_wrapped
                _enabled = None

            # We now need to build the wrapper, but there are a couple of
            # different cases we need to consider.

            if instance is None:
                if isclass(wrapped):
                    # In this case the decorator was applied to a class
                    # type but optional keyword arguments were not supplied
                    # for initialising an instance of the class to be used
                    # as the decorator wrapper function.
                    #
                    #     @decorator
                    #     class mydecoratorclass:
                    #         def __init__(self, arg=None):
                    #             self.arg = arg
                    #         def __call__(self, wrapped, instance,
                    #                 args, kwargs):
                    #             return wrapped(*args, **kwargs)
                    #
                    #     @mydecoratorclass
                    #     def function():
                    #         pass
                    #
                    # We still need to create an instance of the class to
                    # be used as the decorator wrapper function, but no
                    # arguments are pass.

                    target_wrapper = wrapped()

                else:
                    # In this case the decorator was applied to a normal
                    # function, or possibly a static method of a class.
                    #
                    #     @decorator
                    #     def mydecoratorfuntion(wrapped, instance,
                    #             args, kwargs):
                    #         return wrapped(*args, **kwargs)
                    #
                    #     @mydecoratorfunction
                    #     def function():
                    #         pass
                    #
                    # That normal function becomes the decorator wrapper
                    # function.

                    target_wrapper = wrapper

            else:
                if isclass(instance):
                    # In this case the decorator was applied to a class
                    # method.
                    #
                    #     class myclass:
                    #         @decorator
                    #         @classmethod
                    #         def decoratorclassmethod(cls, wrapped,
                    #                 instance, args, kwargs):
                    #             return wrapped(*args, **kwargs)
                    #
                    #     instance = myclass()
                    #
                    #     @instance.decoratorclassmethod
                    #     def function():
                    #         pass
                    #
                    # This one is a bit strange because binding was actually
                    # performed on the wrapper created by our decorator
                    # factory. We need to apply that binding to the decorator
                    # wrapper function that the decorator factory
                    # was applied to.

                    target_wrapper = wrapper.__get__(None, instance)

                else:
                    # In this case the decorator was applied to an instance
                    # method.
                    #
                    #     class myclass:
                    #         @decorator
                    #         def decoratorclassmethod(self, wrapped,
                    #                 instance, args, kwargs):
                    #             return wrapped(*args, **kwargs)
                    #
                    #     instance = myclass()
                    #
                    #     @instance.decoratorclassmethod
                    #     def function():
                    #         pass
                    #
                    # This one is a bit strange because binding was actually
                    # performed on the wrapper created by our decorator
                    # factory. We need to apply that binding to the decorator
                    # wrapper function that the decorator factory
                    # was applied to.

                    target_wrapper = wrapper.__get__(instance, type(instance))

            # Finally build the wrapper itself and return it.

            return _build(target_wrapped, target_wrapper, _enabled, adapter)

        # We first return our magic function wrapper here so we can
        # determine in what context the decorator factory was used. In
        # other words, it is itself a universal decorator. The decorator
        # function is used as the adapter so that linters see a signature
        # corresponding to the decorator and not the wrapper it is being
        # applied to.

        return _build(wrapper, _wrapper, adapter=decorator)

    else:
        # The wrapper still has not been provided, so we are just
        # collecting the optional keyword arguments. Return the
        # decorator again wrapped in a partial using the collected
        # arguments.

        return partial(decorator, enabled=enabled, adapter=adapter, proxy=proxy)


# Descriptor decorator for automatically binding state to a wrapper.
# When applied to a method decorated with function_wrapper or decorator,
# it intercepts descriptor access so that when the method is accessed
# via an instance, the instance is automatically attached to the
# resulting wrapper using __self_setattr__. This is useful when you
# want the wrapper to carry a reference back to the object that owns
# the wrapper factory method, making state accessible through the
# decorated function without manual setup.


class _StateBindingWrapper(BaseObjectProxy):
    """A descriptor decorator that binds the owner instance to wrappers
    produced by a wrapper factory method.

    When applied on top of a method decorated with ``function_wrapper``
    or ``decorator``, it intercepts the descriptor ``__get__`` so that
    when the method is accessed through an instance, the owner instance
    is automatically attached as an attribute on the resulting wrapper
    via ``__self_setattr__``.

    The *name* keyword argument controls the attribute name under which
    the owner instance is stored on the wrapper (default ``"state"``).
    """

    def __init__(self, *, name="state"):
        # Initialise the proxy with a placeholder. The actual wrapper
        # factory is set later when __call__ is invoked as a decorator.

        super().__init__(None)

        self._self_name = name

    def __call__(self, wrapper_factory):
        # Called when used as a decorator, receiving the wrapper factory
        # (the function_wrapper or decorator-decorated method) to wrap.

        self.__wrapped__ = wrapper_factory

        return self

    def __get__(self, instance, owner):
        # When accessed via the class rather than an instance, return
        # the descriptor itself unchanged. This preserves introspection
        # of the underlying wrapper factory via the object proxy.

        if instance is None:
            return self

        # Bind the underlying wrapper factory to the instance so it
        # receives the correct self/cls when called.

        wrapper_factory = self.__wrapped__.__get__(instance, owner)

        name = self._self_name

        def _bind(func):
            wrapper = wrapper_factory(func)

            # Attach the owner instance to the wrapper so it can be
            # accessed as an attribute of the decorated function.

            wrapper.__self_setattr__(name, instance)

            return wrapper

        return _bind


bind_state_to_wrapper = _StateBindingWrapper
