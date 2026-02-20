import inspect
import sys

from .__wrapt__ import FunctionWrapper

# Helper functions for applying wrappers to existing functions.


def resolve_path(target, name):
    """
    Resolves the dotted path supplied as `name` to an attribute on a target
    object. The `target` can be a module, class, or instance of a class. If the
    `target` argument is a string, it is assumed to be the name of a module,
    which will be imported if necessary and then used as the target object.
    Returns a tuple containing the parent object holding the attribute lookup
    resolved to, the attribute name (path prefix removed if present), and the
    original attribute value.
    """

    if isinstance(target, str):
        __import__(target)
        target = sys.modules[target]

    parent = target

    path = name.split(".")
    attribute = path[0]

    # We can't just always use getattr() because in doing
    # that on a class it will cause binding to occur which
    # will complicate things later and cause some things not
    # to work. For the case of a class we therefore access
    # the __dict__ directly. To cope though with the wrong
    # class being given to us, or a method being moved into
    # a base class, we need to walk the class hierarchy to
    # work out exactly which __dict__ the method was defined
    # in, as accessing it from __dict__ will fail if it was
    # not actually on the class given. Fallback to using
    # getattr() if we can't find it. If it truly doesn't
    # exist, then that will fail.

    def lookup_attribute(parent, attribute):
        if inspect.isclass(parent):
            for cls in inspect.getmro(parent):
                if attribute in vars(cls):
                    return vars(cls)[attribute]
            else:
                return getattr(parent, attribute)
        else:
            return getattr(parent, attribute)

    original = lookup_attribute(parent, attribute)

    for attribute in path[1:]:
        parent = original
        original = lookup_attribute(parent, attribute)

    return (parent, attribute, original)


def apply_patch(parent, attribute, replacement):
    """
    Convenience function for applying a patch to an attribute. Currently this
    maps to the standard setattr() function, but in the future may be extended
    to support more complex patching strategies.
    """

    setattr(parent, attribute, replacement)


def wrap_object(target, name, factory, args=(), kwargs={}):
    """
    Wraps an object which is the attribute of a target object with a wrapper
    object created by the `factory` function. The `target` can be a module,
    class, or instance of a class. In the special case of `target` being a
    string, it is assumed to be the name of a module, with the module being
    imported if necessary and then used as the target object. The `name` is a
    string representing the dotted path to the attribute. The `factory` function
    should accept the original object and may accept additional positional and
    keyword arguments which will be set by unpacking input arguments using
    `*args` and `**kwargs` calling conventions. The factory function should
    return a new object that will replace the original object.
    """

    (parent, attribute, original) = resolve_path(target, name)
    wrapper = factory(original, *args, **kwargs)
    apply_patch(parent, attribute, wrapper)

    return wrapper


# Function for applying a proxy object to an attribute of a class
# instance. The wrapper works by defining an attribute of the same name
# on the class which is a descriptor and which intercepts access to the
# instance attribute. Note that this cannot be used on attributes which
# are themselves defined by a property object.


class AttributeWrapper:

    def __init__(self, attribute, factory, args, kwargs):
        self.attribute = attribute
        self.factory = factory
        self.args = args
        self.kwargs = kwargs

    def __get__(self, instance, owner):
        value = instance.__dict__[self.attribute]
        return self.factory(value, *self.args, **self.kwargs)

    def __set__(self, instance, value):
        instance.__dict__[self.attribute] = value

    def __delete__(self, instance):
        del instance.__dict__[self.attribute]


def wrap_object_attribute(module, name, factory, args=(), kwargs={}):
    """
    Wraps an object which is the attribute of a class instance with a wrapper
    object created by the `factory` function. It does this by patching the
    class, not the instance, with a descriptor that intercepts access to the
    instance attribute. The `module` can be a module, class, or instance of a
    class. In the special case of `module` being a string, it is assumed to be
    the name of a module, with the module being imported if necessary and then
    used as the target object. The `name` is a string representing the dotted
    path to the attribute. The `factory` function should accept the original
    object and may accept additional positional and keyword arguments which will
    be set by unpacking input arguments using `*args` and `**kwargs` calling
    conventions. The factory function should return a new object that will
    replace the original object.
    """

    path, attribute = name.rsplit(".", 1)
    parent = resolve_path(module, path)[2]
    wrapper = AttributeWrapper(attribute, factory, args, kwargs)
    apply_patch(parent, attribute, wrapper)
    return wrapper


# Functions for creating a simple decorator using a FunctionWrapper,
# plus short cut functions for applying wrappers to functions. These are
# for use when doing monkey patching. For a more featured way of
# creating decorators see the decorator decorator instead.


def function_wrapper(wrapper):
    """
    Creates a decorator for wrapping a function with a `wrapper` function.
    The decorator which is returned may also be applied to any other callable
    objects such as lambda functions, methods, classmethods, and staticmethods,
    or objects which implement the `__call__()` method. The `wrapper` function
    should accept the `wrapped` function, `instance`, `args`, and `kwargs`,
    arguments and return the result of calling the wrapped function or some
    other appropriate value.
    """

    def _wrapper(wrapped, instance, args, kwargs):
        target_wrapped = args[0]
        if instance is None:
            target_wrapper = wrapper
        elif inspect.isclass(instance):
            target_wrapper = wrapper.__get__(None, instance)
        else:
            target_wrapper = wrapper.__get__(instance, type(instance))
        return FunctionWrapper(target_wrapped, target_wrapper)

    return FunctionWrapper(wrapper, _wrapper)


def wrap_function_wrapper(target, name, wrapper):
    """
    Wraps a function which is the attribute of a target object with a `wrapper`
    function. The `target` can be a module, class, or instance of a class. In
    the special case of `target` being a string, it is assumed to be the name
    of a module, with the module being imported if necessary. The `name` is a
    string representing the dotted path to the attribute. The `wrapper` function
    should accept the `wrapped` function, `instance`, `args`, and `kwargs`
    arguments, and would return the result of calling the wrapped attribute or
    some other appropriate value.
    """

    return wrap_object(target, name, FunctionWrapper, (wrapper,))


def patch_function_wrapper(target, name, enabled=None):
    """
    Creates a decorator which can be applied to a wrapper function, where the
    wrapper function will be used to wrap a function which is the attribute of
    a target object. The `target` can be a module, class, or instance of a class.
    In the special case of `target` being a string, it is assumed to be the name
    of a module, with the module being imported if necessary. The `name` is a
    string representing the dotted path to the attribute. The `enabled`
    argument can be a boolean or a callable that returns a boolean. When a
    callable is provided, it will be called each time the wrapper is invoked to
    determine if the wrapper function should be executed or whether the wrapped
    function should be called directly. If `enabled` is not provided, the
    wrapper is enabled by default.
    """

    def _wrapper(wrapper):
        return wrap_object(target, name, FunctionWrapper, (wrapper, enabled))

    return _wrapper


def transient_function_wrapper(target, name):
    """Creates a decorator that patches a target function with a wrapper
    function, but only for the duration of the call that the decorator was
    applied to. The `target` can be a module, class, or instance of a class.
    In the special case of `target` being a string, it is assumed to be the name
    of a module, with the module being imported if necessary. The `name` is a
    string representing the dotted path to the attribute.
    """

    def _decorator(wrapper):
        def _wrapper(wrapped, instance, args, kwargs):
            target_wrapped = args[0]
            if instance is None:
                target_wrapper = wrapper
            elif inspect.isclass(instance):
                target_wrapper = wrapper.__get__(None, instance)
            else:
                target_wrapper = wrapper.__get__(instance, type(instance))

            def _execute(wrapped, instance, args, kwargs):
                (parent, attribute, original) = resolve_path(target, name)
                replacement = FunctionWrapper(original, target_wrapper)
                setattr(parent, attribute, replacement)
                try:
                    return wrapped(*args, **kwargs)
                finally:
                    setattr(parent, attribute, original)

            return FunctionWrapper(target_wrapped, _execute)

        return FunctionWrapper(wrapper, _wrapper)

    return _decorator
