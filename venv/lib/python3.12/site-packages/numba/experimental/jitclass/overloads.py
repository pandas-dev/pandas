"""
Overloads for ClassInstanceType for built-in functions that call dunder methods
on an object.
"""
from functools import wraps
import inspect
import operator

from numba.core.extending import overload
from numba.core.types import ClassInstanceType


def _get_args(n_args):
    assert n_args in (1, 2)
    return list("xy")[:n_args]


def class_instance_overload(target):
    """
    Decorator to add an overload for target that applies when the first argument
    is a ClassInstanceType.
    """
    def decorator(func):
        @wraps(func)
        def wrapped(*args, **kwargs):
            if not isinstance(args[0], ClassInstanceType):
                return
            return func(*args, **kwargs)

        if target is not complex:
            # complex ctor needs special treatment as it uses kwargs
            params = list(inspect.signature(wrapped).parameters)
            assert params == _get_args(len(params))
        return overload(target)(wrapped)

    return decorator


def extract_template(template, name):
    """
    Extract a code-generated function from a string template.
    """
    namespace = {}
    exec(template, namespace)
    return namespace[name]


def register_simple_overload(func, *attrs, n_args=1,):
    """
    Register an overload for func that checks for methods __attr__ for each
    attr in attrs.
    """
    # Use a template to set the signature correctly.
    arg_names = _get_args(n_args)
    template = f"""
def func({','.join(arg_names)}):
    pass
"""

    @wraps(extract_template(template, "func"))
    def overload_func(*args, **kwargs):
        options = [
            try_call_method(args[0], f"__{attr}__", n_args)
            for attr in attrs
        ]
        return take_first(*options)

    return class_instance_overload(func)(overload_func)


def try_call_method(cls_type, method, n_args=1):
    """
    If method is defined for cls_type, return a callable that calls this method.
    If not, return None.
    """
    if method in cls_type.jit_methods:
        arg_names = _get_args(n_args)
        template = f"""
def func({','.join(arg_names)}):
    return {arg_names[0]}.{method}({','.join(arg_names[1:])})
"""
        return extract_template(template, "func")


def try_call_complex_method(cls_type, method):
    """ __complex__ needs special treatment as the argument names are kwargs
    and therefore specific in name and default value.
    """
    if method in cls_type.jit_methods:
        template = f"""
def func(real=0, imag=0):
    return real.{method}()
"""
        return extract_template(template, "func")


def take_first(*options):
    """
    Take the first non-None option.
    """
    assert all(o is None or inspect.isfunction(o) for o in options), options
    for o in options:
        if o is not None:
            return o


@class_instance_overload(bool)
def class_bool(x):
    using_bool_impl = try_call_method(x, "__bool__")

    if '__len__' in x.jit_methods:
        def using_len_impl(x):
            return bool(len(x))
    else:
        using_len_impl = None

    always_true_impl = lambda x: True

    return take_first(using_bool_impl, using_len_impl, always_true_impl)


@class_instance_overload(complex)
def class_complex(real=0, imag=0):
    return take_first(
        try_call_complex_method(real, "__complex__"),
        lambda real=0, imag=0: complex(float(real))
    )


@class_instance_overload(operator.contains)
def class_contains(x, y):
    # https://docs.python.org/3/reference/expressions.html#membership-test-operations
    return try_call_method(x, "__contains__", 2)
    # TODO: use __iter__ if defined.


@class_instance_overload(float)
def class_float(x):
    options = [try_call_method(x, "__float__")]

    if (
        "__index__" in x.jit_methods
    ):
        options.append(lambda x: float(x.__index__()))

    return take_first(*options)


@class_instance_overload(int)
def class_int(x):
    options = [try_call_method(x, "__int__")]

    options.append(try_call_method(x, "__index__"))

    return take_first(*options)


@class_instance_overload(str)
def class_str(x):
    return take_first(
        try_call_method(x, "__str__"),
        lambda x: repr(x),
    )


@class_instance_overload(operator.ne)
def class_ne(x, y):
    # This doesn't use register_reflected_overload like the other operators
    # because it falls back to inverting __eq__ rather than reflecting its
    # arguments (as per the definition of the Python data model).
    return take_first(
        try_call_method(x, "__ne__", 2),
        lambda x, y: not (x == y),
    )


def register_reflected_overload(func, meth_forward, meth_reflected):
    def class_lt(x, y):
        normal_impl = try_call_method(x, f"__{meth_forward}__", 2)

        if f"__{meth_reflected}__" in y.jit_methods:
            def reflected_impl(x, y):
                return y > x
        else:
            reflected_impl = None

        return take_first(normal_impl, reflected_impl)

    class_instance_overload(func)(class_lt)


register_simple_overload(abs, "abs")
register_simple_overload(len, "len")
register_simple_overload(hash, "hash")

# Comparison operators.
register_reflected_overload(operator.ge, "ge", "le")
register_reflected_overload(operator.gt, "gt", "lt")
register_reflected_overload(operator.le, "le", "ge")
register_reflected_overload(operator.lt, "lt", "gt")

# Note that eq is missing support for fallback to `x is y`, but `is` and
# `operator.is` are presently unsupported in general.
register_reflected_overload(operator.eq, "eq", "eq")

# Arithmetic operators.
register_simple_overload(operator.add, "add", n_args=2)
register_simple_overload(operator.floordiv, "floordiv", n_args=2)
register_simple_overload(operator.lshift, "lshift", n_args=2)
register_simple_overload(operator.mul, "mul", n_args=2)
register_simple_overload(operator.mod, "mod", n_args=2)
register_simple_overload(operator.neg, "neg")
register_simple_overload(operator.pos, "pos")
register_simple_overload(operator.invert, "invert")
register_simple_overload(operator.pow, "pow", n_args=2)
register_simple_overload(operator.rshift, "rshift", n_args=2)
register_simple_overload(operator.sub, "sub", n_args=2)
register_simple_overload(operator.truediv, "truediv", n_args=2)

# Inplace arithmetic operators.
register_simple_overload(operator.iadd, "iadd", "add", n_args=2)
register_simple_overload(operator.ifloordiv, "ifloordiv", "floordiv", n_args=2)
register_simple_overload(operator.ilshift, "ilshift", "lshift", n_args=2)
register_simple_overload(operator.imul, "imul", "mul", n_args=2)
register_simple_overload(operator.imod, "imod", "mod", n_args=2)
register_simple_overload(operator.ipow, "ipow", "pow", n_args=2)
register_simple_overload(operator.irshift, "irshift", "rshift", n_args=2)
register_simple_overload(operator.isub, "isub", "sub", n_args=2)
register_simple_overload(operator.itruediv, "itruediv", "truediv", n_args=2)

# Logical operators.
register_simple_overload(operator.and_, "and", n_args=2)
register_simple_overload(operator.or_, "or", n_args=2)
register_simple_overload(operator.xor, "xor", n_args=2)

# Inplace logical operators.
register_simple_overload(operator.iand, "iand", "and", n_args=2)
register_simple_overload(operator.ior, "ior", "or", n_args=2)
register_simple_overload(operator.ixor, "ixor", "xor", n_args=2)
