import inspect
import typing as py_typing

from numba.core.typing.typeof import typeof
from numba.core import errors, types


class AsNumbaTypeRegistry:
    """
    A registry for Python types. It stores a lookup table for simple cases
    (e.g. ``int``) and a list of functions for more complicated cases (e.g.
    generics like ``List[int]``).

    Python types are used in Python type annotations, and in instance checks.
    Therefore, this registry supports determining the Numba type of Python type
    annotations at compile time, along with determining the type of classinfo
    arguments to ``isinstance()``.

    This registry is not used dynamically on instances at runtime; to check the
    type of an object at runtime, use ``numba.typeof``.
    """

    def __init__(self):
        self.lookup = {
            type(example): typeof(example)
            for example in [
                0,
                0.0,
                complex(0),
                "numba",
                True,
                None,
            ]
        }

        self.functions = [self._builtin_infer, self._numba_type_infer]

    def _numba_type_infer(self, py_type):
        if isinstance(py_type, types.Type):
            return py_type

    def _builtin_infer(self, py_type):
        if not isinstance(py_type, py_typing._GenericAlias):
            return

        if getattr(py_type, "__origin__", None) is py_typing.Union:
            if len(py_type.__args__) != 2:
                raise errors.TypingError(
                    "Cannot type Union of more than two types")

            (arg_1_py, arg_2_py) = py_type.__args__

            if arg_2_py is type(None): # noqa: E721
                return types.Optional(self.infer(arg_1_py))
            elif arg_1_py is type(None): # noqa: E721
                return types.Optional(self.infer(arg_2_py))
            else:
                raise errors.TypingError(
                    "Cannot type Union that is not an Optional "
                    f"(neither type type {arg_2_py} is not NoneType")

        if getattr(py_type, "__origin__", None) is list:
            (element_py,) = py_type.__args__
            return types.ListType(self.infer(element_py))

        if getattr(py_type, "__origin__", None) is dict:
            key_py, value_py = py_type.__args__
            return types.DictType(self.infer(key_py), self.infer(value_py))

        if getattr(py_type, "__origin__", None) is set:
            (element_py,) = py_type.__args__
            return types.Set(self.infer(element_py))

        if getattr(py_type, "__origin__", None) is tuple:
            tys = tuple(map(self.infer, py_type.__args__))
            return types.BaseTuple.from_types(tys)

    def register(self, func_or_py_type, numba_type=None):
        """
        Add support for new Python types (e.g. user-defined JitClasses) to the
        registry. For a simple pair of a Python type and a Numba type, this can
        be called as a function ``register(py_type, numba_type)``. If more
        complex logic is required (e.g. for generic types), ``register`` can be
        used as a decorator for a function that takes a Python type as input
        and returns a Numba type or ``None``.
        """
        if numba_type is not None:
            # register used with a specific (py_type, numba_type) pair.
            assert isinstance(numba_type, types.Type)
            self.lookup[func_or_py_type] = numba_type
        else:
            # register used as a decorator.
            assert inspect.isfunction(func_or_py_type)
            self.functions.append(func_or_py_type)

    def try_infer(self, py_type):
        """
        Try to determine the Numba type of a given Python type. We first
        consider the lookup dictionary. If ``py_type`` is not there, we iterate
        through the registered functions until one returns a Numba type.  If
        type inference fails, return ``None``.
        """
        result = self.lookup.get(py_type, None)

        for func in self.functions:
            if result is not None:
                break
            result = func(py_type)

        if result is not None and not isinstance(result, types.Type):
            raise errors.TypingError(
                f"as_numba_type should return a Numba type, got {result}"
            )
        return result

    def infer(self, py_type):
        result = self.try_infer(py_type)
        if result is None:
            raise errors.TypingError(
                f"Cannot infer Numba type of Python type {py_type}"
            )
        return result

    def __call__(self, py_type):
        return self.infer(py_type)


as_numba_type = AsNumbaTypeRegistry()
