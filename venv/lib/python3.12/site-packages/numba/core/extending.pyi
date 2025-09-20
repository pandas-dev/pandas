import collections
import sys
import weakref

from numba.core.serialize import ReduceMixin


def type_callable(func):
    ...


def overload(func, jit_options={}, strict=True, inline='never',
             prefer_literal=False, **kwargs):
    ...


def register_jitable(*args, **kwargs):
    ...


def overload_attribute(typ, attr, **kwargs):
    ...


def _overload_method_common(typ, attr, **kwargs):
    ...


def overload_method(typ, attr, **kwargs):
    ...


def overload_classmethod(typ, attr, **kwargs):
    ...


def make_attribute_wrapper(typeclass, struct_attr, python_attr):
    ...


class _Intrinsic(ReduceMixin, Generic[P, R]):
    _memo: weakref.WeakValueDictionary
    _recent: collections.deque

    def __init__(
        self,
        name,
        defn: Callable[Concatenate[object, P], R],
        prefer_literal=False,
        **kwargs
    ):
        ...

    @property
    def _uuid(self):
        ...

    def _set_uuid(self, u):
        ...

    def _register(self):
        ...

    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> R:
        ...

    def __repr__(self):
        ...

    def __deepcopy__(self, memo):
        ...

    def _reduce_states(self):
        ...

    @classmethod
    def _rebuild(cls, uuid, name, defn):
        ...


if sys.version_info >= (3, 10):
    from typing import Callable, Concatenate, Generic, ParamSpec, TypeVar

    # Type of the parameters of the function being converted to an intrinsic,
    # excluding the typing context parameter.
    P = ParamSpec("P")
    # Return type of the function being converted to an intrinsic
    R = TypeVar("R")

    def intrinsic(*args, **kwargs) -> Callable[
        [Callable[Concatenate[object, P], R]],
        Callable[P, R]
    ]:
        ...
else:
    def intrinsic(*args, **kwargs):
        ...


def get_cython_function_address(module_name, function_name):
    ...


def include_path():
    ...


def sentry_literal_args(pysig, literal_args, args, kwargs):
    ...


class SentryLiteralArgs(collections.namedtuple(
        '_SentryLiteralArgs', ['literal_args'])):
    def for_function(self, func):
        ...

    def for_pysig(self, pysig):
        ...


class BoundLiteralArgs(collections.namedtuple(
        'BoundLiteralArgs', ['pysig', 'literal_args'])):
    def bind(self, *args, **kwargs):
        ...


def is_jitted(function):
    ...
