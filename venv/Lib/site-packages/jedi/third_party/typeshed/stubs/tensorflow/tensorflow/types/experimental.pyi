import abc
from typing import Any, Generic, TypeVar, overload
from typing_extensions import ParamSpec

import tensorflow as tf
from tensorflow._aliases import ContainerGeneric

_P = ParamSpec("_P")
_R_co = TypeVar("_R_co", covariant=True)

class Callable(Generic[_P, _R_co], metaclass=abc.ABCMeta):
    def __call__(self, *args: _P.args, **kwargs: _P.kwargs) -> _R_co: ...

class ConcreteFunction(Callable[_P, _R_co], metaclass=abc.ABCMeta):
    def __call__(self, *args: _P.args, **kwargs: _P.kwargs) -> _R_co: ...

class PolymorphicFunction(Callable[_P, _R_co], metaclass=abc.ABCMeta):
    @overload
    @abc.abstractmethod
    def get_concrete_function(self, *args: _P.args, **kwargs: _P.kwargs) -> ConcreteFunction[_P, _R_co]: ...
    @overload
    @abc.abstractmethod
    def get_concrete_function(
        self, *args: ContainerGeneric[tf.TypeSpec[Any]], **kwargs: ContainerGeneric[tf.TypeSpec[Any]]
    ) -> ConcreteFunction[_P, _R_co]: ...
    def experimental_get_compiler_ir(self, *args, **kwargs): ...

GenericFunction = PolymorphicFunction

def __getattr__(name: str): ...  # incomplete module
