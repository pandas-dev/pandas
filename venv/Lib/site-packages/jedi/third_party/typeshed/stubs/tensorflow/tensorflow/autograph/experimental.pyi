from collections.abc import Callable, Iterable
from enum import Enum
from typing import TypeVar, overload
from typing_extensions import ParamSpec

import tensorflow as tf
from tensorflow._aliases import Integer

_Param = ParamSpec("_Param")
_RetType = TypeVar("_RetType")

class Feature(Enum):
    ALL = "ALL"
    ASSERT_STATEMENTS = "ASSERT_STATEMENTS"
    AUTO_CONTROL_DEPS = "AUTO_CONTROL_DEPS"
    BUILTIN_FUNCTIONS = "BUILTIN_FUNCTIONS"
    EQUALITY_OPERATORS = "EQUALITY_OPERATORS"
    LISTS = "LISTS"
    NAME_SCOPES = "NAME_SCOPES"

@overload
def do_not_convert(func: Callable[_Param, _RetType]) -> Callable[_Param, _RetType]: ...
@overload
def do_not_convert(func: None = None) -> Callable[[Callable[_Param, _RetType]], Callable[_Param, _RetType]]: ...
def set_loop_options(
    parallel_iterations: Integer = ...,
    swap_memory: bool = ...,
    maximum_iterations: Integer = ...,
    shape_invariants: Iterable[tuple[tf.Tensor, tf.TensorShape]] = ...,
) -> None: ...
