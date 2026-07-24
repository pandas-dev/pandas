from _typeshed import Incomplete
from abc import ABCMeta
from builtins import bool as _bool
from typing import Any

import numpy as np
from tensorflow._aliases import DTypeLike
from tensorflow.python.framework.dtypes import HandleData

class _DTypeMeta(ABCMeta): ...

class DType(metaclass=_DTypeMeta):
    __slots__ = ["_handle_data"]
    def __init__(self, type_enum: int, handle_data: HandleData | None = None) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def as_numpy_dtype(self) -> type[np.number[Any]]: ...
    @property
    def is_numpy_compatible(self) -> _bool: ...
    @property
    def is_bool(self) -> _bool: ...
    @property
    def is_floating(self) -> _bool: ...
    @property
    def is_integer(self) -> _bool: ...
    @property
    def is_quantized(self) -> _bool: ...
    @property
    def is_unsigned(self) -> _bool: ...
    def __getattr__(self, name: str) -> Incomplete: ...

bool: DType
complex128: DType
complex64: DType
bfloat16: DType
float16: DType
half: DType
float32: DType
float64: DType
double: DType
int8: DType
int16: DType
int32: DType
int64: DType
uint8: DType
uint16: DType
uint32: DType
uint64: DType
qint8: DType
qint16: DType
qint32: DType
quint8: DType
quint16: DType
string: DType

def as_dtype(type_value: DTypeLike) -> DType: ...
def __getattr__(name: str): ...  # incomplete module
