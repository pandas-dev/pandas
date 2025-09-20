from __future__ import annotations

import sys
from typing import Any, Protocol, TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeVar, runtime_checkable
else:
    from typing_extensions import TypeAliasType, TypeVar, runtime_checkable

import numpy as np

from optype._utils import set_module

__all__ = ["DType", "HasDType", "ToDType"]


def __dir__() -> list[str]:
    return __all__


###

ST = TypeVar("ST", bound=np.generic, default=Any)
DT_co = TypeVar("DT_co", bound=np.dtype[Any], default=np.dtype[Any], covariant=True)


@runtime_checkable
@set_module("optype.numpy")
class HasDType(Protocol[DT_co]):
    """HasDType[DT: np.dtype[np.generic] = np.dtype[np.generic]]

    Runtime checkable protocol for objects (or types) that have a `dtype`
    attribute (or property), such as `numpy.ndarray` instances, or
    `numpy.generic` "scalar" instances.

    Anything that implements this interface can be used with the `numpy.dtype`
    constructor, i.e. its constructor is compatible with a signature that
    looks something like `(HasDType[DT: numpy.DType], ...) -> DT`.
    """

    @property
    def dtype(self, /) -> DT_co: ...


DType: TypeAlias = np.dtype[ST]
"""Alias for `numpy.dtype[T: numpy.generic = np.generic]`."""


ToDType = TypeAliasType(
    "ToDType",
    type[ST] | np.dtype[ST] | HasDType[np.dtype[ST]],
    type_params=(ST,),
)
"""
Accepts values that that will result in a `np.dtype[ST]` instance when passed to the
`np.dtype()` constructor.

Note that this does not accept builtin types like `float`, strings like `'f8'`, or any
other valid value that isn't directly expressible in terms a numpy scalar type.

But even though `ToDType` accepts anything with a `.dtype` attribute, passing a
`np.ndarray` instance to `np.dtype()` will raise a `TypeError`. But currently it's
impossible to exclude just `np.ndarray` from this `ToDType` type alias.
"""
