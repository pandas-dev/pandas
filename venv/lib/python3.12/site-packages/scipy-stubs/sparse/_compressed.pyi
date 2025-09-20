import abc
from typing import Any, Generic, Literal, Self, SupportsIndex
from typing_extensions import TypeVar, override

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from ._data import _data_matrix, _minmax_mixin
from ._index import IndexMixin

__all__: list[str] = []

_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool_, default=Any, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=onp.AtLeast1D, default=onp.AtLeast0D[Any], covariant=True)

###

class _cs_matrix(
    _data_matrix[_ScalarT_co, _ShapeT_co],
    _minmax_mixin[_ScalarT_co, _ShapeT_co],
    IndexMixin[_ScalarT_co, _ShapeT_co],
    Generic[_ScalarT_co, _ShapeT_co],
):
    data: onp.ArrayND[_ScalarT_co]
    indices: onp.Array1D[np.int32]
    indptr: onp.Array1D[np.int32]

    @property
    @override
    @abc.abstractmethod
    def format(self, /) -> Literal["bsr", "csc", "csr"]: ...

    #
    @property
    def has_canonical_format(self, /) -> bool: ...
    @has_canonical_format.setter
    def has_canonical_format(self, val: bool, /) -> None: ...

    #
    @property
    def has_sorted_indices(self, /) -> bool: ...
    @has_sorted_indices.setter
    def has_sorted_indices(self, val: bool, /) -> None: ...

    #
    def __init__(
        self,
        /,
        arg1: onp.ToComplexND,
        shape: tuple[SupportsIndex, *tuple[SupportsIndex, ...]] | None = None,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        *,
        maxprint: int | None = None,
    ) -> None: ...

    #
    def sorted_indices(self, /) -> Self: ...
    def sort_indices(self, /) -> None: ...

    #
    def check_format(self, /, full_check: bool = True) -> None: ...
    def eliminate_zeros(self, /) -> None: ...
    def sum_duplicates(self, /) -> None: ...
    def prune(self, /) -> None: ...
