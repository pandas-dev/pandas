from collections.abc import Hashable, Iterable, Mapping, Sequence
from typing import Any

import numpy as np

from xarray import Variable
from xarray.core.indexes import Index, PandasIndex
from xarray.core.types import Self


class ScalarIndex(Index):
    def __init__(self, value: int):
        self.value = value

    @classmethod
    def from_variables(cls, variables, *, options) -> Self:
        var = next(iter(variables.values()))
        return cls(int(var.values))

    def equals(self, other, *, exclude=None) -> bool:
        return isinstance(other, ScalarIndex) and other.value == self.value


class XYIndex(Index):
    def __init__(self, x: PandasIndex, y: PandasIndex):
        self.x: PandasIndex = x
        self.y: PandasIndex = y

    @classmethod
    def from_variables(cls, variables, *, options):
        return cls(
            x=PandasIndex.from_variables({"x": variables["x"]}, options=options),
            y=PandasIndex.from_variables({"y": variables["y"]}, options=options),
        )

    def create_variables(
        self, variables: Mapping[Any, Variable] | None = None
    ) -> dict[Any, Variable]:
        return self.x.create_variables() | self.y.create_variables()

    def equals(self, other, exclude=None):
        if exclude is None:
            exclude = frozenset()
        x_eq = True if self.x.dim in exclude else self.x.equals(other.x)
        y_eq = True if self.y.dim in exclude else self.y.equals(other.y)
        return x_eq and y_eq

    @classmethod
    def concat(
        cls,
        indexes: Sequence[Self],
        dim: Hashable,
        positions: Iterable[Iterable[int]] | None = None,
    ) -> Self:
        first = next(iter(indexes))
        if dim == "x":
            newx = PandasIndex.concat(
                tuple(i.x for i in indexes), dim=dim, positions=positions
            )
            newy = first.y
        elif dim == "y":
            newx = first.x
            newy = PandasIndex.concat(
                tuple(i.y for i in indexes), dim=dim, positions=positions
            )
        return cls(x=newx, y=newy)

    def isel(self, indexers: Mapping[Any, int | slice | np.ndarray | Variable]) -> Self:
        newx = self.x.isel({"x": indexers.get("x", slice(None))})
        newy = self.y.isel({"y": indexers.get("y", slice(None))})
        assert newx is not None
        assert newy is not None
        return type(self)(newx, newy)
