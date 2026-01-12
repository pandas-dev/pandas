from __future__ import annotations

from collections.abc import (
    Iterable,
    Sequence,
)
from typing import (
    TYPE_CHECKING,
    ClassVar,
    Self,
    cast,
    get_args,
    overload,
)
from uuid import UUID

import numpy as np

from pandas.core.dtypes.dtypes import ExtensionDtype

import pandas as pd
from pandas.core.algorithms import take
from pandas.core.arrays.base import ExtensionArray
from pandas.core.arrays.boolean import BooleanArray
from pandas.core.indexers.utils import check_array_indexer
from pandas.core.ops.common import unpack_zerodim_and_defer

if TYPE_CHECKING:
    import builtins

    from numpy.typing import NDArray

    from pandas._typing import (
        Dtype,
        PositionalIndexer,
        ScalarIndexer,
        SequenceIndexer,
    )

    from pandas.core.arrays.boolean import BooleanArray


UuidLike = UUID | bytes | int | str

# 16 void bytes: 128 bit, every pattern valid, no funky behavior like 0 stripping.
_UuidNumpyDtype = np.dtype("V16")


def _to_uuid(v: UuidLike) -> UUID:
    match v:
        case UUID():
            return v
        case bytes():
            return UUID(bytes=v)
        case int():
            return UUID(int=v)
        case str():
            return UUID(v)
    msg = f"Unknown type for Uuid: {type(v)} is not {get_args(UuidLike)}"
    raise TypeError(msg)


class UuidDtype(ExtensionDtype):
    # ExtensionDtype essential API (3 class attrs and methods)

    name: ClassVar[str] = "uuid"
    type: ClassVar[builtins.type[UUID]] = UUID

    @classmethod
    def construct_array_type(cls) -> builtins.type[UuidExtensionArray]:
        return UuidExtensionArray

    # ExtensionDtype overrides
    kind: ClassVar[str] = _UuidNumpyDtype.kind


class UuidExtensionArray(ExtensionArray):
    # Implementation details and convenience

    _data: NDArray[np.void]

    def __init__(self, values: Iterable[UuidLike], *, copy: bool = False) -> None:
        if isinstance(values, np.ndarray):
            self._data = values.astype(_UuidNumpyDtype, copy=copy)
        else:
            # TODO: more efficient
            self._data = np.array(
                [_to_uuid(x).bytes for x in values], dtype=_UuidNumpyDtype
            )

        if self._data.ndim != 1:
            raise ValueError("Array only supports 1-d arrays")

    # ExtensionArray essential API (11 class attrs and methods)

    dtype: ClassVar[UuidDtype] = UuidDtype()

    @classmethod
    def _from_sequence(
        cls,
        scalars: Iterable[UuidLike],
        *,
        dtype: Dtype | None = None,
        copy: bool = False,
    ) -> Self:
        if dtype is None:
            dtype = UuidDtype()
        return cls(scalars, copy=copy)

    @overload
    def __getitem__(self, index: ScalarIndexer) -> UUID: ...
    @overload
    def __getitem__(self, index: SequenceIndexer) -> Self: ...

    def __getitem__(self, index: PositionalIndexer) -> Self | UUID:
        if isinstance(index, int | np.integer):
            return UUID(bytes=self._data[index].tobytes())
        index = check_array_indexer(self, index)
        return self._simple_new(self._data[index])

    def __len__(self) -> int:
        return len(self._data)


def test_construct() -> None:
    """Tests that we can construct UuidExtensionArray from a list of valid values."""
    from uuid import uuid4

    a = UuidExtensionArray([0, u := uuid4()])
    assert a[0] == UUID(int=0)
    assert a[1] == u


def test_series() -> None:
    """Tests that Series accepts unstructured void dtypes."""
    from uuid import uuid4

    s = pd.Series([u := uuid4()], dtype=UuidDtype(), name="s")
    assert str(u) in str(s)
