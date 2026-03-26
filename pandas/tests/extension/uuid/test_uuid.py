from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    ClassVar,
    Self,
)
from uuid import UUID

import numpy as np

from pandas.core.dtypes.dtypes import ExtensionDtype

import pandas as pd
from pandas.core.arrays.base import ExtensionArray

if TYPE_CHECKING:
    import builtins
    from collections.abc import Iterable

    from numpy.typing import NDArray

    from pandas._typing import (
        Dtype,
        ScalarIndexer,
    )


# 16 void bytes: 128 bit, every pattern valid, no funky behavior like 0 stripping.
_UuidNumpyDtype = np.dtype("V16")


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

    def __init__(self, values: Iterable[UUID], *, copy: bool = False) -> None:
        self._data = np.array([x.bytes for x in values], dtype=_UuidNumpyDtype)

    # Parts of ExtensionArray's essential API required for tests:

    dtype: ClassVar[UuidDtype] = UuidDtype()

    @classmethod
    def _from_sequence(
        cls,
        scalars: Iterable[UUID],
        *,
        dtype: Dtype | None = None,
        copy: bool = False,
    ) -> Self:
        if dtype is None:
            dtype = UuidDtype()
        return cls(scalars, copy=copy)

    def __getitem__(self, index: ScalarIndexer) -> UUID:  # type: ignore[override]
        assert isinstance(index, int | np.integer)
        return UUID(bytes=self._data[index].tobytes())

    def __len__(self) -> int:
        return len(self._data)


def test_construct() -> None:
    """Tests that we can construct UuidExtensionArray from a list of valid values."""
    from uuid import uuid4

    a = UuidExtensionArray([UUID(int=0), u := uuid4()])
    assert a[0].int == 0
    assert a[1] == u


def test_series() -> None:
    """Tests that Series accepts (unstructured) void ExtensionDtypes."""
    from uuid import uuid4

    s = pd.Series([u := uuid4()], dtype=UuidDtype(), name="s")
    assert str(u) in str(s)
