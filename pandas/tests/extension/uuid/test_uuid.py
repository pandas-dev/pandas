from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    ClassVar,
    Self,
)
from uuid import (
    UUID,
    uuid4,
)

import numpy as np
import pytest

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.dtypes import (
    ExtensionDtype,
    register_extension_dtype,
)

import pandas as pd
from pandas.core.arrays.base import ExtensionArray
from pandas.tests.extension import base

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


@register_extension_dtype
class UuidDtype(ExtensionDtype):
    # ExtensionDtype essential API (3 class attrs and methods)

    name: ClassVar[str] = "uuid"
    type: ClassVar[builtins.type[UUID]] = UUID

    @classmethod
    def construct_array_type(cls) -> builtins.type[UuidExtensionArray]:
        return UuidExtensionArray

    # ExtensionDtype overrides
    kind: ClassVar[str] = _UuidNumpyDtype.kind

    @classmethod
    def _get_plot_converter(cls):
        from matplotlib.category import StrCategoryConverter

        class UuidCategoryConverter(StrCategoryConverter):
            @staticmethod
            def convert(value, unit, axis):
                if isinstance(value, UUID):
                    return StrCategoryConverter.convert(value.bytes, unit, axis)
                elif is_list_like(value) and all(isinstance(u, UUID) for u in value):
                    return StrCategoryConverter.convert(
                        [u.bytes for u in value], unit, axis
                    )
                return StrCategoryConverter.convert(value, unit, axis)

            @staticmethod
            def default_units(data, axis):
                if isinstance(data, UUID):
                    return StrCategoryConverter.default_units(data.bytes, axis)
                elif is_list_like(data) and all(isinstance(u, UUID) for u in data):
                    return StrCategoryConverter.default_units(
                        [u.bytes for u in data], axis
                    )
                return StrCategoryConverter.default_units(data, axis)

        return [(cls.type, UuidCategoryConverter)]


class UuidExtensionArray(ExtensionArray):
    # Implementation details and convenience

    _data: NDArray[np.void]

    def __init__(
        self, values: Iterable[UUID] | NDArray[np.void], *, copy: bool = False
    ) -> None:
        if isinstance(values, np.ndarray) and values.dtype == _UuidNumpyDtype:
            self._data = values.copy() if copy else values
        else:
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

    def __getitem__(self, index: ScalarIndexer | slice) -> UUID:  # type: ignore[override]
        if isinstance(index, slice):
            return type(self)(self._data[index])
        assert isinstance(index, int | np.integer)
        return UUID(bytes=self._data[index].tobytes())

    def __len__(self) -> int:
        return len(self._data)

    def copy(self) -> Self:
        return type(self)([UUID(bytes=x.tobytes()) for x in self._data], copy=True)

    def isna(self) -> NDArray[np.bool_]:
        # Array does not support missing values as np.void does not support them
        return np.zeros(len(self._data), dtype=bool)


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


@pytest.fixture
def data():
    return UuidExtensionArray([uuid4() for _ in range(10)])


class TestPlotting(base.BasePlottingTests):
    pass
