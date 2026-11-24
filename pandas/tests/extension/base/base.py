import pytest

from pandas.core.dtypes.common import (
    is_list_like,
    is_object_dtype,
)


class BaseExtensionTests:
    pass


def skip_if_list_like(*values, dtype=None) -> None:
    """
    Skip the test if any value is list-like.

    pandas reads a list-like value as a collection rather than one scalar,
    GH#23427. Pass ``dtype`` for single-cell setting, which object dtype
    supports, GH#69152.
    """
    if dtype is not None and is_object_dtype(dtype):
        return
    if any(is_list_like(value) for value in values):
        pytest.skip("pandas reads a list-like scalar as a collection of values")
