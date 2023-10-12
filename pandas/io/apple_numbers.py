"""Apple Numbers format"""
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.compat._optional import import_optional_dependency

if TYPE_CHECKING:
    from collections.abc import (
        Hashable,
        Sequence,
    )

    from pandas._typing import (
        FilePath,
        StorageOptions,
        WriteBuffer,
    )

    from pandas.core.api import DataFrame


def read_apple_numbers(
    path: FilePath,
    columns: Sequence[Hashable] | None = None,
) -> DataFrame:
    """Load an Apple Numbers document from the file path."""

    import_optional_dependency("numbers_parser")

    from numbers_parser import (
        Document,
        NumbersError,
    )

    try:
        Document(path)
    except NumbersError:
        raise ValueError("Numbers failed: e!r")


def to_apple_numbers(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes],
    storage_options: StorageOptions | None = None,
    **kwargs: Any,
) -> None:
    """
    Write a DataFrame to an Apple Numbers document.
    """
