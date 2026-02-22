from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

from narwhals.dependencies import is_narwhals_series

if TYPE_CHECKING:
    from typing_extensions import Never, TypeAlias

    # NOTE: These aliases are created to facilitate autocompletion.
    # Feel free to extend them as you please when adding new features.
    # See: https://github.com/narwhals-dev/narwhals/pull/2983#discussion_r2337548736
    ObjectName: TypeAlias = Literal["inputs", "Series", "DataFrames"]
    SeriesDetail: TypeAlias = Literal[
        "dtype mismatch",
        "exact value mismatch",
        "implementation mismatch",
        "length mismatch",
        "name mismatch",
        "nested value mismatch",
        "null value mismatch",
        "values not within tolerance",
    ]
    DataFramesDetail: TypeAlias = Literal[
        "columns are not in the same order",
        "dtypes do not match",
        "height (row count) mismatch",
        "implementation mismatch",
        "in left, but not in right",
        "in right, but not in left",
        "value mismatch for column",
    ]


def raise_assertion_error(
    objects: ObjectName,
    detail: str,
    left: Any,
    right: Any,
    *,
    cause: Exception | None = None,
) -> Never:
    """Raise a detailed assertion error."""
    __tracebackhide__ = True

    trailing_left = "\n" if is_narwhals_series(left) else " "
    trailing_right = "\n" if is_narwhals_series(right) else " "

    msg = (
        f"{objects} are different ({detail})\n"
        f"[left]:{trailing_left}{left}\n"
        f"[right]:{trailing_right}{right}"
    )
    raise AssertionError(msg) from cause


def raise_series_assertion_error(
    detail: SeriesDetail, left: Any, right: Any, *, cause: Exception | None = None
) -> Never:
    raise_assertion_error("Series", detail, left, right, cause=cause)


def raise_frame_assertion_error(
    detail: DataFramesDetail,
    left: Any,
    right: Any,
    *,
    detail_prefix: str = "",
    detail_suffix: str = "",
    cause: Exception | None = None,
) -> Never:
    raise_assertion_error(
        "DataFrames", f"{detail_prefix}{detail}{detail_suffix}", left, right, cause=cause
    )
