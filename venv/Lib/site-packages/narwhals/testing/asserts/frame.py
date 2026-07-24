from __future__ import annotations

from typing import TYPE_CHECKING, Any

from narwhals._utils import Implementation, qualified_type_name
from narwhals.dataframe import DataFrame, LazyFrame
from narwhals.dependencies import is_narwhals_dataframe, is_narwhals_lazyframe
from narwhals.testing.asserts.series import assert_series_equal
from narwhals.testing.asserts.utils import (
    raise_assertion_error,
    raise_frame_assertion_error,
)

if TYPE_CHECKING:
    from narwhals._typing import Arrow, IntoBackend, Pandas, Polars
    from narwhals.typing import DataFrameT, LazyFrameT

GUARANTEES_ROW_ORDER = {
    Implementation.PANDAS,
    Implementation.MODIN,
    Implementation.CUDF,
    Implementation.PYARROW,
    Implementation.POLARS,
    Implementation.DASK,
}


def assert_frame_equal(
    left: DataFrameT | LazyFrameT,
    right: DataFrameT | LazyFrameT,
    *,
    check_row_order: bool = True,
    check_column_order: bool = True,
    check_dtypes: bool = True,
    check_exact: bool = False,
    rel_tol: float = 1e-5,
    abs_tol: float = 1e-8,
    categorical_as_str: bool = False,
    backend: IntoBackend[Polars | Pandas | Arrow] | None = None,
) -> None:
    """Assert that the left and right frames are equal.

    Raises a detailed `AssertionError` if the frames differ.
    This function is intended for use in unit tests.

    Warning:
        1. In the case of backends that do not guarantee the row order, such as DuckDB,
            Ibis, PySpark, and SQLFrame, `check_row_order` argument is ignored and the
            comparands are sorted by all the columns regardless.
        2. In the case of lazy backends a [`collect(...)`](lazyframe.md#narwhals.dataframe.LazyFrame.collect)
            operation is triggered.

    Arguments:
        left: The first DataFrame or LazyFrame to compare.
        right: The second DataFrame or LazyFrame to compare.
        check_row_order: Requires row order to match.

            This flag is ignored for backends that do not guarantee row order such as
            DuckDB, Ibis, PySpark, SQLFrame.
        check_column_order: Requires column order to match.
        check_dtypes: Requires data types to match.
        check_exact: Requires float values to match exactly. If set to `False`, values are
            considered equal when within tolerance of each other (see `rel_tol` and `abs_tol`).

            Only affects columns with a Float data type.
        rel_tol: Relative tolerance for inexact checking. Fraction of values in `right`.
        abs_tol: Absolute tolerance for inexact checking.
        categorical_as_str: Cast categorical columns to string before comparing.
            Enabling this helps compare columns that do not share the same string cache.
        backend: Allows to specify which eager backend to collect to.
            Check out [`narwhals.LazyFrame.collect`](lazyframe.md#narwhals.dataframe.LazyFrame.collect)
            for more information.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> from narwhals.testing import assert_frame_equal
        >>>
        >>> left_native = pl.LazyFrame({"a": [1, 2, 3]})
        >>> right_native = pl.LazyFrame({"a": [1, 5, 3]})
        >>> left = nw.from_native(left_native)
        >>> right = nw.from_native(right_native)
        >>> assert_frame_equal(left, right)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        AssertionError: DataFrames are different (value mismatch for column "a")
        [left]:
        ┌─────────────────┐
        | Narwhals Series |
        |-----------------|
        |shape: (3,)      |
        |Series: 'a' [i64]|
        |[                |
        |        1        |
        |        2        |
        |        3        |
        |]                |
        └─────────────────┘
        [right]:
        ┌─────────────────┐
        | Narwhals Series |
        |-----------------|
        |shape: (3,)      |
        |Series: 'a' [i64]|
        |[                |
        |        1        |
        |        5        |
        |        3        |
        |]                |
        └─────────────────┘
    """
    __tracebackhide__ = True

    if any(
        not (is_narwhals_dataframe(obj) or is_narwhals_lazyframe(obj))
        for obj in (left, right)
    ):
        msg = (
            "Expected `narwhals.DataFrame` or `narwhals.LazyFrame` instance, found:\n"
            f"[left]: {qualified_type_name(type(left))}\n"
            f"[right]: {qualified_type_name(type(right))}\n\n"
            "Hint: Use `nw.from_native(obj, allow_series=False)` to convert each native "
            "object into a `narwhals.DataFrame` or `narwhals.LazyFrame` first."
        )
        raise TypeError(msg)

    left_impl, right_impl = left.implementation, right.implementation
    if left_impl != right_impl:
        raise_frame_assertion_error("implementation mismatch", left_impl, right_impl)

    left_eager, right_eager = _check_correct_input_type(left, right, backend=backend)

    _assert_dataframe_equal(
        left=left_eager,
        right=right_eager,
        impl=left_impl,
        check_row_order=check_row_order,
        check_column_order=check_column_order,
        check_dtypes=check_dtypes,
        check_exact=check_exact,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        categorical_as_str=categorical_as_str,
    )


def _check_correct_input_type(  # noqa: RET503
    left: DataFrameT | LazyFrameT,
    right: DataFrameT | LazyFrameT,
    backend: IntoBackend[Polars | Pandas | Arrow] | None,
) -> tuple[DataFrame[Any], DataFrame[Any]]:
    # Adapted from https://github.com/pola-rs/polars/blob/afdbf3056d1228cf493901e45f536b0905cec8ea/py-polars/src/polars/testing/asserts/frame.py#L15-L17
    if isinstance(left, DataFrame) and isinstance(right, DataFrame):
        return left, right

    if isinstance(left, LazyFrame) and isinstance(right, LazyFrame):
        return left.collect(backend), right.collect(backend)

    raise_assertion_error(
        "inputs",
        "unexpected input types",
        left=type(left).__name__,
        right=type(right).__name__,
    )


def _assert_dataframe_equal(
    left: DataFrameT,
    right: DataFrameT,
    impl: Implementation,
    *,
    check_row_order: bool,
    check_column_order: bool,
    check_dtypes: bool,
    check_exact: bool,
    rel_tol: float,
    abs_tol: float,
    categorical_as_str: bool,
) -> None:
    # Adapted from https://github.com/pola-rs/polars/blob/afdbf3056d1228cf493901e45f536b0905cec8ea/crates/polars-testing/src/asserts/utils.rs#L829
    # NOTE: Here `impl` comes from the original dataframe, not the `.collect`-ed one, and
    # it's used to distinguish between backends that do and do not guarantee row order.
    _check_schema_equal(
        left, right, check_dtypes=check_dtypes, check_column_order=check_column_order
    )

    left_len, right_len = len(left), len(right)
    if left_len != right_len:
        raise_frame_assertion_error("height (row count) mismatch", left_len, right_len)

    if left_len == 0:  # Return early due to same schema but no values
        return

    left_schema = left.schema
    if (not check_row_order) or (impl not in GUARANTEES_ROW_ORDER):
        # !NOTE: Sort by all the non-nested dtypes columns.
        # See: https://github.com/narwhals-dev/narwhals/issues/2939
        # !WARNING: This might lead to wrong results if there are duplicate values in the
        # sorting columns as the final order might still be non fully deterministic.
        sort_by = [name for name, dtype in left_schema.items() if not dtype.is_nested()]

        if not sort_by:
            # If only nested dtypes are available, then we raise an exception.
            msg = "`check_row_order=False` is not supported (yet) with only nested data type."
            raise NotImplementedError(msg)

        left = left.sort(sort_by)
        right = right.sort(sort_by)

    for col_name in left_schema.names():
        _series_left = left.get_column(col_name)
        _series_right = right.get_column(col_name)
        try:
            assert_series_equal(
                _series_left,
                _series_right,
                check_dtypes=False,
                check_names=False,
                check_order=True,
                check_exact=check_exact,
                rel_tol=rel_tol,
                abs_tol=abs_tol,
                categorical_as_str=categorical_as_str,
            )
        except AssertionError:
            raise_frame_assertion_error(
                detail="value mismatch for column",
                left=_series_left,
                right=_series_right,
                detail_suffix=f' "{col_name}"',
            )


def _check_schema_equal(
    left: DataFrameT, right: DataFrameT, *, check_dtypes: bool, check_column_order: bool
) -> None:
    """Compares DataFrame schema based on specified criteria.

    Adapted from https://github.com/pola-rs/polars/blob/afdbf3056d1228cf493901e45f536b0905cec8ea/crates/polars-testing/src/asserts/utils.rs#L667-L698
    """
    lschema, rschema = left.schema, right.schema

    # Fast path for equal DataFrames
    if lschema == rschema:
        return

    lnames, rnames = lschema.names(), rschema.names()
    lset, rset = set(lnames), set(rnames)

    if left_not_in_right := sorted(lset.difference(rset)):
        raise_frame_assertion_error(
            detail="in left, but not in right",
            left=lset,
            right=rset,
            detail_prefix=f"{left_not_in_right} ",
        )
    if right_not_in_left := sorted(rset.difference(lset)):
        raise_frame_assertion_error(
            detail="in right, but not in left",
            left=lset,
            right=rset,
            detail_prefix=f"{right_not_in_left} ",
        )

    if check_column_order and lnames != rnames:
        raise_frame_assertion_error(
            detail="columns are not in the same order", left=lnames, right=rnames
        )

    if check_dtypes:
        rdtypes = (
            rschema.dtypes()
            if check_column_order
            else [rschema[col_name] for col_name in lnames]
        )

        if (ldtypes := lschema.dtypes()) != rdtypes:
            raise_frame_assertion_error(
                detail="dtypes do not match", left=ldtypes, right=rdtypes
            )
