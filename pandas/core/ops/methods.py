"""
Functions to generate methods and pin them to the appropriate classes.
"""
from __future__ import annotations

import operator
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Hashable,
    Literal,
    cast,
)

import numpy as np

from pandas._libs import lib
from pandas._typing import (
    ArrayLike,
    Axis,
    AxisInt,
    Level,
)
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import (
    is_array_like,
    is_list_like,
)

from pandas.core import (
    algorithms,
    roperator,
)
from pandas.core.ops.array_ops import (
    get_array_op,
    maybe_prepare_scalar_for_op,
)
from pandas.core.ops.common import get_op_result_name
from pandas.core.ops.docstrings import make_flex_doc

if TYPE_CHECKING:
    from pandas import (
        DataFrame,
        Series,
    )
    from pandas.core.internals import (
        ArrayManager,
        BlockManager,
    )


class FrameOps:
    _get_axis_number: Callable[[Any], int]
    _mgr: BlockManager | ArrayManager

    def _cmp_method(self, other, op):
        axis: Literal[1] = 1  # only relevant for Series other case

        self, other = self._align_for_op(other, axis, flex=False, level=None)

        # See GH#4537 for discussion of scalar op behavior
        new_data = self._dispatch_frame_op(other, op, axis=axis)
        return self._construct_result(new_data)

    def _arith_method(self, other, op):
        if self._should_reindex_frame_op(other, op, 1, None, None):
            return self._arith_method_with_reindex(other, op)

        axis: Literal[1] = 1  # only relevant for Series other case
        other = maybe_prepare_scalar_for_op(
            other,
            # error: "FrameOps" has no attribute "shape"
            (self.shape[axis],),  # type: ignore[attr-defined]
        )

        self, other = self._align_for_op(other, axis, flex=True, level=None)

        new_data = self._dispatch_frame_op(other, op, axis=axis)
        return self._construct_result(new_data)

    _logical_method = _arith_method

    def _arith_method_with_reindex(self, right: DataFrame, op) -> DataFrame:
        """
        For DataFrame-with-DataFrame operations that require reindexing,
        operate only on shared columns, then reindex.

        Parameters
        ----------
        right : DataFrame
        op : binary operator

        Returns
        -------
        DataFrame
        """
        left = cast("DataFrame", self)

        # GH#31623, only operate on shared columns
        cols, lcols, rcols = left.columns.join(
            right.columns, how="inner", level=None, return_indexers=True
        )

        new_left = left.iloc[:, lcols]
        new_right = right.iloc[:, rcols]
        result = op(new_left, new_right)

        # Do the join on the columns instead of using left._align_for_op
        #  to avoid constructing two potentially large/sparse DataFrames
        join_columns, _, _ = left.columns.join(
            right.columns, how="outer", level=None, return_indexers=True
        )

        if result.columns.has_duplicates:
            # Avoid reindexing with a duplicate axis.
            # https://github.com/pandas-dev/pandas/issues/35194
            indexer, _ = result.columns.get_indexer_non_unique(join_columns)
            indexer = algorithms.unique1d(indexer)
            result = result._reindex_with_indexers(
                {1: [join_columns, indexer]}, allow_dups=True
            )
        else:
            result = result.reindex(join_columns, axis=1)

        return result

    def _should_reindex_frame_op(self, right, op, axis: int, fill_value, level) -> bool:
        """
        Check if this is an operation between DataFrames that will need to reindex.
        """
        if op is operator.pow or op is roperator.rpow:
            # GH#32685 pow has special semantics for operating with null values
            return False

        if not isinstance(right, FrameOps):
            return False

        if fill_value is None and level is None and axis == 1:
            # TODO: any other cases we should handle here?

            # Intersection is always unique so we have to check the unique columns
            # error: "FrameOps" has no attribute "columns"
            left_uniques = self.columns.unique()  # type: ignore[attr-defined]
            # error: "FrameOps" has no attribute "columns"
            right_uniques = right.columns.unique()  # type: ignore[attr-defined]
            cols = left_uniques.intersection(right_uniques)
            if len(cols) and not (
                len(cols) == len(left_uniques) and len(cols) == len(right_uniques)
            ):
                # TODO: is there a shortcut available when len(cols) == 0?
                return True

        return False

    def _align_for_op(
        self, other, axis, flex: bool | None = False, level: Level = None
    ):
        """
        Convert rhs to meet lhs dims if input is list, tuple or np.ndarray.

        Parameters
        ----------
        left : DataFrame
        right : Any
        axis : int, str, or None
        flex : bool or None, default False
            Whether this is a flex op, in which case we reindex.
            None indicates not to check for alignment.
        level : int or level name, default None

        Returns
        -------
        left : DataFrame
        right : Any
        """
        self = cast("DataFrame", self)
        left, right = self, other

        def to_series(right):
            msg = (
                "Unable to coerce to Series, "
                "length must be {req_len}: given {given_len}"
            )

            # pass dtype to avoid doing inference, which would break consistency
            #  with Index/Series ops
            dtype = None
            if getattr(right, "dtype", None) == object:
                # can't pass right.dtype unconditionally as that would break on e.g.
                #  datetime64[h] ndarray
                dtype = object

            if axis is not None and left._get_axis_number(axis) == 0:
                if len(left.index) != len(right):
                    raise ValueError(
                        msg.format(req_len=len(left.index), given_len=len(right))
                    )
                right = left._constructor_sliced(right, index=left.index, dtype=dtype)
            else:
                if len(left.columns) != len(right):
                    raise ValueError(
                        msg.format(req_len=len(left.columns), given_len=len(right))
                    )
                right = left._constructor_sliced(right, index=left.columns, dtype=dtype)
            return right

        if isinstance(right, np.ndarray):
            if right.ndim == 1:
                right = to_series(right)

            elif right.ndim == 2:
                # We need to pass dtype=right.dtype to retain object dtype
                #  otherwise we lose consistency with Index and array ops
                dtype = None
                if right.dtype == object:
                    # can't pass right.dtype unconditionally as that would break on e.g.
                    #  datetime64[h] ndarray
                    dtype = object

                if right.shape == left.shape:
                    right = left._constructor(
                        right, index=left.index, columns=left.columns, dtype=dtype
                    )

                elif right.shape[0] == left.shape[0] and right.shape[1] == 1:
                    # Broadcast across columns
                    right = np.broadcast_to(right, left.shape)
                    right = left._constructor(
                        right, index=left.index, columns=left.columns, dtype=dtype
                    )

                elif right.shape[1] == left.shape[1] and right.shape[0] == 1:
                    # Broadcast along rows
                    right = to_series(right[0, :])

                else:
                    raise ValueError(
                        "Unable to coerce to DataFrame, shape "
                        f"must be {left.shape}: given {right.shape}"
                    )

            elif right.ndim > 2:
                raise ValueError(
                    "Unable to coerce to Series/DataFrame, "
                    f"dimension must be <= 2: {right.shape}"
                )

        elif is_list_like(right) and not isinstance(right, (SeriesOps, FrameOps)):
            # GH#36702. Raise when attempting arithmetic with list of array-like.
            if any(is_array_like(el) for el in right):
                raise ValueError(
                    f"Unable to coerce list of {type(right[0])} to Series/DataFrame"
                )
            # GH#17901
            right = to_series(right)

        if flex is not None and isinstance(right, FrameOps):
            rframe = cast("DataFrame", right)
            if not left._indexed_same(rframe):
                if flex:
                    left, right = left.align(
                        rframe, join="outer", level=level, copy=False
                    )
                else:
                    raise ValueError(
                        "Can only compare identically-labeled (both index and columns) "
                        "DataFrame objects"
                    )
        elif isinstance(right, SeriesOps):
            right = cast("Series", right)

            # axis=1 is default for DataFrame-with-Series op
            axis = left._get_axis_number(axis) if axis is not None else 1

            if not flex:
                if not left.axes[axis].equals(right.index):
                    raise ValueError(
                        "Operands are not aligned. Do "
                        "`left, right = left.align(right, axis=1, copy=False)` "
                        "before operating."
                    )

            left, right = left.align(
                # error: Argument 1 to "align" of "DataFrame" has incompatible
                # type "Series"; expected "DataFrame"
                right,  # type: ignore[arg-type]
                join="outer",
                axis=axis,
                level=level,
                copy=False,
            )
            right = left._maybe_align_series_as_frame(right, axis)

        return left, right

    def _maybe_align_series_as_frame(self, series: Series, axis: AxisInt):
        """
        If the Series operand is not EA-dtype, we can broadcast to 2D and operate
        blockwise.
        """
        rvalues = series._values
        if not isinstance(rvalues, np.ndarray):
            # TODO(EA2D): no need to special-case with 2D EAs
            if rvalues.dtype in ("datetime64[ns]", "timedelta64[ns]"):
                # We can losslessly+cheaply cast to ndarray
                rvalues = np.asarray(rvalues)
            else:
                return series

        if axis == 0:
            rvalues = rvalues.reshape(-1, 1)
        else:
            rvalues = rvalues.reshape(1, -1)

        # error: "FrameOps" has no attribute "shape"
        rvalues = np.broadcast_to(rvalues, self.shape)  # type: ignore[attr-defined]
        # pass dtype to avoid doing inference
        # error: "FrameOps" has no attribute "_constructor"
        return self._constructor(  # type: ignore[attr-defined]
            rvalues,
            # error: "FrameOps" has no attribute "index"
            index=self.index,  # type: ignore[attr-defined]
            # error: "FrameOps" has no attribute "columns"
            columns=self.columns,  # type: ignore[attr-defined]
            dtype=rvalues.dtype,
        )

    def _dispatch_frame_op(self, right, func: Callable, axis: AxisInt | None = None):
        """
        Evaluate the frame operation func(left, right) by evaluating
        column-by-column, dispatching to the Series implementation.

        Parameters
        ----------
        right : scalar, Series, or DataFrame
        func : arithmetic or comparison operator
        axis : {None, 0, 1}

        Returns
        -------
        DataFrame
        """
        # Get the appropriate array-op to apply to each column/block's values.
        array_op = get_array_op(func)

        right = lib.item_from_zerodim(right)
        if not is_list_like(right):
            # i.e. scalar, faster than checking np.ndim(right) == 0
            with np.errstate(all="ignore"):
                bm = self._mgr.apply(array_op, right=right)
            # error: "FrameOps" has no attribute "_constructor"
            return self._constructor(bm)  # type: ignore[attr-defined]

        elif isinstance(right, FrameOps):
            # error: "FrameOps" has no attribute "index"
            assert self.index.equals(right.index)  # type: ignore[attr-defined]
            # error: "FrameOps" has no attribute "columns"
            assert self.columns.equals(right.columns)  # type: ignore[attr-defined]
            # TODO: The previous assertion `assert right._indexed_same(self)`
            #  fails in cases with empty columns reached via
            #  _arith_method_with_reindex

            # TODO operate_blockwise expects a manager of the same type
            with np.errstate(all="ignore"):
                bm = self._mgr.operate_blockwise(
                    # error: Argument 1 to "operate_blockwise" of "ArrayManager" has
                    # incompatible type "Union[ArrayManager, BlockManager]"; expected
                    # "ArrayManager"
                    # error: Argument 1 to "operate_blockwise" of "BlockManager" has
                    # incompatible type "Union[ArrayManager, BlockManager]"; expected
                    # "BlockManager"
                    right._mgr,  # type: ignore[arg-type]
                    array_op,
                )
            # error: "FrameOps" has no attribute "_constructor"
            return self._constructor(bm)  # type: ignore[attr-defined]

        elif isinstance(right, SeriesOps) and axis == 1:
            # axis=1 means we want to operate row-by-row
            # error: "FrameOps" has no attribute "columns"
            assert right.index.equals(self.columns)  # type: ignore[attr-defined]

            right = right._values
            # maybe_align_as_frame ensures we do not have an ndarray here
            assert not isinstance(right, np.ndarray)

            # error: "FrameOps" has no attribute "_iter_column_arrays"
            col_arrays = self._iter_column_arrays()  # type: ignore[attr-defined]
            with np.errstate(all="ignore"):
                arrays = [
                    array_op(_left, _right) for _left, _right in zip(col_arrays, right)
                ]

        elif isinstance(right, SeriesOps):
            # error: "FrameOps" has no attribute "index"
            assert right.index.equals(self.index)  # type: ignore[attr-defined]
            right = right._values

            # error: "FrameOps" has no attribute "_iter_column_arrays"
            col_arrays = self._iter_column_arrays()  # type: ignore[attr-defined]
            with np.errstate(all="ignore"):
                arrays = [array_op(left, right) for left in col_arrays]

        else:
            raise NotImplementedError(right)

        # error: "Type[FrameOps]" has no attribute "_from_arrays"
        return type(self)._from_arrays(  # type: ignore[attr-defined]
            arrays,
            # error: "FrameOps" has no attribute "columns"
            self.columns,  # type: ignore[attr-defined]
            # error: "FrameOps" has no attribute "index"
            self.index,  # type: ignore[attr-defined]
            verify_integrity=False,
        )

    def _combine_frame(self, other: FrameOps, func, fill_value=None):
        # at this point we have `self._indexed_same(other)`

        if fill_value is None:
            # since _arith_op may be called in a loop, avoid function call
            #  overhead if possible by doing this check once
            _arith_op = func

        else:
            from pandas.core.ops import fill_binop

            def _arith_op(left, right):
                # for the mixed_type case where we iterate over columns,
                # _arith_op(left, right) is equivalent to
                # left._binop(right, func, fill_value=fill_value)
                left, right = fill_binop(left, right, fill_value)
                return func(left, right)

        new_data = self._dispatch_frame_op(other, _arith_op)
        return new_data

    def _construct_result(self, result) -> DataFrame:
        """
        Wrap the result of an arithmetic, comparison, or logical operation.

        Parameters
        ----------
        result : DataFrame

        Returns
        -------
        DataFrame
        """
        # error: "FrameOps" has no attribute "_constructor"
        out = self._constructor(result, copy=False)  # type: ignore[attr-defined]
        out = out.__finalize__(self)
        # Pin columns instead of passing to constructor for compat with
        #  non-unique columns case
        # error: "FrameOps" has no attribute "columns"
        out.columns = self.columns  # type: ignore[attr-defined]
        # error: "FrameOps" has no attribute "index"
        out.index = self.index  # type: ignore[attr-defined]
        return out

    def __divmod__(self, other) -> tuple[DataFrame, DataFrame]:
        # Naive implementation, room for optimization
        div = self // other
        mod = self - div * other
        return div, mod

    def __rdivmod__(self, other) -> tuple[DataFrame, DataFrame]:
        # Naive implementation, room for optimization
        div = other // self
        mod = other - div * self
        return div, mod

    def _flex_arith_method(
        self, other, op, *, axis: Axis = "columns", level=None, fill_value=None
    ):
        axis = self._get_axis_number(axis) if axis is not None else 1

        if self._should_reindex_frame_op(other, op, axis, fill_value, level):
            return self._arith_method_with_reindex(other, op)

        if isinstance(other, SeriesOps) and fill_value is not None:
            # TODO: We could allow this in cases where we end up going
            #  through the DataFrame path
            raise NotImplementedError(f"fill_value {fill_value} not supported.")

        other = maybe_prepare_scalar_for_op(
            other,
            # error: "FrameOps" has no attribute "shape"
            self.shape,  # type: ignore[attr-defined]
        )
        self, other = self._align_for_op(other, axis, flex=True, level=level)

        if isinstance(other, FrameOps):
            # Another DataFrame
            new_data = self._combine_frame(other, op, fill_value)

        elif isinstance(other, SeriesOps):
            new_data = self._dispatch_frame_op(other, op, axis=axis)
        else:
            # in this case we always have `np.ndim(other) == 0`
            if fill_value is not None:
                # error: "FrameOps" has no attribute "fillna"
                self = self.fillna(fill_value)  # type: ignore[attr-defined]

            new_data = self._dispatch_frame_op(other, op)

        return self._construct_result(new_data)

    def _flex_cmp_method(self, other, op, *, axis: Axis = "columns", level=None):
        axis = self._get_axis_number(axis) if axis is not None else 1

        self, other = self._align_for_op(other, axis, flex=True, level=level)

        new_data = self._dispatch_frame_op(other, op, axis=axis)
        return self._construct_result(new_data)

    @Appender(make_flex_doc("eq", "dataframe"))
    def eq(self, other, axis: Axis = "columns", level=None):
        return self._flex_cmp_method(other, operator.eq, axis=axis, level=level)

    @Appender(make_flex_doc("ne", "dataframe"))
    def ne(self, other, axis: Axis = "columns", level=None):
        return self._flex_cmp_method(other, operator.ne, axis=axis, level=level)

    @Appender(make_flex_doc("le", "dataframe"))
    def le(self, other, axis: Axis = "columns", level=None):
        return self._flex_cmp_method(other, operator.le, axis=axis, level=level)

    @Appender(make_flex_doc("lt", "dataframe"))
    def lt(self, other, axis: Axis = "columns", level=None):
        return self._flex_cmp_method(other, operator.lt, axis=axis, level=level)

    @Appender(make_flex_doc("ge", "dataframe"))
    def ge(self, other, axis: Axis = "columns", level=None):
        return self._flex_cmp_method(other, operator.ge, axis=axis, level=level)

    @Appender(make_flex_doc("gt", "dataframe"))
    def gt(self, other, axis: Axis = "columns", level=None):
        return self._flex_cmp_method(other, operator.gt, axis=axis, level=level)

    @Appender(make_flex_doc("add", "dataframe"))
    def add(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.add, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("radd", "dataframe"))
    def radd(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.radd, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("sub", "dataframe"))
    def sub(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.sub, level=level, fill_value=fill_value, axis=axis
        )

    subtract = sub

    @Appender(make_flex_doc("rsub", "dataframe"))
    def rsub(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.rsub, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("mul", "dataframe"))
    def mul(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.mul, level=level, fill_value=fill_value, axis=axis
        )

    multiply = mul

    @Appender(make_flex_doc("rmul", "dataframe"))
    def rmul(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.rmul, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("truediv", "dataframe"))
    def truediv(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.truediv, level=level, fill_value=fill_value, axis=axis
        )

    div = truediv
    divide = truediv

    @Appender(make_flex_doc("rtruediv", "dataframe"))
    def rtruediv(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.rtruediv, level=level, fill_value=fill_value, axis=axis
        )

    rdiv = rtruediv

    @Appender(make_flex_doc("floordiv", "dataframe"))
    def floordiv(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.floordiv, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rfloordiv", "dataframe"))
    def rfloordiv(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.rfloordiv, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("mod", "dataframe"))
    def mod(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.mod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rmod", "dataframe"))
    def rmod(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.rmod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("pow", "dataframe"))
    def pow(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, operator.pow, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rpow", "dataframe"))
    def rpow(self, other, axis: Axis = "columns", level=None, fill_value=None):
        return self._flex_arith_method(
            other, roperator.rpow, level=level, fill_value=fill_value, axis=axis
        )


class SeriesOps:
    _get_axis_number: Callable[[Any], int]
    _values: ArrayLike

    def _binop(self, other: SeriesOps, func, level=None, fill_value=None):
        """
        Perform generic binary operation with optional fill value.

        Parameters
        ----------
        other : Series
        func : binary operator
        fill_value : float or object
            Value to substitute for NA/null values. If both Series are NA in a
            location, the result will be NA regardless of the passed fill value.
        level : int or level name, default None
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.

        Returns
        -------
        Series
        """
        if not isinstance(other, SeriesOps):
            raise AssertionError("Other operand must be Series")

        this = self

        # error: "SeriesOps" has no attribute "index"
        if not self.index.equals(other.index):  # type: ignore[attr-defined]
            # error: "SeriesOps" has no attribute "align"
            this, other = self.align(  # type: ignore[attr-defined]
                other, level=level, join="outer", copy=False
            )

        from pandas.core.ops import fill_binop

        this_vals, other_vals = fill_binop(this._values, other._values, fill_value)

        with np.errstate(all="ignore"):
            result = func(this_vals, other_vals)

        name = get_op_result_name(self, other)
        return this._construct_result(result, name)

    def _construct_result(
        self, result: ArrayLike | tuple[ArrayLike, ArrayLike], name: Hashable
    ) -> Series | tuple[Series, Series]:
        """
        Construct an appropriately-labelled Series from the result of an op.

        Parameters
        ----------
        result : ndarray or ExtensionArray
        name : Label

        Returns
        -------
        Series
            In the case of __divmod__ or __rdivmod__, a 2-tuple of Series.
        """
        if isinstance(result, tuple):
            # produced by divmod or rdivmod

            res1 = self._construct_result(result[0], name=name)
            res2 = self._construct_result(result[1], name=name)

            # GH#33427 assertions to keep mypy happy
            assert isinstance(res1, SeriesOps)
            assert isinstance(res2, SeriesOps)
            return (res1, res2)

        # TODO: result should always be ArrayLike, but this fails for some
        #  JSONArray tests
        dtype = getattr(result, "dtype", None)
        # error: "SeriesOps" has no attribute "_constructor"
        out = self._constructor(  # type: ignore[attr-defined]
            result,
            # error: "SeriesOps" has no attribute "index"
            index=self.index,  # type: ignore[attr-defined]
            dtype=dtype,
        )
        out = out.__finalize__(self)

        # Set the result's name after __finalize__ is called because __finalize__
        #  would set it back to self.name
        out.name = name
        return out

    def _flex_method(self, other, op, *, level=None, fill_value=None, axis: Axis = 0):
        if axis is not None:
            self._get_axis_number(axis)

        res_name = get_op_result_name(self, other)

        if isinstance(other, SeriesOps):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, list, tuple)):
            # error: Argument 1 to "len" has incompatible type "SeriesOps";
            # expected "Sized"
            if len(other) != len(self):  # type: ignore[arg-type]
                raise ValueError("Lengths must be equal")
            # error: "SeriesOps" has no attribute "index"
            other = self._constructor(other, self.index)  # type: ignore[attr-defined]
            result = self._binop(other, op, level=level, fill_value=fill_value)
            result.name = res_name
            return result
        else:
            if fill_value is not None:
                # error: "SeriesOps" has no attribute "fillna"
                self = self.fillna(fill_value)  # type: ignore[attr-defined]

            return op(self, other)

    @Appender(make_flex_doc("eq", "series"))
    def eq(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.eq, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("ne", "series"))
    def ne(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.ne, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("le", "series"))
    def le(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.le, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("lt", "series"))
    def lt(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.lt, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("ge", "series"))
    def ge(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.ge, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("gt", "series"))
    def gt(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.gt, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("add", "series"))
    def add(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.add, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("radd", "series"))
    def radd(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.radd, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("sub", "series"))
    def sub(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.sub, level=level, fill_value=fill_value, axis=axis
        )

    subtract = sub

    @Appender(make_flex_doc("rsub", "series"))
    def rsub(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rsub, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("mul", "series"))
    def mul(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.mul, level=level, fill_value=fill_value, axis=axis
        )

    multiply = mul

    @Appender(make_flex_doc("rmul", "series"))
    def rmul(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rmul, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("truediv", "series"))
    def truediv(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.truediv, level=level, fill_value=fill_value, axis=axis
        )

    div = truediv
    divide = truediv

    @Appender(make_flex_doc("rtruediv", "series"))
    def rtruediv(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rtruediv, level=level, fill_value=fill_value, axis=axis
        )

    rdiv = rtruediv

    @Appender(make_flex_doc("floordiv", "series"))
    def floordiv(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.floordiv, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rfloordiv", "series"))
    def rfloordiv(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rfloordiv, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("mod", "series"))
    def mod(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.mod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rmod", "series"))
    def rmod(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rmod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("pow", "series"))
    def pow(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, operator.pow, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rpow", "series"))
    def rpow(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rpow, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("divmod", "series"))
    def divmod(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, divmod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(make_flex_doc("rdivmod", "series"))
    def rdivmod(self, other, level=None, fill_value=None, axis: Axis = 0):
        return self._flex_method(
            other, roperator.rdivmod, level=level, fill_value=fill_value, axis=axis
        )
