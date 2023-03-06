"""
Functions to generate methods and pin them to the appropriate classes.
"""
from __future__ import annotations

import operator
from typing import cast

import numpy as np

from pandas._typing import Axis
from pandas.util._decorators import Appender

from pandas.core import roperator
from pandas.core.ops.array_ops import maybe_prepare_scalar_for_op
from pandas.core.ops.common import get_op_result_name
from pandas.core.ops.docstrings import make_flex_doc


class FrameOps:
    def _flex_arith_method(
        self, other, op, *, axis: Axis = "columns", level=None, fill_value=None
    ):
        axis = self._get_axis_number(axis) if axis is not None else 1
        axis = cast(int, axis)

        if self._should_reindex_frame_op(other, op, axis, fill_value, level):
            return self._arith_method_with_reindex(other, op)

        if isinstance(other, SeriesOps) and fill_value is not None:
            # TODO: We could allow this in cases where we end up going
            #  through the DataFrame path
            raise NotImplementedError(f"fill_value {fill_value} not supported.")

        other = maybe_prepare_scalar_for_op(other, self.shape)
        self, other = self._align_for_op(other, axis, flex=True, level=level)

        if isinstance(other, FrameOps):
            # Another DataFrame
            new_data = self._combine_frame(other, op, fill_value)

        elif isinstance(other, SeriesOps):
            new_data = self._dispatch_frame_op(other, op, axis=axis)
        else:
            # in this case we always have `np.ndim(other) == 0`
            if fill_value is not None:
                self = self.fillna(fill_value)

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
    def _flex_method(self, other, op, *, level=None, fill_value=None, axis: Axis = 0):
        if axis is not None:
            self._get_axis_number(axis)

        res_name = get_op_result_name(self, other)

        if isinstance(other, SeriesOps):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, list, tuple)):
            if len(other) != len(self):
                raise ValueError("Lengths must be equal")
            other = self._constructor(other, self.index)
            result = self._binop(other, op, level=level, fill_value=fill_value)
            result.name = res_name
            return result
        else:
            if fill_value is not None:
                self = self.fillna(fill_value)

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
