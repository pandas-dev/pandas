import operator

import pytest

import pandas as pd
from pandas.core import ops

from .base import BaseExtensionTests


class BaseOpsUtil(BaseExtensionTests):

    def get_op_from_name(self, op_name):
        short_opname = op_name.strip('_')
        try:
            op = getattr(operator, short_opname)
        except AttributeError:
            # Assume it is the reverse operator
            rop = getattr(operator, short_opname[1:])
            op = lambda x, y: rop(y, x)

        return op

    def check_opname(self, s, op_name, other, exc=Exception):
        op = self.get_op_from_name(op_name)

        self._check_op(s, op, other, op_name, exc)

    def _check_op(self, s, op, other, op_name, exc=NotImplementedError):
        if exc is None:
            result = op(s, other)
            expected = s.combine(other, op)
            self.assert_series_equal(result, expected)
        else:
            with pytest.raises(exc):
                op(s, other)

    def _check_divmod_op(self, s, op, other, exc=Exception):
        # divmod has multiple return values, so check separatly
        if exc is None:
            result_div, result_mod = op(s, other)
            if op is divmod:
                expected_div, expected_mod = s // other, s % other
            else:
                expected_div, expected_mod = other // s, other % s
            self.assert_series_equal(result_div, expected_div)
            self.assert_series_equal(result_mod, expected_mod)
        else:
            with pytest.raises(exc):
                divmod(s, other)


class BaseArithmeticOpsTests(BaseOpsUtil):
    """Various Series and DataFrame arithmetic ops methods.

    Subclasses supporting various ops should set the class variables
    to indicate that they support ops of that kind

    * series_scalar_exc = TypeError
    * frame_scalar_exc = TypeError
    * series_array_exc = TypeError
    * divmod_exc = TypeError
    """
    series_scalar_exc = TypeError
    frame_scalar_exc = TypeError
    series_array_exc = TypeError
    divmod_exc = TypeError

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # series & scalar
        op_name = all_arithmetic_operators
        s = pd.Series(data)
        self.check_opname(s, op_name, s.iloc[0], exc=self.series_scalar_exc)

    @pytest.mark.xfail(run=False, reason="_reduce needs implementation")
    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators):
        # frame & scalar
        op_name = all_arithmetic_operators
        df = pd.DataFrame({'A': data})
        self.check_opname(df, op_name, data[0], exc=self.frame_scalar_exc)

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        # ndarray & other series
        op_name = all_arithmetic_operators
        s = pd.Series(data)
        self.check_opname(s, op_name, pd.Series([s.iloc[0]] * len(s)),
                          exc=self.series_array_exc)

    def test_divmod(self, data):
        s = pd.Series(data)
        self._check_divmod_op(s, divmod, 1, exc=self.divmod_exc)
        self._check_divmod_op(1, ops.rdivmod, s, exc=self.divmod_exc)

    def test_divmod_series_array(self, data, data_for_twos):
        s = pd.Series(data)
        self._check_divmod_op(s, divmod, data)

        other = data_for_twos
        self._check_divmod_op(other, ops.rdivmod, s)

        other = pd.Series(other)
        self._check_divmod_op(other, ops.rdivmod, s)

    def test_add_series_with_extension_array(self, data):
        s = pd.Series(data)
        result = s + data
        expected = pd.Series(data + data)
        self.assert_series_equal(result, expected)

    def test_error(self, data, all_arithmetic_operators):
        # invalid ops
        op_name = all_arithmetic_operators
        with pytest.raises(AttributeError):
            getattr(data, op_name)

    def test_direct_arith_with_series_returns_not_implemented(self, data):
        # EAs should return NotImplemented for ops with Series.
        # Pandas takes care of unboxing the series and calling the EA's op.
        other = pd.Series(data)
        if hasattr(data, '__add__'):
            result = data.__add__(other)
            assert result is NotImplemented
        else:
            raise pytest.skip(
                "{} does not implement add".format(data.__class__.__name__)
            )


class BaseComparisonOpsTests(BaseOpsUtil):
    """Various Series and DataFrame comparison ops methods."""

    def _compare_other(self, s, data, op_name, other):
        op = self.get_op_from_name(op_name)
        if op_name == '__eq__':
            assert getattr(data, op_name)(other) is NotImplemented
            assert not op(s, other).all()
        elif op_name == '__ne__':
            assert getattr(data, op_name)(other) is NotImplemented
            assert op(s, other).all()

        else:

            # array
            assert getattr(data, op_name)(other) is NotImplemented

            # series
            s = pd.Series(data)
            with pytest.raises(TypeError):
                op(s, other)

    def test_compare_scalar(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        self._compare_other(s, data, op_name, 0)

    def test_compare_array(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        other = pd.Series([data[0]] * len(data))
        self._compare_other(s, data, op_name, other)

    def test_direct_arith_with_series_returns_not_implemented(self, data):
        # EAs should return NotImplemented for ops with Series.
        # Pandas takes care of unboxing the series and calling the EA's op.
        other = pd.Series(data)
        if hasattr(data, '__eq__'):
            result = data.__eq__(other)
            assert result is NotImplemented
        else:
            raise pytest.skip(
                "{} does not implement __eq__".format(data.__class__.__name__)
            )
