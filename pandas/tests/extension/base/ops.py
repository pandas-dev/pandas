import pytest

import operator

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

    def check_opname(self, s, op_name, other, exc=NotImplementedError):
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

    def _check_divmod_op(self, s, op, other, exc=NotImplementedError):
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
    """Various Series and DataFrame arithmetic ops methods."""

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # series & scalar
        op_name = all_arithmetic_operators
        s = pd.Series(data)
        self.check_opname(s, op_name, s.iloc[0], exc=TypeError)

    @pytest.mark.xfail(run=False, reason="_reduce needs implementation")
    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators):
        # frame & scalar
        op_name = all_arithmetic_operators
        df = pd.DataFrame({'A': data})
        self.check_opname(df, op_name, data[0], exc=TypeError)

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        # ndarray & other series
        op_name = all_arithmetic_operators
        s = pd.Series(data)
        self.check_opname(s, op_name, pd.Series([s.iloc[0]] * len(s)),
                          exc=TypeError)

    def test_divmod(self, data):
        s = pd.Series(data)
        self._check_divmod_op(s, divmod, 1, exc=TypeError)
        self._check_divmod_op(1, ops.rdivmod, s, exc=TypeError)

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
