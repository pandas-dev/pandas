import pytest
import numpy as np
import pandas as pd
from .base import BaseExtensionTests


class BaseOpsTests(BaseExtensionTests):
    """Various Series and DataFrame ops methos."""

    def compare(self, s, op, other, exc=NotImplementedError):

        with pytest.raises(exc):
            getattr(s, op)(other)

    def test_arith_scalar(self, data, all_arithmetic_operators):
        # scalar
        op = all_arithmetic_operators
        s = pd.Series(data)
        self.compare(s, op, 1, exc=TypeError)

    def test_arith_array(self, data, all_arithmetic_operators):
        # ndarray & other series
        op = all_arithmetic_operators
        s = pd.Series(data)
        self.compare(s, op, np.ones(len(s), dtype=s.dtype.type), exc=TypeError)

    def test_compare_scalar(self, data, all_compare_operators):
        op = all_compare_operators

        s = pd.Series(data)

        if op in '__eq__':
            assert getattr(data, op)(0) is NotImplemented
            assert not getattr(s, op)(0).all()
        elif op in '__ne__':
            assert getattr(data, op)(0) is NotImplemented
            assert getattr(s, op)(0).all()

        else:

            # array
            getattr(data, op)(0) is NotImplementedError

            # series
            s = pd.Series(data)
            with pytest.raises(TypeError):
                getattr(s, op)(0)

    def test_error(self, data, all_arithmetic_operators):

        # invalid ops
        op = all_arithmetic_operators
        with pytest.raises(AttributeError):
            getattr(data, op)
