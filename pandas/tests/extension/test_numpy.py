import numpy as np
import pytest

import pandas as pd
from pandas import compat
from pandas.core.arrays.numpy_ import PandasArray, PandasDtype
import pandas.util.testing as tm

from . import base


@pytest.fixture
def dtype():
    return PandasDtype(np.dtype('float'))


@pytest.fixture
def allow_in_pandas(monkeypatch):
    """
    A monkeypatch to tells pandas to let us in.

    By default, passing a PandasArray to an index / series / frame
    constructor will unbox that PandasArray to an ndarray, and treat
    it as a non-EA column. We don't want people using EAs without
    reason.

    The mechanism for this is a check against ABCPandasArray
    in each constructor.

    But, for testing, we need to allow them in pandas. So we patch
    the _typ of PandasArray, so that we evade the ABCPandasArray
    check.
    """
    with monkeypatch.context() as m:
        m.setattr(PandasArray, '_typ', 'extension')
        yield


@pytest.fixture
def data(allow_in_pandas, dtype):
    return PandasArray(np.arange(1, 101, dtype=dtype._dtype))


@pytest.fixture
def data_missing(allow_in_pandas):
    return PandasArray(np.array([np.nan, 1.0]))


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def na_cmp():
    def cmp(a, b):
        return np.isnan(a) and np.isnan(b)
    return cmp


@pytest.fixture
def data_for_sorting(allow_in_pandas):
    """Length-3 array with a known sort order.

    This should be three items [B, C, A] with
    A < B < C
    """
    return PandasArray(
        np.array([1, 2, 0])
    )


@pytest.fixture
def data_missing_for_sorting(allow_in_pandas):
    """Length-3 array with a known sort order.

    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    return PandasArray(
        np.array([1, np.nan, 0])
    )


@pytest.fixture
def data_for_grouping(allow_in_pandas):
    """Data for factorization, grouping, and unique tests.

    Expected to be like [B, B, NA, NA, A, A, B, C]

    Where A < B < C and NA is missing
    """
    a, b, c = np.arange(3)
    return PandasArray(np.array(
        [b, b, np.nan, np.nan, a, a, b, c]
    ))


class BaseNumPyTests(object):
    pass


class TestCasting(BaseNumPyTests, base.BaseCastingTests):
    pass


class TestConstructors(BaseNumPyTests, base.BaseConstructorsTests):
    @pytest.mark.skip(reason="We don't register our dtype")
    # We don't want to register. This test should probably be split in two.
    def test_from_dtype(self, data):
        pass


class TestDtype(BaseNumPyTests, base.BaseDtypeTests):

    @pytest.mark.skip(reason="Incorrect expected.")
    # we unsurprisingly clash with a NumPy name.
    def test_check_dtype(self, data):
        pass


class TestGetitem(BaseNumPyTests, base.BaseGetitemTests):
    pass


class TestGroupby(BaseNumPyTests, base.BaseGroupbyTests):
    pass


class TestInterface(BaseNumPyTests, base.BaseInterfaceTests):
    pass


class TestMethods(BaseNumPyTests, base.BaseMethodsTests):

    @pytest.mark.skip(reason="TODO: remove?")
    def test_value_counts(self, all_data, dropna):
        pass

    @pytest.mark.skip(reason="Incorrect expected")
    # We have a bool dtype, so the result is an ExtensionArray
    # but expected is not
    def test_combine_le(self, data_repeated):
        super(TestMethods, self).test_combine_le(data_repeated)


class TestArithmetics(BaseNumPyTests, base.BaseArithmeticOpsTests):
    divmod_exc = None
    series_scalar_exc = None
    frame_scalar_exc = None
    series_array_exc = None

    def test_divmod_series_array(self, data):
        s = pd.Series(data)
        self._check_divmod_op(s, divmod, data, exc=None)

    @pytest.mark.skip("We implement ops")
    def test_error(self, data, all_arithmetic_operators):
        pass

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        if (compat.PY2 and
                all_arithmetic_operators in {'__div__', '__rdiv__'}):
            raise pytest.skip(
                "Matching NumPy int / int -> float behavior."
            )
        super(TestArithmetics, self).test_arith_series_with_scalar(
            data, all_arithmetic_operators
        )

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        if (compat.PY2 and
                all_arithmetic_operators in {'__div__', '__rdiv__'}):
            raise pytest.skip(
                "Matching NumPy int / int -> float behavior."
            )
        super(TestArithmetics, self).test_arith_series_with_array(
            data, all_arithmetic_operators
        )


class TestPrinting(BaseNumPyTests, base.BasePrintingTests):
    pass


class TestNumericReduce(BaseNumPyTests, base.BaseNumericReduceTests):

    def check_reduce(self, s, op_name, skipna):
        result = getattr(s, op_name)(skipna=skipna)
        # avoid coercing int -> float. Just cast to the actual numpy type.
        expected = getattr(s.astype(s.dtype._dtype), op_name)(skipna=skipna)
        tm.assert_almost_equal(result, expected)


class TestBooleanReduce(BaseNumPyTests, base.BaseBooleanReduceTests):
    pass


class TestMising(BaseNumPyTests, base.BaseMissingTests):
    pass


class TestReshaping(BaseNumPyTests, base.BaseReshapingTests):

    @pytest.mark.skip("Incorrect parent test")
    # not actually a mixed concat, since we concat int and int.
    def test_concat_mixed_dtypes(self, data):
        super(TestReshaping, self).test_concat_mixed_dtypes(data)


class TestSetitem(BaseNumPyTests, base.BaseSetitemTests):
    pass


class TestParsing(BaseNumPyTests, base.BaseParsingTests):
    pass
