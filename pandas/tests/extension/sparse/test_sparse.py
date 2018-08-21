import pytest
import pandas as pd
import numpy as np

from pandas.core.sparse.dtype import SparseDtype
from pandas import SparseArray
from pandas.tests.extension import base
import pandas.util.testing as tm


def make_data():
    data = np.random.uniform(size=100)
    data[2::3] = np.nan
    return data


@pytest.fixture
def dtype():
    return SparseDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    res = SparseArray(make_data())
    return res


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return SparseArray([np.nan, 1.0])


@pytest.fixture
def data_repeated():
    """Return different versions of data for count times"""
    def gen(count):
        for _ in range(count):
            yield SparseArray(make_data())
    yield gen


@pytest.fixture
def data_for_sorting():
    return SparseArray([2, 3, 1])


@pytest.fixture
def data_missing_for_sorting():
    return SparseArray([2, np.nan, 1])


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def na_cmp():
    return lambda left, right: pd.isna(left) and pd.isna(right)


@pytest.fixture
def data_for_grouping():
    return SparseArray([1, 1, np.nan, np.nan, 2, 2, 1, 3])


class TestDtype(base.BaseDtypeTests):

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is SparseArray


class TestInterface(base.BaseInterfaceTests):
    def test_no_values_attribute(self, data):
        pytest.skip("Welp")


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):

    def test_concat_mixed_dtypes(self, data):
        # https://github.com/pandas-dev/pandas/issues/20762
        # This should be the same, aside from concat([sparse, float])
        df1 = pd.DataFrame({'A': data[:3]})
        df2 = pd.DataFrame({"A": [1, 2, 3]})
        df3 = pd.DataFrame({"A": ['a', 'b', 'c']}).astype('category')
        dfs = [df1, df2, df3]

        # dataframes
        result = pd.concat(dfs)
        expected = pd.concat([x.apply(lambda s: np.asarray(s).astype(object))
                              for x in dfs])
        self.assert_frame_equal(result, expected)
        #
        # # series
        # result = pd.concat([x['A'] for x in dfs])
        # expected = pd.concat([x['A'].astype(object) for x in dfs])
        # self.assert_series_equal(result, expected)
        #
        # # simple test for just EA and one other
        # result = pd.concat([df1, df2])
        # # We can preserve float dtype here.
        # # XXX the different behavior between frame and series is bad.
        # # fix this.
        # expected = pd.concat([df1.astype(float), df2.astype(float)])
        # self.assert_frame_equal(result, expected)
        #
        # result = pd.concat([df1['A'], df2['A']])
        # expected = pd.concat([df1['A'].astype(float),
        #                       df2['A'].astype(float)])
        # self.assert_series_equal(result, expected)


class TestGetitem(base.BaseGetitemTests):

    def test_get(self, data):
        s = pd.Series(data, index=[2 * i for i in range(len(data))])
        assert np.isnan(s.get(4)) and np.isnan(s.iloc[2])
        assert s.get(2) == s.iloc[1]


# Skipping TestSetitem, since we don't implement it.

class TestMissing(base.BaseMissingTests):
    @pytest.mark.skip(reason="Unsupported")
    def test_fillna_limit_pad(self):
        pass

    @pytest.mark.skip(reason="Unsupported")
    def test_fillna_limit_backfill(self):
        pass

    @pytest.mark.skip(reason="Unsupported")
    def test_fillna_series_method(self):
        pass

    @pytest.mark.skip(reason="Unsupported")
    def test_fillna_series(self):
        # this one looks doable.
        pass

    def test_fillna_frame(self, data_missing):
        # Have to override to specify that fill_value will change.
        fill_value = data_missing[1]

        result = pd.DataFrame({
            "A": data_missing,
            "B": [1, 2]
        }).fillna(fill_value)

        if pd.isna(data_missing.fill_value):
            dtype = SparseDtype(data_missing.dtype, fill_value)
        else:
            dtype = data_missing.dtype

        expected = pd.DataFrame({
            "A": data_missing._from_sequence([fill_value, fill_value],
                                             dtype=dtype),
            "B": [1, 2],
        })

        self.assert_frame_equal(result, expected)


class TestMethods(base.BaseMethodsTests):

    def test_combine_le(self, data_repeated):
        # We return a Series[SparseArray].__le__ returns a
        # Series[Sparse[bool]]
        # rather than Series[bool]
        orig_data1, orig_data2 = data_repeated(2)
        s1 = pd.Series(orig_data1)
        s2 = pd.Series(orig_data2)
        result = s1.combine(s2, lambda x1, x2: x1 <= x2)
        expected = pd.Series(pd.SparseArray([
            a <= b for (a, b) in
            zip(list(orig_data1), list(orig_data2))
        ], fill_value=False))
        self.assert_series_equal(result, expected)

        val = s1.iloc[0]
        result = s1.combine(val, lambda x1, x2: x1 <= x2)
        expected = pd.Series(pd.SparseArray([
            a <= val for a in list(orig_data1)
        ], fill_value=False))
        self.assert_series_equal(result, expected)


class TestCasting(base.BaseCastingTests):
    pass


class TestArithmeticOps(base.BaseArithmeticOpsTests):
    series_scalar_exc = None
    frame_scalar_exc = None
    divmod_exc = None
    series_array_exc = None

    def test_error(self, data, all_arithmetic_operators):
        # not sure
        pass

    @pytest.mark.xfail(reason="TODO", strict=True)
    def test_divmod(self, data):
        super().test_divmod(data)

    @pytest.mark.xfail(reson="what is this test doing?", strict=True)
    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        super(TestArithmeticOps, self).test_arith_series_with_array(
            data, all_arithmetic_operators
        )


class TestComparisonOps(base.BaseComparisonOpsTests):

    def _compare_other(self, s, data, op_name, other):
        op = self.get_op_from_name(op_name)

        # array
        result = pd.Series(op(data, other))
        assert result.dtype == 'Sparse[bool]'

        with np.errstate(all='ignore'):
            expected = pd.Series(
                pd.SparseArray(op(np.asarray(data), np.asarray(other)),
                               fill_value=result.values.fill_value)
            )

        tm.assert_series_equal(result, expected)

        # series
        s = pd.Series(data)
        result = op(s, other)
        tm.assert_series_equal(result, expected)

    @pytest.mark.skip(reason="segfault")
    def test_compare_array(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        other = [0] * len(data)
        self._compare_other(s, data, op_name, other)


def test_slice():
    import pandas.util.testing as tm

    arr = pd.SparseArray([1, None, 2])
    result = arr[:]
    tm.assert_sp_array_equal(arr, result)
