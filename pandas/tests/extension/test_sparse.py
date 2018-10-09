import pytest
import pandas as pd
import numpy as np

from pandas.core.sparse.dtype import SparseDtype
from pandas import SparseArray
from pandas.tests.extension import base
import pandas.util.testing as tm


def make_data(fill_value):
    if np.isnan(fill_value):
        data = np.random.uniform(size=100)
    else:
        data = np.random.randint(0, 100, size=100)

    data[2::3] = fill_value
    return data


@pytest.fixture
def dtype():
    return SparseDtype()


@pytest.fixture(params=[0, np.nan])
def data(request):
    """Length-100 PeriodArray for semantics test."""
    res = SparseArray(make_data(request.param),
                      fill_value=request.param)
    return res


@pytest.fixture(params=[0, np.nan])
def data_missing(request):
    """Length 2 array with [NA, Valid]"""
    return SparseArray([np.nan, 1], fill_value=request.param)


@pytest.fixture(params=[0, np.nan])
def data_repeated(request):
    """Return different versions of data for count times"""
    def gen(count):
        for _ in range(count):
            yield SparseArray(make_data(request.param),
                              fill_value=request.param)
    yield gen


@pytest.fixture(params=[0, np.nan])
def data_for_sorting(request):
    return SparseArray([2, 3, 1], fill_value=request.param)


@pytest.fixture(params=[0, np.nan])
def data_missing_for_sorting(request):
    return SparseArray([2, np.nan, 1], fill_value=request.param)


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def na_cmp():
    return lambda left, right: pd.isna(left) and pd.isna(right)


@pytest.fixture(params=[0, np.nan])
def data_for_grouping(request):
    return SparseArray([1, 1, np.nan, np.nan, 2, 2, 1, 3],
                       fill_value=request.param)


class BaseSparseTests(object):

    def _check_unsupported(self, data):
        if data.dtype == SparseDtype(int, 0):
            pytest.skip("Can't store nan in int array.")


class TestDtype(BaseSparseTests, base.BaseDtypeTests):

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is SparseArray


class TestInterface(BaseSparseTests, base.BaseInterfaceTests):
    def test_no_values_attribute(self, data):
        pytest.skip("We have values")


class TestConstructors(BaseSparseTests, base.BaseConstructorsTests):
    pass


class TestReshaping(BaseSparseTests, base.BaseReshapingTests):

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

    def test_concat_columns(self, data, na_value):
        self._check_unsupported(data)
        super(TestReshaping, self).test_concat_columns(data, na_value)

    def test_align(self, data, na_value):
        self._check_unsupported(data)
        super(TestReshaping, self).test_align(data, na_value)

    def test_align_frame(self, data, na_value):
        self._check_unsupported(data)
        super(TestReshaping, self).test_align_frame(data, na_value)

    def test_align_series_frame(self, data, na_value):
        self._check_unsupported(data)
        super(TestReshaping, self).test_align_series_frame(data, na_value)

    def test_merge(self, data, na_value):
        self._check_unsupported(data)
        super(TestReshaping, self).test_merge(data, na_value)


class TestGetitem(BaseSparseTests, base.BaseGetitemTests):

    def test_get(self, data):
        s = pd.Series(data, index=[2 * i for i in range(len(data))])
        if np.isnan(s.values.fill_value):
            assert np.isnan(s.get(4)) and np.isnan(s.iloc[2])
        else:
            assert s.get(4) == s.iloc[2]
        assert s.get(2) == s.iloc[1]

    def test_reindex(self, data, na_value):
        self._check_unsupported(data)
        super(TestGetitem, self).test_reindex(data, na_value)


# Skipping TestSetitem, since we don't implement it.

class TestMissing(BaseSparseTests, base.BaseMissingTests):
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


class TestMethods(BaseSparseTests, base.BaseMethodsTests):

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


class TestCasting(BaseSparseTests, base.BaseCastingTests):
    pass


class TestArithmeticOps(BaseSparseTests, base.BaseArithmeticOpsTests):
    series_scalar_exc = None
    frame_scalar_exc = None
    divmod_exc = None
    series_array_exc = None

    def _skip_if_different_combine(self, data):
        if data.fill_value == 0:
            # arith ops call on dtype.fill_value so that the sparsity
            # is maintained. Combine can't be called on a dtype in
            # general, so we can't make the expected. This is tested elsewhere
            raise pytest.skip("Incorrected expected from Series.combine")

    def test_error(self, data, all_arithmetic_operators):
        pass

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        self._skip_if_different_combine(data)
        super(TestArithmeticOps, self).test_arith_series_with_scalar(
            data,
            all_arithmetic_operators
        )

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        self._skip_if_different_combine(data)
        super(TestArithmeticOps, self).test_arith_series_with_array(
            data,
            all_arithmetic_operators
        )


class TestComparisonOps(BaseSparseTests, base.BaseComparisonOpsTests):

    def _compare_other(self, s, data, op_name, other):
        op = self.get_op_from_name(op_name)

        # array
        result = pd.Series(op(data, other))
        # hard to test the fill value, since we don't know what expected
        # is in general.
        # Rely on tests in `tests/sparse` to validate that.
        assert isinstance(result.dtype, SparseDtype)
        assert result.dtype.subtype == np.dtype('bool')

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
