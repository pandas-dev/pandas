from collections import OrderedDict
from collections.abc import Iterator
from datetime import (
    datetime,
    timedelta,
)

from dateutil.tz import tzoffset
import numpy as np
from numpy import ma
import pytest

from pandas._libs import (
    iNaT,
    lib,
)
from pandas.compat import HAS_PYARROW
from pandas.errors import (
    IntCastingNaNError,
    Pandas4Warning,
)

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import (
    IntegerArray,
    IntervalArray,
)
from pandas.core.internals.blocks import NumpyBlock


class TestSeriesConstructors:
    def test_from_ints_with_non_nano_dt64_dtype(self, index_or_series):
        values = np.arange(10)

        res = index_or_series(values, dtype="M8[s]")
        expected = index_or_series(values.astype("M8[s]"))
        tm.assert_equal(res, expected)

        res = index_or_series(list(values), dtype="M8[s]")
        tm.assert_equal(res, expected)

    def test_from_na_value_and_interval_of_datetime_dtype(self):
        # GH#41805
        ser = pd.Series([None], dtype="interval[datetime64[ns]]")
        assert ser.isna().all()
        assert ser.dtype == "interval[datetime64[ns], right]"

    def test_infer_with_date_and_datetime(self):
        # GH#49341 pre-2.0 we inferred datetime-and-date to datetime64, which
        #  was inconsistent with Index behavior
        ts = pd.Timestamp(2016, 1, 1)
        vals = [ts.to_pydatetime(), ts.date()]

        ser = pd.Series(vals)
        expected = pd.Series(vals, dtype=object)
        tm.assert_series_equal(ser, expected)

        idx = pd.Index(vals)
        expected = pd.Index(vals, dtype=object)
        tm.assert_index_equal(idx, expected)

    def test_unparsable_strings_with_dt64_dtype(self):
        # pre-2.0 these would be silently ignored and come back with object dtype
        vals = ["aa"]
        msg = "^Unknown datetime string format, unable to parse: aa$"
        with pytest.raises(ValueError, match=msg):
            pd.Series(vals, dtype="datetime64[ns]")

        with pytest.raises(ValueError, match=msg):
            pd.Series(np.array(vals, dtype=object), dtype="datetime64[ns]")

    def test_invalid_dtype_conversion_datetime_to_timedelta(self):
        # GH#60728
        vals = pd.Series([pd.NaT, pd.Timestamp(2025, 1, 1)], dtype="datetime64[ns]")
        msg = r"^Cannot cast DatetimeArray to dtype timedelta64\[ns\]$"
        with pytest.raises(TypeError, match=msg):
            pd.Series(vals, dtype="timedelta64[ns]")

    @pytest.mark.parametrize(
        "constructor",
        [
            # NOTE: some overlap with test_constructor_empty but that test does not
            # test for None or an empty generator.
            # test_constructor_pass_none tests None but only with the index also
            # passed.
            (lambda idx: pd.Series(index=idx)),
            (lambda idx: pd.Series(None, index=idx)),
            (lambda idx: pd.Series({}, index=idx)),
            (lambda idx: pd.Series((), index=idx)),
            (lambda idx: pd.Series([], index=idx)),
            (lambda idx: pd.Series((_ for _ in []), index=idx)),
            (lambda idx: pd.Series(data=None, index=idx)),
            (lambda idx: pd.Series(data={}, index=idx)),
            (lambda idx: pd.Series(data=(), index=idx)),
            (lambda idx: pd.Series(data=[], index=idx)),
            (lambda idx: pd.Series(data=(_ for _ in []), index=idx)),
        ],
    )
    @pytest.mark.parametrize("empty_index", [None, []])
    def test_empty_constructor(self, constructor, empty_index):
        # GH 49573 (addition of empty_index parameter)
        expected = pd.Series(index=empty_index)
        result = constructor(empty_index)

        assert result.dtype == object
        assert len(result.index) == 0
        tm.assert_series_equal(result, expected, check_index_type=True)

    def test_invalid_dtype(self):
        # GH15520
        msg = "not understood"
        invalid_list = [pd.Timestamp, "Timestamp", list]
        for dtype in invalid_list:
            with pytest.raises(TypeError, match=msg):
                pd.Series([], name="time", dtype=dtype)

    def test_invalid_compound_dtype(self):
        # GH#13296
        c_dtype = np.dtype([("a", "i8"), ("b", "f4")])
        cdt_arr = np.array([(1, 0.4), (256, -13)], dtype=c_dtype)

        with pytest.raises(ValueError, match="Use DataFrame instead"):
            pd.Series(cdt_arr, index=["A", "B"])

    def test_scalar_conversion(self):
        # Pass in scalar is disabled
        scalar = pd.Series(0.5)
        assert not isinstance(scalar, float)

    def test_scalar_extension_dtype(self, ea_scalar_and_dtype):
        # GH 28401

        ea_scalar, ea_dtype = ea_scalar_and_dtype

        ser = pd.Series(ea_scalar, index=range(3))
        expected = pd.Series([ea_scalar] * 3, dtype=ea_dtype)

        assert ser.dtype == ea_dtype
        tm.assert_series_equal(ser, expected)

    def test_constructor(self, datetime_series, using_infer_string):
        empty_series = pd.Series()
        assert datetime_series.index._is_all_dates

        # Pass in Series
        derived = pd.Series(datetime_series)
        assert derived.index._is_all_dates

        tm.assert_index_equal(derived.index, datetime_series.index)
        # Ensure new index is not created
        assert id(datetime_series.index) == id(derived.index)

        # Mixed type Series
        mixed = pd.Series(["hello", np.nan], index=[0, 1])
        assert mixed.dtype == np.object_ if not using_infer_string else "str"
        assert np.isnan(mixed[1])

        assert not empty_series.index._is_all_dates
        assert not pd.Series().index._is_all_dates

        # exception raised is of type ValueError GH35744
        with pytest.raises(
            ValueError,
            match=r"Data must be 1-dimensional, got ndarray of shape \(3, 3\) instead",
        ):
            pd.Series(
                np.random.default_rng(2).standard_normal((3, 3)), index=np.arange(3)
            )

        mixed.name = "Series"
        rs = pd.Series(mixed).name
        xp = "Series"
        assert rs == xp

        # raise on MultiIndex GH4187
        m = pd.MultiIndex.from_arrays([[1, 2], [3, 4]])
        msg = "initializing a Series from a MultiIndex is not supported"
        with pytest.raises(NotImplementedError, match=msg):
            pd.Series(m)

    def test_constructor_index_ndim_gt_1_raises(self):
        # GH#18579
        df = pd.DataFrame([[1, 2], [3, 4], [5, 6]], index=[3, 6, 9])
        with pytest.raises(ValueError, match="Index data must be 1-dimensional"):
            pd.Series([1, 3, 2], index=df)

    @pytest.mark.parametrize("input_class", [list, dict, OrderedDict])
    def test_constructor_empty(self, input_class, using_infer_string):
        empty = pd.Series()
        empty2 = pd.Series(input_class())

        # these are Index() and RangeIndex() which don't compare type equal
        # but are just .equals
        tm.assert_series_equal(empty, empty2, check_index_type=False)

        # With explicit dtype:
        empty = pd.Series(dtype="float64")
        empty2 = pd.Series(input_class(), dtype="float64")
        tm.assert_series_equal(empty, empty2, check_index_type=False)

        # GH 18515 : with dtype=category:
        empty = pd.Series(dtype="category")
        empty2 = pd.Series(input_class(), dtype="category")
        tm.assert_series_equal(empty, empty2, check_index_type=False)

        if input_class is not list:
            # With index:
            empty = pd.Series(index=range(10))
            empty2 = pd.Series(input_class(), index=range(10))
            tm.assert_series_equal(empty, empty2)

            # With index and dtype float64:
            empty = pd.Series(np.nan, index=range(10))
            empty2 = pd.Series(input_class(), index=range(10), dtype="float64")
            tm.assert_series_equal(empty, empty2)

            # GH 19853 : with empty string, index and dtype str
            empty = pd.Series("", dtype=str, index=range(3))
            if using_infer_string:
                empty2 = pd.Series("", index=range(3), dtype="str")
            else:
                empty2 = pd.Series("", index=range(3))
            tm.assert_series_equal(empty, empty2)

    @pytest.mark.parametrize("input_arg", [np.nan, float("nan")])
    def test_constructor_nan(self, input_arg):
        empty = pd.Series(dtype="float64", index=range(10))
        empty2 = pd.Series(input_arg, index=range(10))

        tm.assert_series_equal(empty, empty2, check_index_type=False)

    @pytest.mark.parametrize(
        "dtype",
        ["f8", "i8", "M8[ns]", "m8[ns]", "category", "object", "datetime64[ns, UTC]"],
    )
    @pytest.mark.parametrize("index", [None, pd.Index([])])
    def test_constructor_dtype_only(self, dtype, index):
        # GH-20865
        result = pd.Series(dtype=dtype, index=index)
        assert result.dtype == dtype
        assert len(result) == 0

    def test_constructor_no_data_index_order(self):
        result = pd.Series(index=["b", "a", "c"])
        assert result.index.tolist() == ["b", "a", "c"]

    def test_constructor_no_data_string_type(self):
        # GH 22477
        result = pd.Series(index=[1], dtype=str)
        assert np.isnan(result.iloc[0])

    @pytest.mark.parametrize("item", ["entry", "ѐ", 13])
    def test_constructor_string_element_string_type(self, item):
        # GH 22477
        result = pd.Series(item, index=[1], dtype=str)
        assert result.iloc[0] == str(item)

    def test_constructor_dtype_str_na_values(self, string_dtype):
        # https://github.com/pandas-dev/pandas/issues/21083
        ser = pd.Series(["x", None], dtype=string_dtype)
        result = ser.isna()
        expected = pd.Series([False, True])
        tm.assert_series_equal(result, expected)
        assert ser.iloc[1] is None

        ser = pd.Series(["x", np.nan], dtype=string_dtype)
        assert np.isnan(ser.iloc[1])

    def test_constructor_series(self):
        index1 = ["d", "b", "a", "c"]
        index2 = sorted(index1)
        s1 = pd.Series([4, 7, -5, 3], index=index1)
        s2 = pd.Series(s1, index=index2)

        tm.assert_series_equal(s2, s1.sort_index())

    def test_constructor_iterable(self):
        # GH 21987
        class Iter:
            def __iter__(self) -> Iterator:
                yield from range(10)

        expected = pd.Series(list(range(10)), dtype="int64")
        result = pd.Series(Iter(), dtype="int64")
        tm.assert_series_equal(result, expected)

    def test_constructor_sequence(self):
        # GH 21987
        expected = pd.Series(list(range(10)), dtype="int64")
        result = pd.Series(range(10), dtype="int64")
        tm.assert_series_equal(result, expected)

    def test_constructor_single_str(self):
        # GH 21987
        expected = pd.Series(["abc"])
        result = pd.Series("abc")
        tm.assert_series_equal(result, expected)

    def test_constructor_list_like(self):
        # make sure that we are coercing different
        # list-likes to standard dtypes and not
        # platform specific
        expected = pd.Series([1, 2, 3], dtype="int64")
        for obj in [[1, 2, 3], (1, 2, 3), np.array([1, 2, 3], dtype="int64")]:
            result = pd.Series(obj, index=[0, 1, 2])
            tm.assert_series_equal(result, expected)

    def test_constructor_boolean_index(self):
        # GH#18579
        s1 = pd.Series([1, 2, 3], index=[4, 5, 6])

        index = s1 == 2
        result = pd.Series([1, 3, 2], index=index)
        expected = pd.Series([1, 3, 2], index=[False, True, False])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["bool", "int32", "int64", "float64"])
    def test_constructor_index_dtype(self, dtype):
        # GH 17088

        s = pd.Series(pd.Index([0, 2, 4]), dtype=dtype)
        assert s.dtype == dtype

    @pytest.mark.parametrize(
        "input_vals",
        [
            [1, 2],
            ["1", "2"],
            list(pd.date_range("1/1/2011", periods=2, freq="h")),
            list(pd.date_range("1/1/2011", periods=2, freq="h", tz="US/Eastern")),
            [pd.Interval(left=0, right=5)],
        ],
    )
    def test_constructor_list_str(self, input_vals, string_dtype):
        # GH 16605
        # Ensure that data elements from a list are converted to strings
        # when dtype is str, 'str', or 'U'
        result = pd.Series(input_vals, dtype=string_dtype)
        expected = pd.Series(input_vals).astype(string_dtype)
        tm.assert_series_equal(result, expected)

    def test_constructor_list_str_na(self, string_dtype):
        result = pd.Series([1.0, 2.0, np.nan], dtype=string_dtype)
        expected = pd.Series(["1.0", "2.0", np.nan], dtype=object)
        tm.assert_series_equal(result, expected)
        assert np.isnan(result[2])

    def test_constructor_generator(self):
        gen = (i for i in range(10))

        result = pd.Series(gen)
        exp = pd.Series(range(10))
        tm.assert_series_equal(result, exp)

        # same but with non-default index
        gen = (i for i in range(10))
        result = pd.Series(gen, index=range(10, 20))
        exp.index = range(10, 20)
        tm.assert_series_equal(result, exp)

    def test_constructor_map(self):
        # GH8909
        m = (x for x in range(10))

        result = pd.Series(m)
        exp = pd.Series(range(10))
        tm.assert_series_equal(result, exp)

        # same but with non-default index
        m = (x for x in range(10))
        result = pd.Series(m, index=range(10, 20))
        exp.index = range(10, 20)
        tm.assert_series_equal(result, exp)

    def test_constructor_categorical(self):
        msg = "Constructing a Categorical with a dtype and values containing"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            cat = pd.Categorical([0, 1, 2, 0, 1, 2], ["a", "b", "c"])
        res = pd.Series(cat)
        tm.assert_categorical_equal(res.values, cat)

        # can cast to a new dtype
        result = pd.Series(pd.Categorical([1, 2, 3]), dtype="int64")
        expected = pd.Series([1, 2, 3], dtype="int64")
        tm.assert_series_equal(result, expected)

    def test_construct_from_categorical_with_dtype(self):
        # GH12574
        ser = pd.Series(pd.Categorical([1, 2, 3]), dtype="category")
        assert isinstance(ser.dtype, CategoricalDtype)

    def test_construct_intlist_values_category_dtype(self):
        ser = pd.Series([1, 2, 3], dtype="category")
        assert isinstance(ser.dtype, CategoricalDtype)

    def test_constructor_categorical_with_coercion(self):
        factor = pd.Categorical(["a", "b", "b", "a", "a", "c", "c", "c"])
        # test basic creation / coercion of categoricals
        s = pd.Series(factor, name="A")
        assert s.dtype == "category"
        assert len(s) == len(factor)

        # in a frame
        df = pd.DataFrame({"A": factor})
        result = df["A"]
        tm.assert_series_equal(result, s)
        result = df.iloc[:, 0]
        tm.assert_series_equal(result, s)
        assert len(df) == len(factor)

        df = pd.DataFrame({"A": s})
        result = df["A"]
        tm.assert_series_equal(result, s)
        assert len(df) == len(factor)

        # multiples
        df = pd.DataFrame({"A": s, "B": s, "C": 1})
        result1 = df["A"]
        result2 = df["B"]
        tm.assert_series_equal(result1, s)
        tm.assert_series_equal(result2, s, check_names=False)
        assert result2.name == "B"
        assert len(df) == len(factor)

    def test_constructor_categorical_with_coercion2(self):
        # GH8623
        x = pd.DataFrame(
            [[1, "John P. Doe"], [2, "Jane Dove"], [1, "John P. Doe"]],
            columns=["person_id", "person_name"],
        )
        x["person_name"] = pd.Categorical(x.person_name)  # doing this breaks transform

        expected = x.iloc[0].person_name
        result = x.person_name.iloc[0]
        assert result == expected

        result = x.person_name[0]
        assert result == expected

        result = x.person_name.loc[0]
        assert result == expected

    def test_constructor_series_to_categorical(self):
        # see GH#16524: test conversion of Series to Categorical
        series = pd.Series(["a", "b", "c"])

        result = pd.Series(series, dtype="category")
        expected = pd.Series(["a", "b", "c"], dtype="category")

        tm.assert_series_equal(result, expected)

    def test_constructor_categorical_dtype(self):
        result = pd.Series(
            ["a", "b"], dtype=CategoricalDtype(["a", "b", "c"], ordered=True)
        )
        assert isinstance(result.dtype, CategoricalDtype)
        tm.assert_index_equal(result.cat.categories, pd.Index(["a", "b", "c"]))
        assert result.cat.ordered

        result = pd.Series(["a", "b"], dtype=CategoricalDtype(["b", "a"]))
        assert isinstance(result.dtype, CategoricalDtype)
        tm.assert_index_equal(result.cat.categories, pd.Index(["b", "a"]))
        assert result.cat.ordered is False

        # GH 19565 - Check broadcasting of scalar with Categorical dtype
        result = pd.Series(
            "a", index=[0, 1], dtype=CategoricalDtype(["a", "b"], ordered=True)
        )
        expected = pd.Series(
            ["a", "a"], index=[0, 1], dtype=CategoricalDtype(["a", "b"], ordered=True)
        )
        tm.assert_series_equal(result, expected)

    def test_constructor_categorical_string(self):
        # GH 26336: the string 'category' maintains existing CategoricalDtype
        cdt = CategoricalDtype(categories=list("dabc"), ordered=True)
        expected = pd.Series(list("abcabc"), dtype=cdt)

        # Series(Categorical, dtype='category') keeps existing dtype
        cat = pd.Categorical(list("abcabc"), dtype=cdt)
        result = pd.Series(cat, dtype="category")
        tm.assert_series_equal(result, expected)

        # Series(Series[Categorical], dtype='category') keeps existing dtype
        result = pd.Series(result, dtype="category")
        tm.assert_series_equal(result, expected)

    def test_categorical_sideeffects_free(self):
        # Passing a categorical to a Series and then changing values in either
        # the series or the categorical should not change the values in the
        # other one, IF you specify copy!
        cat = pd.Categorical(["a", "b", "c", "a"])
        s = pd.Series(cat, copy=True)
        assert s.cat is not cat
        s = s.cat.rename_categories([1, 2, 3])
        exp_s = np.array([1, 2, 3, 1], dtype=np.int64)
        exp_cat = np.array(["a", "b", "c", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(s.__array__(), exp_s)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # setting
        s[0] = 2
        exp_s2 = np.array([2, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s2)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # however, copy is False by default
        # so this WILL change values
        cat = pd.Categorical(["a", "b", "c", "a"])
        s = pd.Series(cat, copy=False)
        assert s._values is cat
        s = s.cat.rename_categories([1, 2, 3])
        assert s._values is not cat
        exp_s = np.array([1, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s)

        s[0] = 2
        exp_s2 = np.array([2, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s2)

    def test_unordered_compare_equal(self):
        left = pd.Series(["a", "b", None], dtype=CategoricalDtype(["a", "b"]))
        right = pd.Series(pd.Categorical(["a", "b", np.nan], categories=["a", "b"]))
        tm.assert_series_equal(left, right)

    def test_constructor_maskedarray(self):
        data = ma.masked_all((3,), dtype=float)
        result = pd.Series(data)
        expected = pd.Series([np.nan, np.nan, np.nan])
        tm.assert_series_equal(result, expected)

        data[0] = 0.0
        data[2] = 2.0
        index = ["a", "b", "c"]
        result = pd.Series(data, index=index)
        expected = pd.Series([0.0, np.nan, 2.0], index=index)
        tm.assert_series_equal(result, expected)

        data[1] = 1.0
        result = pd.Series(data, index=index)
        expected = pd.Series([0.0, 1.0, 2.0], index=index)
        tm.assert_series_equal(result, expected)

        data = ma.masked_all((3,), dtype=int)
        result = pd.Series(data)
        expected = pd.Series([np.nan, np.nan, np.nan], dtype=float)
        tm.assert_series_equal(result, expected)

        data[0] = 0
        data[2] = 2
        index = ["a", "b", "c"]
        result = pd.Series(data, index=index)
        expected = pd.Series([0, np.nan, 2], index=index, dtype=float)
        tm.assert_series_equal(result, expected)

        data[1] = 1
        result = pd.Series(data, index=index)
        expected = pd.Series([0, 1, 2], index=index, dtype=int)
        tm.assert_series_equal(result, expected)

        data = ma.masked_all((3,), dtype=bool)
        result = pd.Series(data)
        expected = pd.Series([np.nan, np.nan, np.nan], dtype=object)
        tm.assert_series_equal(result, expected)

        data[0] = True
        data[2] = False
        index = ["a", "b", "c"]
        result = pd.Series(data, index=index)
        expected = pd.Series([True, np.nan, False], index=index, dtype=object)
        tm.assert_series_equal(result, expected)

        data[1] = True
        result = pd.Series(data, index=index)
        expected = pd.Series([True, True, False], index=index, dtype=bool)
        tm.assert_series_equal(result, expected)

        data = ma.masked_all((3,), dtype="M8[ns]")
        result = pd.Series(data)
        expected = pd.Series([iNaT, iNaT, iNaT], dtype="M8[ns]")
        tm.assert_series_equal(result, expected)

        data[0] = datetime(2001, 1, 1)
        data[2] = datetime(2001, 1, 3)
        index = ["a", "b", "c"]
        result = pd.Series(data, index=index)
        expected = pd.Series(
            [datetime(2001, 1, 1), iNaT, datetime(2001, 1, 3)],
            index=index,
            dtype="M8[ns]",
        )
        tm.assert_series_equal(result, expected)

        data[1] = datetime(2001, 1, 2)
        result = pd.Series(data, index=index)
        expected = pd.Series(
            [datetime(2001, 1, 1), datetime(2001, 1, 2), datetime(2001, 1, 3)],
            index=index,
            dtype="M8[ns]",
        )
        tm.assert_series_equal(result, expected)

    def test_constructor_maskedarray_hardened(self):
        # Check numpy masked arrays with hard masks -- from GH24574
        data = ma.masked_all((3,), dtype=float).harden_mask()
        result = pd.Series(data)
        expected = pd.Series([np.nan, np.nan, np.nan])
        tm.assert_series_equal(result, expected)

    def test_series_ctor_plus_datetimeindex(self):
        rng = pd.date_range("20090415", "20090519", freq="B")
        data = dict.fromkeys(rng, 1)

        result = pd.Series(data, index=rng)
        assert result.index.is_(rng)

    def test_constructor_default_index(self):
        s = pd.Series([0, 1, 2])
        tm.assert_index_equal(s.index, pd.Index(range(3)), exact=True)

    @pytest.mark.parametrize(
        "input",
        [
            [1, 2, 3],
            (1, 2, 3),
            list(range(3)),
            pd.Categorical(["a", "b", "a"]),
            (i for i in range(3)),
            (x for x in range(3)),
        ],
    )
    def test_constructor_index_mismatch(self, input):
        # GH 19342
        # test that construction of a Series with an index of different length
        # raises an error
        msg = r"Length of values \(3\) does not match length of index \(4\)"
        with pytest.raises(ValueError, match=msg):
            pd.Series(input, index=np.arange(4))

    def test_constructor_numpy_scalar(self):
        # GH 19342
        # construction with a numpy scalar
        # should not raise
        result = pd.Series(np.array(100), index=np.arange(4), dtype="int64")
        expected = pd.Series(100, index=np.arange(4), dtype="int64")
        tm.assert_series_equal(result, expected)

    def test_constructor_broadcast_list(self):
        # GH 19342
        # construction with single-element container and index
        # should raise
        msg = r"Length of values \(1\) does not match length of index \(3\)"
        with pytest.raises(ValueError, match=msg):
            pd.Series(["foo"], index=["a", "b", "c"])

    def test_constructor_corner(self):
        df = pd.DataFrame(range(5), index=pd.date_range("2020-01-01", periods=5))
        objs = [df, df]
        s = pd.Series(objs, index=[0, 1])
        assert isinstance(s, pd.Series)

    def test_constructor_sanitize(self):
        s = pd.Series(np.array([1.0, 1.0, 8.0]), dtype="i8")
        assert s.dtype == np.dtype("i8")

        msg = r"Cannot convert non-finite values \(NA or inf\) to integer"
        with pytest.raises(IntCastingNaNError, match=msg):
            pd.Series(np.array([1.0, 1.0, np.nan]), copy=True, dtype="i8")

    def test_constructor_copy(self):
        # GH15125
        # test dtype parameter has no side effects on copy=True
        for data in [[1.0], np.array([1.0])]:
            x = pd.Series(data)
            y = pd.Series(x, copy=True, dtype=float)

            # copy=True maintains original data in Series
            tm.assert_series_equal(x, y)

            # changes to origin of copy does not affect the copy
            x[0] = 2.0
            assert not x.equals(y)
            assert x[0] == 2.0
            assert y[0] == 1.0

    @pytest.mark.parametrize(
        "index",
        [
            pd.date_range("20170101", periods=3, tz="US/Eastern"),
            pd.date_range("20170101", periods=3),
            pd.timedelta_range("1 day", periods=3),
            pd.period_range("2012Q1", periods=3, freq="Q"),
            pd.Index(list("abc")),
            pd.Index([1, 2, 3]),
            pd.RangeIndex(0, 3),
        ],
        ids=lambda x: type(x).__name__,
    )
    def test_constructor_limit_copies(self, index):
        # GH 17449
        # limit copies of input
        s = pd.Series(index)

        # we make 1 copy; this is just a smoke test here
        assert s._mgr.blocks[0].values is not index

    def test_constructor_shallow_copy(self):
        # constructing a Series from Series with copy=False should still
        # give a "shallow" copy (share data, not attributes)
        # https://github.com/pandas-dev/pandas/issues/49523
        s = pd.Series([1, 2, 3])
        s_orig = s.copy()
        s2 = pd.Series(s)
        assert s2._mgr is not s._mgr
        # Overwriting index of s2 doesn't change s
        s2.index = ["a", "b", "c"]
        tm.assert_series_equal(s, s_orig)

    def test_constructor_pass_none(self):
        s = pd.Series(None, index=range(5))
        assert s.dtype == np.float64

        s = pd.Series(None, index=range(5), dtype=object)
        assert s.dtype == np.object_

        # GH 7431
        # inference on the index
        s = pd.Series(index=np.array([None]))
        expected = pd.Series(index=pd.Index([None]))
        tm.assert_series_equal(s, expected)

    def test_constructor_pass_nan_nat(self):
        # GH 13467
        exp = pd.Series([np.nan, np.nan], dtype=np.float64)
        assert exp.dtype == np.float64
        tm.assert_series_equal(pd.Series([np.nan, np.nan]), exp)
        tm.assert_series_equal(pd.Series(np.array([np.nan, np.nan])), exp)

        exp = pd.Series([pd.NaT, pd.NaT])
        assert exp.dtype == "datetime64[s]"
        tm.assert_series_equal(pd.Series([pd.NaT, pd.NaT]), exp)
        tm.assert_series_equal(pd.Series(np.array([pd.NaT, pd.NaT])), exp)

        tm.assert_series_equal(pd.Series([pd.NaT, np.nan]), exp)
        tm.assert_series_equal(pd.Series(np.array([pd.NaT, np.nan])), exp)

        tm.assert_series_equal(pd.Series([np.nan, pd.NaT]), exp)
        tm.assert_series_equal(pd.Series(np.array([np.nan, pd.NaT])), exp)

    def test_constructor_cast(self):
        msg = "could not convert string to float"
        with pytest.raises(ValueError, match=msg):
            pd.Series(["a", "b", "c"], dtype=float)

    def test_constructor_signed_int_overflow_raises(self):
        # GH#41734 disallow silent overflow, enforced in 2.0
        msg = "The elements provided in the data cannot all be casted to the dtype"
        with pytest.raises(OverflowError, match=msg):
            pd.Series([1, 200, 923442], dtype="int8")

        with pytest.raises(OverflowError, match=msg):
            pd.Series([1, 200, 923442], dtype="uint8")

    @pytest.mark.parametrize("reverse", [False, True])
    def test_constructor_int_above_float64_range(self, reverse):
        # GH#66519 an int too big for any numeric dtype infers object dtype
        #  wherever it sits. A list is inferred with convert_numeric=True and
        #  used to raise OverflowError in either order; an object array stops at
        #  the first integer, so it only raised with the big value first.
        huge = 10**400
        data = [1, huge] if reverse else [huge, 1]
        arr = np.array(data, dtype=object)

        assert pd.Series(data).dtype == object
        assert pd.Index(data).dtype == object
        assert pd.Series(arr).dtype == object
        assert pd.Index(arr).dtype == object

    @pytest.mark.parametrize("reverse", [False, True])
    def test_constructor_uint64_int64_conflict_after_none(self, reverse):
        # GH#66519 a None before either value used to hide the conflict,
        #  giving a float64 that rounds 2**64 - 1 up to 2**64
        data = [-1, 2**64 - 1] if reverse else [2**64 - 1, -1]

        result = pd.Series([None, *data])
        assert result.dtype == object
        assert result[1] == data[0]

    @pytest.mark.parametrize(
        "values",
        [
            np.array([1], dtype=np.uint16),
            np.array([1], dtype=np.uint32),
            np.array([1], dtype=np.uint64),
            [np.uint16(1)],
            [np.uint32(1)],
            [np.uint64(1)],
        ],
    )
    def test_constructor_numpy_uints(self, values):
        # GH#47294
        value = values[0]
        result = pd.Series(values)

        assert result[0].dtype == value.dtype
        assert result[0] == value

    def test_constructor_unsigned_dtype_overflow(self, any_unsigned_int_numpy_dtype):
        # see gh-15832
        msg = (
            f"The elements provided in the data cannot "
            f"all be casted to the dtype {any_unsigned_int_numpy_dtype}"
        )
        with pytest.raises(OverflowError, match=msg):
            pd.Series([-1], dtype=any_unsigned_int_numpy_dtype)

    def test_constructor_floating_data_int_dtype(self, frame_or_series):
        # GH#40110
        arr = np.random.default_rng(2).standard_normal(2)

        # Long-standing behavior (for Series, new in 2.0 for DataFrame)
        #  has been to ignore the dtype on these;
        #  not clear if this is what we want long-term
        # expected = frame_or_series(arr)

        # GH#49599 as of 2.0 we raise instead of silently retaining float dtype
        msg = "Trying to coerce float values to integer"
        with pytest.raises(ValueError, match=msg):
            frame_or_series(arr, dtype="i8")

        with pytest.raises(ValueError, match=msg):
            frame_or_series(list(arr), dtype="i8")

        # pre-2.0, when we had NaNs, we silently ignored the integer dtype
        arr[0] = np.nan
        # expected = frame_or_series(arr)

        msg = r"Cannot convert non-finite values \(NA or inf\) to integer"
        with pytest.raises(IntCastingNaNError, match=msg):
            frame_or_series(arr, dtype="i8")

        exc = IntCastingNaNError
        if frame_or_series is pd.Series:
            # TODO: try to align these
            exc = ValueError
            msg = "cannot convert float NaN to integer"
        with pytest.raises(exc, match=msg):
            # same behavior if we pass list instead of the ndarray
            frame_or_series(list(arr), dtype="i8")

        # float array that can be losslessly cast to integers
        arr = np.array([1.0, 2.0], dtype="float64")
        expected = frame_or_series(arr.astype("i8"))

        obj = frame_or_series(arr, dtype="i8")
        tm.assert_equal(obj, expected)

        obj = frame_or_series(list(arr), dtype="i8")
        tm.assert_equal(obj, expected)

    def test_constructor_coerce_float_fail(self, any_int_numpy_dtype):
        # see gh-15832
        # Updated: make sure we treat this list the same as we would treat
        #  the equivalent ndarray
        # GH#49599 pre-2.0 we silently retained float dtype, in 2.0 we raise
        vals = [1, 2, 3.5]

        msg = "Trying to coerce float values to integer"
        with pytest.raises(ValueError, match=msg):
            pd.Series(vals, dtype=any_int_numpy_dtype)
        with pytest.raises(ValueError, match=msg):
            pd.Series(np.array(vals), dtype=any_int_numpy_dtype)

    def test_constructor_coerce_float_valid(self, float_numpy_dtype):
        s = pd.Series([1, 2, 3.5], dtype=float_numpy_dtype)
        expected = pd.Series([1, 2, 3.5]).astype(float_numpy_dtype)
        tm.assert_series_equal(s, expected)

    def test_constructor_invalid_coerce_ints_with_float_nan(self, any_int_numpy_dtype):
        # GH 22585
        # Updated: make sure we treat this list the same as we would treat the
        # equivalent ndarray
        vals = [1, 2, np.nan]
        # pre-2.0 this would return with a float dtype, in 2.0 we raise

        msg = "cannot convert float NaN to integer"
        with pytest.raises(ValueError, match=msg):
            pd.Series(vals, dtype=any_int_numpy_dtype)
        msg = r"Cannot convert non-finite values \(NA or inf\) to integer"
        with pytest.raises(IntCastingNaNError, match=msg):
            pd.Series(np.array(vals), dtype=any_int_numpy_dtype)

    def test_constructor_dtype_no_cast(self):
        # see gh-1572
        s = pd.Series([1, 2, 3])
        s2 = pd.Series(s, dtype=np.int64)

        s2[1] = 5
        assert s[1] == 2

    def test_constructor_datelike_coercion(self):
        # GH 9477
        # incorrectly inferring on dateimelike looking when object dtype is
        # specified
        s = pd.Series([pd.Timestamp("20130101"), "NOV"], dtype=object)
        assert s.iloc[0] == pd.Timestamp("20130101")
        assert s.iloc[1] == "NOV"
        assert s.dtype == object

    def test_constructor_datelike_coercion2(self):
        # the dtype was being reset on the slicing and re-inferred to datetime
        # even thought the blocks are mixed
        belly = "216 3T19".split()
        wing1 = "2T15 4H19".split()
        wing2 = "416 4T20".split()
        mat = pd.to_datetime("2016-01-22 2019-09-07".split())
        df = pd.DataFrame({"wing1": wing1, "wing2": wing2, "mat": mat}, index=belly)

        result = df.loc["3T19"]
        assert result.dtype == object
        result = df.loc["216"]
        assert result.dtype == object

    def test_constructor_mixed_int_and_timestamp(self, frame_or_series):
        # specifically Timestamp with nanos, not datetimes
        objs = [pd.Timestamp(9), 10, pd.NaT._value]
        result = frame_or_series(objs, dtype="M8[ns]")

        expected = frame_or_series([pd.Timestamp(9), pd.Timestamp(10), pd.NaT])
        tm.assert_equal(result, expected)

    def test_constructor_datetimes_with_nulls(self):
        # gh-15869
        for arr in [
            np.array([None, None, None, None, datetime(2011, 1, 1), None]),
            np.array([None, None, datetime(2011, 1, 1), None]),
        ]:
            result = pd.Series(arr)
            assert result.dtype == "M8[us]"

    def test_constructor_dtype_datetime64(self):
        s = pd.Series(iNaT, dtype="M8[ns]", index=range(5))
        assert pd.isna(s).all()

        # in theory this should be all nulls, but since
        # we are not specifying a dtype is ambiguous
        s = pd.Series(iNaT, index=range(5))
        assert not pd.isna(s).all()

        s = pd.Series(np.nan, dtype="M8[ns]", index=range(5))
        assert pd.isna(s).all()

        s = pd.Series([datetime(2001, 1, 2, 0, 0), iNaT], dtype="M8[ns]")
        assert pd.isna(s[1])
        assert s.dtype == "M8[ns]"

        s = pd.Series([datetime(2001, 1, 2, 0, 0), np.nan], dtype="M8[ns]")
        assert pd.isna(s[1])
        assert s.dtype == "M8[ns]"

    def test_constructor_dtype_datetime64_10(self):
        # GH3416
        pydates = [datetime(2013, 1, 1), datetime(2013, 1, 2), datetime(2013, 1, 3)]
        dates = [np.datetime64(x) for x in pydates]

        ser = pd.Series(dates)
        assert ser.dtype == "M8[us]"

        ser.iloc[0] = np.nan
        assert ser.dtype == "M8[us]"

        # GH3414 related
        expected = pd.Series(pydates, dtype="datetime64[ms]")

        result = pd.Series(pd.Series(dates).astype(np.int64) / 1000, dtype="M8[ms]")
        tm.assert_series_equal(result, expected)

        result = pd.Series(dates, dtype="datetime64[ms]")
        tm.assert_series_equal(result, expected)

        expected = pd.Series(
            [pd.NaT, datetime(2013, 1, 2), datetime(2013, 1, 3)], dtype="datetime64[ns]"
        )
        result = pd.Series([np.nan, *dates[1:]], dtype="datetime64[ns]")
        tm.assert_series_equal(result, expected)

    def test_constructor_dtype_datetime64_11(self):
        pydates = [datetime(2013, 1, 1), datetime(2013, 1, 2), datetime(2013, 1, 3)]
        dates = [np.datetime64(x) for x in pydates]

        dts = pd.Series(dates, dtype="datetime64[ns]")

        # valid astype
        dts.astype("int64")

        # invalid casting
        msg = r"Converting from datetime64\[ns\] to int32 is not supported"
        with pytest.raises(TypeError, match=msg):
            dts.astype("int32")

        # ints are ok
        # we test with np.int64 to get similar results on
        # windows / 32-bit platforms
        result = pd.Series(dts, dtype=np.int64)
        expected = pd.Series(dts.astype(np.int64))
        tm.assert_series_equal(result, expected)

    def test_constructor_dtype_datetime64_9(self):
        # invalid dates can be help as object
        result = pd.Series([datetime(2, 1, 1)])
        assert result[0] == datetime(2, 1, 1, 0, 0)

        result = pd.Series([datetime(3000, 1, 1)])
        assert result[0] == datetime(3000, 1, 1, 0, 0)

    def test_constructor_dtype_datetime64_8(self):
        # don't mix types
        result = pd.Series([pd.Timestamp("20130101"), 1], index=["a", "b"])
        assert result["a"] == pd.Timestamp("20130101")
        assert result["b"] == 1

    def test_constructor_dtype_datetime64_7(self):
        # GH6529
        # coerce datetime64 non-ns properly
        dates = pd.date_range("01-Jan-2015", "01-Dec-2015", freq="ME")
        values2 = dates.view(np.ndarray).astype("datetime64[ns]")
        expected = pd.Series(values2, index=dates)

        for unit in ["s", "D", "ms", "us", "ns"]:
            dtype = np.dtype(f"M8[{unit}]")
            values1 = dates.view(np.ndarray).astype(dtype)
            result = pd.Series(values1, dates)
            if unit == "D":
                # for unit="D" we cast to nearest-supported reso, i.e. "s"
                dtype = np.dtype("M8[s]")
            assert result.dtype == dtype
            tm.assert_series_equal(result, expected.astype(dtype))

        # GH 13876
        # coerce to non-ns to object properly
        expected = pd.Series(values2, index=dates, dtype=object)
        for dtype in ["s", "D", "ms", "us", "ns"]:
            values1 = dates.view(np.ndarray).astype(f"M8[{dtype}]")
            result = pd.Series(values1, index=dates, dtype=object)
            tm.assert_series_equal(result, expected)

        # leave datetime.date alone
        dates2 = np.array([d.date() for d in dates.to_pydatetime()], dtype=object)
        series1 = pd.Series(dates2, dates)
        tm.assert_numpy_array_equal(series1.values, dates2)
        assert series1.dtype == object

    def test_constructor_dtype_datetime64_6(self):
        # as of 2.0, these no longer infer datetime64 based on the strings,
        #  matching the Index behavior

        ser = pd.Series([None, pd.NaT, "2013-08-05 15:30:00.000001"])
        assert ser.dtype == object

        ser = pd.Series([np.nan, pd.NaT, "2013-08-05 15:30:00.000001"])
        assert ser.dtype == object

        ser = pd.Series([pd.NaT, None, "2013-08-05 15:30:00.000001"])
        assert ser.dtype == object

        ser = pd.Series([pd.NaT, np.nan, "2013-08-05 15:30:00.000001"])
        assert ser.dtype == object

    def test_constructor_dtype_datetime64_5(self):
        # tz-aware (UTC and other tz's)
        # GH 8411
        dr = pd.date_range("20130101", periods=3)
        assert pd.Series(dr).iloc[0].tz is None
        dr = pd.date_range("20130101", periods=3, tz="UTC")
        assert str(pd.Series(dr).iloc[0].tz) == "UTC"
        dr = pd.date_range("20130101", periods=3, tz="US/Eastern")
        assert str(pd.Series(dr).iloc[0].tz) == "US/Eastern"

    def test_constructor_dtype_datetime64_4(self):
        # non-convertible
        ser = pd.Series([1479596223000, -1479590, pd.NaT])
        assert ser.dtype == "object"
        assert ser[2] is pd.NaT
        assert "NaT" in str(ser)

    def test_constructor_dtype_datetime64_3(self):
        # if we passed a NaT it remains
        ser = pd.Series([datetime(2010, 1, 1), datetime(2, 1, 1), pd.NaT])
        assert ser.dtype == "M8[us]"
        assert ser[2] is pd.NaT
        assert "NaT" in str(ser)

    def test_constructor_dtype_datetime64_2(self):
        # if we passed a nan it remains
        ser = pd.Series([datetime(2010, 1, 1), datetime(2, 1, 1), np.nan])
        assert ser.dtype == "M8[us]"
        assert ser[2] is pd.NaT
        assert "NaT" in str(ser)

    def test_constructor_with_datetime_tz(self):
        # 8260
        # support datetime64 with tz

        dr = pd.date_range("20130101", periods=3, tz="US/Eastern", unit="ns")
        s = pd.Series(dr)
        assert s.dtype.name == "datetime64[ns, US/Eastern]"
        assert s.dtype == "datetime64[ns, US/Eastern]"
        assert isinstance(s.dtype, pd.DatetimeTZDtype)
        assert "datetime64[ns, US/Eastern]" in str(s)

        # export
        depr_msg = "Series.values returning an ndarray that drops timezone information"
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            result = s.values
        assert isinstance(result, np.ndarray)
        assert result.dtype == "datetime64[ns]"

        exp = pd.DatetimeIndex(result)
        exp = exp.tz_localize("UTC").tz_convert(tz=s.dt.tz)
        tm.assert_index_equal(dr, exp, check_freq=False)

        # indexing
        result = s.iloc[0]
        assert result == pd.Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern")
        result = s[0]
        assert result == pd.Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern")

        result = s[pd.Series([True, True, False], index=s.index)]
        tm.assert_series_equal(result, s[0:2])

        result = s.iloc[0:1]
        tm.assert_series_equal(result, pd.Series(dr[0:1]))

        # concat
        result = pd.concat([s.iloc[0:1], s.iloc[1:]])
        tm.assert_series_equal(result, s)

        # short str
        assert "datetime64[ns, US/Eastern]" in str(s)

        # formatting with NaT
        result = s.shift()
        assert "datetime64[ns, US/Eastern]" in str(result)
        assert "NaT" in str(result)

        result = pd.DatetimeIndex(s, freq="infer")
        tm.assert_index_equal(result, dr)

    def test_constructor_with_datetime_tz5(self):
        # long str
        ser = pd.Series(
            pd.date_range("20130101", periods=1000, tz="US/Eastern", unit="ns")
        )
        assert "datetime64[ns, US/Eastern]" in str(ser)

    def test_constructor_with_datetime_tz4(self):
        # inference
        ser = pd.Series(
            [
                pd.Timestamp("2013-01-01 13:00:00-0800", tz="US/Pacific").as_unit("s"),
                pd.Timestamp("2013-01-02 14:00:00-0800", tz="US/Pacific").as_unit("s"),
            ]
        )
        assert ser.dtype == "datetime64[s, US/Pacific]"
        assert lib.infer_dtype(ser, skipna=True) == "datetime64"

    def test_constructor_with_datetime_tz3(self):
        ser = pd.Series(
            [
                pd.Timestamp("2013-01-01 13:00:00-0800", tz="US/Pacific"),
                pd.Timestamp("2013-01-02 14:00:00-0800", tz="US/Eastern"),
            ]
        )
        assert ser.dtype == "object"
        assert lib.infer_dtype(ser, skipna=True) == "datetime"

    def test_constructor_with_datetime_tz2(self):
        # with all NaT
        ser = pd.Series(pd.NaT, index=[0, 1], dtype="datetime64[ns, US/Eastern]")
        dti = pd.DatetimeIndex(["NaT", "NaT"], tz="US/Eastern").as_unit("ns")
        expected = pd.Series(dti)
        tm.assert_series_equal(ser, expected)

    def test_constructor_no_partial_datetime_casting(self):
        # GH#40111
        vals = [
            "nan",
            pd.Timestamp("1990-01-01"),
            "2015-03-14T16:15:14.123-08:00",
            "2019-03-04T21:56:32.620-07:00",
            None,
        ]
        ser = pd.Series(vals)
        assert all(ser[i] is vals[i] for i in range(len(vals)))

    @pytest.mark.parametrize("arg", ["2013-01-01 00:00:00", pd.NaT, np.nan, None])
    def test_constructor_with_naive_string_and_datetimetz_dtype(self, arg):
        # GH 17415: With naive string
        result = pd.Series([arg], dtype="datetime64[ns, CET]")
        expected = pd.Series([pd.Timestamp(arg)], dtype="M8[ns]").dt.tz_localize("CET")
        tm.assert_series_equal(result, expected)

    def test_constructor_datetime64_bigendian(self):
        # GH#30976
        ms = np.datetime64(1, "ms")
        arr = np.array([np.datetime64(1, "ms")], dtype=">M8[ms]")

        result = pd.Series(arr)
        expected = pd.Series([pd.Timestamp(ms)]).astype("M8[ms]")
        assert expected.dtype == "M8[ms]"
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("interval_constructor", [pd.IntervalIndex, IntervalArray])
    def test_construction_interval(self, interval_constructor):
        # construction from interval & array of intervals
        intervals = interval_constructor.from_breaks(np.arange(3), closed="right")
        result = pd.Series(intervals)
        expected_subtype = np.dtype(np.intp)
        assert result.dtype == f"interval[{expected_subtype}, right]"
        msg = "Series.values returning an object-dtype ndarray for IntervalDtype"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            tm.assert_index_equal(pd.Index(result.values), pd.Index(intervals))

    @pytest.mark.parametrize(
        "data_constructor", [list, np.array], ids=["list", "ndarray[object]"]
    )
    def test_constructor_infer_interval(self, data_constructor):
        # GH 23563: consistent closed results in interval dtype
        data = [pd.Interval(0, 1), pd.Interval(0, 2), None]
        result = pd.Series(data_constructor(data))
        expected = pd.Series(IntervalArray(data))
        assert result.dtype == "interval[float64, right]"
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "data_constructor", [list, np.array], ids=["list", "ndarray[object]"]
    )
    def test_constructor_interval_mixed_closed(self, data_constructor):
        # GH 23563: mixed closed results in object dtype (not interval dtype)
        data = [pd.Interval(0, 1, closed="both"), pd.Interval(0, 2, closed="neither")]
        result = pd.Series(data_constructor(data))
        assert result.dtype == object
        assert result.tolist() == data

    def test_construction_consistency(self):
        # make sure that we are not re-localizing upon construction
        # GH 14928
        ser = pd.Series(pd.date_range("20130101", periods=3, tz="US/Eastern"))

        result = pd.Series(ser, dtype=ser.dtype)
        tm.assert_series_equal(result, ser)

        result = pd.Series(ser.dt.tz_convert("UTC"), dtype=ser.dtype)
        tm.assert_series_equal(result, ser)

        depr_msg = "Series.values returning an ndarray that drops timezone information"

        # Pre-2.0 dt64 values were treated as utc, which was inconsistent
        #  with DatetimeIndex, which treats them as wall times, see GH#33401
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            result = pd.Series(ser.values, dtype=ser.dtype)
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            expected = pd.Series(ser.values).dt.tz_localize(ser.dtype.tz)
        tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            # one suggested alternative to the deprecated (changed in 2.0) usage
            middle = pd.Series(ser.values).dt.tz_localize("UTC")
            result = middle.dt.tz_convert(ser.dtype.tz)
        tm.assert_series_equal(result, ser)

        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            # the other suggested alternative to the deprecated usage
            result = pd.Series(ser.values.view("int64"), dtype=ser.dtype)
        tm.assert_series_equal(result, ser)

    @pytest.mark.parametrize(
        "data_constructor", [list, np.array], ids=["list", "ndarray[object]"]
    )
    def test_constructor_infer_period(self, data_constructor):
        data = [pd.Period("2000", "D"), pd.Period("2001", "D"), None]
        result = pd.Series(data_constructor(data))
        expected = pd.Series(pd.PeriodIndex(data))
        tm.assert_series_equal(result, expected)
        assert result.dtype == "Period[D]"

    def test_construct_from_ints_including_iNaT_scalar_period_dtype(self):
        # GH#64158 the integers are read as calendar years, so they have to be
        #  years Period(int, freq) accepts
        msg = "Passing integer data"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            series = pd.Series([1000, 2000, 3000, pd._libs.iNaT], dtype="period[D]")

        val = series[3]
        assert pd.isna(val)

        series[2] = val
        assert pd.isna(series[2])

    def test_constructor_period_incompatible_frequency(self):
        data = [pd.Period("2000", "D"), pd.Period("2001", "Y")]
        result = pd.Series(data)
        assert result.dtype == object
        assert result.tolist() == data

    def test_constructor_periodindex(self):
        # GH7932
        # converting a PeriodIndex when put in a Series

        pi = pd.period_range("20130101", periods=5, freq="D")
        s = pd.Series(pi)
        assert s.dtype == "Period[D]"
        expected = pd.Series(pi.astype(object))
        assert expected.dtype == object

    def test_constructor_dict(self):
        d = {"a": 0.0, "b": 1.0, "c": 2.0}

        result = pd.Series(d)
        expected = pd.Series(d, index=sorted(d.keys()))
        tm.assert_series_equal(result, expected)

        result = pd.Series(d, index=["b", "c", "d", "a"])
        expected = pd.Series([1, 2, np.nan, 0], index=["b", "c", "d", "a"])
        tm.assert_series_equal(result, expected)

        pidx = pd.period_range("2020-01-01", periods=10, freq="D")
        d = {pidx[0]: 0, pidx[1]: 1}
        result = pd.Series(d, index=pidx)
        expected = pd.Series(np.nan, pidx, dtype=np.float64)
        expected.iloc[0] = 0
        expected.iloc[1] = 1
        tm.assert_series_equal(result, expected)

    def test_constructor_dict_list_value_explicit_dtype(self):
        # GH 18625
        d = {"a": [[2], [3], [4]]}
        result = pd.Series(d, index=["a"], dtype="object")
        expected = pd.Series(d, index=["a"])
        tm.assert_series_equal(result, expected)

    def test_constructor_dict_order(self):
        # GH19018
        # initialization ordering: by insertion order
        d = {"b": 1, "a": 0, "c": 2}
        result = pd.Series(d)
        expected = pd.Series([1, 0, 2], index=list("bac"))
        tm.assert_series_equal(result, expected)

    def test_constructor_dict_extension(self, ea_scalar_and_dtype):
        ea_scalar, ea_dtype = ea_scalar_and_dtype
        d = {"a": ea_scalar}
        result = pd.Series(d, index=["a"])
        expected = pd.Series(ea_scalar, index=["a"], dtype=ea_dtype)

        assert result.dtype == ea_dtype

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("value", [2, np.nan, None, float("nan")])
    def test_constructor_dict_nan_key(self, value):
        # GH 18480
        d = {1: "a", value: "b", float("nan"): "c", 4: "d"}
        result = pd.Series(d).sort_values()
        expected = pd.Series(["a", "b", "c", "d"], index=[1, value, np.nan, 4])
        tm.assert_series_equal(result, expected)

        # MultiIndex:
        d = {(1, 1): "a", (2, np.nan): "b", (3, value): "c"}
        result = pd.Series(d).sort_values()
        expected = pd.Series(
            ["a", "b", "c"], index=pd.Index([(1, 1), (2, np.nan), (3, value)])
        )
        tm.assert_series_equal(result, expected)

    def test_constructor_dict_datetime64_index(self):
        # GH 9456

        dates_as_str = ["1984-02-19", "1988-11-06", "1989-12-03", "1990-03-15"]
        values = [42544017.198965244, 1234565, 40512335.181958228, -1]

        def create_data(constructor):
            return dict(
                zip((constructor(x) for x in dates_as_str), values, strict=True)
            )

        data_datetime64 = create_data(np.datetime64)
        data_datetime = create_data(lambda x: datetime.strptime(x, "%Y-%m-%d"))
        data_Timestamp = create_data(pd.Timestamp)

        expected = pd.Series(values, (pd.Timestamp(x) for x in dates_as_str))

        result_datetime64 = pd.Series(data_datetime64)
        result_datetime = pd.Series(data_datetime)
        result_Timestamp = pd.Series(data_Timestamp)

        tm.assert_series_equal(
            result_datetime64, expected.set_axis(expected.index.as_unit("s"))
        )
        tm.assert_series_equal(result_datetime, expected)
        tm.assert_series_equal(result_Timestamp, expected)

    def test_constructor_dict_tuple_indexer(self):
        # GH 12948
        data = {(1, 1, None): -1.0}
        result = pd.Series(data)
        expected = pd.Series(
            -1.0,
            index=pd.MultiIndex(levels=[[1], [1], [np.nan]], codes=[[0], [0], [-1]]),
        )
        tm.assert_series_equal(result, expected)

    def test_constructor_mapping(self, non_dict_mapping_subclass):
        # GH 29788
        ndm = non_dict_mapping_subclass({3: "three"})
        result = pd.Series(ndm)
        expected = pd.Series(["three"], index=[3])

        tm.assert_series_equal(result, expected)

    def test_constructor_list_of_tuples(self):
        data = [(1, 1), (2, 2), (2, 3)]
        s = pd.Series(data)
        assert list(s) == data

    def test_constructor_tuple_of_tuples(self):
        data = ((1, 1), (2, 2), (2, 3))
        s = pd.Series(data)
        assert tuple(s) == data

    @pytest.mark.parametrize(
        "data, expected_values, expected_index",
        [
            ({(1, 2): 3, (None, 5): 6}, [3, 6], [(1, 2), (None, 5)]),
            ({(1,): 3, (4, 5): 6}, [3, 6], [(1, None), (4, 5)]),
        ],
    )
    def test_constructor_dict_of_tuples(self, data, expected_values, expected_index):
        # GH 60695
        result = pd.Series(data).sort_values()
        expected = pd.Series(
            expected_values, index=pd.MultiIndex.from_tuples(expected_index)
        )
        tm.assert_series_equal(result, expected)

    # https://github.com/pandas-dev/pandas/issues/22698
    def test_fromDict(self, using_infer_string):
        data = {"a": 0, "b": 1, "c": 2, "d": 3}

        series = pd.Series(data)
        tm.assert_is_sorted(series.index)

        data = {"a": 0, "b": "1", "c": "2", "d": datetime(2011, 1, 1)}
        series = pd.Series(data)
        assert series.dtype == np.object_

        data = {"a": 0, "b": "1", "c": "2", "d": "3"}
        series = pd.Series(data)
        assert series.dtype == np.object_ if not using_infer_string else "str"

        data = {"a": "0", "b": "1"}
        series = pd.Series(data, dtype=float)
        assert series.dtype == np.float64

    def test_fromValue(self, datetime_series, using_infer_string):
        nans = pd.Series(np.nan, index=datetime_series.index, dtype=np.float64)
        assert nans.dtype == np.float64
        assert len(nans) == len(datetime_series)

        strings = pd.Series("foo", index=datetime_series.index)
        assert strings.dtype == np.object_ if not using_infer_string else "str"
        assert len(strings) == len(datetime_series)

        d = datetime(2011, 1, 1)
        dates = pd.Series(d, index=datetime_series.index)
        assert dates.dtype == "M8[us]"
        assert len(dates) == len(datetime_series)

        # GH12336
        # Test construction of categorical series from value
        categorical = pd.Series(0, index=datetime_series.index, dtype="category")
        expected = pd.Series(0, index=datetime_series.index).astype("category")
        assert categorical.dtype == "category"
        assert len(categorical) == len(datetime_series)
        tm.assert_series_equal(categorical, expected)

    def test_constructor_dtype_timedelta64(self):
        # basic
        td = pd.Series([timedelta(days=i) for i in range(3)])
        assert td.dtype == "timedelta64[us]"

        td = pd.Series([timedelta(days=1)])
        assert td.dtype == "timedelta64[us]"

        td = pd.Series([timedelta(days=1), timedelta(days=2), np.timedelta64(1, "s")])

        assert td.dtype == "timedelta64[us]"

        # mixed with NaT
        td = pd.Series([timedelta(days=1), pd.NaT], dtype="m8[ns]")
        assert td.dtype == "timedelta64[ns]"

        td = pd.Series([timedelta(days=1), np.nan], dtype="m8[ns]")
        assert td.dtype == "timedelta64[ns]"

        td = pd.Series([np.timedelta64(300000000, "ns"), pd.NaT], dtype="m8[ns]")
        assert td.dtype == "timedelta64[ns]"

        # improved inference
        # GH5689
        td = pd.Series([np.timedelta64(300000000, "ns"), pd.NaT])
        assert td.dtype == "timedelta64[ns]"

        # because iNaT is int, not coerced to timedelta
        td = pd.Series([np.timedelta64(300000000, "ns"), iNaT])
        assert td.dtype == "object"

        td = pd.Series([np.timedelta64(300000000, "ns"), np.nan])
        assert td.dtype == "timedelta64[ns]"

        td = pd.Series([pd.NaT, np.timedelta64(300000000, "ns")])
        assert td.dtype == "timedelta64[ns]"

        td = pd.Series([np.timedelta64(1, "s")])
        assert td.dtype == "timedelta64[s]"

        # valid astype
        td.astype("int64")

        # invalid casting
        msg = r"Converting from timedelta64\[s\] to int32 is not supported"
        with pytest.raises(TypeError, match=msg):
            td.astype("int32")

        # this is an invalid casting
        msg = "|".join(
            [
                "Could not convert object to NumPy timedelta",
                "Could not convert 'foo' to NumPy timedelta",
            ]
        )
        with pytest.raises(ValueError, match=msg):
            pd.Series([timedelta(days=1), "foo"], dtype="m8[ns]")

        # leave as object here
        td = pd.Series([timedelta(days=i) for i in range(3)] + ["foo"])
        assert td.dtype == "object"

        # as of 2.0, these no longer infer timedelta64 based on the strings,
        #  matching Index behavior
        ser = pd.Series([None, pd.NaT, "1 Day"])
        assert ser.dtype == object

        ser = pd.Series([np.nan, pd.NaT, "1 Day"])
        assert ser.dtype == object

        ser = pd.Series([pd.NaT, None, "1 Day"])
        assert ser.dtype == object

        ser = pd.Series([pd.NaT, np.nan, "1 Day"])
        assert ser.dtype == object

    # GH 16406
    def test_constructor_mixed_tz(self):
        s = pd.Series(
            [pd.Timestamp("20130101"), pd.Timestamp("20130101", tz="US/Eastern")]
        )
        expected = pd.Series(
            [pd.Timestamp("20130101"), pd.Timestamp("20130101", tz="US/Eastern")],
            dtype="object",
        )
        tm.assert_series_equal(s, expected)

    def test_NaT_scalar(self):
        series = pd.Series([0, 1000, 2000, iNaT], dtype="M8[ns]")

        val = series[3]
        assert pd.isna(val)

        series[2] = val
        assert pd.isna(series[2])

    def test_NaT_cast(self):
        # GH10747
        result = pd.Series([np.nan]).astype("M8[ns]")
        expected = pd.Series([pd.NaT], dtype="M8[ns]")
        tm.assert_series_equal(result, expected)

    def test_constructor_name_hashable(self):
        for n in [777, 777.0, "name", datetime(2001, 11, 11), (1,), "\u05d0"]:
            for data in [[1, 2, 3], np.ones(3), {"a": 0, "b": 1}]:
                s = pd.Series(data, name=n)
                assert s.name == n

    def test_constructor_name_unhashable(self):
        msg = r"Series\.name must be a hashable type"
        for n in [["name_list"], np.ones(2), {1: 2}]:
            for data in [["name_list"], np.ones(2), {1: 2}]:
                with pytest.raises(TypeError, match=msg):
                    pd.Series(data, name=n)

    def test_auto_conversion(self):
        series = pd.Series(list(pd.date_range("1/1/2000", periods=10, unit="ns")))
        assert series.dtype == "M8[ns]"

    def test_convert_non_ns(self):
        # convert from a numpy array of non-ns timedelta64
        arr = np.array([1, 2, 3], dtype="timedelta64[s]")
        ser = pd.Series(arr)
        assert ser.dtype == arr.dtype

        tdi = pd.timedelta_range("00:00:01", periods=3, freq="s").as_unit("s")
        expected = pd.Series(tdi)
        assert expected.dtype == arr.dtype
        tm.assert_series_equal(ser, expected)

        # convert from a numpy array of non-ns datetime64
        arr = np.array(
            ["2013-01-01", "2013-01-02", "2013-01-03"], dtype="datetime64[D]"
        )
        ser = pd.Series(arr)
        expected = pd.Series(
            pd.date_range("20130101", periods=3, freq="D"), dtype="M8[s]"
        )
        assert expected.dtype == "M8[s]"
        tm.assert_series_equal(ser, expected)

        arr = np.array(
            ["2013-01-01 00:00:01", "2013-01-01 00:00:02", "2013-01-01 00:00:03"],
            dtype="datetime64[s]",
        )
        ser = pd.Series(arr)
        expected = pd.Series(
            pd.date_range("20130101 00:00:01", periods=3, freq="s"), dtype="M8[s]"
        )
        assert expected.dtype == "M8[s]"
        tm.assert_series_equal(ser, expected)

    @pytest.mark.parametrize(
        "index",
        [
            pd.date_range("1/1/2000", periods=10),
            pd.timedelta_range("1 day", periods=10),
            pd.period_range("2000-Q1", periods=10, freq="Q"),
        ],
        ids=lambda x: type(x).__name__,
    )
    def test_constructor_cant_cast_datetimelike(self, index):
        # floats are not ok
        # strip Index to convert PeriodIndex -> Period
        # We don't care whether the error message says
        # PeriodIndex or PeriodArray
        msg = f"Cannot cast {type(index).__name__.rstrip('Index')}.*? to "

        with pytest.raises(TypeError, match=msg):
            pd.Series(index, dtype=float)

        # ints are ok
        # we test with np.int64 to get similar results on
        # windows / 32-bit platforms
        result = pd.Series(index, dtype=np.int64)
        expected = pd.Series(index.astype(np.int64))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "index",
        [
            pd.date_range("1/1/2000", periods=10),
            pd.timedelta_range("1 day", periods=10),
            pd.period_range("2000-Q1", periods=10, freq="Q"),
        ],
        ids=lambda x: type(x).__name__,
    )
    def test_constructor_cast_object(self, index):
        s = pd.Series(index, dtype=object)
        exp = pd.Series(index).astype(object)
        tm.assert_series_equal(s, exp)

        s = pd.Series(pd.Index(index, dtype=object), dtype=object)
        exp = pd.Series(index).astype(object)
        tm.assert_series_equal(s, exp)

        s = pd.Series(index.astype(object), dtype=object)
        exp = pd.Series(index).astype(object)
        tm.assert_series_equal(s, exp)

    @pytest.mark.parametrize("dtype", [np.datetime64, np.timedelta64])
    def test_constructor_generic_timestamp_no_frequency(self, dtype, request):
        # see gh-15524, gh-15987
        msg = "dtype has no unit. Please pass in"

        if np.dtype(dtype).name not in ["timedelta64", "datetime64"]:
            mark = pytest.mark.xfail(reason="GH#33890 Is assigned ns unit")
            request.applymarker(mark)

        with pytest.raises(ValueError, match=msg):
            pd.Series([], dtype=dtype)

    @pytest.mark.parametrize("unit", ["ps", "as", "fs", "Y", "M", "W", "D", "h", "m"])
    @pytest.mark.parametrize("kind", ["m", "M"])
    def test_constructor_generic_timestamp_bad_frequency(self, kind, unit):
        # see gh-15524, gh-15987
        # as of 2.0 we raise on any non-supported unit rather than silently
        #  cast to nanos; previously we only raised for frequencies higher
        #  than ns
        dtype = f"{kind}8[{unit}]"

        msg = "dtype=.* is not supported. Supported resolutions are"
        with pytest.raises(TypeError, match=msg):
            pd.Series([], dtype=dtype)

        with pytest.raises(TypeError, match=msg):
            # pre-2.0 the DataFrame cast raised but the Series case did not
            pd.DataFrame([[0]], dtype=dtype)

    @pytest.mark.parametrize("dtype", [None, "uint8", "category"])
    def test_constructor_range_dtype(self, dtype):
        # GH 16804
        expected = pd.Series([0, 1, 2, 3, 4], dtype=dtype or "int64")
        result = pd.Series(range(5), dtype=dtype)
        tm.assert_series_equal(result, expected)

    def test_constructor_range_overflows(self):
        # GH#30173 range objects that overflow int64
        rng = range(2**63, 2**63 + 4)
        ser = pd.Series(rng)
        expected = pd.Series(list(rng))
        tm.assert_series_equal(ser, expected)
        assert list(ser) == list(rng)
        assert ser.dtype == np.uint64

        rng2 = range(2**63 + 4, 2**63, -1)
        ser2 = pd.Series(rng2)
        expected2 = pd.Series(list(rng2))
        tm.assert_series_equal(ser2, expected2)
        assert list(ser2) == list(rng2)
        assert ser2.dtype == np.uint64

        rng3 = range(-(2**63), -(2**63) - 4, -1)
        ser3 = pd.Series(rng3)
        expected3 = pd.Series(list(rng3))
        tm.assert_series_equal(ser3, expected3)
        assert list(ser3) == list(rng3)
        assert ser3.dtype == object

        rng4 = range(2**73, 2**73 + 4)
        ser4 = pd.Series(rng4)
        expected4 = pd.Series(list(rng4))
        tm.assert_series_equal(ser4, expected4)
        assert list(ser4) == list(rng4)
        assert ser4.dtype == object

    def test_constructor_tz_mixed_data(self):
        # GH 13051
        dt_list = [
            pd.Timestamp("2016-05-01 02:03:37"),
            pd.Timestamp("2016-04-30 19:03:37-0700", tz="US/Pacific"),
        ]
        result = pd.Series(dt_list)
        expected = pd.Series(dt_list, dtype=object)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("pydt", [True, False])
    def test_constructor_data_aware_dtype_naive(self, tz_aware_fixture, pydt):
        # GH#25843, GH#41555, GH#33401
        tz = tz_aware_fixture
        ts = pd.Timestamp("2019", tz=tz)
        if pydt:
            ts = ts.to_pydatetime()

        msg = (
            "Cannot convert timezone-aware data to timezone-naive dtype. "
            r"Use pd.Series\(values\).dt.tz_localize\(None\) instead."
        )
        with pytest.raises(ValueError, match=msg):
            pd.Series([ts], dtype="datetime64[ns]")

        with pytest.raises(ValueError, match=msg):
            pd.Series(np.array([ts], dtype=object), dtype="datetime64[ns]")

        with pytest.raises(ValueError, match=msg):
            pd.Series({0: ts}, dtype="datetime64[ns]")

        msg = "Cannot unbox tzaware Timestamp to tznaive dtype"
        with pytest.raises(TypeError, match=msg):
            pd.Series(ts, index=[0], dtype="datetime64[ns]")

    def test_constructor_datetime64(self):
        rng = pd.date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        dates = np.asarray(rng)

        series = pd.Series(dates)
        assert np.issubdtype(series.dtype, np.dtype("M8[ns]"))

    def test_constructor_datetimelike_scalar_to_string_dtype(
        self, nullable_string_dtype
    ):
        # https://github.com/pandas-dev/pandas/pull/33846
        result = pd.Series("M", index=[1, 2, 3], dtype=nullable_string_dtype)
        expected = pd.Series(
            ["M", "M", "M"], index=[1, 2, 3], dtype=nullable_string_dtype
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("box", [lambda x: x, np.datetime64])
    def test_constructor_sparse_datetime64(self, box):
        # https://github.com/pandas-dev/pandas/issues/35762
        values = [box("2012-01-01"), box("2013-01-01")]
        dtype = pd.SparseDtype("datetime64[ns]")
        result = pd.Series(values, dtype=dtype)
        arr = pd.arrays.SparseArray(values, dtype=dtype)
        expected = pd.Series(arr)
        tm.assert_series_equal(result, expected)

    def test_construction_from_ordered_collection(self):
        # https://github.com/pandas-dev/pandas/issues/36044
        result = pd.Series({"a": 1, "b": 2}.keys())
        expected = pd.Series(["a", "b"])
        tm.assert_series_equal(result, expected)

        result = pd.Series({"a": 1, "b": 2}.values())
        expected = pd.Series([1, 2])
        tm.assert_series_equal(result, expected)

    def test_construction_from_large_int_scalar_no_overflow(self):
        # https://github.com/pandas-dev/pandas/issues/36291
        n = 1_000_000_000_000_000_000_000
        result = pd.Series(n, index=[0])
        expected = pd.Series(n)
        tm.assert_series_equal(result, expected)

    def test_constructor_list_of_periods_infers_period_dtype(self):
        series = pd.Series(list(pd.period_range("2000-01-01", periods=10, freq="D")))
        assert series.dtype == "Period[D]"

        series = pd.Series(
            [pd.Period("2011-01-01", freq="D"), pd.Period("2011-02-01", freq="D")]
        )
        assert series.dtype == "Period[D]"

    def test_constructor_subclass_dict(self, dict_subclass):
        data = dict_subclass((x, 10.0 * x) for x in range(10))
        series = pd.Series(data)
        expected = pd.Series(dict(data.items()))
        tm.assert_series_equal(series, expected)

    def test_constructor_ordereddict(self):
        # GH3283
        data = OrderedDict(
            (f"col{i}", np.random.default_rng(2).random()) for i in range(12)
        )

        series = pd.Series(data)
        expected = pd.Series(list(data.values()), list(data.keys()))
        tm.assert_series_equal(series, expected)

        # Test with subclass
        class A(OrderedDict):
            pass

        series = pd.Series(A(data))
        tm.assert_series_equal(series, expected)

    @pytest.mark.parametrize(
        "data, expected_index_multi",
        [
            ({("a", "a"): 0.0, ("b", "a"): 1.0, ("b", "c"): 2.0}, True),
            ({("a",): 0.0, ("a", "b"): 1.0}, True),
            ({"z": 111.0, ("a", "a"): 0.0, ("b", "a"): 1.0, ("b", "c"): 2.0}, False),
        ],
    )
    def test_constructor_dict_multiindex(self, data, expected_index_multi):
        # GH#60695
        result = pd.Series(data)

        if expected_index_multi:
            expected = pd.Series(
                list(data.values()),
                index=pd.MultiIndex.from_tuples(list(data.keys())),
            )
            tm.assert_series_equal(result, expected)
        else:
            expected = pd.Series(
                list(data.values()),
                index=pd.Index(list(data.keys())),
            )
            tm.assert_series_equal(result, expected)

    def test_constructor_dict_multiindex_reindex_flat(self):
        # construction involves reindexing with a MultiIndex corner case
        data = {("i", "i"): 0, ("i", "j"): 1, ("j", "i"): 2, "j": np.nan}
        expected = pd.Series(data)

        result = pd.Series(expected[:-1].to_dict(), index=expected.index)
        tm.assert_series_equal(result, expected)

    def test_constructor_dict_timedelta_index(self):
        # GH #12169 : Resample category data with timedelta index
        # construct Series from dict as data and TimedeltaIndex as index
        # will result NaN in result Series data
        expected = pd.Series(
            data=["A", "B", "C"], index=pd.to_timedelta([0, 10, 20], unit="s")
        )

        result = pd.Series(
            data={
                pd.to_timedelta(0, unit="s"): "A",
                pd.to_timedelta(10, unit="s"): "B",
                pd.to_timedelta(20, unit="s"): "C",
            },
            index=pd.to_timedelta([0, 10, 20], unit="s"),
        )
        tm.assert_series_equal(result, expected)

    def test_constructor_infer_index_tz(self):
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [
            datetime(2012, 5, 11, 11, tzinfo=tzinfo),
            datetime(2012, 5, 11, 12, tzinfo=tzinfo),
        ]
        series = pd.Series(data=values, index=index)

        assert series.index.tz == tzinfo

        # it works! GH#2443
        repr(series.index[0])

    def test_constructor_with_pandas_dtype(self):
        # going through 2D->1D path
        vals = [(1,), (2,), (3,)]
        ser = pd.Series(vals)
        dtype = ser.array.dtype  # NumpyEADtype
        ser2 = pd.Series(vals, dtype=dtype)
        tm.assert_series_equal(ser, ser2)

    def test_constructor_int_dtype_missing_values(self):
        # GH#43017
        result = pd.Series(index=[0], dtype="int64")
        expected = pd.Series(np.nan, index=[0], dtype="float64")
        tm.assert_series_equal(result, expected)

    def test_constructor_bool_dtype_missing_values(self):
        # GH#43018
        result = pd.Series(index=[0], dtype="bool")
        expected = pd.Series(True, index=[0], dtype="bool")
        tm.assert_series_equal(result, expected)

    def test_constructor_int64_dtype(self, any_int_dtype):
        # GH#44923
        result = pd.Series(["0", "1", "2"], dtype=any_int_dtype)
        expected = pd.Series([0, 1, 2], dtype=any_int_dtype)
        tm.assert_series_equal(result, expected)

    def test_constructor_raise_on_lossy_conversion_of_strings(self):
        # GH#44923
        with pytest.raises(OverflowError, match="The elements provided in the data"):
            pd.Series(["128"], dtype="int8")

    def test_constructor_dtype_timedelta_alternative_construct(self):
        # GH#35465
        result = pd.Series([1000000, 200000, 3000000], dtype="timedelta64[us]")
        expected = pd.Series(pd.to_timedelta([1000000, 200000, 3000000], unit="us"))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "data", [[1.5, 2.5, 90.0], [1, 2.5, 90], [1.5, pd.NaT, 90]]
    )
    def test_constructor_dtype_timedelta_float_honors_unit(self, data):
        # GH#63499 float and mixed int/float data should be interpreted in
        #  the dtype's unit, matching the integer case, instead of as
        #  nanoseconds (which silently truncated these to zero)
        result = pd.Series(data, dtype="timedelta64[s]")
        expected = pd.Series(pd.to_timedelta(data, unit="s").as_unit("s"))
        tm.assert_series_equal(result, expected)

    def test_constructor_dtype_timedelta_ns_s_astype_int64(self):
        # GH#35465
        result = pd.Series([1000000, 200000, 3000000], dtype="timedelta64[ns]").astype(
            "int64"
        )
        expected = pd.Series([1000000, 200000, 3000000], dtype="timedelta64[s]").astype(
            "int64"
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("func", [pd.Series, pd.DataFrame, pd.Index, pd.array])
    def test_constructor_mismatched_null_nullable_dtype(
        self, func, any_numeric_ea_dtype
    ):
        # GH#44514
        # the non-numeric NA is rejected before any cast is attempted
        msg = "'values' contains non-numeric NA"

        for null in [*tm.NP_NAT_OBJECTS, pd.NaT]:
            with pytest.raises(TypeError, match=msg):
                func([null, 1.0, 3.0], dtype=any_numeric_ea_dtype)

    def test_series_constructor_ea_int_from_bool(self):
        # GH#42137
        result = pd.Series([True, False, True, pd.NA], dtype="Int64")
        expected = pd.Series([1, 0, 1, pd.NA], dtype="Int64")
        tm.assert_series_equal(result, expected)

        result = pd.Series([True, False, True], dtype="Int64")
        expected = pd.Series([1, 0, 1], dtype="Int64")
        tm.assert_series_equal(result, expected)

    def test_series_constructor_ea_int_from_string_bool(self):
        # GH#42137
        with pytest.raises(ValueError, match="invalid literal"):
            pd.Series(["True", "False", "True", pd.NA], dtype="Int64")

    @pytest.mark.parametrize("val", [1, 1.0])
    def test_series_constructor_overflow_uint_ea(self, val):
        # GH#38798
        max_val = np.iinfo(np.uint64).max - 1
        result = pd.Series([max_val, val], dtype="UInt64")
        expected = pd.Series(np.array([max_val, 1], dtype="uint64"), dtype="UInt64")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("val", [1, 1.0])
    def test_series_constructor_overflow_uint_ea_with_na(self, val):
        # GH#38798
        max_val = np.iinfo(np.uint64).max - 1
        result = pd.Series([max_val, val, pd.NA], dtype="UInt64")
        expected = pd.Series(
            IntegerArray(
                np.array([max_val, 1, 0], dtype="uint64"),
                np.array([0, 0, 1], dtype=np.bool_),
            )
        )
        tm.assert_series_equal(result, expected)

    def test_series_constructor_overflow_uint_with_nan(self):
        # GH#38798
        max_val = np.iinfo(np.uint64).max - 1
        result = pd.Series([max_val, pd.NA], dtype="UInt64")
        expected = pd.Series(
            IntegerArray(
                np.array([max_val, 1], dtype="uint64"),
                np.array([0, 1], dtype=np.bool_),
            )
        )
        tm.assert_series_equal(result, expected)

    def test_series_constructor_ea_all_na(self):
        # GH#38798
        result = pd.Series([pd.NA, pd.NA], dtype="UInt64")
        expected = pd.Series(
            IntegerArray(
                np.array([1, 1], dtype="uint64"),
                np.array([1, 1], dtype=np.bool_),
            )
        )
        tm.assert_series_equal(result, expected)

    def test_series_from_index_dtype_equal_does_not_copy(self):
        # GH#52008
        idx = pd.Index([1, 2, 3])
        expected = idx.copy(deep=True)
        ser = pd.Series(idx, dtype="int64")
        ser.iloc[0] = 100
        tm.assert_index_equal(idx, expected)

    def test_series_string_inference(self):
        # GH#54430
        with pd.option_context("future.infer_string", True):
            ser = pd.Series(["a", "b"])
        dtype = pd.StringDtype("pyarrow" if HAS_PYARROW else "python", na_value=np.nan)
        expected = pd.Series(["a", "b"], dtype=dtype)
        tm.assert_series_equal(ser, expected)

        expected = pd.Series(["a", 1], dtype="object")
        ser = pd.Series(["a", 1])
        tm.assert_series_equal(ser, expected)

    @pytest.mark.parametrize("na_value", [None, np.nan, pd.NA])
    def test_series_string_with_na_inference(self, na_value):
        # GH#54430
        with pd.option_context("future.infer_string", True):
            ser = pd.Series(["a", na_value])
        dtype = pd.StringDtype("pyarrow" if HAS_PYARROW else "python", na_value=np.nan)
        expected = pd.Series(["a", None], dtype=dtype)
        tm.assert_series_equal(ser, expected)

    def test_series_string_inference_scalar(self):
        # GH#54430
        with pd.option_context("future.infer_string", True):
            ser = pd.Series("a", index=[1])
        dtype = pd.StringDtype("pyarrow" if HAS_PYARROW else "python", na_value=np.nan)
        expected = pd.Series("a", index=[1], dtype=dtype)
        tm.assert_series_equal(ser, expected)

    def test_series_string_inference_array_string_dtype(self):
        # GH#54496
        with pd.option_context("future.infer_string", True):
            ser = pd.Series(np.array(["a", "b"]))
        dtype = pd.StringDtype("pyarrow" if HAS_PYARROW else "python", na_value=np.nan)
        expected = pd.Series(["a", "b"], dtype=dtype)
        tm.assert_series_equal(ser, expected)

    def test_series_string_inference_storage_definition(self):
        # https://github.com/pandas-dev/pandas/issues/54793
        # but after PDEP-14 (string dtype), it was decided to keep dtype="string"
        # returning the NA string dtype, so expected is changed from
        # "string[pyarrow_numpy]" to "string[python]"
        expected = pd.Series(
            ["a", "b"], dtype="string[pyarrow]" if HAS_PYARROW else "string[python]"
        )
        with pd.option_context("future.infer_string", True):
            result = pd.Series(["a", "b"], dtype="string")
        tm.assert_series_equal(result, expected)

        expected = pd.Series(["a", "b"], dtype=pd.StringDtype(na_value=np.nan))
        with pd.option_context("future.infer_string", True):
            result = pd.Series(["a", "b"], dtype="str")
        tm.assert_series_equal(result, expected)

    def test_series_constructor_infer_string_scalar(self):
        # GH#55537
        with pd.option_context("future.infer_string", True):
            ser = pd.Series("a", index=[1, 2], dtype="string[python]")
        expected = pd.Series(["a", "a"], index=[1, 2], dtype="string[python]")
        tm.assert_series_equal(ser, expected)
        assert ser.dtype.storage == "python"

    def test_series_string_inference_na_first(self):
        # GH#55655
        with pd.option_context("future.infer_string", True):
            result = pd.Series([pd.NA, "b"])
        dtype = pd.StringDtype("pyarrow" if HAS_PYARROW else "python", na_value=np.nan)
        expected = pd.Series([None, "b"], dtype=dtype)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("klass", [pd.Series, pd.Index])
    def test_inference_on_pandas_objects(self, klass):
        # GH#56012
        obj = klass([pd.Timestamp("2019-12-31")], dtype=object)
        # This doesn't do inference
        result = pd.Series(obj)
        assert result.dtype == np.object_


class TestSeriesConstructorIndexCoercion:
    def test_series_constructor_datetimelike_index_coercion(self):
        idx = pd.date_range("2020-01-01", periods=5)
        ser = pd.Series(
            np.random.default_rng(2).standard_normal(len(idx)), idx.astype(object)
        )
        # as of 2.0, we no longer silently cast the object-dtype index
        #  to DatetimeIndex GH#39307, GH#23598
        assert not isinstance(ser.index, pd.DatetimeIndex)

    @pytest.mark.parametrize("container", [None, np.array, pd.Series, pd.Index])
    @pytest.mark.parametrize("data", [1.0, range(4)])
    def test_series_constructor_infer_multiindex(self, container, data):
        indexes = [["a", "a", "b", "b"], ["x", "y", "x", "y"]]
        if container is not None:
            indexes = [container(ind) for ind in indexes]

        multi = pd.Series(data, index=indexes)
        assert isinstance(multi.index, pd.MultiIndex)

    # TODO: make this not cast to object in pandas 3.0
    @pytest.mark.parametrize(
        "data",
        [
            ["a", "b", "c"],
            ["a", "b", np.nan],
        ],
    )
    def test_np_string_array_object_cast(self, data):
        from numpy.dtypes import StringDType

        arr = np.array(data, dtype=StringDType())
        res = pd.Series(arr)
        assert res.dtype == np.object_

        if data[-1] is np.nan:
            # as of GH#62522 the comparison op for `res==data` casts data
            #  using sanitize_array, which casts to 'str' dtype, which does not
            #  consider string 'nan' to be equal to np.nan,
            #  (which apparently numpy does?  weird.)
            assert (res.iloc[:-1] == data[:-1]).all()
            assert res.iloc[-1] == "nan"
        else:
            assert (res == data).all()


class TestSeriesConstructorInternals:
    def test_constructor_no_pandas_array(self):
        ser = pd.Series([1, 2, 3])
        result = pd.Series(ser.array)
        tm.assert_series_equal(ser, result)
        assert isinstance(result._mgr.blocks[0], NumpyBlock)
        assert result._mgr.blocks[0].is_numeric

    def test_from_array(self):
        result = pd.Series(pd.array(["1h", "2h"], dtype="timedelta64[ns]"))
        assert result._mgr.blocks[0].is_extension is False

        result = pd.Series(pd.array(["2015"], dtype="datetime64[ns]"))
        assert result._mgr.blocks[0].is_extension is False

    def test_from_list_dtype(self):
        result = pd.Series(["1h", "2h"], dtype="timedelta64[ns]")
        assert result._mgr.blocks[0].is_extension is False

        result = pd.Series(["2015"], dtype="datetime64[ns]")
        assert result._mgr.blocks[0].is_extension is False


def test_constructor(rand_series_with_duplicate_datetimeindex):
    dups = rand_series_with_duplicate_datetimeindex
    assert isinstance(dups, pd.Series)
    assert isinstance(dups.index, pd.DatetimeIndex)


@pytest.mark.parametrize(
    "input_dict,expected",
    [
        ({0: 0}, np.array([[0]], dtype=np.int64)),
        ({"a": "a"}, np.array([["a"]], dtype=object)),
        ({1: 1}, np.array([[1]], dtype=np.int64)),
    ],
)
def test_numpy_array(input_dict, expected):
    result = np.array([pd.Series(input_dict)])
    tm.assert_numpy_array_equal(result, expected)


def test_index_ordered_dict_keys():
    # GH 22077

    param_index = OrderedDict(
        [
            ((("a", "b"), ("c", "d")), 1),
            ((("a", None), ("c", "d")), 2),
        ]
    )
    series = pd.Series([1, 2], index=param_index.keys())
    expected = pd.Series(
        [1, 2],
        index=pd.MultiIndex.from_tuples(
            [(("a", "b"), ("c", "d")), (("a", None), ("c", "d"))]
        ),
    )
    tm.assert_series_equal(series, expected)


@pytest.mark.parametrize(
    "input_list",
    [
        [1, complex("nan"), 2],
        [1 + 1j, complex("nan"), 2 + 2j],
    ],
)
def test_series_with_complex_nan(input_list):
    # GH#53627
    ser = pd.Series(input_list)
    result = pd.Series(ser.array)
    assert ser.dtype == "complex128"
    tm.assert_series_equal(ser, result)


def test_dict_keys_rangeindex():
    result = pd.Series({0: 1, 1: 2})
    expected = pd.Series([1, 2], index=pd.RangeIndex(2))
    tm.assert_series_equal(result, expected, check_index_type=True)


def test_constructor_from_series_with_incompatible_dtype_raises():
    # GH#59060
    # Series(series_with_mixed_types, dtype=int) should raise consistently
    # whether creating directly or from an existing Series
    ser = pd.Series([1, 2, "x", 4, 5])
    with pytest.raises(ValueError, match="invalid literal"):
        pd.Series(ser, dtype=int)


def test_constructor_preserves_byteorder():
    # GH#43042 non-native byteorder (and the values) must be preserved
    arr = np.array([0, 256, 2**40], dtype=">i8")
    expected = [0, 256, 2**40]

    # inferred from a big-endian ndarray
    result = pd.Series(arr)
    assert result.dtype == np.dtype(">i8")
    assert result.tolist() == expected

    # explicit big-endian dtype from a python list
    result = pd.Series([0, 256, 2**40], dtype=">i8")
    assert result.dtype == np.dtype(">i8")
    assert result.tolist() == expected
