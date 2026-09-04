import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestAstype:
    @pytest.mark.parametrize("cls", [pd.Categorical, pd.CategoricalIndex])
    @pytest.mark.parametrize("values", [[1, np.nan], [pd.Timestamp("2000"), pd.NaT]])
    def test_astype_nan_to_int(self, cls, values):
        # GH#28406
        obj = cls(values)

        msg = "Cannot (cast|convert)"
        with pytest.raises((ValueError, TypeError), match=msg):
            obj.astype(int)

    @pytest.mark.parametrize(
        "expected",
        [
            pd.array(["2019", "2020"], dtype="datetime64[ns, UTC]"),
            pd.array([0, 0], dtype="timedelta64[ns]"),
            pd.array([pd.Period("2019"), pd.Period("2020")], dtype="period[Y-DEC]"),
            pd.array([pd.Interval(0, 1), pd.Interval(1, 2)], dtype="interval"),
            pd.array([1, pd.NA], dtype="Int64"),
        ],
    )
    def test_astype_category_to_extension_dtype(self, expected):
        # GH#28668
        result = expected.astype("category").astype(expected.dtype)

        tm.assert_extension_array_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype, expected",
        [
            (
                "datetime64[ns]",
                np.array(["2015-01-01T00:00:00.000000000"], dtype="datetime64[ns]"),
            ),
            (
                "datetime64[ns, MET]",
                pd.DatetimeIndex(
                    [pd.Timestamp("2015-01-01 00:00:00+0100", tz="MET")]
                ).array,
            ),
        ],
    )
    def test_astype_to_datetime64(self, dtype, expected):
        # GH#28448
        result = pd.Categorical(["2015-01-01"]).astype(dtype)
        assert result == expected

    def test_astype_str_int_categories_to_nullable_int(self):
        # GH#39616
        dtype = pd.CategoricalDtype([str(i) for i in range(5)])
        codes = np.random.default_rng(2).integers(5, size=20)
        arr = pd.Categorical.from_codes(codes, dtype=dtype)

        res = arr.astype("Int64")
        expected = pd.array(codes, dtype="Int64")
        tm.assert_extension_array_equal(res, expected)

    def test_astype_str_int_categories_to_nullable_float(self):
        # GH#39616
        dtype = pd.CategoricalDtype([str(i / 2) for i in range(5)])
        codes = np.random.default_rng(2).integers(5, size=20)
        arr = pd.Categorical.from_codes(codes, dtype=dtype)

        res = arr.astype("Float64")
        expected = pd.array(codes, dtype="Float64") / 2
        tm.assert_extension_array_equal(res, expected)

    def test_astype(self, ordered):
        # string
        cat = pd.Categorical(list("abbaaccc"), ordered=ordered)
        result = cat.astype(object)
        expected = np.array(cat)
        tm.assert_numpy_array_equal(result, expected)

        msg = r"Cannot cast (object|str) dtype to float64"
        with pytest.raises(ValueError, match=msg):
            cat.astype(float)

        # numeric
        cat = pd.Categorical([0, 1, 2, 2, 1, 0, 1, 0, 2], ordered=ordered)
        result = cat.astype(object)
        expected = np.array(cat, dtype=object)
        tm.assert_numpy_array_equal(result, expected)

        result = cat.astype(int)
        expected = np.array(cat, dtype="int")
        tm.assert_numpy_array_equal(result, expected)

        result = cat.astype(float)
        expected = np.array(cat, dtype=float)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("dtype_ordered", [True, False])
    def test_astype_category(self, dtype_ordered, ordered):
        # GH#10696/GH#18593
        data = list("abcaacbab")
        cat = pd.Categorical(data, categories=list("bac"), ordered=ordered)

        # standard categories
        dtype = pd.CategoricalDtype(ordered=dtype_ordered)
        result = cat.astype(dtype)
        expected = pd.Categorical(
            data, categories=cat.categories, ordered=dtype_ordered
        )
        tm.assert_categorical_equal(result, expected)

        # non-standard categories
        dtype = pd.CategoricalDtype(list("adc"), dtype_ordered)
        msg = "Constructing a Categorical with a dtype and values containing"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = cat.astype(dtype)
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            expected = pd.Categorical(data, dtype=dtype)
        tm.assert_categorical_equal(result, expected)

        if dtype_ordered is False:
            # dtype='category' can't specify ordered, so only test once
            result = cat.astype("category")
            expected = cat
            tm.assert_categorical_equal(result, expected)

    def test_astype_category_copy_false_nocopy_codes(self):
        # GH#62000
        cat = pd.Categorical([3, 2, 4, 1])
        new = cat.astype("category", copy=False)
        assert tm.shares_memory(new.codes, cat.codes)
        new = cat.astype("category", copy=True)
        assert not tm.shares_memory(new.codes, cat.codes)

    def test_astype_object_datetime_categories(self):
        # GH#40754
        cat = pd.Categorical(pd.to_datetime(["2021-03-27", pd.NaT]))
        result = cat.astype(object)
        expected = np.array(
            [pd.Timestamp("2021-03-27 00:00:00"), pd.NaT], dtype="object"
        )
        tm.assert_numpy_array_equal(result, expected)

    def test_astype_object_timestamp_categories(self):
        # GH#18024
        cat = pd.Categorical([pd.Timestamp("2014-01-01")])
        result = cat.astype(object)
        expected = np.array([pd.Timestamp("2014-01-01 00:00:00")], dtype="object")
        tm.assert_numpy_array_equal(result, expected)

    def test_astype_category_readonly_mask_values(self):
        # GH#53658
        arr = pd.array([0, 1, 2], dtype="Int64")
        arr._mask.flags["WRITEABLE"] = False
        result = arr.astype("category")
        expected = pd.array([0, 1, 2], dtype="Int64").astype("category")
        tm.assert_extension_array_equal(result, expected)
