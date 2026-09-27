import inspect
import pydoc

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestSeriesMisc:
    def test_tab_completion(self):
        # GH 9910
        s = pd.Series(list("abcd"))
        # Series of str values should have .str but not .dt/.cat in __dir__
        assert "str" in dir(s)
        assert "dt" not in dir(s)
        assert "cat" not in dir(s)

    def test_tab_completion_dt(self):
        # similarly for .dt
        s = pd.Series(pd.date_range("1/1/2015", periods=5))
        assert "dt" in dir(s)
        assert "str" not in dir(s)
        assert "cat" not in dir(s)

    def test_tab_completion_cat(self):
        # Similarly for .cat, but with the twist that str and dt should be
        # there if the categories are of that type first cat and str.
        s = pd.Series(list("abbcd"), dtype="category")
        assert "cat" in dir(s)
        assert "str" in dir(s)  # as it is a string categorical
        assert "dt" not in dir(s)

    def test_tab_completion_cat_str(self):
        # similar to cat and str
        s = pd.Series(pd.date_range("1/1/2015", periods=5)).astype("category")
        assert "cat" in dir(s)
        assert "str" not in dir(s)
        assert "dt" in dir(s)  # as it is a datetime categorical

    def test_tab_completion_with_categorical(self):
        # test the tab completion display
        ok_for_cat = [
            "categories",
            "codes",
            "ordered",
            "set_categories",
            "add_categories",
            "remove_categories",
            "rename_categories",
            "reorder_categories",
            "remove_unused_categories",
            "as_ordered",
            "as_unordered",
        ]

        s = pd.Series(list("aabbcde")).astype("category")
        results = sorted({r for r in s.cat.__dir__() if not r.startswith("_")})
        tm.assert_almost_equal(results, sorted(set(ok_for_cat)))

    @pytest.mark.parametrize(
        "index",
        [
            pd.Index(list("ab") * 5, dtype="category"),
            pd.Index([str(i) for i in range(10)]),
            pd.Index(["foo", "bar", "baz"] * 2),
            pd.date_range("2020-01-01", periods=10),
            pd.period_range("2020-01-01", periods=10, freq="D"),
            pd.timedelta_range("1 day", periods=10),
            pd.Index(np.arange(10), dtype=np.uint64),
            pd.Index(np.arange(10), dtype=np.int64),
            pd.Index(np.arange(10), dtype=np.float64),
            pd.Index([True, False]),
            pd.Index([f"a{i}" for i in range(101)]),
            pd.MultiIndex.from_tuples(zip("ABCD", "EFGH", strict=True)),
            pd.MultiIndex.from_tuples(zip([0, 1, 2, 3], "EFGH", strict=True)),
        ],
    )
    def test_index_tab_completion(self, index):
        # dir contains string-like values of the Index.
        s = pd.Series(index=index, dtype=object)
        dir_s = dir(s)
        for i, x in enumerate(s.index.unique(level=0)):
            if i < 100:
                assert not isinstance(x, str) or not x.isidentifier() or x in dir_s
            else:
                assert x not in dir_s

    @pytest.mark.parametrize("ser", [pd.Series(dtype=object), pd.Series([1])])
    def test_not_hashable(self, ser):
        msg = "unhashable type: 'Series'"
        with pytest.raises(TypeError, match=msg):
            hash(ser)

    def test_contains(self, datetime_series):
        tm.assert_contains_all(datetime_series.index, datetime_series)

    def test_axis_alias(self):
        s = pd.Series([1, 2, np.nan])
        tm.assert_series_equal(s.dropna(axis="rows"), s.dropna(axis="index"))
        assert s.dropna().sum(axis="rows") == 3
        assert s._get_axis_number("rows") == 0
        assert s._get_axis_name("rows") == "index"

    def test_class_axis(self):
        # https://github.com/pandas-dev/pandas/issues/18147
        # no exception and no empty docstring
        assert pydoc.getdoc(pd.Series.index)

    def test_ndarray_compat(self):
        # test numpy compat with Series as sub-class of NDFrame
        tsdf = pd.DataFrame(
            np.random.default_rng(2).standard_normal((1000, 3)),
            columns=["A", "B", "C"],
            index=pd.date_range("1/1/2000", periods=1000),
        )

        def f(x):
            return x[x.idxmax()]

        result = tsdf.apply(f)
        expected = tsdf.max()
        tm.assert_series_equal(result, expected)

    def test_ndarray_compat_like_func(self):
        # using an ndarray like function
        s = pd.Series(np.random.default_rng(2).standard_normal(10))
        result = pd.Series(np.ones_like(s))
        expected = pd.Series(1, index=range(10), dtype="float64")
        tm.assert_series_equal(result, expected)

    def test_empty_method(self):
        s_empty = pd.Series(dtype=object)
        assert s_empty.empty

    @pytest.mark.parametrize("dtype", ["int64", object])
    def test_empty_method_full_series(self, dtype):
        full_series = pd.Series(index=[1], dtype=dtype)
        assert not full_series.empty

    @pytest.mark.parametrize("dtype", [None, "Int64"])
    def test_integer_series_size(self, dtype):
        # GH 25580
        s = pd.Series(range(9), dtype=dtype)
        assert s.size == 9

    def test_attrs(self):
        s = pd.Series([0, 1], name="abc")
        assert s.attrs == {}
        s.attrs["version"] = 1
        result = s + 1
        assert result.attrs == {"version": 1}

    def test_inspect_getmembers(self):
        # GH38782
        ser = pd.Series(dtype=object)
        inspect.getmembers(ser)

    def test_unknown_attribute(self):
        # GH#9680
        tdi = pd.timedelta_range(start=0, periods=10, freq="1s")
        ser = pd.Series(np.random.default_rng(2).normal(size=10), index=tdi)
        assert "foo" not in ser.__dict__
        msg = "'Series' object has no attribute 'foo'"
        with pytest.raises(AttributeError, match=msg):
            ser.foo

    @pytest.mark.parametrize("op", ["year", "day", "second", "weekday"])
    def test_datetime_series_no_datelike_attrs(self, op, datetime_series):
        # GH#7206
        msg = f"'Series' object has no attribute '{op}'"
        with pytest.raises(AttributeError, match=msg):
            getattr(datetime_series, op)

    def test_series_datetimelike_attribute_access(self):
        # attribute access should still work!
        ser = pd.Series({"year": 2000, "month": 1, "day": 10})
        assert ser.year == 2000
        assert ser.month == 1
        assert ser.day == 10

    def test_series_datetimelike_attribute_access_invalid(self):
        ser = pd.Series({"year": 2000, "month": 1, "day": 10})
        msg = "'Series' object has no attribute 'weekday'"
        with pytest.raises(AttributeError, match=msg):
            ser.weekday

    @pytest.mark.parametrize(
        "kernel, has_numeric_only",
        [
            ("skew", True),
            ("var", True),
            ("all", False),
            ("prod", True),
            ("any", False),
            ("idxmin", False),
            ("quantile", False),
            ("idxmax", False),
            ("min", True),
            ("sem", True),
            ("mean", True),
            ("nunique", False),
            ("max", True),
            ("sum", True),
            ("count", False),
            ("median", True),
            ("std", True),
            ("rank", True),
            ("pct_change", False),
            ("cummax", False),
            ("shift", False),
            ("diff", False),
            ("cumsum", False),
            ("cummin", False),
            ("cumprod", False),
            ("fillna", False),
            ("ffill", False),
            ("bfill", False),
            ("sample", False),
            ("tail", False),
            ("take", False),
            ("head", False),
            ("cov", False),
            ("corr", False),
        ],
    )
    @pytest.mark.parametrize("dtype", [bool, int, float, object])
    def test_numeric_only(self, kernel, has_numeric_only, dtype):
        # GH#47500
        ser = pd.Series([0, 1, 1], dtype=dtype)
        if kernel == "corrwith":
            args = (ser,)
        elif kernel == "corr":
            args = (ser,)
        elif kernel == "cov":
            args = (ser,)
        elif kernel == "nth":
            args = (0,)
        elif kernel == "fillna":
            args = (True,)
        elif kernel == "fillna":
            args = ("ffill",)
        elif kernel == "take":
            args = ([0],)
        elif kernel == "quantile":
            args = (0.5,)
        else:
            args = ()
        method = getattr(ser, kernel)
        if not has_numeric_only:
            msg = (
                "(got an unexpected keyword argument 'numeric_only'"
                "|too many arguments passed in)"
            )
            with pytest.raises(TypeError, match=msg):
                method(*args, numeric_only=True)
        elif dtype is object:
            msg = f"Series.{kernel} does not allow numeric_only=True with non-numeric"
            with pytest.raises(TypeError, match=msg):
                method(*args, numeric_only=True)
        else:
            result = method(*args, numeric_only=True)
            expected = method(*args, numeric_only=False)
            if isinstance(expected, pd.Series):
                # transformer
                tm.assert_series_equal(result, expected)
            else:
                # reducer
                assert result == expected
