import pydoc

import numpy as np
import pytest

import pandas.util._test_decorators as td
from pandas.util._test_decorators import async_mark

import pandas as pd
from pandas import DataFrame, Index, Series, Timedelta, Timestamp, date_range
import pandas._testing as tm


class TestSeriesMisc:
    def test_getitem_preserve_name(self, datetime_series):
        result = datetime_series[datetime_series > 0]
        assert result.name == datetime_series.name

        result = datetime_series[[0, 2, 4]]
        assert result.name == datetime_series.name

        result = datetime_series[5:10]
        assert result.name == datetime_series.name

    def test_pickle_datetimes(self, datetime_series):
        unp_ts = self._pickle_roundtrip(datetime_series)
        tm.assert_series_equal(unp_ts, datetime_series)

    def test_pickle_strings(self, string_series):
        unp_series = self._pickle_roundtrip(string_series)
        tm.assert_series_equal(unp_series, string_series)

    def _pickle_roundtrip(self, obj):

        with tm.ensure_clean() as path:
            obj.to_pickle(path)
            unpickled = pd.read_pickle(path)
            return unpickled

    def test_tab_completion(self):
        # GH 9910
        s = Series(list("abcd"))
        # Series of str values should have .str but not .dt/.cat in __dir__
        assert "str" in dir(s)
        assert "dt" not in dir(s)
        assert "cat" not in dir(s)

        # similarly for .dt
        s = Series(date_range("1/1/2015", periods=5))
        assert "dt" in dir(s)
        assert "str" not in dir(s)
        assert "cat" not in dir(s)

        # Similarly for .cat, but with the twist that str and dt should be
        # there if the categories are of that type first cat and str.
        s = Series(list("abbcd"), dtype="category")
        assert "cat" in dir(s)
        assert "str" in dir(s)  # as it is a string categorical
        assert "dt" not in dir(s)

        # similar to cat and str
        s = Series(date_range("1/1/2015", periods=5)).astype("category")
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

        def get_dir(s):
            results = [r for r in s.cat.__dir__() if not r.startswith("_")]
            return sorted(set(results))

        s = Series(list("aabbcde")).astype("category")
        results = get_dir(s)
        tm.assert_almost_equal(results, sorted(set(ok_for_cat)))

    @pytest.mark.parametrize(
        "index",
        [
            tm.makeUnicodeIndex(10),
            tm.makeStringIndex(10),
            tm.makeCategoricalIndex(10),
            Index(["foo", "bar", "baz"] * 2),
            tm.makeDateIndex(10),
            tm.makePeriodIndex(10),
            tm.makeTimedeltaIndex(10),
            tm.makeIntIndex(10),
            tm.makeUIntIndex(10),
            tm.makeIntIndex(10),
            tm.makeFloatIndex(10),
            Index([True, False]),
            Index([f"a{i}" for i in range(101)]),
            pd.MultiIndex.from_tuples(zip("ABCD", "EFGH")),
            pd.MultiIndex.from_tuples(zip([0, 1, 2, 3], "EFGH")),
        ],
    )
    def test_index_tab_completion(self, index):
        # dir contains string-like values of the Index.
        s = Series(index=index, dtype=object)
        dir_s = dir(s)
        for i, x in enumerate(s.index.unique(level=0)):
            if i < 100:
                assert not isinstance(x, str) or not x.isidentifier() or x in dir_s
            else:
                assert x not in dir_s

    def test_not_hashable(self):
        s_empty = Series(dtype=object)
        s = Series([1])
        msg = "'Series' objects are mutable, thus they cannot be hashed"
        with pytest.raises(TypeError, match=msg):
            hash(s_empty)
        with pytest.raises(TypeError, match=msg):
            hash(s)

    def test_contains(self, datetime_series):
        tm.assert_contains_all(datetime_series.index, datetime_series)

    def test_iter_datetimes(self, datetime_series):
        for i, val in enumerate(datetime_series):
            assert val == datetime_series[i]

    def test_iter_strings(self, string_series):
        for i, val in enumerate(string_series):
            assert val == string_series[i]

    def test_keys(self, datetime_series):
        assert datetime_series.keys() is datetime_series.index

    def test_values(self, datetime_series):
        tm.assert_almost_equal(
            datetime_series.values, datetime_series, check_dtype=False
        )

    def test_iteritems_datetimes(self, datetime_series):
        for idx, val in datetime_series.iteritems():
            assert val == datetime_series[idx]

    def test_iteritems_strings(self, string_series):
        for idx, val in string_series.iteritems():
            assert val == string_series[idx]

        # assert is lazy (generators don't define reverse, lists do)
        assert not hasattr(string_series.iteritems(), "reverse")

    def test_items_datetimes(self, datetime_series):
        for idx, val in datetime_series.items():
            assert val == datetime_series[idx]

    def test_items_strings(self, string_series):
        for idx, val in string_series.items():
            assert val == string_series[idx]

        # assert is lazy (generators don't define reverse, lists do)
        assert not hasattr(string_series.items(), "reverse")

    def test_raise_on_info(self):
        s = Series(np.random.randn(10))
        msg = "'Series' object has no attribute 'info'"
        with pytest.raises(AttributeError, match=msg):
            s.info()

    def test_axis_alias(self):
        s = Series([1, 2, np.nan])
        tm.assert_series_equal(s.dropna(axis="rows"), s.dropna(axis="index"))
        assert s.dropna().sum("rows") == 3
        assert s._get_axis_number("rows") == 0
        assert s._get_axis_name("rows") == "index"

    def test_class_axis(self):
        # https://github.com/pandas-dev/pandas/issues/18147
        # no exception and no empty docstring
        assert pydoc.getdoc(Series.index)

    def test_item(self):
        s = Series([1])
        result = s.item()
        assert result == 1
        assert result == s.iloc[0]
        assert isinstance(result, int)  # i.e. not np.int64

        ser = Series([0.5], index=[3])
        result = ser.item()
        assert isinstance(result, float)
        assert result == 0.5

        ser = Series([1, 2])
        msg = "can only convert an array of size 1"
        with pytest.raises(ValueError, match=msg):
            ser.item()

        dti = pd.date_range("2016-01-01", periods=2)
        with pytest.raises(ValueError, match=msg):
            dti.item()
        with pytest.raises(ValueError, match=msg):
            Series(dti).item()

        val = dti[:1].item()
        assert isinstance(val, Timestamp)
        val = Series(dti)[:1].item()
        assert isinstance(val, Timestamp)

        tdi = dti - dti
        with pytest.raises(ValueError, match=msg):
            tdi.item()
        with pytest.raises(ValueError, match=msg):
            Series(tdi).item()

        val = tdi[:1].item()
        assert isinstance(val, Timedelta)
        val = Series(tdi)[:1].item()
        assert isinstance(val, Timedelta)

        # Case where ser[0] would not work
        ser = Series(dti, index=[5, 6])
        val = ser[:1].item()
        assert val == dti[0]

    def test_ndarray_compat(self):

        # test numpy compat with Series as sub-class of NDFrame
        tsdf = DataFrame(
            np.random.randn(1000, 3),
            columns=["A", "B", "C"],
            index=date_range("1/1/2000", periods=1000),
        )

        def f(x):
            return x[x.idxmax()]

        result = tsdf.apply(f)
        expected = tsdf.max()
        tm.assert_series_equal(result, expected)

        # using an ndarray like function
        s = Series(np.random.randn(10))
        result = Series(np.ones_like(s))
        expected = Series(1, index=range(10), dtype="float64")
        tm.assert_series_equal(result, expected)

        # ravel
        s = Series(np.random.randn(10))
        tm.assert_almost_equal(s.ravel(order="F"), s.values.ravel(order="F"))

    def test_empty_method(self):
        s_empty = Series(dtype=object)
        assert s_empty.empty

        s2 = Series(index=[1], dtype=object)
        for full_series in [Series([1]), s2]:
            assert not full_series.empty

    @async_mark()
    @td.check_file_leaks
    async def test_tab_complete_warning(self, ip):
        # https://github.com/pandas-dev/pandas/issues/16409
        pytest.importorskip("IPython", minversion="6.0.0")
        from IPython.core.completer import provisionalcompleter

        code = "import pandas as pd; s = Series(dtype=object)"
        await ip.run_code(code)

        # TODO: remove it when Ipython updates
        # GH 33567, jedi version raises Deprecation warning in Ipython
        import jedi

        if jedi.__version__ < "0.17.0":
            warning = tm.assert_produces_warning(None)
        else:
            warning = tm.assert_produces_warning(
                DeprecationWarning, check_stacklevel=False
            )
        with warning:
            with provisionalcompleter("ignore"):
                list(ip.Completer.completions("s.", 1))

    def test_integer_series_size(self):
        # GH 25580
        s = Series(range(9))
        assert s.size == 9
        s = Series(range(9), dtype="Int64")
        assert s.size == 9

    def test_attrs(self):
        s = Series([0, 1], name="abc")
        assert s.attrs == {}
        s.attrs["version"] = 1
        result = s + 1
        assert result.attrs == {"version": 1}

    @pytest.mark.parametrize("allows_duplicate_labels", [True, False, None])
    def test_set_flags(self, allows_duplicate_labels):
        df = Series([1, 2])
        result = df.set_flags(allows_duplicate_labels=allows_duplicate_labels)
        if allows_duplicate_labels is None:
            # We don't update when it's not provided
            assert result.flags.allows_duplicate_labels is True
        else:
            assert result.flags.allows_duplicate_labels is allows_duplicate_labels

        # We made a copy
        assert df is not result
        # We didn't mutate df
        assert df.flags.allows_duplicate_labels is True

        # But we didn't copy data
        result.iloc[0] = 0
        assert df.iloc[0] == 0

        # Now we do copy.
        result = df.set_flags(
            copy=True, allows_duplicate_labels=allows_duplicate_labels
        )
        result.iloc[0] = 10
        assert df.iloc[0] == 0
