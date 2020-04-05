from collections import OrderedDict
import pydoc
import warnings

import numpy as np
import pytest

from pandas.util._test_decorators import async_mark

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    DatetimeIndex,
    Index,
    Series,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    date_range,
    period_range,
    timedelta_range,
)
import pandas._testing as tm
from pandas.core.arrays import PeriodArray

import pandas.io.formats.printing as printing


class TestSeriesMisc:
    def test_scalarop_preserve_name(self, datetime_series):
        result = datetime_series * 2
        assert result.name == datetime_series.name

    def test_copy_name(self, datetime_series):
        result = datetime_series.copy()
        assert result.name == datetime_series.name

    def test_copy_index_name_checking(self, datetime_series):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy

        datetime_series.index.name = None
        assert datetime_series.index.name is None
        assert datetime_series is datetime_series

        cp = datetime_series.copy()
        cp.index.name = "foo"
        printing.pprint_thing(datetime_series.index.name)
        assert datetime_series.index.name is None

    def test_append_preserve_name(self, datetime_series):
        result = datetime_series[:5].append(datetime_series[5:])
        assert result.name == datetime_series.name

    def test_binop_maybe_preserve_name(self, datetime_series):
        # names match, preserve
        result = datetime_series * datetime_series
        assert result.name == datetime_series.name
        result = datetime_series.mul(datetime_series)
        assert result.name == datetime_series.name

        result = datetime_series * datetime_series[:-2]
        assert result.name == datetime_series.name

        # names don't match, don't preserve
        cp = datetime_series.copy()
        cp.name = "something else"
        result = datetime_series + cp
        assert result.name is None
        result = datetime_series.add(cp)
        assert result.name is None

        ops = ["add", "sub", "mul", "div", "truediv", "floordiv", "mod", "pow"]
        ops = ops + ["r" + op for op in ops]
        for op in ops:
            # names match, preserve
            s = datetime_series.copy()
            result = getattr(s, op)(s)
            assert result.name == datetime_series.name

            # names don't match, don't preserve
            cp = datetime_series.copy()
            cp.name = "changed"
            result = getattr(s, op)(cp)
            assert result.name is None

    def test_combine_first_name(self, datetime_series):
        result = datetime_series.combine_first(datetime_series[:5])
        assert result.name == datetime_series.name

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

    def test_sort_index_name(self, datetime_series):
        result = datetime_series.sort_index(ascending=False)
        assert result.name == datetime_series.name

    def test_constructor_dict(self):
        d = {"a": 0.0, "b": 1.0, "c": 2.0}
        result = Series(d)
        expected = Series(d, index=sorted(d.keys()))
        tm.assert_series_equal(result, expected)

        result = Series(d, index=["b", "c", "d", "a"])
        expected = Series([1, 2, np.nan, 0], index=["b", "c", "d", "a"])
        tm.assert_series_equal(result, expected)

    def test_constructor_subclass_dict(self, dict_subclass):
        data = dict_subclass((x, 10.0 * x) for x in range(10))
        series = Series(data)
        expected = Series(dict(data.items()))
        tm.assert_series_equal(series, expected)

    def test_constructor_ordereddict(self):
        # GH3283
        data = OrderedDict((f"col{i}", np.random.random()) for i in range(12))

        series = Series(data)
        expected = Series(list(data.values()), list(data.keys()))
        tm.assert_series_equal(series, expected)

        # Test with subclass
        class A(OrderedDict):
            pass

        series = Series(A(data))
        tm.assert_series_equal(series, expected)

    def test_constructor_dict_multiindex(self):
        d = {("a", "a"): 0.0, ("b", "a"): 1.0, ("b", "c"): 2.0}
        _d = sorted(d.items())
        result = Series(d)
        expected = Series(
            [x[1] for x in _d], index=pd.MultiIndex.from_tuples([x[0] for x in _d])
        )
        tm.assert_series_equal(result, expected)

        d["z"] = 111.0
        _d.insert(0, ("z", d["z"]))
        result = Series(d)
        expected = Series(
            [x[1] for x in _d], index=pd.Index([x[0] for x in _d], tupleize_cols=False)
        )
        result = result.reindex(index=expected.index)
        tm.assert_series_equal(result, expected)

    def test_constructor_dict_timedelta_index(self):
        # GH #12169 : Resample category data with timedelta index
        # construct Series from dict as data and TimedeltaIndex as index
        # will result NaN in result Series data
        expected = Series(
            data=["A", "B", "C"], index=pd.to_timedelta([0, 10, 20], unit="s")
        )

        result = Series(
            data={
                pd.to_timedelta(0, unit="s"): "A",
                pd.to_timedelta(10, unit="s"): "B",
                pd.to_timedelta(20, unit="s"): "C",
            },
            index=pd.to_timedelta([0, 10, 20], unit="s"),
        )
        tm.assert_series_equal(result, expected)

    def test_sparse_accessor_updates_on_inplace(self):
        s = pd.Series([1, 1, 2, 3], dtype="Sparse[int]")
        s.drop([0, 1], inplace=True)
        assert s.sparse.density == 1.0

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
        s = pd.Series(index=index, dtype=object)
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
        # HACK: By doing this in two stages, we avoid 2to3 wrapping the call
        # to .keys() in a list()
        getkeys = datetime_series.keys
        assert getkeys() is datetime_series.index

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

    def test_copy(self):

        for deep in [None, False, True]:
            s = Series(np.arange(10), dtype="float64")

            # default deep is True
            if deep is None:
                s2 = s.copy()
            else:
                s2 = s.copy(deep=deep)

            s2[::2] = np.NaN

            if deep is None or deep is True:
                # Did not modify original Series
                assert np.isnan(s2[0])
                assert not np.isnan(s[0])
            else:
                # we DID modify the original Series
                assert np.isnan(s2[0])
                assert np.isnan(s[0])

    def test_copy_tzaware(self):
        # GH#11794
        # copy of tz-aware
        expected = Series([Timestamp("2012/01/01", tz="UTC")])
        expected2 = Series([Timestamp("1999/01/01", tz="UTC")])

        for deep in [None, False, True]:

            s = Series([Timestamp("2012/01/01", tz="UTC")])

            if deep is None:
                s2 = s.copy()
            else:
                s2 = s.copy(deep=deep)

            s2[0] = pd.Timestamp("1999/01/01", tz="UTC")

            # default deep is True
            if deep is None or deep is True:
                # Did not modify original Series
                tm.assert_series_equal(s2, expected2)
                tm.assert_series_equal(s, expected)
            else:
                # we DID modify the original Series
                tm.assert_series_equal(s2, expected2)
                tm.assert_series_equal(s, expected2)

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

    def test_numpy_unique(self, datetime_series):
        # it works!
        np.unique(datetime_series)

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

    def test_str_accessor_updates_on_inplace(self):
        s = pd.Series(list("abc"))
        s.drop([0], inplace=True)
        assert len(s.str.lower()) == 2

    def test_str_attribute(self):
        # GH9068
        methods = ["strip", "rstrip", "lstrip"]
        s = Series([" jack", "jill ", " jesse ", "frank"])
        for method in methods:
            expected = Series([getattr(str, method)(x) for x in s.values])
            tm.assert_series_equal(getattr(Series.str, method)(s.str), expected)

        # str accessor only valid with string values
        s = Series(range(5))
        with pytest.raises(AttributeError, match="only use .str accessor"):
            s.str.repeat(2)

    def test_empty_method(self):
        s_empty = pd.Series(dtype=object)
        assert s_empty.empty

        s2 = pd.Series(index=[1], dtype=object)
        for full_series in [pd.Series([1]), s2]:
            assert not full_series.empty

    @async_mark()
    async def test_tab_complete_warning(self, ip):
        # https://github.com/pandas-dev/pandas/issues/16409
        pytest.importorskip("IPython", minversion="6.0.0")
        from IPython.core.completer import provisionalcompleter

        code = "import pandas as pd; s = pd.Series()"
        await ip.run_code(code)
        with tm.assert_produces_warning(None):
            with provisionalcompleter("ignore"):
                list(ip.Completer.completions("s.", 1))

    def test_integer_series_size(self):
        # GH 25580
        s = Series(range(9))
        assert s.size == 9
        s = Series(range(9), dtype="Int64")
        assert s.size == 9

    def test_attrs(self):
        s = pd.Series([0, 1], name="abc")
        assert s.attrs == {}
        s.attrs["version"] = 1
        result = s + 1
        assert result.attrs == {"version": 1}


class TestCategoricalSeries:
    @pytest.mark.parametrize(
        "method",
        [
            lambda x: x.cat.set_categories([1, 2, 3]),
            lambda x: x.cat.reorder_categories([2, 3, 1], ordered=True),
            lambda x: x.cat.rename_categories([1, 2, 3]),
            lambda x: x.cat.remove_unused_categories(),
            lambda x: x.cat.remove_categories([2]),
            lambda x: x.cat.add_categories([4]),
            lambda x: x.cat.as_ordered(),
            lambda x: x.cat.as_unordered(),
        ],
    )
    def test_getname_categorical_accessor(self, method):
        # GH 17509
        s = Series([1, 2, 3], name="A").astype("category")
        expected = "A"
        result = method(s).name
        assert result == expected

    def test_cat_accessor(self):
        s = Series(Categorical(["a", "b", np.nan, "a"]))
        tm.assert_index_equal(s.cat.categories, Index(["a", "b"]))
        assert not s.cat.ordered, False

        exp = Categorical(["a", "b", np.nan, "a"], categories=["b", "a"])
        s.cat.set_categories(["b", "a"], inplace=True)
        tm.assert_categorical_equal(s.values, exp)

        res = s.cat.set_categories(["b", "a"])
        tm.assert_categorical_equal(res.values, exp)

        s[:] = "a"
        s = s.cat.remove_unused_categories()
        tm.assert_index_equal(s.cat.categories, Index(["a"]))

    def test_cat_accessor_api(self):
        # GH 9322
        from pandas.core.arrays.categorical import CategoricalAccessor

        assert Series.cat is CategoricalAccessor
        s = Series(list("aabbcde")).astype("category")
        assert isinstance(s.cat, CategoricalAccessor)

        invalid = Series([1])
        with pytest.raises(AttributeError, match="only use .cat accessor"):
            invalid.cat
        assert not hasattr(invalid, "cat")

    def test_cat_accessor_no_new_attributes(self):
        # https://github.com/pandas-dev/pandas/issues/10673
        c = Series(list("aabbcde")).astype("category")
        with pytest.raises(AttributeError, match="You cannot add any new attribute"):
            c.cat.xlabel = "a"

    def test_cat_accessor_updates_on_inplace(self):
        s = Series(list("abc")).astype("category")
        s.drop(0, inplace=True)
        s.cat.remove_unused_categories(inplace=True)
        assert len(s.cat.categories) == 2

    def test_categorical_delegations(self):

        # invalid accessor
        msg = r"Can only use \.cat accessor with a 'category' dtype"
        with pytest.raises(AttributeError, match=msg):
            Series([1, 2, 3]).cat
        with pytest.raises(AttributeError, match=msg):
            Series([1, 2, 3]).cat()
        with pytest.raises(AttributeError, match=msg):
            Series(["a", "b", "c"]).cat
        with pytest.raises(AttributeError, match=msg):
            Series(np.arange(5.0)).cat
        with pytest.raises(AttributeError, match=msg):
            Series([Timestamp("20130101")]).cat

        # Series should delegate calls to '.categories', '.codes', '.ordered'
        # and the methods '.set_categories()' 'drop_unused_categories()' to the
        # categorical
        s = Series(Categorical(["a", "b", "c", "a"], ordered=True))
        exp_categories = Index(["a", "b", "c"])
        tm.assert_index_equal(s.cat.categories, exp_categories)
        s.cat.categories = [1, 2, 3]
        exp_categories = Index([1, 2, 3])
        tm.assert_index_equal(s.cat.categories, exp_categories)

        exp_codes = Series([0, 1, 2, 0], dtype="int8")
        tm.assert_series_equal(s.cat.codes, exp_codes)

        assert s.cat.ordered
        s = s.cat.as_unordered()
        assert not s.cat.ordered
        s.cat.as_ordered(inplace=True)
        assert s.cat.ordered

        # reorder
        s = Series(Categorical(["a", "b", "c", "a"], ordered=True))
        exp_categories = Index(["c", "b", "a"])
        exp_values = np.array(["a", "b", "c", "a"], dtype=np.object_)
        s = s.cat.set_categories(["c", "b", "a"])
        tm.assert_index_equal(s.cat.categories, exp_categories)
        tm.assert_numpy_array_equal(s.values.__array__(), exp_values)
        tm.assert_numpy_array_equal(s.__array__(), exp_values)

        # remove unused categories
        s = Series(Categorical(["a", "b", "b", "a"], categories=["a", "b", "c"]))
        exp_categories = Index(["a", "b"])
        exp_values = np.array(["a", "b", "b", "a"], dtype=np.object_)
        s = s.cat.remove_unused_categories()
        tm.assert_index_equal(s.cat.categories, exp_categories)
        tm.assert_numpy_array_equal(s.values.__array__(), exp_values)
        tm.assert_numpy_array_equal(s.__array__(), exp_values)

        # This method is likely to be confused, so test that it raises an error
        # on wrong inputs:
        msg = "'Series' object has no attribute 'set_categories'"
        with pytest.raises(AttributeError, match=msg):
            s.set_categories([4, 3, 2, 1])

        # right: s.cat.set_categories([4,3,2,1])

        # GH18862 (let Series.cat.rename_categories take callables)
        s = Series(Categorical(["a", "b", "c", "a"], ordered=True))
        result = s.cat.rename_categories(lambda x: x.upper())
        expected = Series(
            Categorical(["A", "B", "C", "A"], categories=["A", "B", "C"], ordered=True)
        )
        tm.assert_series_equal(result, expected)

    def test_dt_accessor_api_for_categorical(self):
        # https://github.com/pandas-dev/pandas/issues/10661
        from pandas.core.indexes.accessors import Properties

        s_dr = Series(date_range("1/1/2015", periods=5, tz="MET"))
        c_dr = s_dr.astype("category")

        s_pr = Series(period_range("1/1/2015", freq="D", periods=5))
        c_pr = s_pr.astype("category")

        s_tdr = Series(timedelta_range("1 days", "10 days"))
        c_tdr = s_tdr.astype("category")

        # only testing field (like .day)
        # and bool (is_month_start)
        get_ops = lambda x: x._datetimelike_ops

        test_data = [
            ("Datetime", get_ops(DatetimeIndex), s_dr, c_dr),
            ("Period", get_ops(PeriodArray), s_pr, c_pr),
            ("Timedelta", get_ops(TimedeltaIndex), s_tdr, c_tdr),
        ]

        assert isinstance(c_dr.dt, Properties)

        special_func_defs = [
            ("strftime", ("%Y-%m-%d",), {}),
            ("tz_convert", ("EST",), {}),
            ("round", ("D",), {}),
            ("floor", ("D",), {}),
            ("ceil", ("D",), {}),
            ("asfreq", ("D",), {}),
            # FIXME: don't leave commented-out
            # ('tz_localize', ("UTC",), {}),
        ]
        _special_func_names = [f[0] for f in special_func_defs]

        # the series is already localized
        _ignore_names = ["tz_localize", "components"]

        for name, attr_names, s, c in test_data:
            func_names = [
                f
                for f in dir(s.dt)
                if not (
                    f.startswith("_")
                    or f in attr_names
                    or f in _special_func_names
                    or f in _ignore_names
                )
            ]

            func_defs = [(f, (), {}) for f in func_names]
            for f_def in special_func_defs:
                if f_def[0] in dir(s.dt):
                    func_defs.append(f_def)

            for func, args, kwargs in func_defs:
                with warnings.catch_warnings():
                    if func == "to_period":
                        # dropping TZ
                        warnings.simplefilter("ignore", UserWarning)
                    res = getattr(c.dt, func)(*args, **kwargs)
                    exp = getattr(s.dt, func)(*args, **kwargs)

                tm.assert_equal(res, exp)

            for attr in attr_names:
                res = getattr(c.dt, attr)
                exp = getattr(s.dt, attr)

            if isinstance(res, DataFrame):
                tm.assert_frame_equal(res, exp)
            elif isinstance(res, Series):
                tm.assert_series_equal(res, exp)
            else:
                tm.assert_almost_equal(res, exp)

        invalid = Series([1, 2, 3]).astype("category")
        msg = "Can only use .dt accessor with datetimelike"

        with pytest.raises(AttributeError, match=msg):
            invalid.dt
        assert not hasattr(invalid, "str")
