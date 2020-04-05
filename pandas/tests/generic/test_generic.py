from copy import copy, deepcopy

import numpy as np
import pytest

from pandas.compat.numpy import _np_version_under1p17

from pandas.core.dtypes.common import is_scalar

import pandas as pd
from pandas import DataFrame, MultiIndex, Series, date_range
import pandas._testing as tm
import pandas.core.common as com

# ----------------------------------------------------------------------
# Generic types test cases


class Generic:
    @property
    def _ndim(self):
        return self._typ._AXIS_LEN

    def _axes(self):
        """ return the axes for my object typ """
        return self._typ._AXIS_ORDERS

    def _construct(self, shape, value=None, dtype=None, **kwargs):
        """
        construct an object for the given shape
        if value is specified use that if its a scalar
        if value is an array, repeat it as needed
        """
        if isinstance(shape, int):
            shape = tuple([shape] * self._ndim)
        if value is not None:
            if is_scalar(value):
                if value == "empty":
                    arr = None
                    dtype = np.float64

                    # remove the info axis
                    kwargs.pop(self._typ._info_axis_name, None)
                else:
                    arr = np.empty(shape, dtype=dtype)
                    arr.fill(value)
            else:
                fshape = np.prod(shape)
                arr = value.ravel()
                new_shape = fshape / arr.shape[0]
                if fshape % arr.shape[0] != 0:
                    raise Exception("invalid value passed in _construct")

                arr = np.repeat(arr, new_shape).reshape(shape)
        else:
            arr = np.random.randn(*shape)
        return self._typ(arr, dtype=dtype, **kwargs)

    def _compare(self, result, expected):
        self._comparator(result, expected)

    def test_rename(self):

        # single axis
        idx = list("ABCD")
        # relabeling values passed into self.rename
        args = [
            str.lower,
            {x: x.lower() for x in idx},
            Series({x: x.lower() for x in idx}),
        ]

        for axis in self._axes():
            kwargs = {axis: idx}
            obj = self._construct(4, **kwargs)

            for arg in args:
                # rename a single axis
                result = obj.rename(**{axis: arg})
                expected = obj.copy()
                setattr(expected, axis, list("abcd"))
                self._compare(result, expected)

        # multiple axes at once

    def test_get_numeric_data(self):

        n = 4
        kwargs = {self._typ._AXIS_NAMES[i]: list(range(n)) for i in range(self._ndim)}

        # get the numeric data
        o = self._construct(n, **kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # non-inclusion
        result = o._get_bool_data()
        expected = self._construct(n, value="empty", **kwargs)
        self._compare(result, expected)

        # get the bool data
        arr = np.array([True, True, False, True])
        o = self._construct(n, value=arr, **kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # _get_numeric_data is includes _get_bool_data, so can't test for
        # non-inclusion

    def test_nonzero(self):

        # GH 4633
        # look at the boolean/nonzero behavior for objects
        obj = self._construct(shape=4)
        msg = f"The truth value of a {self._typ.__name__} is ambiguous"
        with pytest.raises(ValueError, match=msg):
            bool(obj == 0)
        with pytest.raises(ValueError, match=msg):
            bool(obj == 1)
        with pytest.raises(ValueError, match=msg):
            bool(obj)

        obj = self._construct(shape=4, value=1)
        with pytest.raises(ValueError, match=msg):
            bool(obj == 0)
        with pytest.raises(ValueError, match=msg):
            bool(obj == 1)
        with pytest.raises(ValueError, match=msg):
            bool(obj)

        obj = self._construct(shape=4, value=np.nan)
        with pytest.raises(ValueError, match=msg):
            bool(obj == 0)
        with pytest.raises(ValueError, match=msg):
            bool(obj == 1)
        with pytest.raises(ValueError, match=msg):
            bool(obj)

        # empty
        obj = self._construct(shape=0)
        with pytest.raises(ValueError, match=msg):
            bool(obj)

        # invalid behaviors

        obj1 = self._construct(shape=4, value=1)
        obj2 = self._construct(shape=4, value=1)

        with pytest.raises(ValueError, match=msg):
            if obj1:
                pass

        with pytest.raises(ValueError, match=msg):
            obj1 and obj2
        with pytest.raises(ValueError, match=msg):
            obj1 or obj2
        with pytest.raises(ValueError, match=msg):
            not obj1

    def test_downcast(self):
        # test close downcasting

        o = self._construct(shape=4, value=9, dtype=np.int64)
        result = o.copy()
        result._data = o._data.downcast()
        self._compare(result, o)

        o = self._construct(shape=4, value=9.5)
        result = o.copy()
        result._data = o._data.downcast()
        self._compare(result, o)

    def test_constructor_compound_dtypes(self):
        # see gh-5191
        # Compound dtypes should raise NotImplementedError.

        def f(dtype):
            return self._construct(shape=3, value=1, dtype=dtype)

        msg = (
            "compound dtypes are not implemented "
            f"in the {self._typ.__name__} constructor"
        )

        with pytest.raises(NotImplementedError, match=msg):
            f([("A", "datetime64[h]"), ("B", "str"), ("C", "int32")])

        # these work (though results may be unexpected)
        f("int64")
        f("float64")
        f("M8[ns]")

    def check_metadata(self, x, y=None):
        for m in x._metadata:
            v = getattr(x, m, None)
            if y is None:
                assert v is None
            else:
                assert v == getattr(y, m, None)

    def test_metadata_propagation(self):
        # check that the metadata matches up on the resulting ops

        o = self._construct(shape=3)
        o.name = "foo"
        o2 = self._construct(shape=3)
        o2.name = "bar"

        # ----------
        # preserving
        # ----------

        # simple ops with scalars
        for op in ["__add__", "__sub__", "__truediv__", "__mul__"]:
            result = getattr(o, op)(1)
            self.check_metadata(o, result)

        # ops with like
        for op in ["__add__", "__sub__", "__truediv__", "__mul__"]:
            result = getattr(o, op)(o)
            self.check_metadata(o, result)

        # simple boolean
        for op in ["__eq__", "__le__", "__ge__"]:
            v1 = getattr(o, op)(o)
            self.check_metadata(o, v1)
            self.check_metadata(o, v1 & v1)
            self.check_metadata(o, v1 | v1)

        # combine_first
        result = o.combine_first(o2)
        self.check_metadata(o, result)

        # ---------------------------
        # non-preserving (by default)
        # ---------------------------

        # add non-like
        result = o + o2
        self.check_metadata(result)

        # simple boolean
        for op in ["__eq__", "__le__", "__ge__"]:

            # this is a name matching op
            v1 = getattr(o, op)(o)
            v2 = getattr(o, op)(o2)
            self.check_metadata(v2)
            self.check_metadata(v1 & v2)
            self.check_metadata(v1 | v2)

    def test_head_tail(self, indices):
        # GH5370

        o = self._construct(shape=len(indices))

        axis = o._get_axis_name(0)
        setattr(o, axis, indices)

        o.head()

        self._compare(o.head(), o.iloc[:5])
        self._compare(o.tail(), o.iloc[-5:])

        # 0-len
        self._compare(o.head(0), o.iloc[0:0])
        self._compare(o.tail(0), o.iloc[0:0])

        # bounded
        self._compare(o.head(len(o) + 1), o)
        self._compare(o.tail(len(o) + 1), o)

        # neg index
        self._compare(o.head(-3), o.head(len(indices) - 3))
        self._compare(o.tail(-3), o.tail(len(indices) - 3))

    def test_sample(self):
        # Fixes issue: 2419

        o = self._construct(shape=10)

        ###
        # Check behavior of random_state argument
        ###

        # Check for stability when receives seed or random state -- run 10
        # times.
        for test in range(10):
            seed = np.random.randint(0, 100)
            self._compare(
                o.sample(n=4, random_state=seed), o.sample(n=4, random_state=seed)
            )

            self._compare(
                o.sample(frac=0.7, random_state=seed),
                o.sample(frac=0.7, random_state=seed),
            )

            self._compare(
                o.sample(n=4, random_state=np.random.RandomState(test)),
                o.sample(n=4, random_state=np.random.RandomState(test)),
            )

            self._compare(
                o.sample(frac=0.7, random_state=np.random.RandomState(test)),
                o.sample(frac=0.7, random_state=np.random.RandomState(test)),
            )

            self._compare(
                o.sample(
                    frac=2, replace=True, random_state=np.random.RandomState(test)
                ),
                o.sample(
                    frac=2, replace=True, random_state=np.random.RandomState(test)
                ),
            )

            os1, os2 = [], []
            for _ in range(2):
                np.random.seed(test)
                os1.append(o.sample(n=4))
                os2.append(o.sample(frac=0.7))
            self._compare(*os1)
            self._compare(*os2)

        # Check for error when random_state argument invalid.
        with pytest.raises(ValueError):
            o.sample(random_state="astring!")

        ###
        # Check behavior of `frac` and `N`
        ###

        # Giving both frac and N throws error
        with pytest.raises(ValueError):
            o.sample(n=3, frac=0.3)

        # Check that raises right error for negative lengths
        with pytest.raises(ValueError):
            o.sample(n=-3)
        with pytest.raises(ValueError):
            o.sample(frac=-0.3)

        # Make sure float values of `n` give error
        with pytest.raises(ValueError):
            o.sample(n=3.2)

        # Check lengths are right
        assert len(o.sample(n=4) == 4)
        assert len(o.sample(frac=0.34) == 3)
        assert len(o.sample(frac=0.36) == 4)

        ###
        # Check weights
        ###

        # Weight length must be right
        with pytest.raises(ValueError):
            o.sample(n=3, weights=[0, 1])

        with pytest.raises(ValueError):
            bad_weights = [0.5] * 11
            o.sample(n=3, weights=bad_weights)

        with pytest.raises(ValueError):
            bad_weight_series = Series([0, 0, 0.2])
            o.sample(n=4, weights=bad_weight_series)

        # Check won't accept negative weights
        with pytest.raises(ValueError):
            bad_weights = [-0.1] * 10
            o.sample(n=3, weights=bad_weights)

        # Check inf and -inf throw errors:
        with pytest.raises(ValueError):
            weights_with_inf = [0.1] * 10
            weights_with_inf[0] = np.inf
            o.sample(n=3, weights=weights_with_inf)

        with pytest.raises(ValueError):
            weights_with_ninf = [0.1] * 10
            weights_with_ninf[0] = -np.inf
            o.sample(n=3, weights=weights_with_ninf)

        # All zeros raises errors
        zero_weights = [0] * 10
        with pytest.raises(ValueError):
            o.sample(n=3, weights=zero_weights)

        # All missing weights
        nan_weights = [np.nan] * 10
        with pytest.raises(ValueError):
            o.sample(n=3, weights=nan_weights)

        # Check np.nan are replaced by zeros.
        weights_with_nan = [np.nan] * 10
        weights_with_nan[5] = 0.5
        self._compare(o.sample(n=1, axis=0, weights=weights_with_nan), o.iloc[5:6])

        # Check None are also replaced by zeros.
        weights_with_None = [None] * 10
        weights_with_None[5] = 0.5
        self._compare(o.sample(n=1, axis=0, weights=weights_with_None), o.iloc[5:6])

    def test_sample_upsampling_without_replacement(self):
        # GH27451

        df = pd.DataFrame({"A": list("abc")})
        msg = (
            "Replace has to be set to `True` when "
            "upsampling the population `frac` > 1."
        )
        with pytest.raises(ValueError, match=msg):
            df.sample(frac=2, replace=False)

    def test_sample_is_copy(self):
        # GH-27357, GH-30784: ensure the result of sample is an actual copy and
        # doesn't track the parent dataframe / doesn't give SettingWithCopy warnings
        df = pd.DataFrame(np.random.randn(10, 3), columns=["a", "b", "c"])
        df2 = df.sample(3)

        with tm.assert_produces_warning(None):
            df2["d"] = 1

    def test_size_compat(self):
        # GH8846
        # size property should be defined

        o = self._construct(shape=10)
        assert o.size == np.prod(o.shape)
        assert o.size == 10 ** len(o.axes)

    def test_split_compat(self):
        # xref GH8846
        o = self._construct(shape=10)
        assert len(np.array_split(o, 5)) == 5
        assert len(np.array_split(o, 2)) == 2

    # See gh-12301
    def test_stat_unexpected_keyword(self):
        obj = self._construct(5)
        starwars = "Star Wars"
        errmsg = "unexpected keyword"

        with pytest.raises(TypeError, match=errmsg):
            obj.max(epic=starwars)  # stat_function
        with pytest.raises(TypeError, match=errmsg):
            obj.var(epic=starwars)  # stat_function_ddof
        with pytest.raises(TypeError, match=errmsg):
            obj.sum(epic=starwars)  # cum_function
        with pytest.raises(TypeError, match=errmsg):
            obj.any(epic=starwars)  # logical_function

    @pytest.mark.parametrize("func", ["sum", "cumsum", "any", "var"])
    def test_api_compat(self, func):

        # GH 12021
        # compat for __name__, __qualname__

        obj = self._construct(5)
        f = getattr(obj, func)
        assert f.__name__ == func
        assert f.__qualname__.endswith(func)

    def test_stat_non_defaults_args(self):
        obj = self._construct(5)
        out = np.array([0])
        errmsg = "the 'out' parameter is not supported"

        with pytest.raises(ValueError, match=errmsg):
            obj.max(out=out)  # stat_function
        with pytest.raises(ValueError, match=errmsg):
            obj.var(out=out)  # stat_function_ddof
        with pytest.raises(ValueError, match=errmsg):
            obj.sum(out=out)  # cum_function
        with pytest.raises(ValueError, match=errmsg):
            obj.any(out=out)  # logical_function

    def test_truncate_out_of_bounds(self):
        # GH11382

        # small
        shape = [int(2e3)] + ([1] * (self._ndim - 1))
        small = self._construct(shape, dtype="int8", value=1)
        self._compare(small.truncate(), small)
        self._compare(small.truncate(before=0, after=3e3), small)
        self._compare(small.truncate(before=-1, after=2e3), small)

        # big
        shape = [int(2e6)] + ([1] * (self._ndim - 1))
        big = self._construct(shape, dtype="int8", value=1)
        self._compare(big.truncate(), big)
        self._compare(big.truncate(before=0, after=3e6), big)
        self._compare(big.truncate(before=-1, after=2e6), big)

    @pytest.mark.parametrize(
        "func",
        [copy, deepcopy, lambda x: x.copy(deep=False), lambda x: x.copy(deep=True)],
    )
    @pytest.mark.parametrize("shape", [0, 1, 2])
    def test_copy_and_deepcopy(self, shape, func):
        # GH 15444
        obj = self._construct(shape)
        obj_copy = func(obj)
        assert obj_copy is not obj
        self._compare(obj_copy, obj)

    @pytest.mark.parametrize(
        "periods,fill_method,limit,exp",
        [
            (1, "ffill", None, [np.nan, np.nan, np.nan, 1, 1, 1.5, 0, 0]),
            (1, "ffill", 1, [np.nan, np.nan, np.nan, 1, 1, 1.5, 0, np.nan]),
            (1, "bfill", None, [np.nan, 0, 0, 1, 1, 1.5, np.nan, np.nan]),
            (1, "bfill", 1, [np.nan, np.nan, 0, 1, 1, 1.5, np.nan, np.nan]),
            (-1, "ffill", None, [np.nan, np.nan, -0.5, -0.5, -0.6, 0, 0, np.nan]),
            (-1, "ffill", 1, [np.nan, np.nan, -0.5, -0.5, -0.6, 0, np.nan, np.nan]),
            (-1, "bfill", None, [0, 0, -0.5, -0.5, -0.6, np.nan, np.nan, np.nan]),
            (-1, "bfill", 1, [np.nan, 0, -0.5, -0.5, -0.6, np.nan, np.nan, np.nan]),
        ],
    )
    def test_pct_change(self, periods, fill_method, limit, exp):
        vals = [np.nan, np.nan, 1, 2, 4, 10, np.nan, np.nan]
        obj = self._typ(vals)
        func = getattr(obj, "pct_change")
        res = func(periods=periods, fill_method=fill_method, limit=limit)
        if type(obj) is DataFrame:
            tm.assert_frame_equal(res, DataFrame(exp))
        else:
            tm.assert_series_equal(res, Series(exp))


class TestNDFrame:
    # tests that don't fit elsewhere

    def test_sample(sel):
        # Fixes issue: 2419
        # additional specific object based tests

        # A few dataframe test with degenerate weights.
        easy_weight_list = [0] * 10
        easy_weight_list[5] = 1

        df = pd.DataFrame(
            {
                "col1": range(10, 20),
                "col2": range(20, 30),
                "colString": ["a"] * 10,
                "easyweights": easy_weight_list,
            }
        )
        sample1 = df.sample(n=1, weights="easyweights")
        tm.assert_frame_equal(sample1, df.iloc[5:6])

        # Ensure proper error if string given as weight for Series or
        # DataFrame with axis = 1.
        s = Series(range(10))
        with pytest.raises(ValueError):
            s.sample(n=3, weights="weight_column")

        with pytest.raises(ValueError):
            df.sample(n=1, weights="weight_column", axis=1)

        # Check weighting key error
        with pytest.raises(
            KeyError, match="'String passed to weights not a valid column'"
        ):
            df.sample(n=3, weights="not_a_real_column_name")

        # Check that re-normalizes weights that don't sum to one.
        weights_less_than_1 = [0] * 10
        weights_less_than_1[0] = 0.5
        tm.assert_frame_equal(df.sample(n=1, weights=weights_less_than_1), df.iloc[:1])

        ###
        # Test axis argument
        ###

        # Test axis argument
        df = pd.DataFrame({"col1": range(10), "col2": ["a"] * 10})
        second_column_weight = [0, 1]
        tm.assert_frame_equal(
            df.sample(n=1, axis=1, weights=second_column_weight), df[["col2"]]
        )

        # Different axis arg types
        tm.assert_frame_equal(
            df.sample(n=1, axis="columns", weights=second_column_weight), df[["col2"]]
        )

        weight = [0] * 10
        weight[5] = 0.5
        tm.assert_frame_equal(df.sample(n=1, axis="rows", weights=weight), df.iloc[5:6])
        tm.assert_frame_equal(
            df.sample(n=1, axis="index", weights=weight), df.iloc[5:6]
        )

        # Check out of range axis values
        with pytest.raises(ValueError):
            df.sample(n=1, axis=2)

        with pytest.raises(ValueError):
            df.sample(n=1, axis="not_a_name")

        with pytest.raises(ValueError):
            s = pd.Series(range(10))
            s.sample(n=1, axis=1)

        # Test weight length compared to correct axis
        with pytest.raises(ValueError):
            df.sample(n=1, axis=1, weights=[0.5] * 10)

        # Check weights with axis = 1
        easy_weight_list = [0] * 3
        easy_weight_list[2] = 1

        df = pd.DataFrame(
            {"col1": range(10, 20), "col2": range(20, 30), "colString": ["a"] * 10}
        )
        sample1 = df.sample(n=1, axis=1, weights=easy_weight_list)
        tm.assert_frame_equal(sample1, df[["colString"]])

        # Test default axes
        tm.assert_frame_equal(
            df.sample(n=3, random_state=42), df.sample(n=3, axis=0, random_state=42)
        )

        # Test that function aligns weights with frame
        df = DataFrame({"col1": [5, 6, 7], "col2": ["a", "b", "c"]}, index=[9, 5, 3])
        s = Series([1, 0, 0], index=[3, 5, 9])
        tm.assert_frame_equal(df.loc[[3]], df.sample(1, weights=s))

        # Weights have index values to be dropped because not in
        # sampled DataFrame
        s2 = Series([0.001, 0, 10000], index=[3, 5, 10])
        tm.assert_frame_equal(df.loc[[3]], df.sample(1, weights=s2))

        # Weights have empty values to be filed with zeros
        s3 = Series([0.01, 0], index=[3, 5])
        tm.assert_frame_equal(df.loc[[3]], df.sample(1, weights=s3))

        # No overlap in weight and sampled DataFrame indices
        s4 = Series([1, 0], index=[1, 2])
        with pytest.raises(ValueError):
            df.sample(1, weights=s4)

    @pytest.mark.parametrize(
        "func_str,arg",
        [
            ("np.array", [2, 3, 1, 0]),
            pytest.param(
                "np.random.MT19937",
                3,
                marks=pytest.mark.skipif(_np_version_under1p17, reason="NumPy<1.17"),
            ),
            pytest.param(
                "np.random.PCG64",
                11,
                marks=pytest.mark.skipif(_np_version_under1p17, reason="NumPy<1.17"),
            ),
        ],
    )
    def test_sample_random_state(self, func_str, arg):
        # GH32503
        df = pd.DataFrame({"col1": range(10, 20), "col2": range(20, 30)})
        result = df.sample(n=3, random_state=eval(func_str)(arg))
        expected = df.sample(n=3, random_state=com.random_state(eval(func_str)(arg)))
        tm.assert_frame_equal(result, expected)

    def test_squeeze(self):
        # noop
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(), tm.makeObjectSeries()]:
            tm.assert_series_equal(s.squeeze(), s)
        for df in [tm.makeTimeDataFrame()]:
            tm.assert_frame_equal(df.squeeze(), df)

        # squeezing
        df = tm.makeTimeDataFrame().reindex(columns=["A"])
        tm.assert_series_equal(df.squeeze(), df["A"])

        # don't fail with 0 length dimensions GH11229 & GH8999
        empty_series = Series([], name="five", dtype=np.float64)
        empty_frame = DataFrame([empty_series])
        tm.assert_series_equal(empty_series, empty_series.squeeze())
        tm.assert_series_equal(empty_series, empty_frame.squeeze())

        # axis argument
        df = tm.makeTimeDataFrame(nper=1).iloc[:, :1]
        assert df.shape == (1, 1)
        tm.assert_series_equal(df.squeeze(axis=0), df.iloc[0])
        tm.assert_series_equal(df.squeeze(axis="index"), df.iloc[0])
        tm.assert_series_equal(df.squeeze(axis=1), df.iloc[:, 0])
        tm.assert_series_equal(df.squeeze(axis="columns"), df.iloc[:, 0])
        assert df.squeeze() == df.iloc[0, 0]
        msg = "No axis named 2 for object type DataFrame"
        with pytest.raises(ValueError, match=msg):
            df.squeeze(axis=2)
        msg = "No axis named x for object type DataFrame"
        with pytest.raises(ValueError, match=msg):
            df.squeeze(axis="x")

        df = tm.makeTimeDataFrame(3)
        tm.assert_frame_equal(df.squeeze(axis=0), df)

    def test_numpy_squeeze(self):
        s = tm.makeFloatSeries()
        tm.assert_series_equal(np.squeeze(s), s)

        df = tm.makeTimeDataFrame().reindex(columns=["A"])
        tm.assert_series_equal(np.squeeze(df), df["A"])

    def test_transpose(self):
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(), tm.makeObjectSeries()]:
            # calls implementation in pandas/core/base.py
            tm.assert_series_equal(s.transpose(), s)
        for df in [tm.makeTimeDataFrame()]:
            tm.assert_frame_equal(df.transpose().transpose(), df)

    def test_numpy_transpose(self):
        msg = "the 'axes' parameter is not supported"

        s = tm.makeFloatSeries()
        tm.assert_series_equal(np.transpose(s), s)

        with pytest.raises(ValueError, match=msg):
            np.transpose(s, axes=1)

        df = tm.makeTimeDataFrame()
        tm.assert_frame_equal(np.transpose(np.transpose(df)), df)

        with pytest.raises(ValueError, match=msg):
            np.transpose(df, axes=1)

    def test_take(self):
        indices = [1, 5, -2, 6, 3, -1]
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(), tm.makeObjectSeries()]:
            out = s.take(indices)
            expected = Series(
                data=s.values.take(indices), index=s.index.take(indices), dtype=s.dtype
            )
            tm.assert_series_equal(out, expected)
        for df in [tm.makeTimeDataFrame()]:
            out = df.take(indices)
            expected = DataFrame(
                data=df.values.take(indices, axis=0),
                index=df.index.take(indices),
                columns=df.columns,
            )
            tm.assert_frame_equal(out, expected)

    def test_take_invalid_kwargs(self):
        indices = [-3, 2, 0, 1]
        s = tm.makeFloatSeries()
        df = tm.makeTimeDataFrame()

        for obj in (s, df):
            msg = r"take\(\) got an unexpected keyword argument 'foo'"
            with pytest.raises(TypeError, match=msg):
                obj.take(indices, foo=2)

            msg = "the 'out' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                obj.take(indices, out=indices)

            msg = "the 'mode' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                obj.take(indices, mode="clip")

    @pytest.mark.parametrize("is_copy", [True, False])
    def test_depr_take_kwarg_is_copy(self, is_copy):
        # GH 27357
        df = DataFrame({"A": [1, 2, 3]})
        msg = (
            "is_copy is deprecated and will be removed in a future version. "
            "'take' always returns a copy, so there is no need to specify this."
        )
        with tm.assert_produces_warning(FutureWarning) as w:
            df.take([0, 1], is_copy=is_copy)

        assert w[0].message.args[0] == msg

        s = Series([1, 2, 3])
        with tm.assert_produces_warning(FutureWarning):
            s.take([0, 1], is_copy=is_copy)

    def test_equals(self):
        s1 = pd.Series([1, 2, 3], index=[0, 2, 1])
        s2 = s1.copy()
        assert s1.equals(s2)

        s1[1] = 99
        assert not s1.equals(s2)

        # NaNs compare as equal
        s1 = pd.Series([1, np.nan, 3, np.nan], index=[0, 2, 1, 3])
        s2 = s1.copy()
        assert s1.equals(s2)

        s2[0] = 9.9
        assert not s1.equals(s2)

        idx = MultiIndex.from_tuples([(0, "a"), (1, "b"), (2, "c")])
        s1 = Series([1, 2, np.nan], index=idx)
        s2 = s1.copy()
        assert s1.equals(s2)

        # Add object dtype column with nans
        index = np.random.random(10)
        df1 = DataFrame(np.random.random(10), index=index, columns=["floats"])
        df1["text"] = "the sky is so blue. we could use more chocolate.".split()
        df1["start"] = date_range("2000-1-1", periods=10, freq="T")
        df1["end"] = date_range("2000-1-1", periods=10, freq="D")
        df1["diff"] = df1["end"] - df1["start"]
        df1["bool"] = np.arange(10) % 3 == 0
        df1.loc[::2] = np.nan
        df2 = df1.copy()
        assert df1["text"].equals(df2["text"])
        assert df1["start"].equals(df2["start"])
        assert df1["end"].equals(df2["end"])
        assert df1["diff"].equals(df2["diff"])
        assert df1["bool"].equals(df2["bool"])
        assert df1.equals(df2)
        assert not df1.equals(object)

        # different dtype
        different = df1.copy()
        different["floats"] = different["floats"].astype("float32")
        assert not df1.equals(different)

        # different index
        different_index = -index
        different = df2.set_index(different_index)
        assert not df1.equals(different)

        # different columns
        different = df2.copy()
        different.columns = df2.columns[::-1]
        assert not df1.equals(different)

        # DatetimeIndex
        index = pd.date_range("2000-1-1", periods=10, freq="T")
        df1 = df1.set_index(index)
        df2 = df1.copy()
        assert df1.equals(df2)

        # MultiIndex
        df3 = df1.set_index(["text"], append=True)
        df2 = df1.set_index(["text"], append=True)
        assert df3.equals(df2)

        df2 = df1.set_index(["floats"], append=True)
        assert not df3.equals(df2)

        # NaN in index
        df3 = df1.set_index(["floats"], append=True)
        df2 = df1.set_index(["floats"], append=True)
        assert df3.equals(df2)

        # GH 8437
        a = pd.Series([False, np.nan])
        b = pd.Series([False, np.nan])
        c = pd.Series(index=range(2), dtype=object)
        d = c.copy()
        e = c.copy()
        f = c.copy()
        c[:-1] = d[:-1] = e[0] = f[0] = False
        assert a.equals(a)
        assert a.equals(b)
        assert a.equals(c)
        assert a.equals(d)
        assert a.equals(e)
        assert e.equals(f)

    def test_pipe(self):
        df = DataFrame({"A": [1, 2, 3]})
        f = lambda x, y: x ** y
        result = df.pipe(f, 2)
        expected = DataFrame({"A": [1, 4, 9]})
        tm.assert_frame_equal(result, expected)

        result = df.A.pipe(f, 2)
        tm.assert_series_equal(result, expected.A)

    def test_pipe_tuple(self):
        df = DataFrame({"A": [1, 2, 3]})
        f = lambda x, y: y
        result = df.pipe((f, "y"), 0)
        tm.assert_frame_equal(result, df)

        result = df.A.pipe((f, "y"), 0)
        tm.assert_series_equal(result, df.A)

    def test_pipe_tuple_error(self):
        df = DataFrame({"A": [1, 2, 3]})
        f = lambda x, y: y
        with pytest.raises(ValueError):
            df.pipe((f, "y"), x=1, y=0)

        with pytest.raises(ValueError):
            df.A.pipe((f, "y"), x=1, y=0)

    @pytest.mark.parametrize("box", [pd.Series, pd.DataFrame])
    def test_axis_classmethods(self, box):
        obj = box(dtype=object)
        values = (
            list(box._AXIS_NAMES.keys())
            + list(box._AXIS_NUMBERS.keys())
            + list(box._AXIS_ALIASES.keys())
        )
        for v in values:
            assert obj._get_axis_number(v) == box._get_axis_number(v)
            assert obj._get_axis_name(v) == box._get_axis_name(v)
            assert obj._get_block_manager_axis(v) == box._get_block_manager_axis(v)
