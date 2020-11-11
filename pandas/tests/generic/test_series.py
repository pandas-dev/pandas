from operator import methodcaller

import numpy as np
import pytest

import pandas as pd
from pandas import MultiIndex, Series, date_range
import pandas._testing as tm

from .test_generic import Generic


class TestSeries(Generic):
    _typ = Series
    _comparator = lambda self, x, y: tm.assert_series_equal(x, y)

    @pytest.mark.parametrize("func", ["rename_axis", "_set_axis_name"])
    def test_set_axis_name_mi(self, func):
        s = Series(
            [11, 21, 31],
            index=MultiIndex.from_tuples(
                [("A", x) for x in ["a", "B", "c"]], names=["l1", "l2"]
            ),
        )

        result = methodcaller(func, ["L1", "L2"])(s)
        assert s.index.name is None
        assert s.index.names == ["l1", "l2"]
        assert result.index.name is None
        assert result.index.names, ["L1", "L2"]

    def test_set_axis_name_raises(self):
        s = Series([1])
        msg = "No axis named 1 for object type Series"
        with pytest.raises(ValueError, match=msg):
            s._set_axis_name(name="a", axis=1)

    def test_get_bool_data_preserve_dtype(self):
        o = Series([True, False, True])
        result = o._get_bool_data()
        self._compare(result, o)

    def test_nonzero_single_element(self):

        # allow single item via bool method
        s = Series([True])
        assert s.bool()

        s = Series([False])
        assert not s.bool()

        msg = "The truth value of a Series is ambiguous"
        # single item nan to raise
        for s in [Series([np.nan]), Series([pd.NaT]), Series([True]), Series([False])]:
            with pytest.raises(ValueError, match=msg):
                bool(s)

        msg = "bool cannot act on a non-boolean single element Series"
        for s in [Series([np.nan]), Series([pd.NaT])]:
            with pytest.raises(ValueError, match=msg):
                s.bool()

        # multiple bool are still an error
        msg = "The truth value of a Series is ambiguous"
        for s in [Series([True, True]), Series([False, False])]:
            with pytest.raises(ValueError, match=msg):
                bool(s)
            with pytest.raises(ValueError, match=msg):
                s.bool()

        # single non-bool are an error
        for s in [Series([1]), Series([0]), Series(["a"]), Series([0.0])]:
            msg = "The truth value of a Series is ambiguous"
            with pytest.raises(ValueError, match=msg):
                bool(s)
            msg = "bool cannot act on a non-boolean single element Series"
            with pytest.raises(ValueError, match=msg):
                s.bool()

    def test_metadata_propagation_indiv_resample(self):
        # resample
        ts = Series(
            np.random.rand(1000),
            index=date_range("20130101", periods=1000, freq="s"),
            name="foo",
        )
        result = ts.resample("1T").mean()
        self.check_metadata(ts, result)

        result = ts.resample("1T").min()
        self.check_metadata(ts, result)

        result = ts.resample("1T").apply(lambda x: x.sum())
        self.check_metadata(ts, result)

    def test_metadata_propagation_indiv(self):
        # check that the metadata matches up on the resulting ops

        o = Series(range(3), range(3))
        o.name = "foo"
        o2 = Series(range(3), range(3))
        o2.name = "bar"

        result = o.T
        self.check_metadata(o, result)

        _metadata = Series._metadata
        _finalize = Series.__finalize__
        Series._metadata = ["name", "filename"]
        o.filename = "foo"
        o2.filename = "bar"

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == "concat" and name == "filename":
                    value = "+".join(
                        [getattr(o, name) for o in other.objs if getattr(o, name, None)]
                    )
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        Series.__finalize__ = finalize

        result = pd.concat([o, o2])
        assert result.filename == "foo+bar"
        assert result.name is None

        # reset
        Series._metadata = _metadata
        Series.__finalize__ = _finalize  # FIXME: use monkeypatch
