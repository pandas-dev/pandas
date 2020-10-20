import numpy as np
import pytest

from numpy.random import randn

import pandas as pd
from pandas import (
    Series,
    DatetimeIndex,
    concat,
    Index,
    date_range,
    MultiIndex,
    DataFrame,
)

import pandas._testing as tm


@pytest.fixture(params=[True, False])
def sort(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param


class TestSeriesConcat:
    def test_concat_series(self):
        ts = tm.makeTimeSeries()
        ts.name = "foo"

        pieces = [ts[:5], ts[5:15], ts[15:]]

        result = concat(pieces)
        tm.assert_series_equal(result, ts)
        assert result.name == ts.name

        result = concat(pieces, keys=[0, 1, 2])
        expected = ts.copy()

        ts.index = DatetimeIndex(np.array(ts.index.values, dtype="M8[ns]"))

        exp_codes = [np.repeat([0, 1, 2], [len(x) for x in pieces]), np.arange(len(ts))]
        exp_index = MultiIndex(levels=[[0, 1, 2], ts.index], codes=exp_codes)
        expected.index = exp_index
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype", ["float64", "int8", "uint8", "bool", "m8[ns]", "M8[ns]"]
    )
    def test_concat_empty_series_dtypes_match_roundtrips(self, dtype):
        dtype = np.dtype(dtype)

        result = pd.concat([Series(dtype=dtype)])
        assert result.dtype == dtype

        result = pd.concat([Series(dtype=dtype), Series(dtype=dtype)])
        assert result.dtype == dtype

    def test_concat_empty_series_dtypes_roundtrips(self):

        # round-tripping with self & like self
        dtypes = map(np.dtype, ["float64", "int8", "uint8", "bool", "m8[ns]", "M8[ns]"])

        def int_result_type(dtype, dtype2):
            typs = {dtype.kind, dtype2.kind}
            if not len(typs - {"i", "u", "b"}) and (
                dtype.kind == "i" or dtype2.kind == "i"
            ):
                return "i"
            elif not len(typs - {"u", "b"}) and (
                dtype.kind == "u" or dtype2.kind == "u"
            ):
                return "u"
            return None

        def float_result_type(dtype, dtype2):
            typs = {dtype.kind, dtype2.kind}
            if not len(typs - {"f", "i", "u"}) and (
                dtype.kind == "f" or dtype2.kind == "f"
            ):
                return "f"
            return None

        def get_result_type(dtype, dtype2):
            result = float_result_type(dtype, dtype2)
            if result is not None:
                return result
            result = int_result_type(dtype, dtype2)
            if result is not None:
                return result
            return "O"

        for dtype in dtypes:
            for dtype2 in dtypes:
                if dtype == dtype2:
                    continue

                expected = get_result_type(dtype, dtype2)
                result = pd.concat([Series(dtype=dtype), Series(dtype=dtype2)]).dtype
                assert result.kind == expected

    @pytest.mark.parametrize(
        "left,right,expected",
        [
            # booleans
            (np.bool_, np.int32, np.int32),
            (np.bool_, np.float32, np.object_),
            # datetime-like
            ("m8[ns]", np.bool_, np.object_),
            ("m8[ns]", np.int64, np.object_),
            ("M8[ns]", np.bool_, np.object_),
            ("M8[ns]", np.int64, np.object_),
            # categorical
            ("category", "category", "category"),
            ("category", "object", "object"),
        ],
    )
    def test_concat_empty_series_dtypes(self, left, right, expected):
        result = pd.concat([Series(dtype=left), Series(dtype=right)])
        assert result.dtype == expected

    def test_concat_empty_series_dtypes_triple(self):

        assert (
            pd.concat(
                [Series(dtype="M8[ns]"), Series(dtype=np.bool_), Series(dtype=np.int64)]
            ).dtype
            == np.object_
        )

    def test_concat_empty_series_dtype_category_with_array(self):
        # GH#18515
        assert (
            pd.concat(
                [Series(np.array([]), dtype="category"), Series(dtype="float64")]
            ).dtype
            == "float64"
        )

    def test_concat_empty_series_dtypes_sparse(self):
        result = pd.concat(
            [
                Series(dtype="float64").astype("Sparse"),
                Series(dtype="float64").astype("Sparse"),
            ]
        )
        assert result.dtype == "Sparse[float64]"

        result = pd.concat(
            [Series(dtype="float64").astype("Sparse"), Series(dtype="float64")]
        )
        # TODO: release-note: concat sparse dtype
        expected = pd.SparseDtype(np.float64)
        assert result.dtype == expected

        result = pd.concat(
            [Series(dtype="float64").astype("Sparse"), Series(dtype="object")]
        )
        # TODO: release-note: concat sparse dtype
        expected = pd.SparseDtype("object")
        assert result.dtype == expected

    def test_concat_empty_and_non_empty_series_regression(self):
        # GH 18187 regression test
        s1 = pd.Series([1])
        s2 = pd.Series([], dtype=object)
        expected = s1
        result = pd.concat([s1, s2])
        tm.assert_series_equal(result, expected)

    def test_concat_series_axis1(self, sort=sort):
        ts = tm.makeTimeSeries()

        pieces = [ts[:-2], ts[2:], ts[2:-2]]

        result = concat(pieces, axis=1)
        expected = DataFrame(pieces).T
        tm.assert_frame_equal(result, expected)

        result = concat(pieces, keys=["A", "B", "C"], axis=1)
        expected = DataFrame(pieces, index=["A", "B", "C"]).T
        tm.assert_frame_equal(result, expected)

        # preserve series names, #2489
        s = Series(randn(5), name="A")
        s2 = Series(randn(5), name="B")

        result = concat([s, s2], axis=1)
        expected = DataFrame({"A": s, "B": s2})
        tm.assert_frame_equal(result, expected)

        s2.name = None
        result = concat([s, s2], axis=1)
        tm.assert_index_equal(result.columns, Index(["A", 0], dtype="object"))

        # must reindex, #2603
        s = Series(randn(3), index=["c", "a", "b"], name="A")
        s2 = Series(randn(4), index=["d", "a", "b", "c"], name="B")
        result = concat([s, s2], axis=1, sort=sort)
        expected = DataFrame({"A": s, "B": s2})
        tm.assert_frame_equal(result, expected)

    def test_concat_series_axis1_names_applied(self):
        # ensure names argument is not ignored on axis=1, #23490
        s = Series([1, 2, 3])
        s2 = Series([4, 5, 6])
        result = concat([s, s2], axis=1, keys=["a", "b"], names=["A"])
        expected = DataFrame(
            [[1, 4], [2, 5], [3, 6]], columns=pd.Index(["a", "b"], name="A")
        )
        tm.assert_frame_equal(result, expected)

        result = concat([s, s2], axis=1, keys=[("a", 1), ("b", 2)], names=["A", "B"])
        expected = DataFrame(
            [[1, 4], [2, 5], [3, 6]],
            columns=MultiIndex.from_tuples([("a", 1), ("b", 2)], names=["A", "B"]),
        )
        tm.assert_frame_equal(result, expected)

    def test_concat_series_axis1_same_names_ignore_index(self):
        dates = date_range("01-Jan-2013", "01-Jan-2014", freq="MS")[0:-1]
        s1 = Series(randn(len(dates)), index=dates, name="value")
        s2 = Series(randn(len(dates)), index=dates, name="value")

        result = concat([s1, s2], axis=1, ignore_index=True)
        expected = Index([0, 1])

        tm.assert_index_equal(result.columns, expected)

    @pytest.mark.parametrize("s1name,s2name", [(np.int64(190), (43, 0)), (190, (43, 0))])
    def test_concat_series_name_npscalar_tuple(self, s1name, s2name):
        # GH21015
        s1 = pd.Series({"a": 1, "b": 2}, name=s1name)
        s2 = pd.Series({"c": 5, "d": 6}, name=s2name)
        result = pd.concat([s1, s2])
        expected = pd.Series({"a": 1, "b": 2, "c": 5, "d": 6})
        tm.assert_series_equal(result, expected)
