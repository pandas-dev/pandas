import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm


class TestSeriesCombine:
    def test_combine_scalar(self):
        # GH 21248
        # Note - combine() with another Series is tested elsewhere because
        # it is used when testing operators
        s = pd.Series([i * 10 for i in range(5)])
        result = s.combine(3, lambda x, y: x + y)
        expected = pd.Series([i * 10 + 3 for i in range(5)])
        tm.assert_series_equal(result, expected)

        result = s.combine(22, lambda x, y: min(x, y))
        expected = pd.Series([min(i * 10, 22) for i in range(5)])
        tm.assert_series_equal(result, expected)

    def test_update(self):
        s = Series([1.5, np.nan, 3.0, 4.0, np.nan])
        s2 = Series([np.nan, 3.5, np.nan, 5.0])
        s.update(s2)

        expected = Series([1.5, 3.5, 3.0, 5.0, np.nan])
        tm.assert_series_equal(s, expected)

        # GH 3217
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df["c"] = np.nan

        df["c"].update(Series(["foo"], index=[0]))
        expected = DataFrame(
            [[1, np.nan, "foo"], [3, 2.0, np.nan]], columns=["a", "b", "c"]
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "other, dtype, expected",
        [
            # other is int
            ([61, 63], "int32", pd.Series([10, 61, 12], dtype="int32")),
            ([61, 63], "int64", pd.Series([10, 61, 12])),
            ([61, 63], float, pd.Series([10.0, 61.0, 12.0])),
            ([61, 63], object, pd.Series([10, 61, 12], dtype=object)),
            # other is float, but can be cast to int
            ([61.0, 63.0], "int32", pd.Series([10, 61, 12], dtype="int32")),
            ([61.0, 63.0], "int64", pd.Series([10, 61, 12])),
            ([61.0, 63.0], float, pd.Series([10.0, 61.0, 12.0])),
            ([61.0, 63.0], object, pd.Series([10, 61.0, 12], dtype=object)),
            # others is float, cannot be cast to int
            ([61.1, 63.1], "int32", pd.Series([10.0, 61.1, 12.0])),
            ([61.1, 63.1], "int64", pd.Series([10.0, 61.1, 12.0])),
            ([61.1, 63.1], float, pd.Series([10.0, 61.1, 12.0])),
            ([61.1, 63.1], object, pd.Series([10, 61.1, 12], dtype=object)),
            # other is object, cannot be cast
            ([(61,), (63,)], "int32", pd.Series([10, (61,), 12])),
            ([(61,), (63,)], "int64", pd.Series([10, (61,), 12])),
            ([(61,), (63,)], float, pd.Series([10.0, (61,), 12.0])),
            ([(61,), (63,)], object, pd.Series([10, (61,), 12])),
        ],
    )
    def test_update_dtypes(self, other, dtype, expected):

        s = Series([10, 11, 12], dtype=dtype)
        other = Series(other, index=[1, 3])
        s.update(other)

        tm.assert_series_equal(s, expected)

    def test_concat_empty_series_dtypes_roundtrips(self):

        # round-tripping with self & like self
        dtypes = map(np.dtype, ["float64", "int8", "uint8", "bool", "m8[ns]", "M8[ns]"])

        for dtype in dtypes:
            assert pd.concat([Series(dtype=dtype)]).dtype == dtype
            assert pd.concat([Series(dtype=dtype), Series(dtype=dtype)]).dtype == dtype

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

    def test_concat_empty_series_dtypes(self):

        # booleans
        assert (
            pd.concat([Series(dtype=np.bool_), Series(dtype=np.int32)]).dtype
            == np.int32
        )
        assert (
            pd.concat([Series(dtype=np.bool_), Series(dtype=np.float32)]).dtype
            == np.object_
        )

        # datetime-like
        assert (
            pd.concat([Series(dtype="m8[ns]"), Series(dtype=np.bool)]).dtype
            == np.object_
        )
        assert (
            pd.concat([Series(dtype="m8[ns]"), Series(dtype=np.int64)]).dtype
            == np.object_
        )
        assert (
            pd.concat([Series(dtype="M8[ns]"), Series(dtype=np.bool)]).dtype
            == np.object_
        )
        assert (
            pd.concat([Series(dtype="M8[ns]"), Series(dtype=np.int64)]).dtype
            == np.object_
        )
        assert (
            pd.concat(
                [Series(dtype="M8[ns]"), Series(dtype=np.bool_), Series(dtype=np.int64)]
            ).dtype
            == np.object_
        )

        # categorical
        assert (
            pd.concat([Series(dtype="category"), Series(dtype="category")]).dtype
            == "category"
        )
        # GH 18515
        assert (
            pd.concat(
                [Series(np.array([]), dtype="category"), Series(dtype="float64")]
            ).dtype
            == "float64"
        )
        assert (
            pd.concat([Series(dtype="category"), Series(dtype="object")]).dtype
            == "object"
        )

        # sparse
        # TODO: move?
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
