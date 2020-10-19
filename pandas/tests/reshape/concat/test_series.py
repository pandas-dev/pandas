import numpy as np
import pytest

import pandas as pd
from pandas import Series

import pandas._testing as tm


class TestSeriesConcat:
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


def test_concat_empty_and_non_empty_series_regression():
    # GH 18187 regression test
    s1 = pd.Series([1])
    s2 = pd.Series([], dtype=object)
    expected = s1
    result = pd.concat([s1, s2])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("s1name,s2name", [(np.int64(190), (43, 0)), (190, (43, 0))])
def test_concat_series_name_npscalar_tuple(s1name, s2name):
    # GH21015
    s1 = pd.Series({"a": 1, "b": 2}, name=s1name)
    s2 = pd.Series({"c": 5, "d": 6}, name=s2name)
    result = pd.concat([s1, s2])
    expected = pd.Series({"a": 1, "b": 2, "c": 5, "d": 6})
    tm.assert_series_equal(result, expected)
