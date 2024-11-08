import numpy as np
import pytest

import pandas as pd
from pandas import Index
import pandas._testing as tm


class TestGetIndexer:
    @pytest.mark.parametrize(
        "method,expected",
        [
            ("pad", [-1, 0, 1, 1]),
            ("backfill", [0, 0, 1, -1]),
        ],
    )
    def test_get_indexer_strings(self, any_string_dtype, method, expected):
        expected = np.array(expected, dtype=np.intp)
        index = Index(["b", "c"], dtype=any_string_dtype)
        actual = index.get_indexer(["a", "b", "c", "d"], method=method)

        tm.assert_numpy_array_equal(actual, expected)

    def test_get_indexer_strings_raises(self, any_string_dtype):
        index = Index(["b", "c"], dtype=any_string_dtype)

        msg = "|".join(
            [
                "operation 'sub' not supported for dtype 'str",
                r"unsupported operand type\(s\) for -: 'str' and 'str'",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            index.get_indexer(["a", "b", "c", "d"], method="nearest")

        with pytest.raises(TypeError, match=msg):
            index.get_indexer(["a", "b", "c", "d"], method="pad", tolerance=2)

        with pytest.raises(TypeError, match=msg):
            index.get_indexer(
                ["a", "b", "c", "d"], method="pad", tolerance=[2, 2, 2, 2]
            )


class TestGetIndexerNonUnique:
    @pytest.mark.xfail(reason="TODO(infer_string)", strict=False)
    def test_get_indexer_non_unique_nas(self, any_string_dtype, nulls_fixture):
        index = Index(["a", "b", None], dtype=any_string_dtype)
        indexer, missing = index.get_indexer_non_unique([nulls_fixture])

        expected_indexer = np.array([2], dtype=np.intp)
        expected_missing = np.array([], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected_indexer)
        tm.assert_numpy_array_equal(missing, expected_missing)

        # actually non-unique
        index = Index(["a", None, "b", None], dtype=any_string_dtype)
        indexer, missing = index.get_indexer_non_unique([nulls_fixture])

        expected_indexer = np.array([1, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected_indexer)
        tm.assert_numpy_array_equal(missing, expected_missing)


class TestSliceLocs:
    @pytest.mark.parametrize(
        "in_slice,expected",
        [
            # error: Slice index must be an integer or None
            (pd.IndexSlice[::-1], "yxdcb"),
            (pd.IndexSlice["b":"y":-1], ""),  # type: ignore[misc]
            (pd.IndexSlice["b"::-1], "b"),  # type: ignore[misc]
            (pd.IndexSlice[:"b":-1], "yxdcb"),  # type: ignore[misc]
            (pd.IndexSlice[:"y":-1], "y"),  # type: ignore[misc]
            (pd.IndexSlice["y"::-1], "yxdcb"),  # type: ignore[misc]
            (pd.IndexSlice["y"::-4], "yb"),  # type: ignore[misc]
            # absent labels
            (pd.IndexSlice[:"a":-1], "yxdcb"),  # type: ignore[misc]
            (pd.IndexSlice[:"a":-2], "ydb"),  # type: ignore[misc]
            (pd.IndexSlice["z"::-1], "yxdcb"),  # type: ignore[misc]
            (pd.IndexSlice["z"::-3], "yc"),  # type: ignore[misc]
            (pd.IndexSlice["m"::-1], "dcb"),  # type: ignore[misc]
            (pd.IndexSlice[:"m":-1], "yx"),  # type: ignore[misc]
            (pd.IndexSlice["a":"a":-1], ""),  # type: ignore[misc]
            (pd.IndexSlice["z":"z":-1], ""),  # type: ignore[misc]
            (pd.IndexSlice["m":"m":-1], ""),  # type: ignore[misc]
        ],
    )
    def test_slice_locs_negative_step(self, in_slice, expected, any_string_dtype):
        index = Index(list("bcdxy"), dtype=any_string_dtype)

        s_start, s_stop = index.slice_locs(in_slice.start, in_slice.stop, in_slice.step)
        result = index[s_start : s_stop : in_slice.step]
        expected = Index(list(expected), dtype=any_string_dtype)
        tm.assert_index_equal(result, expected)

    def test_slice_locs_negative_step_oob(self, any_string_dtype):
        index = Index(list("bcdxy"), dtype=any_string_dtype)

        result = index[-10:5:1]
        tm.assert_index_equal(result, index)

        result = index[4:-10:-1]
        expected = Index(list("yxdcb"), dtype=any_string_dtype)
        tm.assert_index_equal(result, expected)

    def test_slice_locs_dup(self, any_string_dtype):
        index = Index(["a", "a", "b", "c", "d", "d"], dtype=any_string_dtype)
        assert index.slice_locs("a", "d") == (0, 6)
        assert index.slice_locs(end="d") == (0, 6)
        assert index.slice_locs("a", "c") == (0, 4)
        assert index.slice_locs("b", "d") == (2, 6)

        index2 = index[::-1]
        assert index2.slice_locs("d", "a") == (0, 6)
        assert index2.slice_locs(end="a") == (0, 6)
        assert index2.slice_locs("d", "b") == (0, 4)
        assert index2.slice_locs("c", "a") == (2, 6)
