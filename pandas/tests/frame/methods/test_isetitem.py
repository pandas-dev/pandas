import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestDataFrameSetItem:
    def test_isetitem_ea_df(self):
        # GH#49922
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        rhs = pd.DataFrame([[11, 12], [13, 14]], dtype="Int64")

        df.isetitem([0, 1], rhs)
        expected = pd.DataFrame(
            {
                0: pd.Series([11, 13], dtype="Int64"),
                1: pd.Series([12, 14], dtype="Int64"),
                2: [3, 6],
            }
        )
        tm.assert_frame_equal(df, expected)

    def test_isetitem_ea_df_scalar_indexer(self):
        # GH#49922
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        rhs = pd.DataFrame([[11], [13]], dtype="Int64")

        df.isetitem(2, rhs)
        expected = pd.DataFrame(
            {
                0: [1, 4],
                1: [2, 5],
                2: pd.Series([11, 13], dtype="Int64"),
            }
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("loc", [slice(0, 3, 2), np.array([True, False, True])])
    def test_isetitem_df_slice_or_mask(self, loc):
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
        value = pd.DataFrame({"x": [10, 20, 30], "y": [40, 50, 60]})

        df.isetitem(loc, value)

        expected = pd.DataFrame({"a": [10, 20, 30], "b": [4, 5, 6], "c": [40, 50, 60]})
        tm.assert_frame_equal(df, expected)

    def test_isetitem_dimension_mismatch(self):
        # GH#51701
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6]})
        value = df.copy()
        with pytest.raises(ValueError, match="Got 2 positions but value has 3 columns"):
            df.isetitem([1, 2], value)

        value = df.copy()
        with pytest.raises(ValueError, match="Got 2 positions but value has 1 columns"):
            df.isetitem([1, 2], value[["a"]])

    @pytest.mark.parametrize("loc", [[1, 0], [2, 0, 1], [3, 1, 0]])
    def test_isetitem_unsorted_loc(self, loc):
        # GH#68930 Block.delete assumes increasing locs, so an out-of-order loc
        # that replaces only part of a block used to leave stale blocks behind
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8]})
        value = np.arange(100, 100 + 2 * len(loc)).reshape(len(loc), 2).T

        df.isetitem(loc, value)

        expected = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8]})
        for i, pos in enumerate(loc):
            expected.isetitem(pos, value[:, i])
        tm.assert_frame_equal(df, expected)
        # values read back fine even when corrupt; this is what fails on revert
        df.describe()

    @pytest.mark.parametrize(
        "loc, expected",
        [
            ([0, 2], {"a": [5, 5, 5], "b": [4, 5, 6], "c": [5, 5, 5]}),
            (slice(0, 3, 2), {"a": [5, 5, 5], "b": [4, 5, 6], "c": [5, 5, 5]}),
            (
                np.array([True, False, True]),
                {"a": [5, 5, 5], "b": [4, 5, 6], "c": [5, 5, 5]},
            ),
            (slice(None), {"a": [5, 5, 5], "b": [5, 5, 5], "c": [5, 5, 5]}),
        ],
    )
    def test_isetitem_scalar_fills_every_selected_position(self, loc, expected):
        # GH#71098 a scalar used to be stored as one column for several positions
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})

        df.isetitem(loc, 5)

        tm.assert_frame_equal(df, pd.DataFrame(expected))

    def test_isetitem_non_list_like_object_fills_every_selected_position(self):
        # GH#71098 a non-list-like value that is_scalar rejects, as in df[cols] = obj
        class Thing:
            pass

        thing = Thing()
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})

        df.isetitem([0, 2], thing)

        expected = pd.DataFrame(
            {"a": [thing] * 3, "b": [4, 5, 6], "c": [thing] * 3},
        )
        tm.assert_frame_equal(df, expected)

    def test_isetitem_scalar_does_not_share_one_array(self):
        # GH#71098 a later write to one filled column must not leak into another
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        df.isetitem([0, 1], 5)

        df.iloc[0, 0] = 99

        tm.assert_frame_equal(df, pd.DataFrame({"a": [99, 5, 5], "b": [5, 5, 5]}))

    @pytest.mark.parametrize(
        "value",
        [
            [10, 20, 30],
            np.array([10, 20, 30]),
            pd.Series([10, 20, 30]),
            pd.array([10, 20, 30], dtype="Int64"),
            pd.Series(pd.date_range("2030", periods=3, tz="UTC")),
        ],
    )
    @pytest.mark.parametrize("loc", [[0, 2], slice(0, 3, 2)])
    def test_isetitem_1d_value_several_positions_raises(self, value, loc):
        # GH#71098 a 1-D-only extension array (Int64) used to be stored once for
        #  several positions; the other values already raised and are controls
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
        expected = df.copy()

        with pytest.raises(ValueError, match="Got 2 positions but value has 1 columns"):
            df.isetitem(loc, value)

        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "value, ncols",
        [
            (np.array([7, 8, 9]), 1),
            (np.arange(6).reshape(3, 2), 2),
            (np.arange(12).reshape(3, 4), 4),
        ],
    )
    def test_isetitem_array_dimension_mismatch(self, value, ncols):
        # GH#68445 used to corrupt the manager or drop value columns instead
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
        expected = df.copy()
        msg = f"Got 3 positions but value has {ncols} columns"
        with pytest.raises(ValueError, match=msg):
            df.isetitem([0, 1, 2], value)
        tm.assert_frame_equal(df, expected)

    def test_isetitem_array_matching_width_still_writes(self):
        # GH#68445 positive control for the check above: the positions and the
        #  value's columns line up in order
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})

        # dtype pinned so the written columns match expected on 32-bit builds
        df.isetitem([0, 2], np.arange(6, dtype=np.int64).reshape(3, 2))

        expected = pd.DataFrame({"a": [0, 2, 4], "b": [4, 5, 6], "c": [1, 3, 5]})
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "value", [5, np.array([7, 8, 9]), np.arange(6).reshape(3, 2)]
    )
    @pytest.mark.parametrize(
        "loc", [[], slice(0, 0), np.array([False, False]), np.array([], dtype=np.intp)]
    )
    def test_isetitem_empty_loc_is_a_no_op(self, loc, value):
        # GH#68445 the column-count check must not turn an empty selection,
        #  which writes nothing, into a raise
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        expected = df.copy()

        df.isetitem(loc, value)

        tm.assert_frame_equal(df, expected)
