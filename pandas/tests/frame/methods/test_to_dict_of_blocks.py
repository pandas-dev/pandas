from pandas import (
    DataFrame,
    MultiIndex,
)
import pandas._testing as tm


class TestToDictOfBlocks:
    def test_no_copy_blocks(self, float_frame):
        # GH#9607
        df = DataFrame(float_frame, copy=True)
        column = df.columns[0]

        _last_df = None
        # use the copy=False, change a column
        blocks = df._to_dict_of_blocks()
        for _df in blocks.values():
            _last_df = _df
            if column in _df:
                _df.loc[:, column] = _df[column] + 1
        assert _last_df is not None and not _last_df[column].equals(df[column])


def test_set_change_dtype_slice():
    # GH#8850
    cols = MultiIndex.from_tuples([("1st", "a"), ("2nd", "b"), ("3rd", "c")])
    df = DataFrame([[1.0, 2, 3], [4.0, 5, 6]], columns=cols)
    df["2nd"] = df["2nd"] * 2.0

    blocks = df._to_dict_of_blocks()
    assert sorted(blocks.keys()) == ["float64", "int64"]
    tm.assert_frame_equal(
        blocks["float64"], DataFrame([[1.0, 4.0], [4.0, 10.0]], columns=cols[:2])
    )
    tm.assert_frame_equal(blocks["int64"], DataFrame([[3], [6]], columns=cols[2:]))
