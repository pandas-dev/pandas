import numpy as np

from pandas import DataFrame
from pandas.core.arrays import PandasArray


class TestToDictOfBlocks:
    def test_copy_blocks(self, float_frame):
        # GH#9607
        df = DataFrame(float_frame, copy=True)
        column = df.columns[0]

        # use the default copy=True, change a column
        blocks = df._to_dict_of_blocks(copy=True)
        for dtype, _df in blocks.items():
            if column in _df:
                _df.loc[:, column] = _df[column] + 1

        # make sure we did not change the original DataFrame
        assert not _df[column].equals(df[column])

    def test_no_copy_blocks(self, float_frame):
        # GH#9607
        df = DataFrame(float_frame, copy=True)
        column = df.columns[0]

        # use the copy=False, change a column
        blocks = df._to_dict_of_blocks(copy=False)
        for dtype, _df in blocks.items():
            if column in _df:
                _df.loc[:, column] = _df[column] + 1

        # make sure we did change the original DataFrame
        assert _df[column].equals(df[column])


def test_to_dict_of_blocks_item_cache():
    # Calling to_dict_of_blocks should not poison item_cache
    df = DataFrame({"a": [1, 2, 3, 4], "b": ["a", "b", "c", "d"]})
    df["c"] = PandasArray(np.array([1, 2, None, 3], dtype=object))
    mgr = df._mgr
    assert len(mgr.blocks) == 3  # i.e. not consolidated

    ser = df["b"]  # populations item_cache["b"]

    df._to_dict_of_blocks()

    # Check that the to_dict_of_blocks didnt break link between ser and df
    ser.values[0] = "foo"
    assert df.loc[0, "b"] == "foo"

    assert df["b"] is ser
