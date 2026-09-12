import pytest

import pandas as pd
import pandas._testing as tm


class TestDropLevel:
    def test_droplevel(self, frame_or_series):
        # GH#20342
        cols = pd.MultiIndex.from_tuples(
            [("c", "e"), ("d", "f")], names=["level_1", "level_2"]
        )
        mi = pd.MultiIndex.from_tuples([(1, 2), (5, 6), (9, 10)], names=["a", "b"])
        df = pd.DataFrame([[3, 4], [7, 8], [11, 12]], index=mi, columns=cols)
        if frame_or_series is not pd.DataFrame:
            df = df.iloc[:, 0]

        # test that dropping of a level in index works
        expected = df.reset_index("a", drop=True)
        result = df.droplevel("a", axis="index")
        tm.assert_equal(result, expected)

        if frame_or_series is pd.DataFrame:
            # test that dropping of a level in columns works
            expected = df.copy()
            expected.columns = pd.Index(["c", "d"], name="level_1")
            result = df.droplevel("level_2", axis="columns")
            tm.assert_equal(result, expected)
        else:
            # test that droplevel raises ValueError on axis != 0
            with pytest.raises(ValueError, match="No axis named columns"):
                df.droplevel(1, axis="columns")
