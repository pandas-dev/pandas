"""Test For axis=None."""
from pandas import (
    DataFrame,
)


class TestDataFrameNone:
    """TestDataFrameReductions class."""
    def test_any_all_int(self):
        """Other methods should model this behavior."""
        # DataFrame consists of ints
        df_1 = DataFrame(
            {
                "foo1": [1, 10, 15],
                "foo": [4, 5, 6]
            }
        )

        res = df_1.any(axis=None)
        exp = True
        assert res == exp

        # DataFrame consists of ints
        df_2 = DataFrame(
            {
                "foo1": [2, 4, 1],
                "foo": [1, 1, 0]
            }
        )

        res = df_2.all(axis=None)
        exp = False
        assert res == exp

    def test_min_max_int(self):
        """Testing methods min and max."""
        # DataFrame consists of ints
        df_1 = DataFrame(
            {
                "foo1": [1, 10, 15],
                "foo": [4, 5, 6]
            }
        )

        res = df_1.max(axis=None, numeric_only=True)
        exp = 15
        assert res == exp

        # DataFrame consists of ints
        df_2 = DataFrame(
            {
                "foo1": [10, 99, 13],
                "foo": [2, 5, 7]
            }
        )

        res = df_2.max(axis=None, numeric_only=True)
        exp = 99
        assert res == exp

        # DataFrame consists of ints
        df_3 = DataFrame(
            {
                "foo1": [1, 10, 15],
                "foo": [4, 5, 6]
            }
        )

        res = df_3.min(axis=None, numeric_only=True)
        exp = 1
        assert res == exp

        # DataFrame consists of ints
        df_4 = DataFrame(
            {
                "foo1": [10, 99, 13],
                "foo": [2, 5, 7]
            }
        )

        res = df_4.min(axis=None, numeric_only=True)
        exp = 2
        assert res == exp

    def test_mean_int(self):
        """Testing method mean."""
        # DataFrame consists of ints
        df_1 = DataFrame(
            {
                "foo1": [1, 10, 15],
                "foo": [4, 5, 6]
            }
        )

        res = df_1.mean(axis=None, numeric_only=True)
        exp = 41 / 6
        assert res == exp

        # DataFrame consists of ints
        df_2 = DataFrame(
            {
                "foo1": [10, 99, 13],
                "foo": [2, 5, 7]
            }
        )

        res = df_2.mean(axis=None, numeric_only=True)
        exp = 136 / 6
        assert res == exp

    def test_median_int(self):
        """Testing method median."""
        # DataFrame consists of ints
        df_1 = DataFrame(
            {
                "foo1": [1, 10, 15],
                "foo": [4, 5, 6]
            }
        )

        res = df_1.median(axis=None, numeric_only=True)
        exp = 5.5
        assert res == exp

        # DataFrame consists of ints
        df_2 = DataFrame(
            {
                "foo1": [10, 99, 13],
                "foo": [2, 5, 7],
                "foo2": [1, 6, 5]
            }
        )

        res = df_2.median(axis=None, numeric_only=True)
        exp = 6
        assert res == exp
