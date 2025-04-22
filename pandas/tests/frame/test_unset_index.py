import numpy as np
import pytest

import pandas as pd
from pandas.testing import assert_frame_equal


class TestUnsetIndex:
    # GH: 60869
    @pytest.fixture
    def df(self):
        # Fixture with custom string index
        return pd.DataFrame(
            [("bird", 389.0), ("bird", 24.0), ("mammal", 80.5), ("mammal", np.nan)],
            index=["falcon", "parrot", "lion", "monkey"],
            columns=("class", "max_speed"),
        )

    def test_unset_index_returns_default_index(self, df):
        result = df.unset_index()
        expected = pd.DataFrame(
            {
                "class": ["bird", "bird", "mammal", "mammal"],
                "max_speed": [389.0, 24.0, 80.5, np.nan],
            },
            index=pd.RangeIndex(0, 4),
        )
        assert_frame_equal(result, expected)

    def test_unset_index_inplace(self, df):
        out = df.unset_index(inplace=True)
        assert out is None
        assert isinstance(df.index, pd.RangeIndex)
        assert df.index.equals(pd.RangeIndex(0, 4))

    def test_unset_index_on_default_index(self):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
        original = df.copy()
        result = df.unset_index()
        assert_frame_equal(result, original)
        assert result is not df  # Should return copy

    def test_unset_index_on_non_default_rangeindex(self):
        # RangeIndex with start=2, step=1
        df = pd.DataFrame({"A": [1, 2, 3]}, index=pd.RangeIndex(start=2, stop=5))
        result = df.unset_index()
        assert result.index.equals(pd.RangeIndex(0, 3))

    def test_unset_index_on_rangeindex_with_step_2(self):
        # RangeIndex with non-default step
        df = pd.DataFrame(
            {"A": [1, 2, 3]}, index=pd.RangeIndex(start=0, step=2, stop=6)
        )
        result = df.unset_index()
        assert result.index.equals(pd.RangeIndex(0, 3))

    def test_unset_index_on_empty_dataframe(self):
        df = pd.DataFrame(columns=["A", "B"])
        result = df.unset_index()
        assert_frame_equal(result, df)
        assert result.index.equals(pd.RangeIndex(0, 0))

    def test_unset_index_on_multiindex(self):
        index = pd.MultiIndex.from_tuples(
            [("a", 1), ("a", 2)], names=["letter", "number"]
        )
        df = pd.DataFrame({"data": [10, 20]}, index=index)
        result = df.unset_index()
        assert result.index.equals(pd.RangeIndex(0, 2))

    def test_idempotent_on_default_index(self):
        df = pd.DataFrame({"A": [1, 2]})
        result = df.unset_index().unset_index()
        assert_frame_equal(result, df.unset_index())
