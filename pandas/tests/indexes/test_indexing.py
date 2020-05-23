"""
test_indexing tests the following Index methods:
    __getitem__
    get_loc
    get_value
    __contains__
    take
    where
    get_indexer
    slice_locs
    asof_locs

The corresponding tests.indexes.[index_type].test_indexing files
contain tests for the corresponding methods specific to those Index subclasses.
"""
import numpy as np
import pytest

from pandas import Float64Index, Index, Int64Index, UInt64Index
import pandas._testing as tm


class TestContains:
    @pytest.mark.parametrize(
        "index,val",
        [
            (Index([0, 1, 2]), 2),
            (Index([0, 1, "2"]), "2"),
            (Index([0, 1, 2, np.inf, 4]), 4),
            (Index([0, 1, 2, np.nan, 4]), 4),
            (Index([0, 1, 2, np.inf]), np.inf),
            (Index([0, 1, 2, np.nan]), np.nan),
        ],
    )
    def test_index_contains(self, index, val):
        assert val in index

    @pytest.mark.parametrize(
        "index,val",
        [
            (Index([0, 1, 2]), "2"),
            (Index([0, 1, "2"]), 2),
            (Index([0, 1, 2, np.inf]), 4),
            (Index([0, 1, 2, np.nan]), 4),
            (Index([0, 1, 2, np.inf]), np.nan),
            (Index([0, 1, 2, np.nan]), np.inf),
            # Checking if np.inf in Int64Index should not cause an OverflowError
            # Related to GH 16957
            (Int64Index([0, 1, 2]), np.inf),
            (Int64Index([0, 1, 2]), np.nan),
            (UInt64Index([0, 1, 2]), np.inf),
            (UInt64Index([0, 1, 2]), np.nan),
        ],
    )
    def test_index_not_contains(self, index, val):
        assert val not in index

    @pytest.mark.parametrize(
        "index,val", [(Index([0, 1, "2"]), 0), (Index([0, 1, "2"]), "2")]
    )
    def test_mixed_index_contains(self, index, val):
        # GH#19860
        assert val in index

    @pytest.mark.parametrize(
        "index,val", [(Index([0, 1, "2"]), "1"), (Index([0, 1, "2"]), 2)]
    )
    def test_mixed_index_not_contains(self, index, val):
        # GH#19860
        assert val not in index

    def test_contains_with_float_index(self):
        # GH#22085
        integer_index = Int64Index([0, 1, 2, 3])
        uinteger_index = UInt64Index([0, 1, 2, 3])
        float_index = Float64Index([0.1, 1.1, 2.2, 3.3])

        for index in (integer_index, uinteger_index):
            assert 1.1 not in index
            assert 1.0 in index
            assert 1 in index

        assert 1.1 in float_index
        assert 1.0 not in float_index
        assert 1 not in float_index


@pytest.mark.parametrize(
    "idx", [Index([1, 2, 3]), Index([0.1, 0.2, 0.3]), Index(["a", "b", "c"])]
)
def test_getitem_deprecated_float(idx):
    # https://github.com/pandas-dev/pandas/issues/34191

    with tm.assert_produces_warning(FutureWarning):
        result = idx[1.0]

    expected = idx[1]
    assert result == expected
