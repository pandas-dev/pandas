from decimal import Decimal

import numpy as np
import pytest

from pandas._libs.missing import is_matching_na
from pandas.compat.numpy import np_version_gt2_2

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
    def test_get_indexer_strings(self, method, expected):
        expected = np.array(expected, dtype=np.intp)
        index = Index(["b", "c"], dtype=object)
        actual = index.get_indexer(["a", "b", "c", "d"], method=method)

        tm.assert_numpy_array_equal(actual, expected)

    def test_get_indexer_strings_raises(self):
        index = Index(["b", "c"], dtype=object)

        msg = "|".join(
            [
                "operation 'sub' not supported for dtype 'str'",
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

    def test_get_indexer_with_NA_values(
        self, unique_nulls_fixture, unique_nulls_fixture2
    ):
        # GH#22332
        # check pairwise, that no pair of na values
        # is mangled
        if unique_nulls_fixture is unique_nulls_fixture2:
            return  # skip it, values are not unique
        arr = np.array([unique_nulls_fixture, unique_nulls_fixture2], dtype=object)
        index = Index(arr, dtype=object)
        result = index.get_indexer(
            Index(
                [unique_nulls_fixture, unique_nulls_fixture2, "Unknown"], dtype=object
            )
        )
        expected = np.array([0, 1, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_get_indexer_infer_string_missing_values(self):
        # ensure the passed list is not cast to string but to object so that
        # the None value is matched in the index
        # https://github.com/pandas-dev/pandas/issues/55834
        idx = Index(["a", "b", None], dtype="object")
        result = idx.get_indexer([None, "x"])
        expected = np.array([2, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)


class TestGetIndexerNonUnique:
    def test_get_indexer_non_unique_nas(self, nulls_fixture):
        # even though this isn't non-unique, this should still work
        index = Index(["a", "b", nulls_fixture], dtype=object)
        indexer, missing = index.get_indexer_non_unique([nulls_fixture])

        expected_indexer = np.array([2], dtype=np.intp)
        expected_missing = np.array([], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected_indexer)
        tm.assert_numpy_array_equal(missing, expected_missing)

        # actually non-unique
        index = Index(["a", nulls_fixture, "b", nulls_fixture], dtype=object)
        indexer, missing = index.get_indexer_non_unique([nulls_fixture])

        expected_indexer = np.array([1, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected_indexer)
        tm.assert_numpy_array_equal(missing, expected_missing)

        # matching-but-not-identical nans
        if is_matching_na(nulls_fixture, float("NaN")):
            index = Index(["a", float("NaN"), "b", float("NaN")], dtype=object)
            match_but_not_identical = True
        elif is_matching_na(nulls_fixture, Decimal("NaN")):
            index = Index(["a", Decimal("NaN"), "b", Decimal("NaN")], dtype=object)
            match_but_not_identical = True
        else:
            match_but_not_identical = False

        if match_but_not_identical:
            indexer, missing = index.get_indexer_non_unique([nulls_fixture])

            expected_indexer = np.array([1, 3], dtype=np.intp)
            tm.assert_numpy_array_equal(indexer, expected_indexer)
            tm.assert_numpy_array_equal(missing, expected_missing)

    @pytest.mark.filterwarnings("ignore:elementwise comp:DeprecationWarning")
    def test_get_indexer_non_unique_np_nats(self, np_nat_fixture, np_nat_fixture2):
        expected_missing = np.array([], dtype=np.intp)
        # matching-but-not-identical nats
        if is_matching_na(np_nat_fixture, np_nat_fixture2):
            # ensure nats are different objects
            index = Index(
                np.array(
                    ["2021-10-02", np_nat_fixture.copy(), np_nat_fixture2.copy()],
                    dtype=object,
                ),
                dtype=object,
            )
            # pass as index to prevent target from being casted to DatetimeIndex
            indexer, missing = index.get_indexer_non_unique(
                Index([np_nat_fixture], dtype=object)
            )
            expected_indexer = np.array([1, 2], dtype=np.intp)
            tm.assert_numpy_array_equal(indexer, expected_indexer)
            tm.assert_numpy_array_equal(missing, expected_missing)
        # dt64nat vs td64nat
        else:
            try:
                np_nat_fixture == np_nat_fixture2
            except (TypeError, OverflowError):
                # Numpy will raise on uncomparable types, like
                # np.datetime64('NaT', 'Y') and np.datetime64('NaT', 'ps')
                # https://github.com/numpy/numpy/issues/22762
                return
            index = Index(
                np.array(
                    [
                        "2021-10-02",
                        np_nat_fixture,
                        np_nat_fixture2,
                        np_nat_fixture,
                        np_nat_fixture2,
                    ],
                    dtype=object,
                ),
                dtype=object,
            )
            # pass as index to prevent target from being casted to DatetimeIndex
            indexer, missing = index.get_indexer_non_unique(
                Index([np_nat_fixture], dtype=object)
            )
            expected_indexer = np.array([1, 3], dtype=np.intp)
            tm.assert_numpy_array_equal(indexer, expected_indexer)
            tm.assert_numpy_array_equal(missing, expected_missing)


class TestMixedResolutionDatetime64:
    # GH#50690 - object-dtype Index with mixed-resolution datetime64 values
    # should respect the invariant x == y => hash(x) == hash(y)

    _xfail_np_hash = pytest.mark.xfail(
        not np_version_gt2_2,
        reason="numpy < 2.2 violates hash invariant for mixed-resolution datetime64",
    )

    @pytest.fixture
    def dt_ms(self):
        return np.datetime64(1, "ms")

    @pytest.fixture
    def dt_us(self):
        return np.datetime64(1000, "us")

    @_xfail_np_hash
    def test_contains(self, dt_ms, dt_us):
        # GH#50690
        left = Index([dt_ms], dtype=object)
        right = Index([dt_us], dtype=object)

        assert dt_us in left
        assert dt_ms in right

    @_xfail_np_hash
    def test_get_loc(self, dt_ms, dt_us):
        # GH#50690
        left = Index([dt_ms], dtype=object)
        right = Index([dt_us], dtype=object)

        assert left.get_loc(dt_us) == 0
        assert right.get_loc(dt_ms) == 0

    def test_get_indexer_monotonic(self, dt_ms, dt_us):
        # GH#50690 - monotonic case uses binary search, not hashtable
        left = Index([dt_ms], dtype=object)
        right = Index([dt_us], dtype=object)

        result = left.get_indexer(right)
        expected = np.array([0], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    @_xfail_np_hash
    def test_get_indexer_non_monotonic(self, dt_ms, dt_us):
        # GH#50690 - non-monotonic case uses hashtable; this is where
        # the hash invariant violation caused incorrect results
        sec = np.datetime64("9999-01-01", "s")
        day = np.datetime64("2016-01-01", "D")
        idx = Index([dt_ms, sec, day], dtype=object)

        result = idx.get_indexer(Index([dt_us], dtype=object))
        expected = np.array([0], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    @_xfail_np_hash
    def test_get_indexer_non_unique(self, dt_ms, dt_us):
        # GH#50690
        idx = Index([dt_ms, "foo", dt_ms], dtype=object)

        indexer, missing = idx.get_indexer_non_unique(Index([dt_us], dtype=object))

        expected_indexer = np.array([0, 2], dtype=np.intp)
        expected_missing = np.array([], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected_indexer)
        tm.assert_numpy_array_equal(missing, expected_missing)
