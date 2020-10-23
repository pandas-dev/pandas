import numpy as np
import pytest

from pandas import Float64Index, Int64Index, Series, UInt64Index
import pandas._testing as tm


@pytest.fixture
def index_large():
    # large values used in UInt64Index tests where no compat needed with Int64/Float64
    large = [2 ** 63, 2 ** 63 + 10, 2 ** 63 + 15, 2 ** 63 + 20, 2 ** 63 + 25]
    return UInt64Index(large)


class TestGetLoc:
    def test_get_loc_float64(self):
        idx = Float64Index([0.0, 1.0, 2.0])
        for method in [None, "pad", "backfill", "nearest"]:
            assert idx.get_loc(1, method) == 1
            if method is not None:
                assert idx.get_loc(1, method, tolerance=0) == 1

        for method, loc in [("pad", 1), ("backfill", 2), ("nearest", 1)]:
            assert idx.get_loc(1.1, method) == loc
            assert idx.get_loc(1.1, method, tolerance=0.9) == loc

        with pytest.raises(KeyError, match="^'foo'$"):
            idx.get_loc("foo")
        with pytest.raises(KeyError, match=r"^1\.5$"):
            idx.get_loc(1.5)
        with pytest.raises(KeyError, match=r"^1\.5$"):
            idx.get_loc(1.5, method="pad", tolerance=0.1)
        with pytest.raises(KeyError, match="^True$"):
            idx.get_loc(True)
        with pytest.raises(KeyError, match="^False$"):
            idx.get_loc(False)

        with pytest.raises(ValueError, match="must be numeric"):
            idx.get_loc(1.4, method="nearest", tolerance="foo")

        with pytest.raises(ValueError, match="must contain numeric elements"):
            idx.get_loc(1.4, method="nearest", tolerance=np.array(["foo"]))

        with pytest.raises(
            ValueError, match="tolerance size must match target index size"
        ):
            idx.get_loc(1.4, method="nearest", tolerance=np.array([1, 2]))

    def test_get_loc_na(self):
        idx = Float64Index([np.nan, 1, 2])
        assert idx.get_loc(1) == 1
        assert idx.get_loc(np.nan) == 0

        idx = Float64Index([np.nan, 1, np.nan])
        assert idx.get_loc(1) == 1

        # FIXME: dont leave commented-out
        # representable by slice [0:2:2]
        # pytest.raises(KeyError, idx.slice_locs, np.nan)
        sliced = idx.slice_locs(np.nan)
        assert isinstance(sliced, tuple)
        assert sliced == (0, 3)

        # not representable by slice
        idx = Float64Index([np.nan, 1, np.nan, np.nan])
        assert idx.get_loc(1) == 1
        msg = "'Cannot get left slice bound for non-unique label: nan"
        with pytest.raises(KeyError, match=msg):
            idx.slice_locs(np.nan)

    def test_get_loc_missing_nan(self):
        # GH#8569
        idx = Float64Index([1, 2])
        assert idx.get_loc(1) == 0
        with pytest.raises(KeyError, match=r"^3$"):
            idx.get_loc(3)
        with pytest.raises(KeyError, match="^nan$"):
            idx.get_loc(np.nan)
        with pytest.raises(TypeError, match=r"'\[nan\]' is an invalid key"):
            # listlike/non-hashable raises TypeError
            idx.get_loc([np.nan])


class TestGetIndexer:
    def test_get_indexer_float64(self):
        idx = Float64Index([0.0, 1.0, 2.0])
        tm.assert_numpy_array_equal(
            idx.get_indexer(idx), np.array([0, 1, 2], dtype=np.intp)
        )

        target = [-0.1, 0.5, 1.1]
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, "pad"), np.array([-1, 0, 1], dtype=np.intp)
        )
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, "backfill"), np.array([0, 1, 2], dtype=np.intp)
        )
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, "nearest"), np.array([0, 1, 1], dtype=np.intp)
        )

    def test_get_indexer_nan(self):
        # GH#7820
        result = Float64Index([1, 2, np.nan]).get_indexer([np.nan])
        expected = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_get_indexer_int64(self):
        index = Int64Index(range(0, 20, 2))
        target = Int64Index(np.arange(10))
        indexer = index.get_indexer(target)
        expected = np.array([0, -1, 1, -1, 2, -1, 3, -1, 4, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

        target = Int64Index(np.arange(10))
        indexer = index.get_indexer(target, method="pad")
        expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

        target = Int64Index(np.arange(10))
        indexer = index.get_indexer(target, method="backfill")
        expected = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_uint64(self, index_large):
        target = UInt64Index(np.arange(10).astype("uint64") * 5 + 2 ** 63)
        indexer = index_large.get_indexer(target)
        expected = np.array([0, -1, 1, 2, 3, 4, -1, -1, -1, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

        target = UInt64Index(np.arange(10).astype("uint64") * 5 + 2 ** 63)
        indexer = index_large.get_indexer(target, method="pad")
        expected = np.array([0, 0, 1, 2, 3, 4, 4, 4, 4, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

        target = UInt64Index(np.arange(10).astype("uint64") * 5 + 2 ** 63)
        indexer = index_large.get_indexer(target, method="backfill")
        expected = np.array([0, 1, 1, 2, 3, 4, -1, -1, -1, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)


class TestWhere:
    @pytest.mark.parametrize(
        "index",
        [
            Float64Index(np.arange(5, dtype="float64")),
            Int64Index(range(0, 20, 2)),
            UInt64Index(np.arange(5, dtype="uint64")),
        ],
    )
    @pytest.mark.parametrize("klass", [list, tuple, np.array, Series])
    def test_where(self, klass, index):
        cond = [True] * len(index)
        expected = index
        result = index.where(klass(cond))

        cond = [False] + [True] * (len(index) - 1)
        expected = Float64Index([index._na_value] + index[1:].tolist())
        result = index.where(klass(cond))
        tm.assert_index_equal(result, expected)


class TestTake:
    @pytest.mark.parametrize("klass", [Float64Index, Int64Index, UInt64Index])
    def test_take_preserve_name(self, klass):
        index = klass([1, 2, 3, 4], name="foo")
        taken = index.take([3, 0, 1])
        assert index.name == taken.name

    def test_take_fill_value_float64(self):
        # GH 12631
        idx = Float64Index([1.0, 2.0, 3.0], name="xxx")
        result = idx.take(np.array([1, 0, -1]))
        expected = Float64Index([2.0, 1.0, 3.0], name="xxx")
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = Float64Index([2.0, 1.0, np.nan], name="xxx")
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = Float64Index([2.0, 1.0, 3.0], name="xxx")
        tm.assert_index_equal(result, expected)

        msg = (
            "When allow_fill=True and fill_value is not None, "
            "all indices must be >= -1"
        )
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        msg = "index -5 is out of bounds for (axis 0 with )?size 3"
        with pytest.raises(IndexError, match=msg):
            idx.take(np.array([1, -5]))

    @pytest.mark.parametrize("klass", [Int64Index, UInt64Index])
    def test_take_fill_value_ints(self, klass):
        # see gh-12631
        idx = klass([1, 2, 3], name="xxx")
        result = idx.take(np.array([1, 0, -1]))
        expected = klass([2, 1, 3], name="xxx")
        tm.assert_index_equal(result, expected)

        name = klass.__name__
        msg = f"Unable to fill values because {name} cannot contain NA"

        # fill_value=True
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -1]), fill_value=True)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = klass([2, 1, 3], name="xxx")
        tm.assert_index_equal(result, expected)

        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        msg = "index -5 is out of bounds for (axis 0 with )?size 3"
        with pytest.raises(IndexError, match=msg):
            idx.take(np.array([1, -5]))


class TestContains:
    @pytest.mark.parametrize("klass", [Float64Index, Int64Index, UInt64Index])
    def test_contains_none(self, klass):
        # GH#35788 should return False, not raise TypeError
        index = klass([0, 1, 2, 3, 4])
        assert None not in index

    def test_contains_float64_nans(self):
        index = Float64Index([1.0, 2.0, np.nan])
        assert np.nan in index

    def test_contains_float64_not_nans(self):
        index = Float64Index([1.0, 2.0, np.nan])
        assert 1.0 in index
