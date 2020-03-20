import numpy as np
import pytest

from pandas.core.dtypes.common import ensure_platform_int

import pandas as pd
from pandas import Float64Index, Index, Int64Index, RangeIndex
import pandas._testing as tm

from ..test_numeric import Numeric

# aliases to make some tests easier to read
RI = RangeIndex
I64 = Int64Index
F64 = Float64Index
OI = Index


class TestRangeIndex(Numeric):
    _holder = RangeIndex
    _compat_props = ["shape", "ndim", "size"]

    @pytest.fixture(
        params=[
            RangeIndex(start=0, stop=20, step=2, name="foo"),
            RangeIndex(start=18, stop=-1, step=-2, name="bar"),
        ],
        ids=["index_inc", "index_dec"],
    )
    def indices(self, request):
        return request.param

    def create_index(self) -> RangeIndex:
        return RangeIndex(start=0, stop=20, step=2)

    def test_can_hold_identifiers(self):
        idx = self.create_index()
        key = idx[0]
        assert idx._can_hold_identifiers_and_holds_name(key) is False

    def test_too_many_names(self):
        index = self.create_index()
        with pytest.raises(ValueError, match="^Length"):
            index.names = ["roger", "harold"]

    @pytest.mark.parametrize(
        "index, start, stop, step",
        [
            (RangeIndex(5), 0, 5, 1),
            (RangeIndex(0, 5), 0, 5, 1),
            (RangeIndex(5, step=2), 0, 5, 2),
            (RangeIndex(1, 5, 2), 1, 5, 2),
        ],
    )
    def test_start_stop_step_attrs(self, index, start, stop, step):
        # GH 25710
        assert index.start == start
        assert index.stop == stop
        assert index.step == step

    @pytest.mark.parametrize("attr_name", ["_start", "_stop", "_step"])
    def test_deprecated_start_stop_step_attrs(self, attr_name):
        # GH 26581
        idx = self.create_index()
        with tm.assert_produces_warning(FutureWarning):
            getattr(idx, attr_name)

    def test_copy(self):
        i = RangeIndex(5, name="Foo")
        i_copy = i.copy()
        assert i_copy is not i
        assert i_copy.identical(i)
        assert i_copy._range == range(0, 5, 1)
        assert i_copy.name == "Foo"

    def test_repr(self):
        i = RangeIndex(5, name="Foo")
        result = repr(i)
        expected = "RangeIndex(start=0, stop=5, step=1, name='Foo')"
        assert result == expected

        result = eval(result)
        tm.assert_index_equal(result, i, exact=True)

        i = RangeIndex(5, 0, -1)
        result = repr(i)
        expected = "RangeIndex(start=5, stop=0, step=-1)"
        assert result == expected

        result = eval(result)
        tm.assert_index_equal(result, i, exact=True)

    def test_insert(self):

        idx = RangeIndex(5, name="Foo")
        result = idx[1:4]

        # test 0th element
        tm.assert_index_equal(idx[0:4], result.insert(0, idx[0]))

        # GH 18295 (test missing)
        expected = Float64Index([0, np.nan, 1, 2, 3, 4])
        for na in (np.nan, pd.NaT, None):
            result = RangeIndex(5).insert(1, na)
            tm.assert_index_equal(result, expected)

    def test_delete(self):

        idx = RangeIndex(5, name="Foo")
        expected = idx[1:].astype(int)
        result = idx.delete(0)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

        expected = idx[:-1].astype(int)
        result = idx.delete(-1)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

        with pytest.raises((IndexError, ValueError)):
            # either depending on numpy version
            result = idx.delete(len(idx))

    def test_view(self):
        i = RangeIndex(0, name="Foo")
        i_view = i.view()
        assert i_view.name == "Foo"

        i_view = i.view("i8")
        tm.assert_numpy_array_equal(i.values, i_view)

        i_view = i.view(RangeIndex)
        tm.assert_index_equal(i, i_view)

    def test_dtype(self):
        index = self.create_index()
        assert index.dtype == np.int64

    def test_cached_data(self):
        # GH 26565, GH26617
        # Calling RangeIndex._data caches an int64 array of the same length at
        # self._cached_data. This test checks whether _cached_data has been set
        idx = RangeIndex(0, 100, 10)

        assert idx._cached_data is None

        repr(idx)
        assert idx._cached_data is None

        str(idx)
        assert idx._cached_data is None

        idx.get_loc(20)
        assert idx._cached_data is None

        90 in idx
        assert idx._cached_data is None

        91 in idx
        assert idx._cached_data is None

        idx.all()
        assert idx._cached_data is None

        idx.any()
        assert idx._cached_data is None

        df = pd.DataFrame({"a": range(10)}, index=idx)

        df.loc[50]
        assert idx._cached_data is None

        with pytest.raises(KeyError, match="51"):
            df.loc[51]
        assert idx._cached_data is None

        df.loc[10:50]
        assert idx._cached_data is None

        df.iloc[5:10]
        assert idx._cached_data is None

        # actually calling idx._data
        assert isinstance(idx._data, np.ndarray)
        assert isinstance(idx._cached_data, np.ndarray)

    def test_is_monotonic(self):
        index = RangeIndex(0, 20, 2)
        assert index.is_monotonic is True
        assert index.is_monotonic_increasing is True
        assert index.is_monotonic_decreasing is False
        assert index._is_strictly_monotonic_increasing is True
        assert index._is_strictly_monotonic_decreasing is False

        index = RangeIndex(4, 0, -1)
        assert index.is_monotonic is False
        assert index._is_strictly_monotonic_increasing is False
        assert index.is_monotonic_decreasing is True
        assert index._is_strictly_monotonic_decreasing is True

        index = RangeIndex(1, 2)
        assert index.is_monotonic is True
        assert index.is_monotonic_increasing is True
        assert index.is_monotonic_decreasing is True
        assert index._is_strictly_monotonic_increasing is True
        assert index._is_strictly_monotonic_decreasing is True

        index = RangeIndex(2, 1)
        assert index.is_monotonic is True
        assert index.is_monotonic_increasing is True
        assert index.is_monotonic_decreasing is True
        assert index._is_strictly_monotonic_increasing is True
        assert index._is_strictly_monotonic_decreasing is True

        index = RangeIndex(1, 1)
        assert index.is_monotonic is True
        assert index.is_monotonic_increasing is True
        assert index.is_monotonic_decreasing is True
        assert index._is_strictly_monotonic_increasing is True
        assert index._is_strictly_monotonic_decreasing is True

    def test_equals_range(self):
        equiv_pairs = [
            (RangeIndex(0, 9, 2), RangeIndex(0, 10, 2)),
            (RangeIndex(0), RangeIndex(1, -1, 3)),
            (RangeIndex(1, 2, 3), RangeIndex(1, 3, 4)),
            (RangeIndex(0, -9, -2), RangeIndex(0, -10, -2)),
        ]
        for left, right in equiv_pairs:
            assert left.equals(right)
            assert right.equals(left)

    def test_logical_compat(self):
        idx = self.create_index()
        assert idx.all() == idx.values.all()
        assert idx.any() == idx.values.any()

    def test_identical(self):
        index = self.create_index()
        i = Index(index.copy())
        assert i.identical(index)

        # we don't allow object dtype for RangeIndex
        if isinstance(index, RangeIndex):
            return

        same_values_different_type = Index(i, dtype=object)
        assert not i.identical(same_values_different_type)

        i = index.copy(dtype=object)
        i = i.rename("foo")
        same_values = Index(i, dtype=object)
        assert same_values.identical(index.copy(dtype=object))

        assert not i.identical(index)
        assert Index(same_values, name="foo", dtype=object).identical(i)

        assert not index.copy(dtype=object).identical(index.copy(dtype="int64"))

    def test_get_indexer(self):
        index = self.create_index()
        target = RangeIndex(10)
        indexer = index.get_indexer(target)
        expected = np.array([0, -1, 1, -1, 2, -1, 3, -1, 4, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_pad(self):
        index = self.create_index()
        target = RangeIndex(10)
        indexer = index.get_indexer(target, method="pad")
        expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_backfill(self):
        index = self.create_index()
        target = RangeIndex(10)
        indexer = index.get_indexer(target, method="backfill")
        expected = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_limit(self):
        # GH 28631
        idx = RangeIndex(4)
        target = RangeIndex(6)
        result = idx.get_indexer(target, method="pad", limit=1)
        expected = np.array([0, 1, 2, 3, 3, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("stop", [0, -1, -2])
    def test_get_indexer_decreasing(self, stop):
        # GH 28678
        index = RangeIndex(7, stop, -3)
        result = index.get_indexer(range(9))
        expected = np.array([-1, 2, -1, -1, 1, -1, -1, 0, -1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_nbytes(self):

        # memory savings vs int index
        i = RangeIndex(0, 1000)
        assert i.nbytes < i._int64index.nbytes / 10

        # constant memory usage
        i2 = RangeIndex(0, 10)
        assert i.nbytes == i2.nbytes

    def test_cant_or_shouldnt_cast(self):
        # can't
        with pytest.raises(TypeError):
            RangeIndex("foo", "bar", "baz")

        # shouldn't
        with pytest.raises(TypeError):
            RangeIndex("0", "1", "2")

    def test_view_index(self):
        index = self.create_index()
        index.view(Index)

    def test_prevent_casting(self):
        index = self.create_index()
        result = index.astype("O")
        assert result.dtype == np.object_

    def test_take_preserve_name(self):
        index = RangeIndex(1, 5, name="foo")
        taken = index.take([3, 0, 1])
        assert index.name == taken.name

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.RangeIndex(1, 4, name="xxx")
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.Int64Index([2, 1, 3], name="xxx")
        tm.assert_index_equal(result, expected)

        # fill_value
        msg = "Unable to fill values because RangeIndex cannot contain NA"
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -1]), fill_value=True)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = pd.Int64Index([2, 1, 3], name="xxx")
        tm.assert_index_equal(result, expected)

        msg = "Unable to fill values because RangeIndex cannot contain NA"
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with pytest.raises(IndexError):
            idx.take(np.array([1, -5]))

    def test_print_unicode_columns(self):
        df = pd.DataFrame({"\u05d0": [1, 2, 3], "\u05d1": [4, 5, 6], "c": [7, 8, 9]})
        repr(df.columns)  # should not raise UnicodeDecodeError

    def test_repr_roundtrip(self):
        index = self.create_index()
        tm.assert_index_equal(eval(repr(index)), index)

    def test_slice_keep_name(self):
        idx = RangeIndex(1, 2, name="asdf")
        assert idx.name == idx[1:].name

    def test_explicit_conversions(self):

        # GH 8608
        # add/sub are overridden explicitly for Float/Int Index
        idx = RangeIndex(5)

        # float conversions
        arr = np.arange(5, dtype="int64") * 3.2
        expected = Float64Index(arr)
        fidx = idx * 3.2
        tm.assert_index_equal(fidx, expected)
        fidx = 3.2 * idx
        tm.assert_index_equal(fidx, expected)

        # interops with numpy arrays
        expected = Float64Index(arr)
        a = np.zeros(5, dtype="float64")
        result = fidx - a
        tm.assert_index_equal(result, expected)

        expected = Float64Index(-arr)
        a = np.zeros(5, dtype="float64")
        result = a - fidx
        tm.assert_index_equal(result, expected)

    def test_has_duplicates(self, indices):
        assert indices.is_unique
        assert not indices.has_duplicates

    def test_extended_gcd(self):
        index = self.create_index()
        result = index._extended_gcd(6, 10)
        assert result[0] == result[1] * 6 + result[2] * 10
        assert 2 == result[0]

        result = index._extended_gcd(10, 6)
        assert 2 == result[1] * 10 + result[2] * 6
        assert 2 == result[0]

    def test_min_fitting_element(self):
        result = RangeIndex(0, 20, 2)._min_fitting_element(1)
        assert 2 == result

        result = RangeIndex(1, 6)._min_fitting_element(1)
        assert 1 == result

        result = RangeIndex(18, -2, -2)._min_fitting_element(1)
        assert 2 == result

        result = RangeIndex(5, 0, -1)._min_fitting_element(1)
        assert 1 == result

        big_num = 500000000000000000000000

        result = RangeIndex(5, big_num * 2, 1)._min_fitting_element(big_num)
        assert big_num == result

    def test_max_fitting_element(self):
        result = RangeIndex(0, 20, 2)._max_fitting_element(17)
        assert 16 == result

        result = RangeIndex(1, 6)._max_fitting_element(4)
        assert 4 == result

        result = RangeIndex(18, -2, -2)._max_fitting_element(17)
        assert 16 == result

        result = RangeIndex(5, 0, -1)._max_fitting_element(4)
        assert 4 == result

        big_num = 500000000000000000000000

        result = RangeIndex(5, big_num * 2, 1)._max_fitting_element(big_num)
        assert big_num == result

    def test_pickle_compat_construction(self):
        # RangeIndex() is a valid constructor
        pass

    def test_slice_specialised(self):
        index = self.create_index()
        index.name = "foo"

        # scalar indexing
        res = index[1]
        expected = 2
        assert res == expected

        res = index[-1]
        expected = 18
        assert res == expected

        # slicing
        # slice value completion
        index_slice = index[:]
        expected = index
        tm.assert_index_equal(index_slice, expected)

        # positive slice values
        index_slice = index[7:10:2]
        expected = Index(np.array([14, 18]), name="foo")
        tm.assert_index_equal(index_slice, expected)

        # negative slice values
        index_slice = index[-1:-5:-2]
        expected = Index(np.array([18, 14]), name="foo")
        tm.assert_index_equal(index_slice, expected)

        # stop overshoot
        index_slice = index[2:100:4]
        expected = Index(np.array([4, 12]), name="foo")
        tm.assert_index_equal(index_slice, expected)

        # reverse
        index_slice = index[::-1]
        expected = Index(index.values[::-1], name="foo")
        tm.assert_index_equal(index_slice, expected)

        index_slice = index[-8::-1]
        expected = Index(np.array([4, 2, 0]), name="foo")
        tm.assert_index_equal(index_slice, expected)

        index_slice = index[-40::-1]
        expected = Index(np.array([], dtype=np.int64), name="foo")
        tm.assert_index_equal(index_slice, expected)

        index_slice = index[40::-1]
        expected = Index(index.values[40::-1], name="foo")
        tm.assert_index_equal(index_slice, expected)

        index_slice = index[10::-1]
        expected = Index(index.values[::-1], name="foo")
        tm.assert_index_equal(index_slice, expected)

    @pytest.mark.parametrize("step", set(range(-5, 6)) - {0})
    def test_len_specialised(self, step):
        # make sure that our len is the same as np.arange calc
        start, stop = (0, 5) if step > 0 else (5, 0)

        arr = np.arange(start, stop, step)
        index = RangeIndex(start, stop, step)
        assert len(index) == len(arr)

        index = RangeIndex(stop, start, step)
        assert len(index) == 0

    @pytest.fixture(
        params=[
            ([RI(1, 12, 5)], RI(1, 12, 5)),
            ([RI(0, 6, 4)], RI(0, 6, 4)),
            ([RI(1, 3), RI(3, 7)], RI(1, 7)),
            ([RI(1, 5, 2), RI(5, 6)], RI(1, 6, 2)),
            ([RI(1, 3, 2), RI(4, 7, 3)], RI(1, 7, 3)),
            ([RI(-4, 3, 2), RI(4, 7, 2)], RI(-4, 7, 2)),
            ([RI(-4, -8), RI(-8, -12)], RI(0, 0)),
            ([RI(-4, -8), RI(3, -4)], RI(0, 0)),
            ([RI(-4, -8), RI(3, 5)], RI(3, 5)),
            ([RI(-4, -2), RI(3, 5)], I64([-4, -3, 3, 4])),
            ([RI(-2), RI(3, 5)], RI(3, 5)),
            ([RI(2), RI(2)], I64([0, 1, 0, 1])),
            ([RI(2), RI(2, 5), RI(5, 8, 4)], RI(0, 6)),
            ([RI(2), RI(3, 5), RI(5, 8, 4)], I64([0, 1, 3, 4, 5])),
            ([RI(-2, 2), RI(2, 5), RI(5, 8, 4)], RI(-2, 6)),
            ([RI(3), I64([-1, 3, 15])], I64([0, 1, 2, -1, 3, 15])),
            ([RI(3), F64([-1, 3.1, 15.0])], F64([0, 1, 2, -1, 3.1, 15.0])),
            ([RI(3), OI(["a", None, 14])], OI([0, 1, 2, "a", None, 14])),
            ([RI(3, 1), OI(["a", None, 14])], OI(["a", None, 14])),
        ]
    )
    def appends(self, request):
        """Inputs and expected outputs for RangeIndex.append test"""
        return request.param

    def test_append(self, appends):
        # GH16212

        indices, expected = appends

        result = indices[0].append(indices[1:])
        tm.assert_index_equal(result, expected, exact=True)

        if len(indices) == 2:
            # Append single item rather than list
            result2 = indices[0].append(indices[1])
            tm.assert_index_equal(result2, expected, exact=True)

    def test_engineless_lookup(self):
        # GH 16685
        # Standard lookup on RangeIndex should not require the engine to be
        # created
        idx = RangeIndex(2, 10, 3)

        assert idx.get_loc(5) == 1
        tm.assert_numpy_array_equal(
            idx.get_indexer([2, 8]), ensure_platform_int(np.array([0, 2]))
        )
        with pytest.raises(KeyError, match="3"):
            idx.get_loc(3)

        assert "_engine" not in idx._cache

        # Different types of scalars can be excluded immediately, no need to
        #  use the _engine
        with pytest.raises(KeyError, match="'a'"):
            idx.get_loc("a")

        assert "_engine" not in idx._cache
