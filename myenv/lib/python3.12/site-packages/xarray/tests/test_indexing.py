from __future__ import annotations

import itertools
from typing import Any

import numpy as np
import pandas as pd
import pytest

from xarray import DataArray, Dataset, Variable
from xarray.core import indexing, nputils
from xarray.core.indexes import PandasIndex, PandasMultiIndex
from xarray.core.types import T_Xarray
from xarray.tests import (
    IndexerMaker,
    ReturnItem,
    assert_array_equal,
    assert_identical,
    raise_if_dask_computes,
    requires_dask,
)

B = IndexerMaker(indexing.BasicIndexer)


class TestIndexCallable:
    def test_getitem(self):
        def getter(key):
            return key * 2

        indexer = indexing.IndexCallable(getter)
        assert indexer[3] == 6
        assert indexer[0] == 0
        assert indexer[-1] == -2

    def test_setitem(self):
        def getter(key):
            return key * 2

        def setter(key, value):
            raise NotImplementedError("Setter not implemented")

        indexer = indexing.IndexCallable(getter, setter)
        with pytest.raises(NotImplementedError):
            indexer[3] = 6


class TestIndexers:
    def set_to_zero(self, x, i):
        x = x.copy()
        x[i] = 0
        return x

    def test_expanded_indexer(self) -> None:
        x = np.random.randn(10, 11, 12, 13, 14)
        y = np.arange(5)
        arr = ReturnItem()
        for i in [
            arr[:],
            arr[...],
            arr[0, :, 10],
            arr[..., 10],
            arr[:5, ..., 0],
            arr[..., 0, :],
            arr[y],
            arr[y, y],
            arr[..., y, y],
            arr[..., 0, 1, 2, 3, 4],
        ]:
            j = indexing.expanded_indexer(i, x.ndim)
            assert_array_equal(x[i], x[j])
            assert_array_equal(self.set_to_zero(x, i), self.set_to_zero(x, j))
        with pytest.raises(IndexError, match=r"too many indices"):
            indexing.expanded_indexer(arr[1, 2, 3], 2)

    def test_stacked_multiindex_min_max(self) -> None:
        data = np.random.randn(3, 23, 4)
        da = DataArray(
            data,
            name="value",
            dims=["replicate", "rsample", "exp"],
            coords=dict(
                replicate=[0, 1, 2], exp=["a", "b", "c", "d"], rsample=list(range(23))
            ),
        )
        da2 = da.stack(sample=("replicate", "rsample"))
        s = da2.sample
        assert_array_equal(da2.loc["a", s.max()], data[2, 22, 0])
        assert_array_equal(da2.loc["b", s.min()], data[0, 0, 1])

    def test_group_indexers_by_index(self) -> None:
        mindex = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two"))
        data = DataArray(
            np.zeros((4, 2, 2)), coords={"x": mindex, "y": [1, 2]}, dims=("x", "y", "z")
        )
        data.coords["y2"] = ("y", [2.0, 3.0])

        grouped_indexers = indexing.group_indexers_by_index(
            data, {"z": 0, "one": "a", "two": 1, "y": 0}, {}
        )

        for idx, indexers in grouped_indexers:
            if idx is None:
                assert indexers == {"z": 0}
            elif idx.equals(data.xindexes["x"]):
                assert indexers == {"one": "a", "two": 1}
            elif idx.equals(data.xindexes["y"]):
                assert indexers == {"y": 0}
        assert len(grouped_indexers) == 3

        with pytest.raises(KeyError, match=r"no index found for coordinate 'y2'"):
            indexing.group_indexers_by_index(data, {"y2": 2.0}, {})
        with pytest.raises(
            KeyError, match=r"'w' is not a valid dimension or coordinate"
        ):
            indexing.group_indexers_by_index(data, {"w": "a"}, {})
        with pytest.raises(ValueError, match=r"cannot supply.*"):
            indexing.group_indexers_by_index(data, {"z": 1}, {"method": "nearest"})

    def test_map_index_queries(self) -> None:
        def create_sel_results(
            x_indexer,
            x_index,
            other_vars,
            drop_coords,
            drop_indexes,
            rename_dims,
        ):
            dim_indexers = {"x": x_indexer}
            index_vars = x_index.create_variables()
            indexes = {k: x_index for k in index_vars}
            variables = {}
            variables.update(index_vars)
            variables.update(other_vars)

            return indexing.IndexSelResult(
                dim_indexers=dim_indexers,
                indexes=indexes,
                variables=variables,
                drop_coords=drop_coords,
                drop_indexes=drop_indexes,
                rename_dims=rename_dims,
            )

        def test_indexer(
            data: T_Xarray,
            x: Any,
            expected: indexing.IndexSelResult,
        ) -> None:
            results = indexing.map_index_queries(data, {"x": x})

            assert results.dim_indexers.keys() == expected.dim_indexers.keys()
            assert_array_equal(results.dim_indexers["x"], expected.dim_indexers["x"])

            assert results.indexes.keys() == expected.indexes.keys()
            for k in results.indexes:
                assert results.indexes[k].equals(expected.indexes[k])

            assert results.variables.keys() == expected.variables.keys()
            for k in results.variables:
                assert_array_equal(results.variables[k], expected.variables[k])

            assert set(results.drop_coords) == set(expected.drop_coords)
            assert set(results.drop_indexes) == set(expected.drop_indexes)
            assert results.rename_dims == expected.rename_dims

        data = Dataset({"x": ("x", [1, 2, 3])})
        mindex = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2], [-1, -2]], names=("one", "two", "three")
        )
        mdata = DataArray(range(8), [("x", mindex)])

        test_indexer(data, 1, indexing.IndexSelResult({"x": 0}))
        test_indexer(data, np.int32(1), indexing.IndexSelResult({"x": 0}))
        test_indexer(data, Variable([], 1), indexing.IndexSelResult({"x": 0}))
        test_indexer(mdata, ("a", 1, -1), indexing.IndexSelResult({"x": 0}))

        expected = create_sel_results(
            [True, True, False, False, False, False, False, False],
            PandasIndex(pd.Index([-1, -2]), "three"),
            {"one": Variable((), "a"), "two": Variable((), 1)},
            ["x"],
            ["one", "two"],
            {"x": "three"},
        )
        test_indexer(mdata, ("a", 1), expected)

        expected = create_sel_results(
            slice(0, 4, None),
            PandasMultiIndex(
                pd.MultiIndex.from_product([[1, 2], [-1, -2]], names=("two", "three")),
                "x",
            ),
            {"one": Variable((), "a")},
            [],
            ["one"],
            {},
        )
        test_indexer(mdata, "a", expected)

        expected = create_sel_results(
            [True, True, True, True, False, False, False, False],
            PandasMultiIndex(
                pd.MultiIndex.from_product([[1, 2], [-1, -2]], names=("two", "three")),
                "x",
            ),
            {"one": Variable((), "a")},
            [],
            ["one"],
            {},
        )
        test_indexer(mdata, ("a",), expected)

        test_indexer(
            mdata, [("a", 1, -1), ("b", 2, -2)], indexing.IndexSelResult({"x": [0, 7]})
        )
        test_indexer(
            mdata, slice("a", "b"), indexing.IndexSelResult({"x": slice(0, 8, None)})
        )
        test_indexer(
            mdata,
            slice(("a", 1), ("b", 1)),
            indexing.IndexSelResult({"x": slice(0, 6, None)}),
        )
        test_indexer(
            mdata,
            {"one": "a", "two": 1, "three": -1},
            indexing.IndexSelResult({"x": 0}),
        )

        expected = create_sel_results(
            [True, True, False, False, False, False, False, False],
            PandasIndex(pd.Index([-1, -2]), "three"),
            {"one": Variable((), "a"), "two": Variable((), 1)},
            ["x"],
            ["one", "two"],
            {"x": "three"},
        )
        test_indexer(mdata, {"one": "a", "two": 1}, expected)

        expected = create_sel_results(
            [True, False, True, False, False, False, False, False],
            PandasIndex(pd.Index([1, 2]), "two"),
            {"one": Variable((), "a"), "three": Variable((), -1)},
            ["x"],
            ["one", "three"],
            {"x": "two"},
        )
        test_indexer(mdata, {"one": "a", "three": -1}, expected)

        expected = create_sel_results(
            [True, True, True, True, False, False, False, False],
            PandasMultiIndex(
                pd.MultiIndex.from_product([[1, 2], [-1, -2]], names=("two", "three")),
                "x",
            ),
            {"one": Variable((), "a")},
            [],
            ["one"],
            {},
        )
        test_indexer(mdata, {"one": "a"}, expected)

    def test_read_only_view(self) -> None:
        arr = DataArray(
            np.random.rand(3, 3),
            coords={"x": np.arange(3), "y": np.arange(3)},
            dims=("x", "y"),
        )  # Create a 2D DataArray
        arr = arr.expand_dims({"z": 3}, -1)  # New dimension 'z'
        arr["z"] = np.arange(3)  # New coords to dimension 'z'
        with pytest.raises(ValueError, match="Do you want to .copy()"):
            arr.loc[0, 0, 0] = 999


class TestLazyArray:
    def test_slice_slice(self) -> None:
        arr = ReturnItem()
        for size in [100, 99]:
            # We test even/odd size cases
            x = np.arange(size)
            slices = [
                arr[:3],
                arr[:4],
                arr[2:4],
                arr[:1],
                arr[:-1],
                arr[5:-1],
                arr[-5:-1],
                arr[::-1],
                arr[5::-1],
                arr[:3:-1],
                arr[:30:-1],
                arr[10:4:],
                arr[::4],
                arr[4:4:4],
                arr[:4:-4],
                arr[::-2],
            ]
            for i in slices:
                for j in slices:
                    expected = x[i][j]
                    new_slice = indexing.slice_slice(i, j, size=size)
                    actual = x[new_slice]
                    assert_array_equal(expected, actual)

    def test_lazily_indexed_array(self) -> None:
        original = np.random.rand(10, 20, 30)
        x = indexing.NumpyIndexingAdapter(original)
        v = Variable(["i", "j", "k"], original)
        lazy = indexing.LazilyIndexedArray(x)
        v_lazy = Variable(["i", "j", "k"], lazy)
        arr = ReturnItem()
        # test orthogonally applied indexers
        indexers = [arr[:], 0, -2, arr[:3], [0, 1, 2, 3], [0], np.arange(10) < 5]
        for i in indexers:
            for j in indexers:
                for k in indexers:
                    if isinstance(j, np.ndarray) and j.dtype.kind == "b":
                        j = np.arange(20) < 5
                    if isinstance(k, np.ndarray) and k.dtype.kind == "b":
                        k = np.arange(30) < 5
                    expected = np.asarray(v[i, j, k])
                    for actual in [
                        v_lazy[i, j, k],
                        v_lazy[:, j, k][i],
                        v_lazy[:, :, k][:, j][i],
                    ]:
                        assert expected.shape == actual.shape
                        assert_array_equal(expected, actual)
                        assert isinstance(actual._data, indexing.LazilyIndexedArray)
                        assert isinstance(v_lazy._data, indexing.LazilyIndexedArray)

                        # make sure actual.key is appropriate type
                        if all(
                            isinstance(k, (int, slice)) for k in v_lazy._data.key.tuple
                        ):
                            assert isinstance(v_lazy._data.key, indexing.BasicIndexer)
                        else:
                            assert isinstance(v_lazy._data.key, indexing.OuterIndexer)

        # test sequentially applied indexers
        indexers = [
            (3, 2),
            (arr[:], 0),
            (arr[:2], -1),
            (arr[:4], [0]),
            ([4, 5], 0),
            ([0, 1, 2], [0, 1]),
            ([0, 3, 5], arr[:2]),
        ]
        for i, j in indexers:

            expected_b = v[i][j]
            actual = v_lazy[i][j]
            assert expected_b.shape == actual.shape
            assert_array_equal(expected_b, actual)

            # test transpose
            if actual.ndim > 1:
                order = np.random.choice(actual.ndim, actual.ndim)
                order = np.array(actual.dims)
                transposed = actual.transpose(*order)
                assert_array_equal(expected_b.transpose(*order), transposed)
                assert isinstance(
                    actual._data,
                    (
                        indexing.LazilyVectorizedIndexedArray,
                        indexing.LazilyIndexedArray,
                    ),
                )

            assert isinstance(actual._data, indexing.LazilyIndexedArray)
            assert isinstance(actual._data.array, indexing.NumpyIndexingAdapter)

    def test_vectorized_lazily_indexed_array(self) -> None:
        original = np.random.rand(10, 20, 30)
        x = indexing.NumpyIndexingAdapter(original)
        v_eager = Variable(["i", "j", "k"], x)
        lazy = indexing.LazilyIndexedArray(x)
        v_lazy = Variable(["i", "j", "k"], lazy)
        arr = ReturnItem()

        def check_indexing(v_eager, v_lazy, indexers):
            for indexer in indexers:
                actual = v_lazy[indexer]
                expected = v_eager[indexer]
                assert expected.shape == actual.shape
                assert isinstance(
                    actual._data,
                    (
                        indexing.LazilyVectorizedIndexedArray,
                        indexing.LazilyIndexedArray,
                    ),
                )
                assert_array_equal(expected, actual)
                v_eager = expected
                v_lazy = actual

        # test orthogonal indexing
        indexers = [(arr[:], 0, 1), (Variable("i", [0, 1]),)]
        check_indexing(v_eager, v_lazy, indexers)

        # vectorized indexing
        indexers = [
            (Variable("i", [0, 1]), Variable("i", [0, 1]), slice(None)),
            (slice(1, 3, 2), 0),
        ]
        check_indexing(v_eager, v_lazy, indexers)

        indexers = [
            (slice(None, None, 2), 0, slice(None, 10)),
            (Variable("i", [3, 2, 4, 3]), Variable("i", [3, 2, 1, 0])),
            (Variable(["i", "j"], [[0, 1], [1, 2]]),),
        ]
        check_indexing(v_eager, v_lazy, indexers)

        indexers = [
            (Variable("i", [3, 2, 4, 3]), Variable("i", [3, 2, 1, 0])),
            (Variable(["i", "j"], [[0, 1], [1, 2]]),),
        ]
        check_indexing(v_eager, v_lazy, indexers)

    def test_lazily_indexed_array_vindex_setitem(self) -> None:

        lazy = indexing.LazilyIndexedArray(np.random.rand(10, 20, 30))

        # vectorized indexing
        indexer = indexing.VectorizedIndexer(
            (np.array([0, 1]), np.array([0, 1]), slice(None, None, None))
        )
        with pytest.raises(
            NotImplementedError,
            match=r"Lazy item assignment with the vectorized indexer is not yet",
        ):
            lazy.vindex[indexer] = 0

    @pytest.mark.parametrize(
        "indexer_class, key, value",
        [
            (indexing.OuterIndexer, (0, 1, slice(None, None, None)), 10),
            (indexing.BasicIndexer, (0, 1, slice(None, None, None)), 10),
        ],
    )
    def test_lazily_indexed_array_setitem(self, indexer_class, key, value) -> None:
        original = np.random.rand(10, 20, 30)
        x = indexing.NumpyIndexingAdapter(original)
        lazy = indexing.LazilyIndexedArray(x)

        if indexer_class is indexing.BasicIndexer:
            indexer = indexer_class(key)
            lazy[indexer] = value
        elif indexer_class is indexing.OuterIndexer:
            indexer = indexer_class(key)
            lazy.oindex[indexer] = value

        assert_array_equal(original[key], value)


class TestCopyOnWriteArray:
    def test_setitem(self) -> None:
        original = np.arange(10)
        wrapped = indexing.CopyOnWriteArray(original)
        wrapped[B[:]] = 0
        assert_array_equal(original, np.arange(10))
        assert_array_equal(wrapped, np.zeros(10))

    def test_sub_array(self) -> None:
        original = np.arange(10)
        wrapped = indexing.CopyOnWriteArray(original)
        child = wrapped[B[:5]]
        assert isinstance(child, indexing.CopyOnWriteArray)
        child[B[:]] = 0
        assert_array_equal(original, np.arange(10))
        assert_array_equal(wrapped, np.arange(10))
        assert_array_equal(child, np.zeros(5))

    def test_index_scalar(self) -> None:
        # regression test for GH1374
        x = indexing.CopyOnWriteArray(np.array(["foo", "bar"]))
        assert np.array(x[B[0]][B[()]]) == "foo"


class TestMemoryCachedArray:
    def test_wrapper(self) -> None:
        original = indexing.LazilyIndexedArray(np.arange(10))
        wrapped = indexing.MemoryCachedArray(original)
        assert_array_equal(wrapped, np.arange(10))
        assert isinstance(wrapped.array, indexing.NumpyIndexingAdapter)

    def test_sub_array(self) -> None:
        original = indexing.LazilyIndexedArray(np.arange(10))
        wrapped = indexing.MemoryCachedArray(original)
        child = wrapped[B[:5]]
        assert isinstance(child, indexing.MemoryCachedArray)
        assert_array_equal(child, np.arange(5))
        assert isinstance(child.array, indexing.NumpyIndexingAdapter)
        assert isinstance(wrapped.array, indexing.LazilyIndexedArray)

    def test_setitem(self) -> None:
        original = np.arange(10)
        wrapped = indexing.MemoryCachedArray(original)
        wrapped[B[:]] = 0
        assert_array_equal(original, np.zeros(10))

    def test_index_scalar(self) -> None:
        # regression test for GH1374
        x = indexing.MemoryCachedArray(np.array(["foo", "bar"]))
        assert np.array(x[B[0]][B[()]]) == "foo"


def test_base_explicit_indexer() -> None:
    with pytest.raises(TypeError):
        indexing.ExplicitIndexer(())

    class Subclass(indexing.ExplicitIndexer):
        pass

    value = Subclass((1, 2, 3))
    assert value.tuple == (1, 2, 3)
    assert repr(value) == "Subclass((1, 2, 3))"


@pytest.mark.parametrize(
    "indexer_cls",
    [indexing.BasicIndexer, indexing.OuterIndexer, indexing.VectorizedIndexer],
)
def test_invalid_for_all(indexer_cls) -> None:
    with pytest.raises(TypeError):
        indexer_cls(None)
    with pytest.raises(TypeError):
        indexer_cls(([],))
    with pytest.raises(TypeError):
        indexer_cls((None,))
    with pytest.raises(TypeError):
        indexer_cls(("foo",))
    with pytest.raises(TypeError):
        indexer_cls((1.0,))
    with pytest.raises(TypeError):
        indexer_cls((slice("foo"),))
    with pytest.raises(TypeError):
        indexer_cls((np.array(["foo"]),))


def check_integer(indexer_cls):
    value = indexer_cls((1, np.uint64(2))).tuple
    assert all(isinstance(v, int) for v in value)
    assert value == (1, 2)


def check_slice(indexer_cls):
    (value,) = indexer_cls((slice(1, None, np.int64(2)),)).tuple
    assert value == slice(1, None, 2)
    assert isinstance(value.step, int)


def check_array1d(indexer_cls):
    (value,) = indexer_cls((np.arange(3, dtype=np.int32),)).tuple
    assert value.dtype == np.int64
    np.testing.assert_array_equal(value, [0, 1, 2])


def check_array2d(indexer_cls):
    array = np.array([[1, 2], [3, 4]], dtype=np.int64)
    (value,) = indexer_cls((array,)).tuple
    assert value.dtype == np.int64
    np.testing.assert_array_equal(value, array)


def test_basic_indexer() -> None:
    check_integer(indexing.BasicIndexer)
    check_slice(indexing.BasicIndexer)
    with pytest.raises(TypeError):
        check_array1d(indexing.BasicIndexer)
    with pytest.raises(TypeError):
        check_array2d(indexing.BasicIndexer)


def test_outer_indexer() -> None:
    check_integer(indexing.OuterIndexer)
    check_slice(indexing.OuterIndexer)
    check_array1d(indexing.OuterIndexer)
    with pytest.raises(TypeError):
        check_array2d(indexing.OuterIndexer)


def test_vectorized_indexer() -> None:
    with pytest.raises(TypeError):
        check_integer(indexing.VectorizedIndexer)
    check_slice(indexing.VectorizedIndexer)
    check_array1d(indexing.VectorizedIndexer)
    check_array2d(indexing.VectorizedIndexer)
    with pytest.raises(ValueError, match=r"numbers of dimensions"):
        indexing.VectorizedIndexer(
            (np.array(1, dtype=np.int64), np.arange(5, dtype=np.int64))
        )


class Test_vectorized_indexer:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.data = indexing.NumpyIndexingAdapter(np.random.randn(10, 12, 13))
        self.indexers = [
            np.array([[0, 3, 2]]),
            np.array([[0, 3, 3], [4, 6, 7]]),
            slice(2, -2, 2),
            slice(2, -2, 3),
            slice(None),
        ]

    def test_arrayize_vectorized_indexer(self) -> None:
        for i, j, k in itertools.product(self.indexers, repeat=3):
            vindex = indexing.VectorizedIndexer((i, j, k))
            vindex_array = indexing._arrayize_vectorized_indexer(
                vindex, self.data.shape
            )
            np.testing.assert_array_equal(
                self.data.vindex[vindex], self.data.vindex[vindex_array]
            )

        actual = indexing._arrayize_vectorized_indexer(
            indexing.VectorizedIndexer((slice(None),)), shape=(5,)
        )
        np.testing.assert_array_equal(actual.tuple, [np.arange(5)])

        actual = indexing._arrayize_vectorized_indexer(
            indexing.VectorizedIndexer((np.arange(5),) * 3), shape=(8, 10, 12)
        )
        expected = np.stack([np.arange(5)] * 3)
        np.testing.assert_array_equal(np.stack(actual.tuple), expected)

        actual = indexing._arrayize_vectorized_indexer(
            indexing.VectorizedIndexer((np.arange(5), slice(None))), shape=(8, 10)
        )
        a, b = actual.tuple
        np.testing.assert_array_equal(a, np.arange(5)[:, np.newaxis])
        np.testing.assert_array_equal(b, np.arange(10)[np.newaxis, :])

        actual = indexing._arrayize_vectorized_indexer(
            indexing.VectorizedIndexer((slice(None), np.arange(5))), shape=(8, 10)
        )
        a, b = actual.tuple
        np.testing.assert_array_equal(a, np.arange(8)[np.newaxis, :])
        np.testing.assert_array_equal(b, np.arange(5)[:, np.newaxis])


def get_indexers(shape, mode):
    if mode == "vectorized":
        indexed_shape = (3, 4)
        indexer = tuple(np.random.randint(0, s, size=indexed_shape) for s in shape)
        return indexing.VectorizedIndexer(indexer)

    elif mode == "outer":
        indexer = tuple(np.random.randint(0, s, s + 2) for s in shape)
        return indexing.OuterIndexer(indexer)

    elif mode == "outer_scalar":
        indexer = (np.random.randint(0, 3, 4), 0, slice(None, None, 2))
        return indexing.OuterIndexer(indexer[: len(shape)])

    elif mode == "outer_scalar2":
        indexer = (np.random.randint(0, 3, 4), -2, slice(None, None, 2))
        return indexing.OuterIndexer(indexer[: len(shape)])

    elif mode == "outer1vec":
        indexer = [slice(2, -3) for s in shape]
        indexer[1] = np.random.randint(0, shape[1], shape[1] + 2)
        return indexing.OuterIndexer(tuple(indexer))

    elif mode == "basic":  # basic indexer
        indexer = [slice(2, -3) for s in shape]
        indexer[0] = 3
        return indexing.BasicIndexer(tuple(indexer))

    elif mode == "basic1":  # basic indexer
        return indexing.BasicIndexer((3,))

    elif mode == "basic2":  # basic indexer
        indexer = [0, 2, 4]
        return indexing.BasicIndexer(tuple(indexer[: len(shape)]))

    elif mode == "basic3":  # basic indexer
        indexer = [slice(None) for s in shape]
        indexer[0] = slice(-2, 2, -2)
        indexer[1] = slice(1, -1, 2)
        return indexing.BasicIndexer(tuple(indexer[: len(shape)]))


@pytest.mark.parametrize("size", [100, 99])
@pytest.mark.parametrize(
    "sl", [slice(1, -1, 1), slice(None, -1, 2), slice(-1, 1, -1), slice(-1, 1, -2)]
)
def test_decompose_slice(size, sl) -> None:
    x = np.arange(size)
    slice1, slice2 = indexing._decompose_slice(sl, size)
    expected = x[sl]
    actual = x[slice1][slice2]
    assert_array_equal(expected, actual)


@pytest.mark.parametrize("shape", [(10, 5, 8), (10, 3)])
@pytest.mark.parametrize(
    "indexer_mode",
    [
        "vectorized",
        "outer",
        "outer_scalar",
        "outer_scalar2",
        "outer1vec",
        "basic",
        "basic1",
        "basic2",
        "basic3",
    ],
)
@pytest.mark.parametrize(
    "indexing_support",
    [
        indexing.IndexingSupport.BASIC,
        indexing.IndexingSupport.OUTER,
        indexing.IndexingSupport.OUTER_1VECTOR,
        indexing.IndexingSupport.VECTORIZED,
    ],
)
def test_decompose_indexers(shape, indexer_mode, indexing_support) -> None:
    data = np.random.randn(*shape)
    indexer = get_indexers(shape, indexer_mode)

    backend_ind, np_ind = indexing.decompose_indexer(indexer, shape, indexing_support)
    indexing_adapter = indexing.NumpyIndexingAdapter(data)

    # Dispatch to appropriate indexing method
    if indexer_mode.startswith("vectorized"):
        expected = indexing_adapter.vindex[indexer]

    elif indexer_mode.startswith("outer"):
        expected = indexing_adapter.oindex[indexer]

    else:
        expected = indexing_adapter[indexer]  # Basic indexing

    if isinstance(backend_ind, indexing.VectorizedIndexer):
        array = indexing_adapter.vindex[backend_ind]
    elif isinstance(backend_ind, indexing.OuterIndexer):
        array = indexing_adapter.oindex[backend_ind]
    else:
        array = indexing_adapter[backend_ind]

    if len(np_ind.tuple) > 0:
        array_indexing_adapter = indexing.NumpyIndexingAdapter(array)
        if isinstance(np_ind, indexing.VectorizedIndexer):
            array = array_indexing_adapter.vindex[np_ind]
        elif isinstance(np_ind, indexing.OuterIndexer):
            array = array_indexing_adapter.oindex[np_ind]
        else:
            array = array_indexing_adapter[np_ind]
    np.testing.assert_array_equal(expected, array)

    if not all(isinstance(k, indexing.integer_types) for k in np_ind.tuple):
        combined_ind = indexing._combine_indexers(backend_ind, shape, np_ind)
        assert isinstance(combined_ind, indexing.VectorizedIndexer)
        array = indexing_adapter.vindex[combined_ind]
        np.testing.assert_array_equal(expected, array)


def test_implicit_indexing_adapter() -> None:
    array = np.arange(10, dtype=np.int64)
    implicit = indexing.ImplicitToExplicitIndexingAdapter(
        indexing.NumpyIndexingAdapter(array), indexing.BasicIndexer
    )
    np.testing.assert_array_equal(array, np.asarray(implicit))
    np.testing.assert_array_equal(array, implicit[:])


def test_implicit_indexing_adapter_copy_on_write() -> None:
    array = np.arange(10, dtype=np.int64)
    implicit = indexing.ImplicitToExplicitIndexingAdapter(
        indexing.CopyOnWriteArray(array)
    )
    assert isinstance(implicit[:], indexing.ImplicitToExplicitIndexingAdapter)


def test_outer_indexer_consistency_with_broadcast_indexes_vectorized() -> None:
    def nonzero(x):
        if isinstance(x, np.ndarray) and x.dtype.kind == "b":
            x = x.nonzero()[0]
        return x

    original = np.random.rand(10, 20, 30)
    v = Variable(["i", "j", "k"], original)
    arr = ReturnItem()
    # test orthogonally applied indexers
    indexers = [
        arr[:],
        0,
        -2,
        arr[:3],
        np.array([0, 1, 2, 3]),
        np.array([0]),
        np.arange(10) < 5,
    ]
    for i, j, k in itertools.product(indexers, repeat=3):
        if isinstance(j, np.ndarray) and j.dtype.kind == "b":  # match size
            j = np.arange(20) < 4
        if isinstance(k, np.ndarray) and k.dtype.kind == "b":
            k = np.arange(30) < 8

        _, expected, new_order = v._broadcast_indexes_vectorized((i, j, k))
        expected_data = nputils.NumpyVIndexAdapter(v.data)[expected.tuple]
        if new_order:
            old_order = range(len(new_order))
            expected_data = np.moveaxis(expected_data, old_order, new_order)

        outer_index = indexing.OuterIndexer((nonzero(i), nonzero(j), nonzero(k)))
        actual = indexing._outer_to_numpy_indexer(outer_index, v.shape)
        actual_data = v.data[actual]
        np.testing.assert_array_equal(actual_data, expected_data)


def test_create_mask_outer_indexer() -> None:
    indexer = indexing.OuterIndexer((np.array([0, -1, 2]),))
    expected = np.array([False, True, False])
    actual = indexing.create_mask(indexer, (5,))
    np.testing.assert_array_equal(expected, actual)

    indexer = indexing.OuterIndexer((1, slice(2), np.array([0, -1, 2])))
    expected = np.array(2 * [[False, True, False]])
    actual = indexing.create_mask(indexer, (5, 5, 5))
    np.testing.assert_array_equal(expected, actual)


def test_create_mask_vectorized_indexer() -> None:
    indexer = indexing.VectorizedIndexer((np.array([0, -1, 2]), np.array([0, 1, -1])))
    expected = np.array([False, True, True])
    actual = indexing.create_mask(indexer, (5,))
    np.testing.assert_array_equal(expected, actual)

    indexer = indexing.VectorizedIndexer(
        (np.array([0, -1, 2]), slice(None), np.array([0, 1, -1]))
    )
    expected = np.array([[False, True, True]] * 2).T
    actual = indexing.create_mask(indexer, (5, 2))
    np.testing.assert_array_equal(expected, actual)


def test_create_mask_basic_indexer() -> None:
    indexer = indexing.BasicIndexer((-1,))
    actual = indexing.create_mask(indexer, (3,))
    np.testing.assert_array_equal(True, actual)

    indexer = indexing.BasicIndexer((0,))
    actual = indexing.create_mask(indexer, (3,))
    np.testing.assert_array_equal(False, actual)


def test_create_mask_dask() -> None:
    da = pytest.importorskip("dask.array")

    indexer = indexing.OuterIndexer((1, slice(2), np.array([0, -1, 2])))
    expected = np.array(2 * [[False, True, False]])
    actual = indexing.create_mask(
        indexer, (5, 5, 5), da.empty((2, 3), chunks=((1, 1), (2, 1)))
    )
    assert actual.chunks == ((1, 1), (2, 1))
    np.testing.assert_array_equal(expected, actual)

    indexer_vec = indexing.VectorizedIndexer(
        (np.array([0, -1, 2]), slice(None), np.array([0, 1, -1]))
    )
    expected = np.array([[False, True, True]] * 2).T
    actual = indexing.create_mask(
        indexer_vec, (5, 2), da.empty((3, 2), chunks=((3,), (2,)))
    )
    assert isinstance(actual, da.Array)
    np.testing.assert_array_equal(expected, actual)

    with pytest.raises(ValueError):
        indexing.create_mask(indexer_vec, (5, 2), da.empty((5,), chunks=(1,)))


def test_create_mask_error() -> None:
    with pytest.raises(TypeError, match=r"unexpected key type"):
        indexing.create_mask((1, 2), (3, 4))  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "indices, expected",
    [
        (np.arange(5), np.arange(5)),
        (np.array([0, -1, -1]), np.array([0, 0, 0])),
        (np.array([-1, 1, -1]), np.array([1, 1, 1])),
        (np.array([-1, -1, 2]), np.array([2, 2, 2])),
        (np.array([-1]), np.array([0])),
        (np.array([0, -1, 1, -1, -1]), np.array([0, 0, 1, 1, 1])),
        (np.array([0, -1, -1, -1, 1]), np.array([0, 0, 0, 0, 1])),
    ],
)
def test_posify_mask_subindexer(indices, expected) -> None:
    actual = indexing._posify_mask_subindexer(indices)
    np.testing.assert_array_equal(expected, actual)


def test_indexing_1d_object_array() -> None:
    items = (np.arange(3), np.arange(6))
    arr = DataArray(np.array(items, dtype=object))

    actual = arr[0]

    expected_data = np.empty((), dtype=object)
    expected_data[()] = items[0]
    expected = DataArray(expected_data)

    assert [actual.data.item()] == [expected.data.item()]


@requires_dask
def test_indexing_dask_array():
    import dask.array

    da = DataArray(
        np.ones(10 * 3 * 3).reshape((10, 3, 3)),
        dims=("time", "x", "y"),
    ).chunk(dict(time=-1, x=1, y=1))
    with raise_if_dask_computes():
        actual = da.isel(time=dask.array.from_array([9], chunks=(1,)))
    expected = da.isel(time=[9])
    assert_identical(actual, expected)


@requires_dask
def test_indexing_dask_array_scalar():
    # GH4276
    import dask.array

    a = dask.array.from_array(np.linspace(0.0, 1.0))
    da = DataArray(a, dims="x")
    x_selector = da.argmax(dim=...)
    with raise_if_dask_computes():
        actual = da.isel(x_selector)
    expected = da.isel(x=-1)
    assert_identical(actual, expected)


@requires_dask
def test_vectorized_indexing_dask_array():
    # https://github.com/pydata/xarray/issues/2511#issuecomment-563330352
    darr = DataArray(data=[0.2, 0.4, 0.6], coords={"z": range(3)}, dims=("z",))
    indexer = DataArray(
        data=np.random.randint(0, 3, 8).reshape(4, 2).astype(int),
        coords={"y": range(4), "x": range(2)},
        dims=("y", "x"),
    )
    with pytest.raises(ValueError, match="Vectorized indexing with Dask arrays"):
        darr[indexer.chunk({"y": 2})]


@requires_dask
def test_advanced_indexing_dask_array():
    # GH4663
    import dask.array as da

    ds = Dataset(
        dict(
            a=("x", da.from_array(np.random.randint(0, 100, 100))),
            b=(("x", "y"), da.random.random((100, 10))),
        )
    )
    expected = ds.b.sel(x=ds.a.compute())
    with raise_if_dask_computes():
        actual = ds.b.sel(x=ds.a)
    assert_identical(expected, actual)

    with raise_if_dask_computes():
        actual = ds.b.sel(x=ds.a.data)
    assert_identical(expected, actual)
