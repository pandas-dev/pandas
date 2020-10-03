import numpy as np

from pandas import Categorical, CategoricalIndex, Index
import pandas._testing as tm


class TestReindex:
    def test_reindex_dtype(self):
        c = CategoricalIndex(["a", "b", "c", "a"])
        res, indexer = c.reindex(["a", "c"])
        tm.assert_index_equal(res, Index(["a", "a", "c"]), exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

        c = CategoricalIndex(["a", "b", "c", "a"])
        res, indexer = c.reindex(Categorical(["a", "c"]))

        exp = CategoricalIndex(["a", "a", "c"], categories=["a", "c"])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

        c = CategoricalIndex(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        res, indexer = c.reindex(["a", "c"])
        exp = Index(["a", "a", "c"], dtype="object")
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

        c = CategoricalIndex(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        res, indexer = c.reindex(Categorical(["a", "c"]))
        exp = CategoricalIndex(["a", "a", "c"], categories=["a", "c"])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

    def test_reindex_duplicate_target(self):
        # See GH25459
        cat = CategoricalIndex(["a", "b", "c"], categories=["a", "b", "c", "d"])
        res, indexer = cat.reindex(["a", "c", "c"])
        exp = Index(["a", "c", "c"], dtype="object")
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 2, 2], dtype=np.intp))

        res, indexer = cat.reindex(
            CategoricalIndex(["a", "c", "c"], categories=["a", "b", "c", "d"])
        )
        exp = CategoricalIndex(["a", "c", "c"], categories=["a", "b", "c", "d"])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 2, 2], dtype=np.intp))

    def test_reindex_empty_index(self):
        # See GH16770
        c = CategoricalIndex([])
        res, indexer = c.reindex(["a", "b"])
        tm.assert_index_equal(res, Index(["a", "b"]), exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([-1, -1], dtype=np.intp))
