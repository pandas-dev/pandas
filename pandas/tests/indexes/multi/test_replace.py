from pandas import MultiIndex
import pandas._testing as tm


class TestMultiIndexReplace:
    def test_multiindex_replace_elementwise(self):
        idx = MultiIndex.from_tuples([("a", 1), ("b", 2), ("a", 3)])

        result = idx.replace("a", "x")

        expected = MultiIndex.from_tuples([("x", 1), ("b", 2), ("x", 3)])

        tm.assert_index_equal(result, expected)

    def test_multiindex_name_preserved(self):
        idx = MultiIndex.from_tuples([("a", 1), ("b", 2)], names=["letter", "number"])

        result = idx.replace("a", "x")

        expected = MultiIndex.from_tuples(
            [("x", 1), ("b", 2)], names=["letter", "number"]
        )

        tm.assert_index_equal(result, expected)
