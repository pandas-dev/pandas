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

    def test_multiindex_replace_scalar(self):
        idx = MultiIndex.from_product(
            [["a", "b"], [1, 2]],
            names=["letter", "number"],
        )

        result = idx.replace("a", "x")

        expected = MultiIndex.from_tuples(
            [("x", 1), ("x", 2), ("b", 1), ("b", 2)],
            names=["letter", "number"],
        )

        tm.assert_index_equal(result, expected)

    def test_multiindex_replace_dict(self):
        idx = MultiIndex.from_product(
            [["a", "b"], [1, 2]],
            names=["letter", "number"],
        )

        result = idx.replace({1: 10})

        expected = MultiIndex.from_tuples(
            [("a", 10), ("a", 2), ("b", 10), ("b", 2)],
            names=["letter", "number"],
        )

        tm.assert_index_equal(result, expected)

    def test_multiindex_replace_preserves_names(self):
        idx = MultiIndex.from_product(
            [["foo", "bar"], [0, 1]],
            names=["level_1", "level_2"],
        )

        result = idx.replace("foo", "baz")

        assert result.names == ["level_1", "level_2"]
