import pytest

from pandas import (
    DatetimeIndex,
    Index,
    to_datetime,
)
import pandas._testing as tm


class TestIndexReplace:
    def test_replace_scalar(self):
        idx = Index([1, 2, 3])

        result = idx.replace(2, 9)
        expected = Index([1, 9, 3])

        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "data, to_replace, value, expected",
        [
            ([1, 2, 3], [1, 3], ["a", "c"], ["a", 2, "c"]),
            ([1, 2, 3], 1, "a", ["a", 2, 3]),
            ([1, None, 2], [1, 2], "a", ["a", None, "a"]),
        ],
    )
    def test_replace_list_and_scalar(self, data, to_replace, value, expected):
        idx = Index(data)
        expected = Index(expected, dtype="object")

        result = idx.replace(to_replace=to_replace, value=value)

        tm.assert_index_equal(result, expected)

    def test_replace_dict(self):
        idx = Index(["a", "b", "c"])

        result = idx.replace({"b": "x"})
        expected = Index(["a", "x", "c"])

        tm.assert_index_equal(result, expected)

    def test_replace_dict_and_value_raises(self):
        idx = Index([1, 2, 3])

        with pytest.raises(ValueError, match="dict-like"):
            idx.replace({1: "a", 3: "c"}, "x")

    def test_name_preserved(self):
        idx = Index([0, 1], name="foo")

        result = idx.replace(1, 0)
        expected = Index([0, 0], name="foo")

        tm.assert_index_equal(result, expected)

    def test_preserves_subclass(self):
        idx = DatetimeIndex(["2020-01-01", "2020-01-02"])

        result = idx.replace("2020", "2021")

        assert isinstance(result, DatetimeIndex)

    def test_preserves_specialized_dtype(self):
        idx = to_datetime(["2020-01-01", "2020-01-02"])

        result = idx.replace("2020-01-01", "2020-01-03")
        expected = to_datetime(["2020-01-03", "2020-01-02"])

        tm.assert_index_equal(result, expected)
        assert isinstance(result, DatetimeIndex)

    def test_doc_example_scalar(self):
        idx = Index(["a", "b", "c", "b"])

        result = idx.replace("b", "x")
        expected = Index(["a", "x", "c", "x"])

        tm.assert_index_equal(result, expected)

    def test_doc_example_dict(self):
        idx = Index(["a", "b", "c", "b"])

        result = idx.replace({"a": "foo", "c": "bar"})
        expected = Index(["foo", "b", "bar", "b"])

        tm.assert_index_equal(result, expected)

    def test_doc_example_regex(self):
        idx = Index(["cat", "dog", "cow"])

        result = idx.replace("c.*", "animal", regex=True)
        expected = Index(["animal", "dog", "animal"])

        tm.assert_index_equal(result, expected)
