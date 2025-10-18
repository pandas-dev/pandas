import pytest

import pandas as pd
import pandas._testing as tm


class TestMultiIndexInsertLevel:
    @pytest.mark.parametrize(
        "position, value, name, expected_tuples, expected_names",
        [
            (
                0,
                "new_value",
                None,
                [("new_value", "A", 1), ("new_value", "B", 2), ("new_value", "C", 3)],
                [None, "level1", "level2"],
            ),
            (
                1,
                "middle",
                None,
                [("A", "middle", 1), ("B", "middle", 2), ("C", "middle", 3)],
                ["level1", None, "level2"],
            ),
            (
                0,
                "new_val",
                "new_level",
                [("new_val", "A", 1), ("new_val", "B", 2), ("new_val", "C", 3)],
                ["new_level", "level1", "level2"],
            ),
            (
                1,
                "middle",
                "custom_name",
                [("A", "middle", 1), ("B", "middle", 2), ("C", "middle", 3)],
                ["level1", "custom_name", "level2"],
            ),
            (
                0,
                "start",
                None,
                [("start", "A", 1), ("start", "B", 2), ("start", "C", 3)],
                [None, "level1", "level2"],
            ),
            (
                2,
                "end",
                None,
                [("A", 1, "end"), ("B", 2, "end"), ("C", 3, "end")],
                ["level1", "level2", None],
            ),
            (
                1,
                100,
                None,
                [("A", 100, 1), ("B", 100, 2), ("C", 100, 3)],
                ["level1", None, "level2"],
            ),
            (
                1,
                1.5,
                None,
                [("A", 1.5, 1), ("B", 1.5, 2), ("C", 1.5, 3)],
                ["level1", None, "level2"],
            ),
            (
                1,
                None,
                None,
                [("A", None, 1), ("B", None, 2), ("C", None, 3)],
                ["level1", None, "level2"],
            ),
            (
                1,
                ["X", "Y", "Z"],
                None,
                [("A", "X", 1), ("B", "Y", 2), ("C", "Z", 3)],
                ["level1", None, "level2"],
            ),
            (
                0,
                "",
                "empty_string",
                [("", "A", 1), ("", "B", 2), ("", "C", 3)],
                ["empty_string", "level1", "level2"],
            ),
            (
                1,
                True,
                None,
                [("A", True, 1), ("B", True, 2), ("C", True, 3)],
                ["level1", None, "level2"],
            ),
        ],
    )
    def test_insert_level_basic(
        self, position, value, name, expected_tuples, expected_names
    ):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        result = simple_idx.insert_level(position, value, name=name)
        expected = pd.MultiIndex.from_tuples(expected_tuples, names=expected_names)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "position, value, expected_error",
        [
            (5, "invalid", "position must be between"),
            (-1, "invalid", "position must be between"),
            (1, ["too", "few"], "Length of values must match"),
            (3, "value", "position must be between"),
        ],
    )
    def test_insert_level_error_cases(self, position, value, expected_error):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        with pytest.raises(ValueError, match=expected_error):
            simple_idx.insert_level(position, value)

    def test_insert_level_preserves_original(self):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        original = simple_idx.copy()
        simple_idx.insert_level(1, "temp")

        tm.assert_index_equal(original, simple_idx)

    def test_insert_level_empty_index(self):
        empty_idx = pd.MultiIndex.from_tuples([], names=["level1", "level2"])

        result = empty_idx.insert_level(0, [])
        expected = pd.MultiIndex.from_tuples([], names=[None, "level1", "level2"])
        tm.assert_index_equal(result, expected)

    def test_insert_level_single_element(self):
        single_idx = pd.MultiIndex.from_tuples([("A", 1)], names=["level1", "level2"])

        result = single_idx.insert_level(1, "middle")
        expected = pd.MultiIndex.from_tuples(
            [("A", "middle", 1)], names=["level1", None, "level2"]
        )
        tm.assert_index_equal(result, expected)
