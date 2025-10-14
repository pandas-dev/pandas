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
        "position, value",
        [
            (0, "start"),
            (2, "end"),
        ],
    )
    def test_insert_level_edge_positions(self, position, value):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        result = simple_idx.insert_level(position, value)
        assert result.nlevels == 3

    @pytest.mark.parametrize(
        "position, value, expected_error",
        [
            (5, "invalid", "position must be between"),
            (-1, "invalid", "position must be between"),
            (1, ["too", "few"], "Length of values must match"),
        ],
    )
    def test_insert_level_error_cases(self, position, value, expected_error):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        with pytest.raises(ValueError, match=expected_error):
            simple_idx.insert_level(position, value)

    @pytest.mark.parametrize(
        "value",
        [100, 1.5, None],
    )
    def test_insert_level_with_different_data_types(self, value):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        result = simple_idx.insert_level(1, value)
        assert result.nlevels == 3

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
        assert result.nlevels == 3
        assert len(result) == 0

    def test_insert_level_with_different_values(self):
        simple_idx = pd.MultiIndex.from_tuples(
            [("A", 1), ("B", 2), ("C", 3)], names=["level1", "level2"]
        )

        new_values = ["X", "Y", "Z"]
        result = simple_idx.insert_level(1, new_values)
        expected = pd.MultiIndex.from_tuples(
            [("A", "X", 1), ("B", "Y", 2), ("C", "Z", 3)],
            names=["level1", None, "level2"],
        )
        tm.assert_index_equal(result, expected)
