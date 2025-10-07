import pytest
import pandas as pd
import numpy as np
import pandas._testing as tm



class TestMultiIndexInsertLevel:

    def setup_method(self):
        self.simple_idx = pd.MultiIndex.from_tuples(
            [('A', 1), ('B', 2), ('C', 3)], names=['level1', 'level2']
        )
        self.empty_idx = pd.MultiIndex.from_tuples([], names=['level1', 'level2'])

    def test_insert_level_basic(self):
        result = self.simple_idx.insert_level(0, 'new_value')
        expected = pd.MultiIndex.from_tuples(
            [('new_value', 'A', 1), ('new_value', 'B', 2), ('new_value', 'C', 3)],
            names=[None, 'level1', 'level2']
        )
        tm.assert_index_equal(result, expected)

        result = self.simple_idx.insert_level(1, 'middle')
        expected = pd.MultiIndex.from_tuples(
            [('A', 'middle', 1), ('B', 'middle', 2), ('C', 'middle', 3)],
            names=['level1', None, 'level2']
        )
        tm.assert_index_equal(result, expected)

    def test_insert_level_with_different_values(self):
        new_values = ['X', 'Y', 'Z']
        result = self.simple_idx.insert_level(1, new_values)
        expected = pd.MultiIndex.from_tuples(
            [('A', 'X', 1), ('B', 'Y', 2), ('C', 'Z', 3)],
            names=['level1', None, 'level2']
        )
        tm.assert_index_equal(result, expected)

    def test_insert_level_with_name(self):
        result = self.simple_idx.insert_level(0, 'new_val', name='new_level')
        assert result.names[0] == 'new_level'

    def test_insert_level_edge_positions(self):
        result_start = self.simple_idx.insert_level(0, 'start')
        assert result_start.nlevels == 3

        result_end = self.simple_idx.insert_level(2, 'end')
        assert result_end.nlevels == 3

    def test_insert_level_error_cases(self):
        with pytest.raises(ValueError, match="position must be between"):
            self.simple_idx.insert_level(5, 'invalid')

        with pytest.raises(ValueError, match="position must be between"):
            self.simple_idx.insert_level(-1, 'invalid')

        with pytest.raises(ValueError, match="Length of values must match"):
            self.simple_idx.insert_level(1, ['too', 'few'])

    def test_insert_level_with_different_data_types(self):
        result_int = self.simple_idx.insert_level(1, 100)

        result_float = self.simple_idx.insert_level(1, 1.5)

        result_none = self.simple_idx.insert_level(1, None)

        assert result_int.nlevels == 3
        assert result_float.nlevels == 3
        assert result_none.nlevels == 3

    def test_insert_level_preserves_original(self):
        original = self.simple_idx.copy()
        result = self.simple_idx.insert_level(1, 'temp')

        tm.assert_index_equal(original, self.simple_idx)

        assert result.nlevels == original.nlevels + 1

        def test_debug_names():
            idx = pd.MultiIndex.from_tuples(
                [('A', 1), ('B', 2), ('C', 3)],
                names=['level1', 'level2']
            )
            print("Original names:", idx.names)

            result = idx.insert_level(0, 'new_value')
            print("Result names:", result.names)

            expected = pd.MultiIndex.from_tuples(
                [('new_value', 'A', 1), ('new_value', 'B', 2), ('new_value', 'C', 3)],
                names=[None, 'level1', 'level2']
            )
            print("Expected names:", expected.names)