import pytest

import pandas as pd


class TestDet:
    def test_det_empty(self):
        df = pd.DataFrame()
        msg = "Empty Dataframe/matrix has no determinant"
        with pytest.raises(ValueError, match=msg):
            df.det()

    def test_det_1x1(self):
        df = pd.DataFrame([[1]])
        result = df.det()
        expected = 1
        assert type(result) == type(expected)
        assert result == expected

    def test_det_2x2(self):
        df = pd.DataFrame([[1, 2], [3, 4]])
        result = df.det()
        expected = -2
        assert type(result) == type(expected)
        assert result == expected

    def test_det_3x3(self):
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = df.det()
        expected = 0
        assert type(result) == type(expected)
        assert result == expected

    def test_det_with_different_primitive_dtypes(self):
        df = pd.DataFrame([[1, 2], [3, 4]], dtype=float)
        result = df.det()
        expected = -2.0
        assert type(result) == type(expected)
        assert result == expected

    def test_det_with_non_square_matrix(self):
        df = pd.DataFrame([[1, 2], [3, 4], [5, 6]])
        msg = "The DataFrame must be square"
        with pytest.raises(ValueError, match=msg):
            df.det()

    def test_det_with_valid_values_but_not_dtype(self):
        df = pd.DataFrame([[1, 2], [3, 4]], dtype=str)
        msg = "The DataFrame must only contain numeric values"
        with pytest.raises(ValueError, match=msg):
            df.det()

    def test_det_with_values_from_mixed_valid_dtypes(self):
        df = pd.DataFrame([[2, 3.0, 0], [1, 2.0, 1], [1, 1.0, 1.0]])
        result = df.det()
        expected = 2.0
        assert type(result) == type(expected)
        assert result == expected

    def test_det_with_bigger_input(self):
        df = pd.DataFrame(
            [
                [3, 7, 6, 6, 1],
                [4, 1, 8, 3, 4],
                [1, 4, 7, 9, 6],
                [1, 3, 2, 5, 0],
                [6, 1, 0, 2, 3],
            ]
        )
        result = df.det()
        expected = -3992
        assert type(result) == type(expected)
        assert result == expected

    def test_det_with_linearly_dependent_rows(self):
        # The determinant of a matrix with linearly dependent rows is 0.
        # The algorithm would have given a 0.0 (because the df has float values)
        # if it tried to calculate the determinant recursively.
        # If it gives a 0, if df has float values,
        # it means it detected the linearly dependent rows
        df = pd.DataFrame([[1.0, 2, 3], [2, 4.0, 6], [3, 6, 9.0]])
        result = df.det()
        expected = 0
        assert type(result) == type(expected)
        assert result == expected
