from unittest import TestCase
from pandas.core.frame import DataFrame


class TestDataFrameValidate(TestCase):
    """Tests for error handling related to data types of method arguments."""
    df = DataFrame({'a': [1, 2], 'b': [3, 4]})

    def test_validate_bool_args(self):
        # Tests for error handling related to boolean arguments.
        invalid_values = [1, "True", [1, 2, 3], 5.0]

        for value in invalid_values:
            with self.assertRaises(ValueError):
                self.df.query('a > b', inplace=value)

            with self.assertRaises(ValueError):
                self.df.eval('a + b', inplace=value)

            with self.assertRaises(ValueError):
                self.df.set_index(keys=['a'], inplace=value)

            with self.assertRaises(ValueError):
                self.df.reset_index(inplace=value)

            with self.assertRaises(ValueError):
                self.df.dropna(inplace=value)

            with self.assertRaises(ValueError):
                self.df.drop_duplicates(inplace=value)

            with self.assertRaises(ValueError):
                self.df.sort_values(by=['a'], inplace=value)
