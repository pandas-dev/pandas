from unittest import TestCase
from pandas.core.series import Series


class TestSeriesValidate(TestCase):
    """Tests for error handling related to data types of method arguments."""
    s = Series([1, 2, 3, 4, 5])

    def test_validate_bool_args(self):
        # Tests for error handling related to boolean arguments.
        invalid_values = [1, "True", [1, 2, 3], 5.0]

        for value in invalid_values:
            with self.assertRaises(ValueError):
                self.s.reset_index(inplace=value)

            with self.assertRaises(ValueError):
                self.s._set_name(name='hello', inplace=value)

            with self.assertRaises(ValueError):
                self.s.sort_values(inplace=value)

            with self.assertRaises(ValueError):
                self.s.sort_index(inplace=value)

            with self.assertRaises(ValueError):
                self.s.sort_index(inplace=value)

            with self.assertRaises(ValueError):
                self.s.rename(inplace=value)

            with self.assertRaises(ValueError):
                self.s.dropna(inplace=value)
