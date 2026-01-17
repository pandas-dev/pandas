from io import StringIO

from pandas import (
    DataFrame,
    Series,
    read_json,
)
from pandas.testing import (
    assert_frame_equal,
    assert_series_equal,
)


class TestOrjson:
    def test_read_json_very_large_integer(self):
        json_str = """
      [
        {
          "composition": "Nb:100",
          "n_atoms": 128000,
          "space_group": 229,
          "time": 1000000000000000000000000
        }
      ]
      """

        result = read_json(StringIO(json_str), engine="orjson")

        expected = DataFrame(
            {
                "composition": "Nb:100",
                "n_atoms": [128000],
                "space_group": [229],
                "time": [1e24],  # orjson parses very large integers as float
            }
        )

        assert_frame_equal(result, expected)

    def test_read_json_very_large_integer_series(self):
        json_str = """
      [
        {
          "value": 1000000000000000000000000
        }
      ]
      """

        result = read_json(
            StringIO(json_str),
            engine="orjson",
            typ="series",
        )

        expected = Series(
            [{"value": 1e24}]
        )  # orjson parses very large integers as float

        assert_series_equal(result, expected)
