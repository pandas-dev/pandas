import pytest

from pandas import DataFrame

from pandas.io.formats.format import DataFrameFormatter
from pandas.io.formats.latex import (
    RegularTableBuilder,
    RowBodyIterator,
    RowHeaderIterator,
    RowStringConverter,
)


class TestTableBuilder:
    @pytest.fixture
    def dataframe(self):
        return DataFrame({"a": [1, 2], "b": ["b1", "b2"]})

    @pytest.fixture
    def table_builder(self, dataframe):
        return RegularTableBuilder(formatter=DataFrameFormatter(dataframe))

    def test_create_row_iterator(self, table_builder):
        iterator = table_builder._create_row_iterator(over="header")
        assert isinstance(iterator, RowHeaderIterator)

    def test_create_body_iterator(self, table_builder):
        iterator = table_builder._create_row_iterator(over="body")
        assert isinstance(iterator, RowBodyIterator)

    def test_create_body_wrong_kwarg_raises(self, table_builder):
        with pytest.raises(ValueError, match="must be either 'header' or 'body'"):
            table_builder._create_row_iterator(over="SOMETHING BAD")


class TestRowStringConverter:
    @pytest.mark.parametrize(
        "row_num, expected",
        [
            (0, r"{} &  Design &  ratio &  xy \\"),
            (1, r"0 &       1 &      4 &  10 \\"),
            (2, r"1 &       2 &      5 &  11 \\"),
        ],
    )
    def test_get_strrow_normal_without_escape(self, row_num, expected):
        df = DataFrame({r"Design": [1, 2, 3], r"ratio": [4, 5, 6], r"xy": [10, 11, 12]})
        row_string_converter = RowStringConverter(
            formatter=DataFrameFormatter(df, escape=True),
        )
        assert row_string_converter.get_strrow(row_num=row_num) == expected

    @pytest.mark.parametrize(
        "row_num, expected",
        [
            (0, r"{} &  Design \# &  ratio, \% &  x\&y \\"),
            (1, r"0 &         1 &         4 &   10 \\"),
            (2, r"1 &         2 &         5 &   11 \\"),
        ],
    )
    def test_get_strrow_normal_with_escape(self, row_num, expected):
        df = DataFrame(
            {r"Design #": [1, 2, 3], r"ratio, %": [4, 5, 6], r"x&y": [10, 11, 12]}
        )
        row_string_converter = RowStringConverter(
            formatter=DataFrameFormatter(df, escape=True),
        )
        assert row_string_converter.get_strrow(row_num=row_num) == expected

    @pytest.mark.parametrize(
        "row_num, expected",
        [
            (0, r"{} & \multicolumn{2}{r}{c1} & \multicolumn{2}{r}{c2} & c3 \\"),
            (1, r"{} &  0 &  1 &  0 &  1 &  0 \\"),
            (2, r"0 &  0 &  5 &  0 &  5 &  0 \\"),
        ],
    )
    def test_get_strrow_multindex_multicolumn(self, row_num, expected):
        df = DataFrame(
            {
                ("c1", 0): {x: x for x in range(5)},
                ("c1", 1): {x: x + 5 for x in range(5)},
                ("c2", 0): {x: x for x in range(5)},
                ("c2", 1): {x: x + 5 for x in range(5)},
                ("c3", 0): {x: x for x in range(5)},
            }
        )

        row_string_converter = RowStringConverter(
            formatter=DataFrameFormatter(df),
            multicolumn=True,
            multicolumn_format="r",
            multirow=True,
        )

        assert row_string_converter.get_strrow(row_num=row_num) == expected
