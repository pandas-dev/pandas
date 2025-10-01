from __future__ import annotations

from io import StringIO

import pytest

from pandas.errors import ParserWarning
import pandas.util._test_decorators as td

from pandas import (
    NA,
    DataFrame,
    DatetimeIndex,
    Series,
    to_datetime,
)
import pandas._testing as tm

from pandas.io.xml import read_xml


@pytest.fixture(params=[pytest.param("lxml", marks=td.skip_if_no("lxml")), "etree"])
def parser(request):
    return request.param


@pytest.fixture(
    params=[None, {"book": ["category", "title", "author", "year", "price"]}]
)
def iterparse(request):
    return request.param


def read_xml_iterparse(data, temp_file, **kwargs):
    with open(temp_file, "w", encoding="utf-8") as f:
        f.write(data)
    return read_xml(temp_file, **kwargs)


xml_types = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <shape>square</shape>
    <degrees>00360</degrees>
    <sides>4.0</sides>
   </row>
  <row>
    <shape>circle</shape>
    <degrees>00360</degrees>
    <sides/>
  </row>
  <row>
    <shape>triangle</shape>
    <degrees>00180</degrees>
    <sides>3.0</sides>
  </row>
</data>"""

xml_dates = """<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <shape>square</shape>
    <degrees>00360</degrees>
    <sides>4.0</sides>
    <date>2020-01-01</date>
   </row>
  <row>
    <shape>circle</shape>
    <degrees>00360</degrees>
    <sides/>
    <date>2021-01-01</date>
  </row>
  <row>
    <shape>triangle</shape>
    <degrees>00180</degrees>
    <sides>3.0</sides>
    <date>2022-01-01</date>
  </row>
</data>"""


# DTYPE


def test_dtype_single_str(parser, temp_file):
    df_result = read_xml(StringIO(xml_types), dtype={"degrees": "str"}, parser=parser)
    df_iter = read_xml_iterparse(
        xml_types,
        parser=parser,
        dtype={"degrees": "str"},
        iterparse={"row": ["shape", "degrees", "sides"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_dtypes_all_str(parser, temp_file):
    df_result = read_xml(StringIO(xml_dates), dtype="string", parser=parser)
    df_iter = read_xml_iterparse(
        xml_dates,
        parser=parser,
        dtype="string",
        iterparse={"row": ["shape", "degrees", "sides", "date"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": ["4.0", None, "3.0"],
            "date": ["2020-01-01", "2021-01-01", "2022-01-01"],
        },
        dtype="string",
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_dtypes_with_names(parser, temp_file):
    df_result = read_xml(
        StringIO(xml_dates),
        names=["Col1", "Col2", "Col3", "Col4"],
        dtype={"Col2": "string", "Col3": "Int64", "Col4": "datetime64[ns]"},
        parser=parser,
    )
    df_iter = read_xml_iterparse(
        xml_dates,
        parser=parser,
        names=["Col1", "Col2", "Col3", "Col4"],
        dtype={"Col2": "string", "Col3": "Int64", "Col4": "datetime64[ns]"},
        iterparse={"row": ["shape", "degrees", "sides", "date"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "Col1": ["square", "circle", "triangle"],
            "Col2": Series(["00360", "00360", "00180"], dtype="string"),
            "Col3": Series([4.0, NA, 3.0], dtype="Int64"),
            "Col4": DatetimeIndex(
                ["2020-01-01", "2021-01-01", "2022-01-01"], dtype="M8[ns]"
            ),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_dtype_nullable_int(parser, temp_file):
    df_result = read_xml(StringIO(xml_types), dtype={"sides": "Int64"}, parser=parser)
    df_iter = read_xml_iterparse(
        xml_types,
        parser=parser,
        dtype={"sides": "Int64"},
        iterparse={"row": ["shape", "degrees", "sides"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": Series([4.0, NA, 3.0], dtype="Int64"),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_dtype_float(parser, temp_file):
    df_result = read_xml(StringIO(xml_types), dtype={"degrees": "float"}, parser=parser)
    df_iter = read_xml_iterparse(
        xml_types,
        parser=parser,
        dtype={"degrees": "float"},
        iterparse={"row": ["shape", "degrees", "sides"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": Series([360, 360, 180]).astype("float"),
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_wrong_dtype(xml_books, parser, iterparse):
    with pytest.raises(
        ValueError, match=('Unable to parse string "Everyday Italian" at position 0')
    ):
        read_xml(
            xml_books, dtype={"title": "Int64"}, parser=parser, iterparse=iterparse
        )


def test_both_dtype_converters(parser, temp_file):
    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    with tm.assert_produces_warning(ParserWarning, match="Both a converter and dtype"):
        df_result = read_xml(
            StringIO(xml_types),
            dtype={"degrees": "str"},
            converters={"degrees": str},
            parser=parser,
        )
        df_iter = read_xml_iterparse(
            xml_types,
            dtype={"degrees": "str"},
            converters={"degrees": str},
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides"]},
            temp_file=temp_file,
        )

        tm.assert_frame_equal(df_result, df_expected)
        tm.assert_frame_equal(df_iter, df_expected)


# CONVERTERS


def test_converters_str(parser, temp_file):
    df_result = read_xml(
        StringIO(xml_types), converters={"degrees": str}, parser=parser
    )
    df_iter = read_xml_iterparse(
        xml_types,
        parser=parser,
        converters={"degrees": str},
        iterparse={"row": ["shape", "degrees", "sides"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_converters_date(parser, temp_file):
    convert_to_datetime = lambda x: to_datetime(x)
    df_result = read_xml(
        StringIO(xml_dates), converters={"date": convert_to_datetime}, parser=parser
    )
    df_iter = read_xml_iterparse(
        xml_dates,
        parser=parser,
        converters={"date": convert_to_datetime},
        iterparse={"row": ["shape", "degrees", "sides", "date"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_wrong_converters_type(xml_books, parser, iterparse):
    with pytest.raises(TypeError, match=("Type converters must be a dict or subclass")):
        read_xml(
            xml_books, converters={"year", str}, parser=parser, iterparse=iterparse
        )


def test_callable_func_converters(xml_books, parser, iterparse):
    with pytest.raises(TypeError, match=("'float' object is not callable")):
        read_xml(
            xml_books, converters={"year": float()}, parser=parser, iterparse=iterparse
        )


def test_callable_str_converters(xml_books, parser, iterparse):
    with pytest.raises(TypeError, match=("'str' object is not callable")):
        read_xml(
            xml_books, converters={"year": "float"}, parser=parser, iterparse=iterparse
        )


# PARSE DATES


def test_parse_dates_column_name(parser, temp_file):
    df_result = read_xml(StringIO(xml_dates), parse_dates=["date"], parser=parser)
    df_iter = read_xml_iterparse(
        xml_dates,
        parser=parser,
        parse_dates=["date"],
        iterparse={"row": ["shape", "degrees", "sides", "date"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_parse_dates_column_index(parser, temp_file):
    df_result = read_xml(StringIO(xml_dates), parse_dates=[3], parser=parser)
    df_iter = read_xml_iterparse(
        xml_dates,
        parser=parser,
        parse_dates=[3],
        iterparse={"row": ["shape", "degrees", "sides", "date"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_parse_dates_true(parser, temp_file):
    df_result = read_xml(StringIO(xml_dates), parse_dates=True, parser=parser)

    df_iter = read_xml_iterparse(
        xml_dates,
        parser=parser,
        parse_dates=True,
        iterparse={"row": ["shape", "degrees", "sides", "date"]},
        temp_file=temp_file,
    )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": ["2020-01-01", "2021-01-01", "2022-01-01"],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)
    tm.assert_frame_equal(df_iter, df_expected)


def test_day_first_parse_dates(parser, temp_file):
    xml = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <shape>square</shape>
    <degrees>00360</degrees>
    <sides>4.0</sides>
    <date>31/12/2020</date>
   </row>
  <row>
    <shape>circle</shape>
    <degrees>00360</degrees>
    <sides/>
    <date>31/12/2021</date>
  </row>
  <row>
    <shape>triangle</shape>
    <degrees>00180</degrees>
    <sides>3.0</sides>
    <date>31/12/2022</date>
  </row>
</data>"""

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-12-31", "2021-12-31", "2022-12-31"]),
        }
    )

    with tm.assert_produces_warning(
        UserWarning, match="Parsing dates in %d/%m/%Y format"
    ):
        df_result = read_xml(StringIO(xml), parse_dates=["date"], parser=parser)
        df_iter = read_xml_iterparse(
            xml,
            parse_dates=["date"],
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
            temp_file=temp_file,
        )

        tm.assert_frame_equal(df_result, df_expected)
        tm.assert_frame_equal(df_iter, df_expected)


def test_wrong_parse_dates_type(xml_books, parser, iterparse):
    with pytest.raises(TypeError, match="Only booleans and lists are accepted"):
        read_xml(xml_books, parse_dates={"date"}, parser=parser, iterparse=iterparse)
