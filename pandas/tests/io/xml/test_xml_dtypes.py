from __future__ import annotations

import pytest

from pandas.errors import ParserWarning
import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
    to_datetime,
)
import pandas._testing as tm

from pandas.io.xml import read_xml


@pytest.fixture(params=[pytest.param("lxml", marks=td.skip_if_no("lxml")), "etree"])
def parser(request):
    return request.param


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


def test_dtype_single_str(parser):
    df_result = read_xml(xml_types, dtype={"degrees": "str"}, parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_dtypes_all_str(parser):
    df_result = read_xml(xml_dates, dtype="string", parser=parser)

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


def test_dtypes_with_names(parser):
    df_result = read_xml(
        xml_dates,
        names=["Col1", "Col2", "Col3", "Col4"],
        dtype={"Col2": "string", "Col3": "Int64", "Col4": "datetime64"},
        parser=parser,
    )

    df_expected = DataFrame(
        {
            "Col1": ["square", "circle", "triangle"],
            "Col2": Series(["00360", "00360", "00180"]).astype("string"),
            "Col3": Series([4.0, float("nan"), 3.0]).astype("Int64"),
            "Col4": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_dtype_nullable_int(parser):
    df_result = read_xml(xml_types, dtype={"sides": "Int64"}, parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": Series([4.0, float("nan"), 3.0]).astype("Int64"),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_dtype_float(parser):
    df_result = read_xml(xml_types, dtype={"degrees": "float"}, parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": Series([360, 360, 180]).astype("float"),
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_wrong_dtype(parser):
    with pytest.raises(
        ValueError, match=('Unable to parse string "square" at position 0')
    ):
        read_xml(xml_types, dtype={"shape": "Int64"}, parser=parser)


def test_both_dtype_converters(parser):
    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    with tm.assert_produces_warning(ParserWarning, match="Both a converter and dtype"):
        df_result = read_xml(
            xml_types,
            dtype={"degrees": "str"},
            converters={"degrees": str},
            parser=parser,
        )

        tm.assert_frame_equal(df_result, df_expected)


# CONVERTERS


def test_converters_str(parser):
    df_result = read_xml(xml_types, converters={"degrees": str}, parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": ["00360", "00360", "00180"],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_converters_date(parser):
    convert_to_datetime = lambda x: to_datetime(x)
    df_result = read_xml(
        xml_dates, converters={"date": convert_to_datetime}, parser=parser
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


def test_wrong_converters_type(parser):
    with pytest.raises(TypeError, match=("Type converters must be a dict or subclass")):
        read_xml(xml_types, converters={"degrees", str}, parser=parser)


def test_callable_func_converters(parser):
    with pytest.raises(TypeError, match=("'float' object is not callable")):
        read_xml(xml_types, converters={"degrees": float()}, parser=parser)


def test_callable_str_converters(parser):
    with pytest.raises(TypeError, match=("'str' object is not callable")):
        read_xml(xml_types, converters={"degrees": "float"}, parser=parser)


# PARSE DATES


def test_parse_dates_column_name(parser):
    df_result = read_xml(xml_dates, parse_dates=["date"], parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_parse_dates_column_index(parser):
    df_result = read_xml(xml_dates, parse_dates=[3], parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_parse_dates_true(parser):
    df_result = read_xml(xml_dates, parse_dates=True, parser=parser)

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": ["2020-01-01", "2021-01-01", "2022-01-01"],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_parse_dates_dictionary(parser):
    xml = """<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <shape>square</shape>
    <degrees>360</degrees>
    <sides>4.0</sides>
    <year>2020</year>
    <month>12</month>
    <day>31</day>
   </row>
  <row>
    <shape>circle</shape>
    <degrees>360</degrees>
    <sides/>
    <year>2021</year>
    <month>12</month>
    <day>31</day>
  </row>
  <row>
    <shape>triangle</shape>
    <degrees>180</degrees>
    <sides>3.0</sides>
    <year>2022</year>
    <month>12</month>
    <day>31</day>
  </row>
</data>"""

    df_result = read_xml(
        xml, parse_dates={"date_end": ["year", "month", "day"]}, parser=parser
    )

    df_expected = DataFrame(
        {
            "date_end": to_datetime(["2020-12-31", "2021-12-31", "2022-12-31"]),
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_result, df_expected)


def test_day_first_parse_dates(parser):
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
        UserWarning, match="Parsing '31/12/2020' in DD/MM/YYYY format"
    ):
        df_result = read_xml(xml, parse_dates=["date"], parser=parser)
        tm.assert_frame_equal(df_result, df_expected)


def test_wrong_parse_dates_type(parser):
    with pytest.raises(
        TypeError, match=("Only booleans, lists, and dictionaries are accepted")
    ):
        read_xml(xml_dates, parse_dates={"date"}, parser=parser)
