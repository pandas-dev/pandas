from __future__ import annotations

import pytest

from pandas.errors import ParserError
import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
    to_datetime,
)
import pandas._testing as tm

from pandas.io.common import get_handle
from pandas.io.xml import read_xml


@pytest.fixture(params=[pytest.param("lxml", marks=td.skip_if_no("lxml")), "etree"])
def parser(request):
    return request.param


@pytest.fixture(params=["rb", "r"])
def mode(request):
    return request.param


geom_df = DataFrame(
    {
        "shape": ["square", "circle", "triangle"],
        "degrees": [360, 360, 180],
        "sides": [4, float("nan"), 3],
    }
)

xml_str = """\
<?xml version='1.0' encoding='utf-8'?>
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

xml_prefix_nmsp = """\
<?xml version='1.0' encoding='utf-8'?>
<doc:data xmlns:doc="http://example.com">
  <doc:row>
    <doc:shape>square</doc:shape>
    <doc:degrees>360</doc:degrees>
    <doc:sides>4.0</doc:sides>
  </doc:row>
  <doc:row>
    <doc:shape>circle</doc:shape>
    <doc:degrees>360</doc:degrees>
    <doc:sides/>
  </doc:row>
  <doc:row>
    <doc:shape>triangle</doc:shape>
    <doc:degrees>180</doc:degrees>
    <doc:sides>3.0</doc:sides>
  </doc:row>
</doc:data>"""

bad_xml = """\
<?xml version='1.0' encoding='utf-8'?>
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
"""

# FILE


def test_file(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_iter = read_xml(
        filename,
        parser=parser,
        iterparse={"book": ["category", "title", "year", "author", "price"]},
    )

    df_expected = DataFrame(
        {
            "category": ["cooking", "children", "web"],
            "title": ["Everyday Italian", "Harry Potter", "Learning XML"],
            "author": ["Giada De Laurentiis", "J K. Rowling", "Erik T. Ray"],
            "year": [2005, 2005, 2003],
            "price": [30.00, 29.99, 39.95],
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_file_xpath_compare(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_xpath = read_xml(filename, parser=parser)
    df_iter = read_xml(
        filename,
        parser=parser,
        iterparse={"book": ["category", "title", "author", "year", "price"]},
    )

    tm.assert_frame_equal(df_xpath, df_iter)


# LARGE FILE


@tm.network
@pytest.mark.slow
def test_large_url_xpath_compare(parser):
    with tm.ensure_clean(filename="cta.xml") as path:
        url = (
            "https://data.cityofchicago.org/api/views/"
            "8pix-ypme/rows.xml?accessType=DOWNLOAD"
        )
        (read_xml(url, xpath=".//row/row", parser=parser).to_xml(path, index=False))

        df_xpath = read_xml(path, parser=parser)
        df_iter = read_xml(
            path,
            parser=parser,
            iterparse={
                "row": [
                    "_id",
                    "_uuid",
                    "_position",
                    "_address",
                    "stop_id",
                    "direction_id",
                    "stop_name",
                    "station_name",
                    "station_descriptive_name",
                    "map_id",
                    "ada",
                    "red",
                    "blue",
                    "g",
                    "brn",
                    "p",
                    "pexp",
                    "y",
                    "pnk",
                    "o",
                    "location",
                ]
            },
        )

    tm.assert_frame_equal(df_xpath, df_iter)


# NAMESPACES


def test_namespace_prefix(parser):
    with tm.ensure_clean(filename="xml_prefix_nmsp.xml") as path:
        with open(path, "w") as f:
            f.write(xml_prefix_nmsp)

        df_iter = read_xml(
            path, parser=parser, iterparse={"row": ["shape", "degrees", "sides"]}
        )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_namespace_prefix_xpath_compare(parser):
    with tm.ensure_clean(filename="xml_prefix_nmsp.xml") as path:
        with open(path, "w") as f:
            f.write(xml_prefix_nmsp)

        df_xpath = read_xml(
            path,
            xpath=".//ns:row",
            namespaces={"ns": "http://example.com"},
            parser=parser,
        )
        df_iter = read_xml(
            path, parser=parser, iterparse={"row": ["shape", "degrees", "sides"]}
        )

        tm.assert_frame_equal(df_xpath, df_iter)


def test_default_namespace_xpath_compare(datapath):
    kml = datapath("io", "data", "xml", "cta_rail_lines.kml")

    df_xpath = read_xml(
        kml, xpath=".//k:Placemark", namespaces={"k": "http://www.opengis.net/kml/2.2"}
    )

    df_iter = read_xml(
        kml,
        iterparse={
            "Placemark": [
                "id",
                "name",
                "Snippet",
                "description",
                "styleUrl",
                "MultiGeometry",
            ]
        },
    )

    tm.assert_frame_equal(df_xpath, df_iter)


# ELEMS_ONLY


def test_elems_only(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")

    df_iter = read_xml(
        filename,
        parser=parser,
        iterparse={"book": ["title", "author", "year", "price"]},
    )

    df_expected = DataFrame(
        {
            "title": ["Everyday Italian", "Harry Potter", "Learning XML"],
            "author": ["Giada De Laurentiis", "J K. Rowling", "Erik T. Ray"],
            "year": [2005, 2005, 2003],
            "price": [30.00, 29.99, 39.95],
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_elems_only_xpath_compare(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_xpath = read_xml(filename, elems_only=True, parser=parser)
    df_iter = read_xml(
        filename,
        parser=parser,
        iterparse={"book": ["title", "author", "year", "price"]},
    )

    tm.assert_frame_equal(df_xpath, df_iter)


# ATTRS_ONLY


def test_attrs_only(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_iter = read_xml(filename, parser=parser, iterparse={"book": ["category"]})
    df_expected = DataFrame({"category": ["cooking", "children", "web"]})

    tm.assert_frame_equal(df_iter, df_expected)


def test_attrs_only_xpath_compare(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_xpath = read_xml(filename, attrs_only=True, parser=parser)
    df_iter = read_xml(filename, parser=parser, iterparse={"book": ["category"]})

    tm.assert_frame_equal(df_xpath, df_iter)


# NAMES


def test_names(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")

    df_iter = read_xml(
        filename,
        parser=parser,
        names=["b_category", "b_title", "b_author", "b_year", "b_price"],
        iterparse={"book": ["category", "title", "author", "year", "price"]},
    )

    df_expected = DataFrame(
        {
            "b_category": ["cooking", "children", "web"],
            "b_title": ["Everyday Italian", "Harry Potter", "Learning XML"],
            "b_author": ["Giada De Laurentiis", "J K. Rowling", "Erik T. Ray"],
            "b_year": [2005, 2005, 2003],
            "b_price": [30.00, 29.99, 39.95],
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_names_xpath_compare(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_xpath = read_xml(
        filename,
        parser=parser,
        names=["b_category", "b_title", "b_author", "b_year", "b_price"],
    )
    df_iter = read_xml(
        filename,
        parser=parser,
        names=["b_category", "b_title", "b_author", "b_year", "b_price"],
        iterparse={"book": ["category", "title", "author", "year", "price"]},
    )

    tm.assert_frame_equal(df_xpath, df_iter)


# DTYPE


def test_dtypes(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")

    df_iter = read_xml(
        filename,
        parser=parser,
        dtype={"year": "Int64", "price": "Float64"},
        iterparse={"book": ["category", "title", "year", "author", "price"]},
    )

    df_expected = DataFrame(
        {
            "category": ["cooking", "children", "web"],
            "title": ["Everyday Italian", "Harry Potter", "Learning XML"],
            "author": ["Giada De Laurentiis", "J K. Rowling", "Erik T. Ray"],
            "year": Series([2005, 2005, 2003]).astype("Int64"),
            "price": Series([30.00, 29.99, 39.95]).astype("Float64"),
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_dtypes_xpath_compare(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")

    df_xpath = read_xml(
        filename, parser=parser, dtype={"year": "Int64", "price": "Float64"}
    )

    df_iter = read_xml(
        filename,
        parser=parser,
        dtype={"year": "Int64", "price": "Float64"},
        iterparse={"book": ["category", "title", "year", "author", "price"]},
    )

    tm.assert_frame_equal(df_xpath, df_iter)


# CONVERTERS


def test_converters(parser):
    convert_to_datetime = lambda x: to_datetime(x)
    with tm.ensure_clean(filename="xml_string.xml") as path:
        with open(path, "w") as f:
            f.write(xml_str)

        df_iter = read_xml(
            path,
            converters={"date": convert_to_datetime},
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
        )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_converters_xpath_compare(parser):
    convert_to_datetime = lambda x: to_datetime(x)
    with tm.ensure_clean(filename="xml_string.xml") as path:
        with open(path, "w") as f:
            f.write(xml_str)

        df_xpath = read_xml(
            path, converters={"date": convert_to_datetime}, parser=parser
        )

        df_iter = read_xml(
            path,
            converters={"date": convert_to_datetime},
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
        )

    tm.assert_frame_equal(df_xpath, df_iter)


# PARSE_DATES


def test_date_parse(parser):
    with tm.ensure_clean(filename="xml_string.xml") as path:
        with open(path, "w") as f:
            f.write(xml_str)

        df_iter = read_xml(
            path,
            parse_dates=["date"],
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
        )

    df_expected = DataFrame(
        {
            "shape": ["square", "circle", "triangle"],
            "degrees": [360, 360, 180],
            "sides": [4.0, float("nan"), 3.0],
            "date": to_datetime(["2020-01-01", "2021-01-01", "2022-01-01"]),
        }
    )

    tm.assert_frame_equal(df_iter, df_expected)


def test_date_parse_xpath_compare(parser):
    with tm.ensure_clean(filename="xml_string.xml") as path:
        with open(path, "w") as f:
            f.write(xml_str)

        df_xpath = read_xml(path, parse_dates=["date"], parser=parser)

        df_iter = read_xml(
            path,
            parse_dates=["date"],
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
        )

    tm.assert_frame_equal(df_xpath, df_iter)


# ENCODING


def test_encoding(datapath, parser):
    filename = datapath("io", "data", "xml", "baby_names.xml")

    df_iter = read_xml(
        filename,
        parser=parser,
        encoding="ISO-8859-1",
        iterparse={"row": ["rank", "malename", "femalename"]},
    )

    df_expected = DataFrame(
        {
            "rank": [1, 2, 3, 4, 5],
            "malename": ["José", "Luis", "Carlos", "Juan", "Jorge"],
            "femalename": ["Sofía", "Valentina", "Isabella", "Camila", "Valeria"],
        }
    )

    tm.assert_frame_equal(df_iter.head(), df_expected)


def test_encoding_xpath_compare(datapath, parser):
    filename = datapath("io", "data", "xml", "baby_names.xml")
    df_xpath = read_xml(filename, parser=parser, encoding="ISO-8859-1")

    df_iter = read_xml(
        filename,
        parser=parser,
        encoding="ISO-8859-1",
        iterparse={"row": ["rank", "malename", "femalename"]},
    )

    tm.assert_frame_equal(df_xpath, df_iter)


# STYLESHEET


@td.skip_if_no("lxml")
def test_stylesheet_xpath_compare(datapath):
    kml = datapath("io", "data", "xml", "cta_rail_lines.kml")
    xsl = datapath("io", "data", "xml", "flatten_doc.xsl")

    df_style = read_xml(
        kml,
        xpath=".//k:Placemark",
        namespaces={"k": "http://www.opengis.net/kml/2.2"},
        stylesheet=xsl,
    )

    df_iter = read_xml(
        kml,
        iterparse={
            "Placemark": [
                "id",
                "name",
                "styleUrl",
                "extrude",
                "altitudeMode",
                "coordinates",
            ]
        },
    )

    tm.assert_frame_equal(df_style, df_iter)


# COMPRESSION


def test_compression_compare(parser, compression_only):
    with tm.ensure_clean() as comp_path, tm.ensure_clean() as ext_path:
        geom_df.to_xml(comp_path, parser=parser, compression=compression_only)

        with get_handle(comp_path, "r", compression=compression_only) as handles:
            with open(ext_path, "w") as f:
                f.write(handles.handle.read())

            df_iter = read_xml(
                ext_path,
                parser=parser,
                iterparse={"row": ["shape", "degrees", "sides"]},
                compression=compression_only,
            )

    tm.assert_frame_equal(geom_df, df_iter)


# STORAGE OPTIONS


@tm.network
@pytest.mark.slow
def test_s3_xpath_compare(parser):
    # Python Software Foundation (2019 IRS-990 RETURN)
    s3_path = "s3://irs-form-990/201923199349319487_public.xml"

    df_xpath = read_xml(
        s3_path,
        xpath=".//irs:Form990PartVIISectionAGrp",
        namespaces={"irs": "http://www.irs.gov/efile"},
        parser=parser,
        storage_options={"anon": True},
    )

    with tm.ensure_clean(filename="irs990.xml") as path:
        with get_handle(s3_path, "rb", is_text=False) as handles:
            with open(path, "wb") as f:
                f.write(handles.handle.read())

        df_iter = read_xml(
            path,
            parser=parser,
            iterparse={
                "Form990PartVIISectionAGrp": [
                    "PersonNm",
                    "TitleTxt",
                    "AverageHoursPerWeekRt",
                    "AverageHoursPerWeekRltdOrgRt",
                    "IndividualTrusteeOrDirectorInd",
                    "OfficerInd",
                    "ReportableCompFromOrgAmt",
                    "ReportableCompFromRltdOrgAmt",
                    "OtherCompensationAmt",
                    "HighestCompensatedEmployeeInd",
                ]
            },
        )

    tm.assert_frame_equal(df_xpath, df_iter)


# PARSER ERROR


def test_string_error(parser):
    with pytest.raises(
        ParserError, match=("iterparse is designed for large XML files")
    ):
        read_xml(
            xml_str,
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
        )


def test_file_like_error(datapath, parser, mode):
    filename = datapath("io", "data", "xml", "books.xml")
    with pytest.raises(
        ParserError, match=("iterparse is designed for large XML files")
    ):
        with open(filename) as f:
            read_xml(
                f,
                parser=parser,
                iterparse={"book": ["category", "title", "year", "author", "price"]},
            )


@tm.network
def test_url_path_error(parser):
    url = "https://www.w3schools.com/xml/books.xml"
    with pytest.raises(
        ParserError, match=("iterparse is designed for large XML files")
    ):
        read_xml(
            url,
            parser=parser,
            iterparse={"row": ["shape", "degrees", "sides", "date"]},
        )


def test_compression_error(parser, compression_only):
    with tm.ensure_clean(filename="geom_xml.zip") as path:
        geom_df.to_xml(path, parser=parser, compression=compression_only)

        with pytest.raises(
            ParserError, match=("iterparse is designed for large XML files")
        ):
            read_xml(
                path,
                parser=parser,
                iterparse={"row": ["shape", "degrees", "sides", "date"]},
                compression=compression_only,
            )


@tm.network
@td.skip_if_no("s3fs")
def test_storage_options_error(parser):
    # Python Software Foundation (2019 IRS-990 RETURN)
    s3 = "s3://irs-form-990/201923199349319487_public.xml"
    with pytest.raises(
        ParserError, match=("iterparse is designed for large XML files")
    ):
        read_xml(
            s3,
            parser=parser,
            iterparse={
                "Form990PartVIISectionAGrp": [
                    "PersonNm",
                    "TitleTxt",
                    "AverageHoursPerWeekRt",
                    "AverageHoursPerWeekRltdOrgRt",
                    "IndividualTrusteeOrDirectorInd",
                    "OfficerInd",
                    "ReportableCompFromOrgAmt",
                    "ReportableCompFromRltdOrgAmt",
                    "OtherCompensationAmt",
                ]
            },
            storage_options={"anon": True},
        )


# OTHER EXCEPTIONS


def test_wrong_dict_type(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    with pytest.raises(TypeError, match="list is not a valid type for iterparse"):
        read_xml(
            filename,
            parser=parser,
            iterparse=["category", "title", "year", "author", "price"],
        )


def test_wrong_dict_value(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    with pytest.raises(
        TypeError, match="<class 'str'> is not a valid type for value in iterparse"
    ):
        read_xml(filename, parser=parser, iterparse={"book": "category"})


def test_bad_xml(datapath, parser):
    with tm.ensure_clean(filename="bad.xml") as path:
        with open(path, "w") as f:
            f.write(bad_xml)

        with pytest.raises(
            SyntaxError, match="Extra content at the end of the document"
        ):
            read_xml(
                path,
                parse_dates=["date"],
                iterparse={"row": ["shape", "degrees", "sides", "date"]},
            )


def test_no_result(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    with pytest.raises(
        ParserError, match="No result from selected items in iterparse."
    ):
        read_xml(
            filename,
            parser=parser,
            iterparse={"node": ["attr1", "elem1", "elem2", "elem3"]},
        )
