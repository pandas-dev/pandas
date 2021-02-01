from io import BytesIO, StringIO

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.xml import read_xml

geom_df = DataFrame(
    {
        "shape": ["square", "circle", "triangle"],
        "degrees": [360, 360, 180],
        "sides": [4, np.nan, 3],
    }
)

planet_df = DataFrame(
    {
        "planet": [
            "Mercury",
            "Venus",
            "Earth",
            "Mars",
            "Jupiter",
            "Saturn",
            "Uranus",
            "Neptune",
        ],
        "type": [
            "terrestrial",
            "terrestrial",
            "terrestrial",
            "terrestrial",
            "gas giant",
            "gas giant",
            "ice giant",
            "ice giant",
        ],
        "location": [
            "inner",
            "inner",
            "inner",
            "inner",
            "outer",
            "outer",
            "outer",
            "outer",
        ],
        "mass": [
            0.330114,
            4.86747,
            5.97237,
            0.641712,
            1898.187,
            568.3174,
            86.8127,
            102.4126,
        ],
    }
)

from_file_expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <index>0</index>
    <category>cooking</category>
    <title>Everyday Italian</title>
    <author>Giada De Laurentiis</author>
    <year>2005</year>
    <price>30.0</price>
  </row>
  <row>
    <index>1</index>
    <category>children</category>
    <title>Harry Potter</title>
    <author>J K. Rowling</author>
    <year>2005</year>
    <price>29.99</price>
  </row>
  <row>
    <index>2</index>
    <category>web</category>
    <title>Learning XML</title>
    <author>Erik T. Ray</author>
    <year>2003</year>
    <price>39.95</price>
  </row>
</data>"""


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_file_output_str_read(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_file = read_xml(filename, parser=parser)

    with tm.ensure_clean("test.xml") as path:
        df_file.to_xml(path, parser=parser)
        with open(path, "rb") as f:
            output = f.read().decode("utf-8").strip()

        # etree and lxml differs on quotes and case in xml declaration
        output = output.replace(
            '<?xml version="1.0" encoding="utf-8"?',
            "<?xml version='1.0' encoding='utf-8'?",
        )

        assert output == from_file_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_file_output_bytes_read(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_file = read_xml(filename, parser=parser)

    with tm.ensure_clean("test.xml") as path:
        df_file.to_xml(path, parser=parser)
        with open(path, "rb") as f:
            output = f.read().decode("utf-8").strip()

        # etree and lxml differs on quotes and case in xml declaration
        output = output.replace(
            '<?xml version="1.0" encoding="utf-8"?',
            "<?xml version='1.0' encoding='utf-8'?",
        )

        assert output == from_file_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_str_output(datapath, parser):
    filename = datapath("io", "data", "xml", "books.xml")
    df_file = read_xml(filename, parser=parser)

    output = df_file.to_xml()

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == from_file_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_index_false(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <category>cooking</category>
    <title>Everyday Italian</title>
    <author>Giada De Laurentiis</author>
    <year>2005</year>
    <price>30.0</price>
  </row>
  <row>
    <category>children</category>
    <title>Harry Potter</title>
    <author>J K. Rowling</author>
    <year>2005</year>
    <price>29.99</price>
  </row>
  <row>
    <category>web</category>
    <title>Learning XML</title>
    <author>Erik T. Ray</author>
    <year>2003</year>
    <price>39.95</price>
  </row>
</data>"""

    filename = datapath("io", "data", "xml", "books.xml")
    df_file = read_xml(filename, parser=parser)

    with tm.ensure_clean("test.xml") as path:
        df_file.to_xml(path, index=False, parser=parser)
        with open(path, "rb") as f:
            output = f.read().decode("utf-8").strip()

        # etree and lxml differs on quotes and case in xml declaration
        output = output.replace(
            '<?xml version="1.0" encoding="utf-8"?',
            "<?xml version='1.0' encoding='utf-8'?",
        )

        assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_index_false_rename_row_root(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<books>
  <book>
    <category>cooking</category>
    <title>Everyday Italian</title>
    <author>Giada De Laurentiis</author>
    <year>2005</year>
    <price>30.0</price>
  </book>
  <book>
    <category>children</category>
    <title>Harry Potter</title>
    <author>J K. Rowling</author>
    <year>2005</year>
    <price>29.99</price>
  </book>
  <book>
    <category>web</category>
    <title>Learning XML</title>
    <author>Erik T. Ray</author>
    <year>2003</year>
    <price>39.95</price>
  </book>
</books>"""

    filename = datapath("io", "data", "xml", "books.xml")
    df_file = read_xml(filename, parser=parser)

    with tm.ensure_clean("test.xml") as path:
        df_file.to_xml(
            path, index=False, root_name="books", row_name="book", parser=parser
        )
        with open(path, "rb") as f:
            output = f.read().decode("utf-8").strip()

        # etree and lxml differs on quotes and case in xml declaration
        output = output.replace(
            '<?xml version="1.0" encoding="utf-8"?',
            "<?xml version='1.0' encoding='utf-8'?",
        )

        assert output == expected


na_expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <index>0</index>
    <shape>square</shape>
    <degrees>360</degrees>
    <sides>4.0</sides>
  </row>
  <row>
    <index>1</index>
    <shape>circle</shape>
    <degrees>360</degrees>
    <sides/>
  </row>
  <row>
    <index>2</index>
    <shape>triangle</shape>
    <degrees>180</degrees>
    <sides>3.0</sides>
  </row>
</data>"""


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_na_elem_output(datapath, parser):
    output = geom_df.to_xml(parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == na_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_na_empty_str_elem_option(datapath, parser):
    output = geom_df.to_xml(na_rep="", parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == na_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_na_empty_elem_option(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <index>0</index>
    <shape>square</shape>
    <degrees>360</degrees>
    <sides>4.0</sides>
  </row>
  <row>
    <index>1</index>
    <shape>circle</shape>
    <degrees>360</degrees>
    <sides>0.0</sides>
  </row>
  <row>
    <index>2</index>
    <shape>triangle</shape>
    <degrees>180</degrees>
    <sides>3.0</sides>
  </row>
</data>"""

    output = geom_df.to_xml(na_rep="0.0", parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_attrs_cols_nan_output(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row index="0" shape="square" degrees="360" sides="4.0"/>
  <row index="1" shape="circle" degrees="360"/>
  <row index="2" shape="triangle" degrees="180" sides="3.0"/>
</data>"""

    output = geom_df.to_xml(attr_cols=["shape", "degrees", "sides"], parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_attrs_cols_prefix(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<doc:data xmlns:doc="http://example.xom">
  <doc:row doc:index="0" doc:shape="square" doc:degrees="360" doc:sides="4.0"/>
  <doc:row doc:index="1" doc:shape="circle" doc:degrees="360"/>
  <doc:row doc:index="2" doc:shape="triangle" doc:degrees="180" doc:sides="3.0"/>
</doc:data>"""

    output = geom_df.to_xml(
        attr_cols=["index", "shape", "degrees", "sides"],
        namespaces={"doc": "http://example.xom"},
        prefix="doc",
        parser=parser,
    )

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_attrs_unknown_column(parser):
    with pytest.raises(KeyError, match=("no valid column")):
        geom_df.to_xml(attr_cols=["shape", "degreees", "sides"], parser=parser)


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_attrs_wrong_type(parser):
    with pytest.raises(TypeError, match=("is not a valid type for attr_cols")):
        geom_df.to_xml(attr_cols='"shape", "degreees", "sides"', parser=parser)


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_elems_cols_nan_output(datapath, parser):
    elems_cols_expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <degrees>360</degrees>
    <sides>4.0</sides>
    <shape>square</shape>
  </row>
  <row>
    <degrees>360</degrees>
    <sides/>
    <shape>circle</shape>
  </row>
  <row>
    <degrees>180</degrees>
    <sides>3.0</sides>
    <shape>triangle</shape>
  </row>
</data>"""

    output = geom_df.to_xml(
        index=False, elem_cols=["degrees", "sides", "shape"], parser=parser
    )

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == elems_cols_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_elems_unknown_column(parser):
    with pytest.raises(KeyError, match=("no valid column")):
        geom_df.to_xml(elem_cols=["shape", "degreees", "sides"], parser=parser)


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_elems_wrong_type(parser):
    with pytest.raises(TypeError, match=("is not a valid type for elem_cols")):
        geom_df.to_xml(elem_cols='"shape", "degreees", "sides"', parser=parser)


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_elems_and_attrs_cols(datapath, parser):
    elems_cols_expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row shape="square">
    <degrees>360</degrees>
    <sides>4.0</sides>
  </row>
  <row shape="circle">
    <degrees>360</degrees>
    <sides/>
  </row>
  <row shape="triangle">
    <degrees>180</degrees>
    <sides>3.0</sides>
  </row>
</data>"""

    output = geom_df.to_xml(
        index=False,
        elem_cols=["degrees", "sides"],
        attr_cols=["shape"],
        parser=parser,
    )

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == elems_cols_expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_hierarchical_columns(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <location>inner</location>
    <type>terrestrial</type>
    <count_mass>4</count_mass>
    <sum_mass>11.811666</sum_mass>
    <mean_mass>2.9529165</mean_mass>
  </row>
  <row>
    <location>outer</location>
    <type>gas giant</type>
    <count_mass>2</count_mass>
    <sum_mass>2466.5044</sum_mass>
    <mean_mass>1233.2522</mean_mass>
  </row>
  <row>
    <location>outer</location>
    <type>ice giant</type>
    <count_mass>2</count_mass>
    <sum_mass>189.2253</sum_mass>
    <mean_mass>94.61265</mean_mass>
  </row>
  <row>
    <location>All</location>
    <type/>
    <count_mass>8</count_mass>
    <sum_mass>2667.541366</sum_mass>
    <mean_mass>333.44267075</mean_mass>
  </row>
</data>"""

    pvt = planet_df.pivot_table(
        index=["location", "type"],
        values="mass",
        aggfunc=["count", "sum", "mean"],
        margins=True,
    )

    output = pvt.to_xml(parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_hierarchical_attrs_columns(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row location="inner" type="terrestrial" count_mass="4" \
sum_mass="11.811666" mean_mass="2.9529165"/>
  <row location="outer" type="gas giant" count_mass="2" \
sum_mass="2466.5044" mean_mass="1233.2522"/>
  <row location="outer" type="ice giant" count_mass="2" \
sum_mass="189.2253" mean_mass="94.61265"/>
  <row location="All" type="" count_mass="8" \
sum_mass="2667.541366" mean_mass="333.44267075"/>
</data>"""

    pvt = planet_df.pivot_table(
        index=["location", "type"],
        values="mass",
        aggfunc=["count", "sum", "mean"],
        margins=True,
    )

    output = pvt.to_xml(attr_cols=list(pvt.reset_index().columns.values), parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_multi_index(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row>
    <location>inner</location>
    <type>terrestrial</type>
    <count>4</count>
    <sum>11.811666</sum>
    <mean>2.9529165</mean>
  </row>
  <row>
    <location>outer</location>
    <type>gas giant</type>
    <count>2</count>
    <sum>2466.5044</sum>
    <mean>1233.2522</mean>
  </row>
  <row>
    <location>outer</location>
    <type>ice giant</type>
    <count>2</count>
    <sum>189.2253</sum>
    <mean>94.61265</mean>
  </row>
</data>"""

    agg = planet_df.groupby(["location", "type"])["mass"].agg(["count", "sum", "mean"])

    output = agg.to_xml(parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_multi_index_attrs_cols(datapath, parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data>
  <row location="inner" type="terrestrial" count="4" sum="11.811666" mean="2.9529165"/>
  <row location="outer" type="gas giant" count="2" sum="2466.5044" mean="1233.2522"/>
  <row location="outer" type="ice giant" count="2" sum="189.2253" mean="94.61265"/>
</data>"""

    agg = planet_df.groupby(["location", "type"])["mass"].agg(["count", "sum", "mean"])

    output = agg.to_xml(attr_cols=list(agg.reset_index().columns.values), parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_default_namespace(parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<data xmlns="http://example.com">
  <row>
    <index>0</index>
    <shape>square</shape>
    <degrees>360</degrees>
    <sides>4.0</sides>
  </row>
  <row>
    <index>1</index>
    <shape>circle</shape>
    <degrees>360</degrees>
    <sides/>
  </row>
  <row>
    <index>2</index>
    <shape>triangle</shape>
    <degrees>180</degrees>
    <sides>3.0</sides>
  </row>
</data>"""

    output = geom_df.to_xml(namespaces={"": "http://example.com"}, parser=parser)

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_namespace_prefix(parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<doc:data xmlns:doc="http://example.com">
  <doc:row>
    <doc:index>0</doc:index>
    <doc:shape>square</doc:shape>
    <doc:degrees>360</doc:degrees>
    <doc:sides>4.0</doc:sides>
  </doc:row>
  <doc:row>
    <doc:index>1</doc:index>
    <doc:shape>circle</doc:shape>
    <doc:degrees>360</doc:degrees>
    <doc:sides/>
  </doc:row>
  <doc:row>
    <doc:index>2</doc:index>
    <doc:shape>triangle</doc:shape>
    <doc:degrees>180</doc:degrees>
    <doc:sides>3.0</doc:sides>
  </doc:row>
</doc:data>"""

    output = geom_df.to_xml(
        namespaces={"doc": "http://example.com"}, prefix="doc", parser=parser
    )

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_missing_prefix_in_nmsp(parser):
    with pytest.raises(KeyError, match=("prefix is not included in namespaces")):

        geom_df.to_xml(
            namespaces={"": "http://example.com"}, prefix="doc", parser=parser
        )


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_namespace_prefix_and_default(parser):
    expected = """\
<?xml version='1.0' encoding='utf-8'?>
<doc:data xmlns="http://example.com" xmlns:doc="http://other.org">
  <doc:row>
    <doc:index>0</doc:index>
    <doc:shape>square</doc:shape>
    <doc:degrees>360</doc:degrees>
    <doc:sides>4.0</doc:sides>
  </doc:row>
  <doc:row>
    <doc:index>1</doc:index>
    <doc:shape>circle</doc:shape>
    <doc:degrees>360</doc:degrees>
    <doc:sides/>
  </doc:row>
  <doc:row>
    <doc:index>2</doc:index>
    <doc:shape>triangle</doc:shape>
    <doc:degrees>180</doc:degrees>
    <doc:sides>3.0</doc:sides>
  </doc:row>
</doc:data>"""

    output = geom_df.to_xml(
        namespaces={"": "http://example.com", "doc": "http://other.org"},
        prefix="doc",
        parser=parser,
    )

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    # etree and lxml differs on order of namespace prefixes
    output = output.replace(
        'xmlns:doc="http://other.org" xmlns="http://example.com"',
        'xmlns="http://example.com" xmlns:doc="http://other.org"',
    )

    assert output == expected


encoding_expected = """\
<?xml version='1.0' encoding='ISO-8859-1'?>
<data>
  <row>
    <index>0</index>
    <rank>1</rank>
    <malename>José</malename>
    <femalename>Sofía</femalename>
  </row>
  <row>
    <index>1</index>
    <rank>2</rank>
    <malename>Luis</malename>
    <femalename>Valentina</femalename>
  </row>
  <row>
    <index>2</index>
    <rank>3</rank>
    <malename>Carlos</malename>
    <femalename>Isabella</femalename>
  </row>
  <row>
    <index>3</index>
    <rank>4</rank>
    <malename>Juan</malename>
    <femalename>Camila</femalename>
  </row>
  <row>
    <index>4</index>
    <rank>5</rank>
    <malename>Jorge</malename>
    <femalename>Valeria</femalename>
  </row>
</data>"""


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_encoding_option_str(datapath, parser):
    filename = datapath("io", "data", "xml", "baby_names.xml")
    df_file = read_xml(filename, parser=parser, encoding="ISO-8859-1").head(5)

    output = df_file.to_xml(encoding="ISO-8859-1")

    # etree and lxml differs on quotes and case in xml declaration
    output = output.replace(
        '<?xml version="1.0" encoding="ISO-8859-1"?',
        "<?xml version='1.0' encoding='ISO-8859-1'?",
    )

    assert output == encoding_expected


@td.skip_if_no("lxml")
def test_correct_encoding_file(datapath):
    filename = datapath("io", "data", "xml", "baby_names.xml")
    df_file = read_xml(filename, encoding="ISO-8859-1")

    with tm.ensure_clean("test.xml") as path:
        df_file.to_xml(path, index=False, encoding="ISO-8859-1")


@pytest.mark.parametrize("parser", ["lxml", "etree"])
@pytest.mark.parametrize("encoding", ["UTF-8", "UTF-16", "ISO-8859-1"])
def test_wrong_encoding_option_lxml(datapath, parser, encoding):
    filename = datapath("io", "data", "xml", "baby_names.xml")
    df_file = read_xml(filename, encoding="ISO-8859-1")

    with tm.ensure_clean("test.xml") as path:
        df_file.to_xml(path, index=False, encoding=encoding, parser=parser)


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_misspelled_encoding(parser):
    with pytest.raises(LookupError, match=("unknown encoding")):
        geom_df.to_xml(parser=parser, encoding="uft-8")


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_xml_declaration_pretty_print(parser):
    expected = """\
<data>
  <row>
    <index>0</index>
    <shape>square</shape>
    <degrees>360</degrees>
    <sides>4.0</sides>
  </row>
  <row>
    <index>1</index>
    <shape>circle</shape>
    <degrees>360</degrees>
    <sides/>
  </row>
  <row>
    <index>2</index>
    <shape>triangle</shape>
    <degrees>180</degrees>
    <sides>3.0</sides>
  </row>
</data>"""

    output = geom_df.to_xml(xml_declaration=False, parser=parser)

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_no_pretty_print_with_decl(parser):
    expected = (
        "<?xml version='1.0' encoding='utf-8'?>\n"
        "<data><row><index>0</index><shape>square</shape>"
        "<degrees>360</degrees><sides>4.0</sides></row><row>"
        "<index>1</index><shape>circle</shape><degrees>360"
        "</degrees><sides/></row><row><index>2</index><shape>"
        "triangle</shape><degrees>180</degrees><sides>3.0</sides>"
        "</row></data>"
    )

    output = geom_df.to_xml(pretty_print=False)

    output = output.replace(
        '<?xml version="1.0" encoding="utf-8"?',
        "<?xml version='1.0' encoding='utf-8'?",
    )

    assert output == expected


@pytest.mark.parametrize("parser", ["lxml", "etree"])
def test_no_pretty_print_no_decl(parser):
    expected = (
        "<data><row><index>0</index><shape>square</shape>"
        "<degrees>360</degrees><sides>4.0</sides></row><row>"
        "<index>1</index><shape>circle</shape><degrees>360"
        "</degrees><sides/></row><row><index>2</index><shape>"
        "triangle</shape><degrees>180</degrees><sides>3.0</sides>"
        "</row></data>"
    )

    output = geom_df.to_xml(xml_declaration=False, pretty_print=False)

    assert output == expected


xsl_expected = """\
<?xml version="1.0" encoding="utf-8"?>
<data>
  <row>
    <field field="index">0</field>
    <field field="shape">square</field>
    <field field="degrees">360</field>
    <field field="sides">4.0</field>
  </row>
  <row>
    <field field="index">1</field>
    <field field="shape">circle</field>
    <field field="degrees">360</field>
    <field field="sides"/>
  </row>
  <row>
    <field field="index">2</field>
    <field field="shape">triangle</field>
    <field field="degrees">180</field>
    <field field="sides">3.0</field>
  </row>
</data>"""


@td.skip_if_no("lxml")
@pytest.mark.parametrize("mode", ["rb", "r"])
def test_stylesheet_file_like(datapath, mode):
    xsl = datapath("io", "data", "xml", "row_field_output.xsl")

    with open(xsl, mode) as f:
        assert geom_df.to_xml(stylesheet=f) == xsl_expected


@td.skip_if_no("lxml")
@pytest.mark.parametrize("mode", ["rb", "r"])
def test_stylesheet_io(datapath, mode):
    xsl = datapath("io", "data", "xml", "row_field_output.xsl")

    with open(xsl, mode) as f:
        xsl_obj = f.read()

    xsl_io = BytesIO(xsl_obj) if isinstance(xsl_obj, bytes) else StringIO(xsl_obj)

    output = geom_df.to_xml(stylesheet=xsl_io)

    assert output == xsl_expected


@td.skip_if_no("lxml")
@pytest.mark.parametrize("mode", ["rb", "r"])
def test_stylesheet_buffered_reader(datapath, mode):
    xsl = datapath("io", "data", "xml", "row_field_output.xsl")

    with open(xsl, mode) as f:
        xsl_obj = f.read()

    output = geom_df.to_xml(stylesheet=xsl_obj)

    assert output == xsl_expected


@td.skip_if_no("lxml")
def test_style_to_csv():
    xsl = """\
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="text" indent="yes" />
    <xsl:strip-space elements="*"/>

    <xsl:param name="delim">,</xsl:param>
    <xsl:template match="/data">
        <xsl:text>,shape,degrees,sides&#xa;</xsl:text>
        <xsl:apply-templates select="row"/>
    </xsl:template>

    <xsl:template match="row">
        <xsl:value-of select="concat(index, $delim, shape, $delim,
                                     degrees, $delim, sides)"/>
         <xsl:text>&#xa;</xsl:text>
    </xsl:template>
</xsl:stylesheet>"""

    out_csv = geom_df.to_csv().strip()
    out_xml = geom_df.to_xml(stylesheet=xsl)

    assert out_csv == out_xml


@td.skip_if_no("lxml")
def test_style_to_string():
    xsl = """\
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="text" indent="yes" />
    <xsl:strip-space elements="*"/>

    <xsl:param name="delim"><xsl:text>               </xsl:text></xsl:param>
    <xsl:template match="/data">
        <xsl:text>      shape  degrees  sides&#xa;</xsl:text>
        <xsl:apply-templates select="row"/>
    </xsl:template>

    <xsl:template match="row">
        <xsl:value-of select="concat(index, ' ',
                                     substring($delim, 1, string-length('triangle')
                                               - string-length(shape) + 1),
                                     shape,
                                     substring($delim, 1, string-length(name(degrees))
                                               - string-length(degrees) + 2),
                                     degrees,
                                     substring($delim, 1, string-length(name(sides))
                                               - string-length(sides) + 2),
                                     sides)"/>
         <xsl:text>&#xa;</xsl:text>
    </xsl:template>
</xsl:stylesheet>"""

    out_str = geom_df.to_string()
    out_xml = geom_df.to_xml(na_rep="NaN", stylesheet=xsl)

    assert out_xml == out_str


@td.skip_if_no("lxml")
def test_style_to_json():
    xsl = """\
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="text" indent="yes" />
    <xsl:strip-space elements="*"/>

    <xsl:param name="quot">"</xsl:param>

    <xsl:template match="/data">
        <xsl:text>{"shape":{</xsl:text>
        <xsl:apply-templates select="descendant::row/shape"/>
        <xsl:text>},"degrees":{</xsl:text>
        <xsl:apply-templates select="descendant::row/degrees"/>
        <xsl:text>},"sides":{</xsl:text>
        <xsl:apply-templates select="descendant::row/sides"/>
        <xsl:text>}}</xsl:text>
    </xsl:template>

    <xsl:template match="shape|degrees|sides">
        <xsl:variable name="val">
            <xsl:if test = ".=''">
                <xsl:value-of select="'null'"/>
            </xsl:if>
            <xsl:if test = "number(text()) = text()">
                <xsl:value-of select="text()"/>
            </xsl:if>
            <xsl:if test = "number(text()) != text()">
                <xsl:value-of select="concat($quot, text(), $quot)"/>
            </xsl:if>
        </xsl:variable>
        <xsl:value-of select="concat($quot, preceding-sibling::index,
                                     $quot,':', $val)"/>
        <xsl:if test="preceding-sibling::index != //row[last()]/index">
            <xsl:text>,</xsl:text>
        </xsl:if>
    </xsl:template>
</xsl:stylesheet>"""

    out_json = geom_df.to_json()
    out_xml = geom_df.to_xml(stylesheet=xsl)

    assert out_json == out_xml


@pytest.mark.skip(
    reason="incorrect <th> tag in <tbody> from to_html() to be skipped until fix"
)
def test_style_to_html():
    xsl = """\
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="xml" omit-xml-declaration="yes" indent="yes" />
    <xsl:strip-space elements="*"/>

    <xsl:template match="/data">
        <table border="1" class="dataframe">
            <thead>
                <tr style="text-align: right;">
                    <th></th>
                    <th>shape</th>
                    <th>degrees</th>
                    <th>sides</th>
                </tr>
            </thead>
            <tbody>
                <xsl:apply-templates select="@*|node()"/>
            </tbody>
        </table>
    </xsl:template>

    <xsl:template match="row">
        <tr>
            <xsl:apply-templates select="@*|node()"/>
        </tr>
    </xsl:template>

    <xsl:template match="row/*">
        <td>
            <xsl:value-of select="text()"/>
        </td>
    </xsl:template>
</xsl:stylesheet>"""

    out_html = geom_df.to_html()
    out_xml = geom_df.to_xml(stylesheet=xsl)

    assert out_html == out_xml
