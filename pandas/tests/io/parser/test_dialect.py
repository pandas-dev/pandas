# -*- coding: utf-8 -*-

"""
Tests that dialects are properly handled during parsing
for all of the parsers defined in parsers.py
"""

import csv

import pytest

from pandas.compat import StringIO
from pandas.errors import ParserWarning

from pandas import DataFrame
import pandas.util.testing as tm


def test_dialect(all_parsers):
    parser = all_parsers
    data = """\
label1,label2,label3
index1,"a,c,e
index2,b,d,f
"""

    dia = csv.excel()
    dia.quoting = csv.QUOTE_NONE

    # Conflicting dialect quoting.
    with tm.assert_produces_warning(ParserWarning):
        df = parser.read_csv(StringIO(data), dialect=dia)

    data = """\
label1,label2,label3
index1,a,c,e
index2,b,d,f
"""
    exp = parser.read_csv(StringIO(data))
    exp.replace("a", "\"a", inplace=True)
    tm.assert_frame_equal(df, exp)


def test_dialect_str(all_parsers):
    dialect_name = "mydialect"
    parser = all_parsers
    data = """\
fruit:vegetable
apple:broccoli
pear:tomato
"""
    exp = DataFrame({
        "fruit": ["apple", "pear"],
        "vegetable": ["broccoli", "tomato"]
    })
    csv.register_dialect(dialect_name, delimiter=":")

    # Conflicting dialect delimiter.
    with tm.assert_produces_warning(ParserWarning):
        df = parser.read_csv(StringIO(data), dialect=dialect_name)

    tm.assert_frame_equal(df, exp)
    csv.unregister_dialect(dialect_name)


def test_invalid_dialect(all_parsers):
    class InvalidDialect(object):
        pass

    data = "a\n1"
    parser = all_parsers
    msg = "Invalid dialect"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), dialect=InvalidDialect)


@pytest.mark.parametrize("delimiter", [",", "."])
def test_dialect_conflict(all_parsers, delimiter):
    data = "a,b\n1,2"
    dialect = "excel"
    parser = all_parsers

    expected = DataFrame({"a": [1], "b": [2]})
    warning_klass = None if delimiter == "," else ParserWarning

    with tm.assert_produces_warning(warning_klass):
        result = parser.read_csv(StringIO(data),
                                 delimiter=delimiter,
                                 dialect=dialect)
        tm.assert_frame_equal(result, expected)
