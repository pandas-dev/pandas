"""
Tests that dialects are properly handled during parsing
for all of the parsers defined in parsers.py
"""

import csv
from io import StringIO

import pytest

from pandas.errors import ParserWarning

import pandas as pd
import pandas._testing as tm

from pandas.io.parsers.readers import _bool_kwargs


@pytest.fixture
def custom_dialect():
    dialect_name = "weird"
    dialect_kwargs = {
        "doublequote": False,
        "escapechar": "~",
        "delimiter": ":",
        "skipinitialspace": False,
        "quotechar": "`",
        "quoting": 3,
    }
    return dialect_name, dialect_kwargs


def test_dialect(all_parsers):
    parser = all_parsers
    data = """\
label1,label2,label3
index1,"a,c,e
index2,b,d,f
"""

    dia = csv.excel()
    dia.quoting = csv.QUOTE_NONE

    if parser.engine == "pyarrow":
        msg = "The 'dialect' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), dialect=dia)
        return

    df = parser.read_csv(StringIO(data), dialect=dia)

    data = """\
label1,label2,label3
index1,a,c,e
index2,b,d,f
"""
    exp = parser.read_csv(StringIO(data))
    exp.replace("a", '"a', inplace=True)
    tm.assert_frame_equal(df, exp)


def test_dialect_str(all_parsers):
    dialect_name = "mydialect"
    parser = all_parsers
    data = """\
fruit:vegetable
apple:broccoli
pear:tomato
"""
    exp = pd.DataFrame(
        {"fruit": ["apple", "pear"], "vegetable": ["broccoli", "tomato"]}
    )

    with tm.with_csv_dialect(dialect_name, delimiter=":"):
        if parser.engine == "pyarrow":
            msg = "The 'dialect' option is not supported with the 'pyarrow' engine"
            with pytest.raises(ValueError, match=msg):
                parser.read_csv(StringIO(data), dialect=dialect_name)
            return

        df = parser.read_csv(StringIO(data), dialect=dialect_name)
        tm.assert_frame_equal(df, exp)


def test_invalid_dialect(all_parsers):
    class InvalidDialect:
        pass

    data = "a\n1"
    parser = all_parsers
    msg = "Invalid dialect"

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), dialect=InvalidDialect)


@pytest.mark.parametrize(
    "arg",
    [None, "doublequote", "escapechar", "skipinitialspace", "quotechar", "quoting"],
)
@pytest.mark.parametrize("value", ["dialect", "default", "other"])
def test_dialect_conflict_except_delimiter(all_parsers, custom_dialect, arg, value):
    # see gh-23761.
    dialect_name, dialect_kwargs = custom_dialect
    parser = all_parsers

    expected = pd.DataFrame({"a": [1], "b": [2]})
    data = "a:b\n1:2"

    warning_klass = None
    kwds = {}

    # arg=None tests when we pass in the dialect without any other arguments.
    if arg is not None:
        if value == "dialect":  # No conflict --> no warning.
            kwds[arg] = dialect_kwargs[arg]
        elif value == "default":  # Default --> no warning.
            from pandas.io.parsers.base_parser import parser_defaults

            kwds[arg] = parser_defaults[arg]
        else:  # Non-default + conflict with dialect --> warning.
            warning_klass = ParserWarning
            # skipinitialspace's dialect value and parser default are both
            # False, so a valid True still conflicts; for the others no valid
            # conflicting value exists, so use a bad one
            kwds[arg] = True if arg == "skipinitialspace" else "blah"

    with tm.with_csv_dialect(dialect_name, **dialect_kwargs):
        if arg in _bool_kwargs and isinstance(kwds.get(arg), str):
            # GH#68318 doublequote has no valid conflicting value left (the
            # dialect says False, the parser default True), so its conflict is
            # only reachable with a non-bool, which is now rejected
            msg = f'For argument "{arg}" expected type bool'
            with pytest.raises(ValueError, match=msg):
                parser.read_csv(StringIO(data), dialect=dialect_name, **kwds)
            return

        if parser.engine == "pyarrow":
            msg = "The 'dialect' option is not supported with the 'pyarrow' engine"
            with pytest.raises(ValueError, match=msg):
                parser.read_csv_check_warnings(
                    # No warning bc we raise
                    None,
                    "Conflicting values for",
                    StringIO(data),
                    dialect=dialect_name,
                    **kwds,
                )
            return

        result = parser.read_csv_check_warnings(
            warning_klass,
            "Conflicting values for",
            StringIO(data),
            dialect=dialect_name,
            **kwds,
        )
        tm.assert_frame_equal(result, expected)


def test_dialect_supplies_non_bool(all_parsers):
    # GH#68318 a dialect object is not required to normalize its attributes the
    # way csv.register_dialect does, so the post-merge check is the only thing
    # standing between a non-bool it supplies and the parser
    parser = all_parsers

    class BadDialect:
        delimiter = ","
        doublequote = "False"
        escapechar = None
        skipinitialspace = False
        quotechar = '"'
        quoting = 0

    if parser.engine == "pyarrow":
        msg = "The 'dialect' option is not supported with the 'pyarrow' engine"
    else:
        msg = 'For argument "doublequote" expected type bool'

    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO('a\n"x""y"\n'), dialect=BadDialect)


@pytest.mark.parametrize(
    "kwargs,warning_klass",
    [
        ({"sep": ","}, None),  # sep is default --> sep_override=True
        ({"sep": "."}, ParserWarning),  # sep isn't default --> sep_override=False
        ({"delimiter": ":"}, None),  # No conflict
        ({"delimiter": None}, None),  # Default arguments --> sep_override=True
        ({"delimiter": ","}, ParserWarning),  # Conflict
        ({"delimiter": "."}, ParserWarning),  # Conflict
    ],
    ids=[
        "sep-override-true",
        "sep-override-false",
        "delimiter-no-conflict",
        "delimiter-default-arg",
        "delimiter-conflict",
        "delimiter-conflict2",
    ],
)
def test_dialect_conflict_delimiter(all_parsers, custom_dialect, kwargs, warning_klass):
    # see gh-23761.
    dialect_name, dialect_kwargs = custom_dialect
    parser = all_parsers

    expected = pd.DataFrame({"a": [1], "b": [2]})
    data = "a:b\n1:2"

    with tm.with_csv_dialect(dialect_name, **dialect_kwargs):
        if parser.engine == "pyarrow":
            msg = "The 'dialect' option is not supported with the 'pyarrow' engine"
            with pytest.raises(ValueError, match=msg):
                parser.read_csv_check_warnings(
                    # no warning bc we raise
                    None,
                    "Conflicting values for 'delimiter'",
                    StringIO(data),
                    dialect=dialect_name,
                    **kwargs,
                )
            return
        result = parser.read_csv_check_warnings(
            warning_klass,
            "Conflicting values for 'delimiter'",
            StringIO(data),
            dialect=dialect_name,
            **kwargs,
        )
        tm.assert_frame_equal(result, expected)
