"""
Tests for the JSON engine adapter: orjson when installed, stdlib otherwise.
"""

from io import StringIO
import re

import numpy as np
import pytest

from pandas._libs import writers

from pandas import (
    DataFrame,
    NaT,
    Series,
    Timestamp,
    read_json,
)
import pandas._testing as tm

from pandas.io.json import _engine
from pandas.io.json._encode import (
    EncodeOptions,
    encode,
)
from pandas.io.json._engine import (
    dumps,
    loads,
)


@pytest.fixture(params=["orjson", "stdlib"])
def engine(request, monkeypatch):
    if request.param == "orjson":
        pytest.importorskip("orjson")
    else:
        monkeypatch.setattr(_engine, "orjson", None)
    return request.param


class TestLoads:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("[NaN, Infinity, -Infinity]", [None, np.inf, -np.inf]),
            ("[1e400, -1e400]", [np.inf, -np.inf]),
            ('{"a": {"b": [NaN]}}', {"a": {"b": [None]}}),
            (
                "[4.9e-324, 0.1, 1.7976931348623157e308]",
                [5e-324, 0.1, 1.7976931348623157e308],
            ),
        ],
    )
    def test_loads_beyond_strict_json(self, engine, text, expected):
        assert loads(text) == expected

    def test_loads_types_exact(self, engine):
        result = loads('[1, 1.0, "1", true, null]')
        assert [type(value) for value in result] == [int, float, str, bool, type(None)]

    @pytest.mark.parametrize("text", ["[1,]", "[,1]", "[]]", "{}\n a", '{"a":b}'])
    def test_loads_invalid(self, engine, text):
        msg = "|".join(["Expecting value", "Extra data"])
        with pytest.raises(ValueError, match=msg):
            loads(text)

    def test_read_json_int_in_range(self, engine):
        text = (
            f'{{"a":[{2**64 - 1}],"b":[{-(2**63)}],"c":["{2**64}"],'
            '"d":[1577836800000000000],"e":[0.1234567890123456789],'
            '"f":[1234567890123456789e1],"g":[-1.5e-19]}'
        )
        result = read_json(StringIO(text))
        expected = DataFrame(
            {
                "a": np.array([2**64 - 1], dtype="uint64"),
                "b": [-(2**63)],
                "c": [float(2**64)],
                "d": [1577836800000000000],
                "e": [0.1234567890123456789],
                "f": [1234567890123456789e1],
                "g": [-1.5e-19],
            }
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "text, msg",
        [
            (f'{{"a":[1,{2**64}]}}', "Value is too big! at position 8"),
            (f'{{"a":[{-(2**63) - 1}]}}', "Value is too small at position 6"),
            (f'{{"a":["{2**64}",{2**64}]}}', "Value is too big! at position 29"),
        ],
    )
    def test_read_json_int_out_of_range(self, engine, text, msg):
        with pytest.raises(ValueError, match=re.escape(msg)):
            read_json(StringIO(text))

    def test_read_json_nan_inf(self, engine):
        data = StringIO('{"a":[NaN,Infinity,-Infinity,1.5],"b":["x",null,NaN,"y"]}')
        result = read_json(data)
        expected = DataFrame(
            {"a": [np.nan, np.inf, -np.inf, 1.5], "b": ["x", None, None, "y"]}
        )
        tm.assert_frame_equal(result, expected)

    def test_read_json_nan_in_strings(self, engine):
        # only the bare literal is a missing value
        data = StringIO('{"a":["NaN","x\\"NaN","Infinity"],"b":[NaN,1.5,2.5]}')
        result = read_json(data, dtype=False)
        expected = DataFrame({"a": ["NaN", 'x"NaN', "Infinity"], "b": [None, 1.5, 2.5]})
        tm.assert_frame_equal(result, expected)

    def test_loads_nan_uses_orjson(self, monkeypatch):
        # NaN is rewritten to null so that orjson reads the text; the stdlib
        # decoder is only used for Infinity, which orjson cannot produce
        pytest.importorskip("orjson")
        calls = []
        stdlib_loads = _engine.json.loads

        def spy(*args, **kwargs):
            calls.append(args[0])
            return stdlib_loads(*args, **kwargs)

        monkeypatch.setattr(_engine.json, "loads", spy)
        assert loads('{"a":[NaN,1],"b":"NaN"}') == {"a": [None, 1], "b": "NaN"}
        assert calls == []
        assert loads("[Infinity,NaN]") == [np.inf, None]
        assert loads("[1e400]") == [np.inf]
        assert len(calls) == 2

    @pytest.mark.parametrize(
        "text, expected",
        [
            (b"[]", (-1, 0, 0, False)),
            (b'[NaN,"NaN",NaN,-Infinity]', (-1, 0, 2, True)),
            (b"[9.9e307,99e306,1e-400,1.5]", (-1, 0, 0, False)),
            # 1e308 and 1.8e308 look alike to the scanner; only the second
            # overflows, so both go to the stdlib decoder
            (b"[1e308]", (-1, 0, 0, True)),
            (b"[1.8e308]", (-1, 0, 0, True)),
            (b"[123456789e300]", (-1, 0, 0, True)),
            (b"[18446744073709551615,-9223372036854775808]", (-1, 0, 0, False)),
            (b'{"NaN":18446744073709551616}', (7, 1, 0, False)),
            (b"[NaN,-9223372036854775809]", (5, -1, 1, False)),
        ],
    )
    def test_json_scan(self, text, expected):
        assert writers.json_scan(text) == expected

    def test_json_nan_to_null(self):
        text = b'{"NaN":[NaN,"x\\"NaN",NaN,1]}'
        result = writers.json_nan_to_null(text, 2)
        assert result == b'{"NaN":[null,"x\\"NaN",null,1]}'
        assert writers.json_nan_to_null(b"", 0) == b""


class TestDumps:
    def test_dumps_basic(self, engine):
        obj = {"a": [1, 2.5, "x/y", True, None, "é"], "b": {"c": []}}
        assert dumps(obj) == '{"a":[1,2.5,"x/y",true,null,"é"],"b":{"c":[]}}'

    def test_dumps_indent(self, engine):
        obj = {"a": [1, {"b": []}], "c": {}}
        expected = (
            '{\n  "a": [\n    1,\n    {\n      "b": []\n    }\n  ],\n  "c": {}\n}'
        )
        assert dumps(obj, indent=2) == expected

    def test_dumps_big_int(self, engine):
        assert dumps([2**70, -(2**70)]) == f"[{2**70},{-(2**70)}]"

    def test_dumps_float_repr(self, engine):
        result = dumps([0.1 + 0.2, 1e16, 1.5, -0.0])
        assert result == "[0.30000000000000004,1e+16,1.5,-0.0]"


class TestJSONZip:
    def test_zip_object(self):
        keys = b'["a","b\\"c","d/e"]'
        values = b'[1,[2,{"x":3}],"1,2"]'
        result = writers.json_zip_object(keys, values)
        assert result == b'{"a":1,"b\\"c":[2,{"x":3}],"d/e":"1,2"}'
        assert writers.json_zip_object(b"[]", b"[]") == b"{}"

    def test_zip_object_length_mismatch(self):
        with pytest.raises(ValueError, match="same length"):
            writers.json_zip_object(b'["a"]', b"[1,2]")

    def test_zip_records(self):
        keys = b'["a","b"]'
        rows = b'[[1,2],["x,y",null],[[1],{"k":"]"}]]'
        result = writers.json_zip_records(keys, rows)
        expected = b'[{"a":1,"b":2},{"a":"x,y","b":null},{"a":[1],"b":{"k":"]"}}]'
        assert result == expected
        assert writers.json_zip_records(b'["a"]', b"[]") == b"[]"

    @pytest.mark.parametrize("rows", [b"[[1]]", b"[[1,2,3]]", b"[1]"])
    def test_zip_records_bad_row(self, rows):
        with pytest.raises(ValueError, match="one value per key"):
            writers.json_zip_records(b'["a","b"]', rows)

    @pytest.mark.parametrize("orient", ["records", "index", "columns"])
    def test_fragments_match_dicts(self, orient):
        # the text-splicing path and the dict-building path give the same text
        pytest.importorskip("orjson")
        df = DataFrame(
            {"a": [1, 2], "b/c": ['x"y', None], "d": [[1, {"k": "]"}], 1.5]},
            index=["i,j", "k"],
        )
        with_fragments = df.to_json(orient=orient)
        options = EncodeOptions(iso_dates=True, date_unit="ms", fragments=False)
        without = dumps(encode(df, orient, options)).replace("/", "\\/")
        assert with_fragments == without


class TestEngineParity:
    """Both engines must produce the same text for the same pandas object."""

    @pytest.fixture
    def frame(self):
        return DataFrame(
            {
                "i": [1, -2, 3],
                "f": [1.5, np.nan, np.inf],
                "s": ["a", "é/b", None],
                "b": [True, False, True],
                "d": [Timestamp("2020-01-01"), NaT, Timestamp("2020-01-03")],
                "o": [[1, {"x": np.nan}], (2, 3), {5: "five"}],
                "big": [2**70, 1, 2],
            }
        )

    @pytest.mark.parametrize(
        "orient", ["split", "records", "index", "columns", "values", "table"]
    )
    @pytest.mark.parametrize("indent", [0, 2])
    def test_to_json_same_output(self, monkeypatch, frame, orient, indent):
        pytest.importorskip("orjson")
        with_orjson = frame.to_json(orient=orient, date_format="iso", indent=indent)
        monkeypatch.setattr(_engine, "orjson", None)
        with_stdlib = frame.to_json(orient=orient, date_format="iso", indent=indent)
        assert with_orjson == with_stdlib

    def test_series_to_json_same_output(self, monkeypatch):
        pytest.importorskip("orjson")
        ser = Series([1.5, np.nan, 2**70, "x"], index=["a", "b", "c", "d"])
        with_orjson = ser.to_json(orient="split")
        monkeypatch.setattr(_engine, "orjson", None)
        with_stdlib = ser.to_json(orient="split")
        assert with_orjson == with_stdlib

    def test_read_json_same_result(self, monkeypatch):
        pytest.importorskip("orjson")
        text = '{"a":[1,NaN,3],"b":[1.5,Infinity,2.5],"c":["x","y",null]}'
        with_orjson = read_json(StringIO(text))
        monkeypatch.setattr(_engine, "orjson", None)
        with_stdlib = read_json(StringIO(text))
        tm.assert_frame_equal(with_orjson, with_stdlib)
