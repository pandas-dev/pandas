"""
Tests column conversion functionality during parsing
for all of the parsers defined in parsers.py
"""

from io import StringIO

from dateutil.parser import parse
import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_converters_type_must_be_dict(all_parsers):
    parser = all_parsers
    data = """index,A,B,C,D
foo,2,3,4,5
"""
    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), converters=0)
        return
    with pytest.raises(TypeError, match="Type converters.+"):
        parser.read_csv(StringIO(data), converters=0)


@pytest.mark.parametrize("column", [3, "D"])
@pytest.mark.parametrize(
    "converter",
    [parse, lambda x: int(x.split("/")[2])],  # Produce integer.
)
def test_converters(all_parsers, column, converter):
    parser = all_parsers
    data = """A,B,C,D
a,1,2,01/01/2009
b,3,4,01/02/2009
c,4,5,01/03/2009
"""
    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), converters={column: converter})
        return

    result = parser.read_csv(StringIO(data), converters={column: converter})

    expected = parser.read_csv(StringIO(data))
    expected["D"] = expected["D"].map(converter)

    tm.assert_frame_equal(result, expected)


def test_converters_no_implicit_conv(all_parsers):
    # see gh-2184
    parser = all_parsers
    data = """000102,1.2,A\n001245,2,B"""

    converters = {0: lambda x: x.strip()}

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), header=None, converters=converters)
        return

    result = parser.read_csv(StringIO(data), header=None, converters=converters)

    # Column 0 should not be casted to numeric and should remain as object.
    expected = pd.DataFrame([["000102", 1.2, "A"], ["001245", 2, "B"]])
    tm.assert_frame_equal(result, expected)


def test_converters_euro_decimal_format(all_parsers):
    # see gh-583
    converters = {}
    parser = all_parsers

    data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,7387
2;121,12;14897,76;DEF;uyt;0,3773
3;878,158;108013,434;GHI;rez;2,7356"""
    converters["Number1"] = converters["Number2"] = converters["Number3"] = lambda x: (
        float(x.replace(",", "."))
    )

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), sep=";", converters=converters)
        return

    result = parser.read_csv(StringIO(data), sep=";", converters=converters)
    expected = pd.DataFrame(
        [
            [1, 1521.1541, 187101.9543, "ABC", "poi", 4.7387],
            [2, 121.12, 14897.76, "DEF", "uyt", 0.3773],
            [3, 878.158, 108013.434, "GHI", "rez", 2.7356],
        ],
        columns=["Id", "Number1", "Number2", "Text1", "Text2", "Number3"],
    )
    tm.assert_frame_equal(result, expected)


def test_converters_corner_with_nans(all_parsers):
    parser = all_parsers
    data = """id,score,days
1,2,12
2,2-5,
3,,14+
4,6-12,2"""

    # Example converters.
    def convert_days(x):
        x = x.strip()

        if not x:
            return np.nan

        is_plus = x.endswith("+")

        if is_plus:
            x = int(x[:-1]) + 1
        else:
            x = int(x)

        return x

    def convert_days_sentinel(x):
        x = x.strip()

        if not x:
            return np.nan

        is_plus = x.endswith("+")

        if is_plus:
            x = int(x[:-1]) + 1
        else:
            x = int(x)

        return x

    def convert_score(x):
        x = x.strip()

        if not x:
            return np.nan

        if x.find("-") > 0:
            val_min, val_max = map(int, x.split("-"))
            val = 0.5 * (val_min + val_max)
        else:
            val = float(x)

        return val

    results = []

    for day_converter in [convert_days, convert_days_sentinel]:
        if parser.engine == "pyarrow":
            msg = "The 'converters' option is not supported with the 'pyarrow' engine"
            with pytest.raises(ValueError, match=msg):
                parser.read_csv(
                    StringIO(data),
                    converters={"score": convert_score, "days": day_converter},
                    na_values=["", None],
                )
            continue

        result = parser.read_csv(
            StringIO(data),
            converters={"score": convert_score, "days": day_converter},
            na_values=["", None],
        )
        assert pd.isna(result["days"][1])
        results.append(result)

    if parser.engine != "pyarrow":
        tm.assert_frame_equal(results[0], results[1])


@pytest.mark.parametrize("conv_f", [lambda x: x, str])
def test_converter_index_col_bug(all_parsers, conv_f):
    # see gh-1835 , GH#40589
    parser = all_parsers
    data = "A;B\n1;2\n3;4"

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data), sep=";", index_col="A", converters={"A": conv_f}
            )
        return

    rs = parser.read_csv(
        StringIO(data), sep=";", index_col="A", converters={"A": conv_f}
    )

    xp = pd.DataFrame({"B": [2, 4]}, index=pd.Index(["1", "3"], name="A"))
    tm.assert_frame_equal(rs, xp)


def test_converter_identity_object(all_parsers):
    # GH#40589
    parser = all_parsers
    data = "A,B\n1,2\n3,4"

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), converters={"A": lambda x: x})
        return

    rs = parser.read_csv(StringIO(data), converters={"A": lambda x: x})

    xp = pd.DataFrame({"A": ["1", "3"], "B": [2, 4]})
    tm.assert_frame_equal(rs, xp)


def test_converter_multi_index(all_parsers):
    # GH 42446
    parser = all_parsers
    data = "A,B,B\nX,Y,Z\n1,2,3"

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data),
                header=list(range(2)),
                converters={
                    ("A", "X"): np.int32,
                    ("B", "Y"): np.int32,
                    ("B", "Z"): np.float32,
                },
            )
        return

    result = parser.read_csv(
        StringIO(data),
        header=list(range(2)),
        converters={
            ("A", "X"): np.int32,
            ("B", "Y"): np.int32,
            ("B", "Z"): np.float32,
        },
    )

    expected = pd.DataFrame(
        {
            ("A", "X"): np.int32([1]),
            ("B", "Y"): np.int32([2]),
            ("B", "Z"): np.float32([3]),
        }
    )

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("keep_default_na", [True, False])
@pytest.mark.parametrize(
    "na_values,values",
    [
        ("CAT", ["1", np.nan, "3"]),
        ([-999, "CAT"], ["1", np.nan, "3"]),
        ({"A": ["CAT"]}, ["1", np.nan, "3"]),
        ({0: ["CAT"]}, ["1", np.nan, "3"]),
        # no entry for the converter column, so it gets the defaults or nothing
        ({"B": ["CAT"]}, ["1", "CAT", "3"]),
    ],
)
def test_converters_na_values(all_parsers, na_values, values, keep_default_na):
    # GH#13302 the c engine skipped na_values entirely on a converter column
    parser = all_parsers
    data = "A,B\nx1,x\nCAT,y\nx3,z"
    converters = {"A": lambda x: x.lstrip("x")}
    kwargs = {
        "converters": converters,
        "na_values": na_values,
        "keep_default_na": keep_default_na,
    }

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), **kwargs)
        return

    result = parser.read_csv(StringIO(data), **kwargs)
    expected = pd.DataFrame({"A": values, "B": ["x", "y", "z"]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("keep_default_na", [True, False])
def test_converters_keep_default_na(all_parsers, keep_default_na):
    # GH#13302
    parser = all_parsers
    data = "A\n1\nNA\n3"
    converters = {"A": lambda x: x}

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data), converters=converters, keep_default_na=keep_default_na
            )
        return

    result = parser.read_csv(
        StringIO(data), converters=converters, keep_default_na=keep_default_na
    )
    values = ["1", np.nan, "3"] if keep_default_na else ["1", "NA", "3"]
    tm.assert_frame_equal(result, pd.DataFrame({"A": values}))


def test_converters_na_filter_false(all_parsers):
    # GH#13302 na_filter=False turns off na_values on a converter column too
    parser = all_parsers
    data = "A\n1\nCAT\n3"
    converters = {"A": lambda x: x}

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data), converters=converters, na_values="CAT", na_filter=False
            )
        return

    result = parser.read_csv(
        StringIO(data), converters=converters, na_values="CAT", na_filter=False
    )
    tm.assert_frame_equal(result, pd.DataFrame({"A": ["1", "CAT", "3"]}))


def test_converters_na_values_numeric_output(all_parsers):
    # GH#13302 na_values are matched against the converter's output, so the
    # column still infers a numeric dtype
    parser = all_parsers
    data = "A\n1\n-999\n3"
    converters = {"A": int}

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), converters=converters, na_values=[-999])
        return

    result = parser.read_csv(StringIO(data), converters=converters, na_values=[-999])
    tm.assert_frame_equal(result, pd.DataFrame({"A": [1.0, np.nan, 3.0]}))


def test_converters_na_values_index_col(all_parsers):
    # GH#13302 a converter column used as the index honors na_values too.
    # The converter is an identity so that this passes on the python engine,
    # which never calls a converter for an index column.
    parser = all_parsers
    data = "A,B\nq,x\nCAT,y"
    converters = {"A": lambda x: x}

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(
                StringIO(data), converters=converters, na_values=["CAT"], index_col=0
            )
        return

    result = parser.read_csv(
        StringIO(data), converters=converters, na_values=["CAT"], index_col=0
    )
    expected = pd.DataFrame({"B": ["x", "y"]}, index=pd.Index(["q", np.nan], name="A"))
    tm.assert_frame_equal(result, expected)


def test_converters_na_values_bool_output(all_parsers):
    # GH#13302 True == 1, so bools a converter returns match a numeric
    # na_values entry, the same cells the read without a converter drops.
    parser = all_parsers
    data = "A\n1\n0\n1"
    converters = {"A": lambda x: bool(int(x))}

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), converters=converters, na_values=[1])
        return

    result = parser.read_csv(StringIO(data), converters=converters, na_values=[1])
    if parser.engine == "python":
        # map_infer converts to a real bool array first, so _infer_types
        # takes the numeric path and the match is lost
        expected = pd.DataFrame({"A": [True, False, True]})
    else:
        expected = pd.DataFrame({"A": [np.nan, False, np.nan]}, dtype=object)
    tm.assert_frame_equal(result, expected)


def test_converters_raising_on_empty_field(all_parsers):
    # GH#13302 na_values apply to the converter's output, so a converter that
    # raises on an empty field still raises; the python engine masks it.
    parser = all_parsers
    data = "a,b\n1,x\n,y\n3,z"

    if parser.engine == "pyarrow":
        msg = "The 'converters' option is not supported with the 'pyarrow' engine"
        with pytest.raises(ValueError, match=msg):
            parser.read_csv(StringIO(data), converters={"a": float})
        return

    if parser.engine == "python":
        # the empty field is masked before the converter sees it, leaving the
        # column object dtype
        result = parser.read_csv(StringIO(data), converters={"a": float})
        expected = pd.DataFrame({"a": [1.0, np.nan, 3.0], "b": ["x", "y", "z"]})
        expected["a"] = expected["a"].astype(object)
        tm.assert_frame_equal(result, expected)
        return

    with pytest.raises(ValueError, match="could not convert string to float"):
        parser.read_csv(StringIO(data), converters={"a": float})
