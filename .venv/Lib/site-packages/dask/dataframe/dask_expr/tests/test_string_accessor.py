from __future__ import annotations

import pytest

from dask.dataframe.dask_expr._collection import DataFrame, from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq
from dask.dataframe.utils import pyarrow_strings_enabled

pd = _backend_library()


@pytest.fixture()
def ser():
    return pd.Series(["a", "b", "1", "aaa", "bbb", "ccc", "ddd", "abcd"])


@pytest.fixture()
def dser(ser):
    import dask.dataframe as dd

    return dd.from_pandas(ser, npartitions=3)


@pytest.mark.parametrize(
    "func, kwargs",
    [
        ("len", {}),
        ("capitalize", {}),
        ("casefold", {}),
        ("contains", {"pat": "a"}),
        ("count", {"pat": "a"}),
        ("endswith", {"pat": "a"}),
        ("extract", {"pat": r"[ab](\d)"}),
        ("extractall", {"pat": r"[ab](\d)"}),
        ("find", {"sub": "a"}),
        ("findall", {"pat": "a"}),
        ("fullmatch", {"pat": "a"}),
        ("get", {"i": 0}),
        ("isalnum", {}),
        ("isalpha", {}),
        ("isdecimal", {}),
        ("isdigit", {}),
        ("islower", {}),
        ("isspace", {}),
        ("istitle", {}),
        ("isupper", {}),
        ("join", {"sep": "-"}),
        ("len", {}),
        ("ljust", {"width": 3}),
        ("lower", {}),
        ("lstrip", {}),
        ("match", {"pat": r"[ab](\d)"}),
        ("normalize", {"form": "NFC"}),
        ("pad", {"width": 3}),
        ("removeprefix", {"prefix": "a"}),
        ("removesuffix", {"suffix": "a"}),
        ("repeat", {"repeats": 2}),
        ("replace", {"pat": "a", "repl": "b"}),
        ("rfind", {"sub": "a"}),
        ("rjust", {"width": 3}),
        ("rstrip", {}),
        ("slice", {"start": 0, "stop": 1}),
        ("slice_replace", {"start": 0, "stop": 1, "repl": "a"}),
        ("startswith", {"pat": "a"}),
        ("strip", {}),
        ("swapcase", {}),
        ("title", {}),
        ("upper", {}),
        ("wrap", {"width": 2}),
        ("zfill", {"width": 2}),
        ("split", {"pat": "a"}),
        ("rsplit", {"pat": "a"}),
        ("cat", {}),
        ("cat", {"others": pd.Series(["a"])}),
    ],
)
def test_string_accessor(ser, dser, func, kwargs):
    if pyarrow_strings_enabled():
        ser = ser.astype("string[pyarrow]")

    assert_eq(getattr(ser.str, func)(**kwargs), getattr(dser.str, func)(**kwargs))

    if func in (
        "contains",
        "endswith",
        "fullmatch",
        "isalnum",
        "isalpha",
        "isdecimal",
        "isdigit",
        "islower",
        "isspace",
        "istitle",
        "isupper",
        "startswith",
        "match",
    ):
        # This returns arrays and doesn't work in dask/dask either
        return

    ser.index = ser.values
    ser = ser.sort_index()
    dser = from_pandas(ser, npartitions=3)
    pdf_result = getattr(ser.index.str, func)(**kwargs)

    if func == "cat" and len(kwargs) > 0:
        # Doesn't work with others on Index
        return
    if isinstance(pdf_result, pd.DataFrame):
        assert_eq(
            getattr(dser.index.str, func)(**kwargs), pdf_result, check_index=False
        )
    else:
        assert_eq(getattr(dser.index.str, func)(**kwargs), pdf_result)


def test_str_accessor_cat(ser, dser):
    sol = ser.str.cat(ser.str.upper(), sep=":")
    assert_eq(dser.str.cat(dser.str.upper(), sep=":"), sol)
    assert_eq(dser.str.cat(ser.str.upper(), sep=":"), sol)
    assert_eq(
        dser.str.cat([dser.str.upper(), ser.str.lower()], sep=":"),
        ser.str.cat([ser.str.upper(), ser.str.lower()], sep=":"),
    )
    assert_eq(dser.str.cat(sep=":"), ser.str.cat(sep=":"))

    for o in ["foo", ["foo"]]:
        with pytest.raises(TypeError):
            dser.str.cat(o)


@pytest.mark.parametrize("index", [None, [0]], ids=["range_index", "other index"])
def test_str_split_(index):
    df = pd.DataFrame({"a": ["a\nb"]}, index=index)
    ddf = from_pandas(df, npartitions=1)

    pd_a = df["a"].str.split("\n", n=1, expand=True)
    dd_a = ddf["a"].str.split("\n", n=1, expand=True)

    assert_eq(dd_a, pd_a)


def test_str_accessor_not_available():
    pdf = pd.DataFrame({"a": [1, 2, 3]})
    df = from_pandas(pdf, npartitions=2)
    # Not available on invalid dtypes
    with pytest.raises(AttributeError, match=".str accessor"):
        df.a.str

    assert "str" not in dir(df.a)


def test_partition():
    df = DataFrame.from_dict({"A": ["A|B", "C|D"]}, npartitions=2)["A"].str.partition(
        "|"
    )
    result = df[1]
    expected = pd.DataFrame.from_dict({"A": ["A|B", "C|D"]})["A"].str.partition("|")[1]
    assert_eq(result, expected)
