from datetime import datetime, timedelta

import numpy as np
import pytest

from pandas._libs import lib

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series, isna, notna
import pandas._testing as tm
import pandas.core.strings as strings


def assert_series_or_index_equal(left, right):
    if isinstance(left, Series):
        tm.assert_series_equal(left, right)
    else:  # Index
        tm.assert_index_equal(left, right)


_any_string_method = [
    ("cat", (), {"sep": ","}),
    ("cat", (Series(list("zyx")),), {"sep": ",", "join": "left"}),
    ("center", (10,), {}),
    ("contains", ("a",), {}),
    ("count", ("a",), {}),
    ("decode", ("UTF-8",), {}),
    ("encode", ("UTF-8",), {}),
    ("endswith", ("a",), {}),
    ("endswith", ("a",), {"na": True}),
    ("endswith", ("a",), {"na": False}),
    ("extract", ("([a-z]*)",), {"expand": False}),
    ("extract", ("([a-z]*)",), {"expand": True}),
    ("extractall", ("([a-z]*)",), {}),
    ("find", ("a",), {}),
    ("findall", ("a",), {}),
    ("get", (0,), {}),
    # because "index" (and "rindex") fail intentionally
    # if the string is not found, search only for empty string
    ("index", ("",), {}),
    ("join", (",",), {}),
    ("ljust", (10,), {}),
    ("match", ("a",), {}),
    ("fullmatch", ("a",), {}),
    ("normalize", ("NFC",), {}),
    ("pad", (10,), {}),
    ("partition", (" ",), {"expand": False}),
    ("partition", (" ",), {"expand": True}),
    ("repeat", (3,), {}),
    ("replace", ("a", "z"), {}),
    ("rfind", ("a",), {}),
    ("rindex", ("",), {}),
    ("rjust", (10,), {}),
    ("rpartition", (" ",), {"expand": False}),
    ("rpartition", (" ",), {"expand": True}),
    ("slice", (0, 1), {}),
    ("slice_replace", (0, 1, "z"), {}),
    ("split", (" ",), {"expand": False}),
    ("split", (" ",), {"expand": True}),
    ("startswith", ("a",), {}),
    ("startswith", ("a",), {"na": True}),
    ("startswith", ("a",), {"na": False}),
    # translating unicode points of "a" to "d"
    ("translate", ({97: 100},), {}),
    ("wrap", (2,), {}),
    ("zfill", (10,), {}),
] + list(
    zip(
        [
            # methods without positional arguments: zip with empty tuple and empty dict
            "capitalize",
            "cat",
            "get_dummies",
            "isalnum",
            "isalpha",
            "isdecimal",
            "isdigit",
            "islower",
            "isnumeric",
            "isspace",
            "istitle",
            "isupper",
            "len",
            "lower",
            "lstrip",
            "partition",
            "rpartition",
            "rsplit",
            "rstrip",
            "slice",
            "slice_replace",
            "split",
            "strip",
            "swapcase",
            "title",
            "upper",
            "casefold",
        ],
        [()] * 100,
        [{}] * 100,
    )
)
ids, _, _ = zip(*_any_string_method)  # use method name as fixture-id


# test that the above list captures all methods of StringMethods
missing_methods = {
    f for f in dir(strings.StringMethods) if not f.startswith("_")
} - set(ids)
assert not missing_methods


@pytest.fixture(params=_any_string_method, ids=ids)
def any_string_method(request):
    """
    Fixture for all public methods of `StringMethods`

    This fixture returns a tuple of the method name and sample arguments
    necessary to call the method.

    Returns
    -------
    method_name : str
        The name of the method in `StringMethods`
    args : tuple
        Sample values for the positional arguments
    kwargs : dict
        Sample values for the keyword arguments

    Examples
    --------
    >>> def test_something(any_string_method):
    ...     s = Series(['a', 'b', np.nan, 'd'])
    ...
    ...     method_name, args, kwargs = any_string_method
    ...     method = getattr(s.str, method_name)
    ...     # will not raise
    ...     method(*args, **kwargs)
    """
    return request.param


# subset of the full set from pandas/conftest.py
_any_allowed_skipna_inferred_dtype = [
    ("string", ["a", np.nan, "c"]),
    ("bytes", [b"a", np.nan, b"c"]),
    ("empty", [np.nan, np.nan, np.nan]),
    ("empty", []),
    ("mixed-integer", ["a", np.nan, 2]),
]
ids, _ = zip(*_any_allowed_skipna_inferred_dtype)  # use inferred type as id


@pytest.fixture(params=_any_allowed_skipna_inferred_dtype, ids=ids)
def any_allowed_skipna_inferred_dtype(request):
    """
    Fixture for all (inferred) dtypes allowed in StringMethods.__init__

    The covered (inferred) types are:
    * 'string'
    * 'empty'
    * 'bytes'
    * 'mixed'
    * 'mixed-integer'

    Returns
    -------
    inferred_dtype : str
        The string for the inferred dtype from _libs.lib.infer_dtype
    values : np.ndarray
        An array of object dtype that will be inferred to have
        `inferred_dtype`

    Examples
    --------
    >>> import pandas._libs.lib as lib
    >>>
    >>> def test_something(any_allowed_skipna_inferred_dtype):
    ...     inferred_dtype, values = any_allowed_skipna_inferred_dtype
    ...     # will pass
    ...     assert lib.infer_dtype(values, skipna=True) == inferred_dtype
    ...
    ...     # constructor for .str-accessor will also pass
    ...     Series(values).str
    """
    inferred_dtype, values = request.param
    values = np.array(values, dtype=object)  # object dtype to avoid casting

    # correctness of inference tested in tests/dtypes/test_inference.py
    return inferred_dtype, values


def test_iter():
    # GH3638
    strs = "google", "wikimedia", "wikipedia", "wikitravel"
    ds = Series(strs)

    with tm.assert_produces_warning(FutureWarning):
        for s in ds.str:
            # iter must yield a Series
            assert isinstance(s, Series)

            # indices of each yielded Series should be equal to the index of
            # the original Series
            tm.assert_index_equal(s.index, ds.index)

            for el in s:
                # each element of the series is either a basestring/str or nan
                assert isinstance(el, str) or isna(el)

    # desired behavior is to iterate until everything would be nan on the
    # next iter so make sure the last element of the iterator was 'l' in
    # this case since 'wikitravel' is the longest string
    assert s.dropna().values.item() == "l"


def test_iter_empty():
    ds = Series([], dtype=object)

    i, s = 100, 1

    with tm.assert_produces_warning(FutureWarning):
        for i, s in enumerate(ds.str):
            pass

    # nothing to iterate over so nothing defined values should remain
    # unchanged
    assert i == 100
    assert s == 1


def test_iter_single_element():
    ds = Series(["a"])

    with tm.assert_produces_warning(FutureWarning):
        for i, s in enumerate(ds.str):
            pass

    assert not i
    tm.assert_series_equal(ds, s)


def test_iter_object_try_string():
    ds = Series(
        [
            slice(None, np.random.randint(10), np.random.randint(10, 20))
            for _ in range(4)
        ]
    )

    i, s = 100, "h"

    with tm.assert_produces_warning(FutureWarning):
        for i, s in enumerate(ds.str):
            pass

    assert i == 100
    assert s == "h"


# test integer/float dtypes (inferred by constructor) and mixed


def test_count():
    values = np.array(["foo", "foofoo", np.nan, "foooofooofommmfoo"], dtype=np.object_)

    result = Series(values).str.count("f[o]+")
    exp = Series([1, 2, np.nan, 4])
    assert isinstance(result, Series)
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = np.array(
        ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
        dtype=object,
    )
    rs = Series(mixed).str.count("a")
    xp = Series([1, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_series_equal(rs, xp)


def test_title():
    values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])

    result = values.str.title()
    exp = Series(["Foo", "Bar", np.nan, "Blah", "Blurg"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["FOO", np.nan, "bar", True, datetime.today(), "blah", None, 1, 2.0])
    mixed = mixed.str.title()
    exp = Series(["Foo", np.nan, "Bar", np.nan, np.nan, "Blah", np.nan, np.nan, np.nan])
    tm.assert_almost_equal(mixed, exp)


def test_lower_upper():
    values = Series(["om", np.nan, "nom", "nom"])

    result = values.str.upper()
    exp = Series(["OM", np.nan, "NOM", "NOM"])
    tm.assert_series_equal(result, exp)

    result = result.str.lower()
    tm.assert_series_equal(result, values)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0])
    mixed = mixed.str.upper()
    rs = Series(mixed).str.lower()
    xp = Series(["a", np.nan, "b", np.nan, np.nan, "foo", np.nan, np.nan, np.nan])
    assert isinstance(rs, Series)
    tm.assert_series_equal(rs, xp)


def test_capitalize():
    values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])
    result = values.str.capitalize()
    exp = Series(["Foo", "Bar", np.nan, "Blah", "Blurg"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["FOO", np.nan, "bar", True, datetime.today(), "blah", None, 1, 2.0])
    mixed = mixed.str.capitalize()
    exp = Series(["Foo", np.nan, "Bar", np.nan, np.nan, "Blah", np.nan, np.nan, np.nan])
    tm.assert_almost_equal(mixed, exp)


def test_swapcase():
    values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])
    result = values.str.swapcase()
    exp = Series(["foo", "bar", np.nan, "bLAH", "BLURG"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["FOO", np.nan, "bar", True, datetime.today(), "Blah", None, 1, 2.0])
    mixed = mixed.str.swapcase()
    exp = Series(["foo", np.nan, "BAR", np.nan, np.nan, "bLAH", np.nan, np.nan, np.nan])
    tm.assert_almost_equal(mixed, exp)


def test_casemethods():
    values = ["aaa", "bbb", "CCC", "Dddd", "eEEE"]
    s = Series(values)
    assert s.str.lower().tolist() == [v.lower() for v in values]
    assert s.str.upper().tolist() == [v.upper() for v in values]
    assert s.str.title().tolist() == [v.title() for v in values]
    assert s.str.capitalize().tolist() == [v.capitalize() for v in values]
    assert s.str.swapcase().tolist() == [v.swapcase() for v in values]


def test_repeat():
    values = Series(["a", "b", np.nan, "c", np.nan, "d"])

    result = values.str.repeat(3)
    exp = Series(["aaa", "bbb", np.nan, "ccc", np.nan, "ddd"])
    tm.assert_series_equal(result, exp)

    result = values.str.repeat([1, 2, 3, 4, 5, 6])
    exp = Series(["a", "bb", np.nan, "cccc", np.nan, "dddddd"])
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0])

    rs = Series(mixed).str.repeat(3)
    xp = Series(
        ["aaa", np.nan, "bbb", np.nan, np.nan, "foofoofoo", np.nan, np.nan, np.nan]
    )
    assert isinstance(rs, Series)
    tm.assert_series_equal(rs, xp)


def test_repeat_with_null():
    # GH: 31632
    values = Series(["a", None], dtype="string")
    result = values.str.repeat([3, 4])
    exp = Series(["aaa", None], dtype="string")
    tm.assert_series_equal(result, exp)

    values = Series(["a", "b"], dtype="string")
    result = values.str.repeat([3, None])
    exp = Series(["aaa", None], dtype="string")
    tm.assert_series_equal(result, exp)


def test_empty_str_methods():
    empty_str = empty = Series(dtype=object)
    empty_int = Series(dtype="int64")
    empty_bool = Series(dtype=bool)
    empty_bytes = Series(dtype=object)

    # GH7241
    # (extract) on empty series

    tm.assert_series_equal(empty_str, empty.str.cat(empty))
    assert "" == empty.str.cat()
    tm.assert_series_equal(empty_str, empty.str.title())
    tm.assert_series_equal(empty_int, empty.str.count("a"))
    tm.assert_series_equal(empty_bool, empty.str.contains("a"))
    tm.assert_series_equal(empty_bool, empty.str.startswith("a"))
    tm.assert_series_equal(empty_bool, empty.str.endswith("a"))
    tm.assert_series_equal(empty_str, empty.str.lower())
    tm.assert_series_equal(empty_str, empty.str.upper())
    tm.assert_series_equal(empty_str, empty.str.replace("a", "b"))
    tm.assert_series_equal(empty_str, empty.str.repeat(3))
    tm.assert_series_equal(empty_bool, empty.str.match("^a"))
    tm.assert_frame_equal(
        DataFrame(columns=[0], dtype=str), empty.str.extract("()", expand=True)
    )
    tm.assert_frame_equal(
        DataFrame(columns=[0, 1], dtype=str), empty.str.extract("()()", expand=True)
    )
    tm.assert_series_equal(empty_str, empty.str.extract("()", expand=False))
    tm.assert_frame_equal(
        DataFrame(columns=[0, 1], dtype=str),
        empty.str.extract("()()", expand=False),
    )
    tm.assert_frame_equal(DataFrame(dtype=str), empty.str.get_dummies())
    tm.assert_series_equal(empty_str, empty_str.str.join(""))
    tm.assert_series_equal(empty_int, empty.str.len())
    tm.assert_series_equal(empty_str, empty_str.str.findall("a"))
    tm.assert_series_equal(empty_int, empty.str.find("a"))
    tm.assert_series_equal(empty_int, empty.str.rfind("a"))
    tm.assert_series_equal(empty_str, empty.str.pad(42))
    tm.assert_series_equal(empty_str, empty.str.center(42))
    tm.assert_series_equal(empty_str, empty.str.split("a"))
    tm.assert_series_equal(empty_str, empty.str.rsplit("a"))
    tm.assert_series_equal(empty_str, empty.str.partition("a", expand=False))
    tm.assert_series_equal(empty_str, empty.str.rpartition("a", expand=False))
    tm.assert_series_equal(empty_str, empty.str.slice(stop=1))
    tm.assert_series_equal(empty_str, empty.str.slice(step=1))
    tm.assert_series_equal(empty_str, empty.str.strip())
    tm.assert_series_equal(empty_str, empty.str.lstrip())
    tm.assert_series_equal(empty_str, empty.str.rstrip())
    tm.assert_series_equal(empty_str, empty.str.wrap(42))
    tm.assert_series_equal(empty_str, empty.str.get(0))
    tm.assert_series_equal(empty_str, empty_bytes.str.decode("ascii"))
    tm.assert_series_equal(empty_bytes, empty.str.encode("ascii"))
    # ismethods should always return boolean (GH 29624)
    tm.assert_series_equal(empty_bool, empty.str.isalnum())
    tm.assert_series_equal(empty_bool, empty.str.isalpha())
    tm.assert_series_equal(empty_bool, empty.str.isdigit())
    tm.assert_series_equal(empty_bool, empty.str.isspace())
    tm.assert_series_equal(empty_bool, empty.str.islower())
    tm.assert_series_equal(empty_bool, empty.str.isupper())
    tm.assert_series_equal(empty_bool, empty.str.istitle())
    tm.assert_series_equal(empty_bool, empty.str.isnumeric())
    tm.assert_series_equal(empty_bool, empty.str.isdecimal())
    tm.assert_series_equal(empty_str, empty.str.capitalize())
    tm.assert_series_equal(empty_str, empty.str.swapcase())
    tm.assert_series_equal(empty_str, empty.str.normalize("NFC"))

    table = str.maketrans("a", "b")
    tm.assert_series_equal(empty_str, empty.str.translate(table))


def test_empty_str_methods_to_frame():
    empty = Series(dtype=str)
    empty_df = DataFrame()
    tm.assert_frame_equal(empty_df, empty.str.partition("a"))
    tm.assert_frame_equal(empty_df, empty.str.rpartition("a"))


def test_ismethods():
    values = ["A", "b", "Xy", "4", "3A", "", "TT", "55", "-", "  "]
    str_s = Series(values)
    alnum_e = [True, True, True, True, True, False, True, True, False, False]
    alpha_e = [True, True, True, False, False, False, True, False, False, False]
    digit_e = [False, False, False, True, False, False, False, True, False, False]

    # TODO: unused
    num_e = [  # noqa
        False,
        False,
        False,
        True,
        False,
        False,
        False,
        True,
        False,
        False,
    ]

    space_e = [False, False, False, False, False, False, False, False, False, True]
    lower_e = [False, True, False, False, False, False, False, False, False, False]
    upper_e = [True, False, False, False, True, False, True, False, False, False]
    title_e = [True, False, True, False, True, False, False, False, False, False]

    tm.assert_series_equal(str_s.str.isalnum(), Series(alnum_e))
    tm.assert_series_equal(str_s.str.isalpha(), Series(alpha_e))
    tm.assert_series_equal(str_s.str.isdigit(), Series(digit_e))
    tm.assert_series_equal(str_s.str.isspace(), Series(space_e))
    tm.assert_series_equal(str_s.str.islower(), Series(lower_e))
    tm.assert_series_equal(str_s.str.isupper(), Series(upper_e))
    tm.assert_series_equal(str_s.str.istitle(), Series(title_e))

    assert str_s.str.isalnum().tolist() == [v.isalnum() for v in values]
    assert str_s.str.isalpha().tolist() == [v.isalpha() for v in values]
    assert str_s.str.isdigit().tolist() == [v.isdigit() for v in values]
    assert str_s.str.isspace().tolist() == [v.isspace() for v in values]
    assert str_s.str.islower().tolist() == [v.islower() for v in values]
    assert str_s.str.isupper().tolist() == [v.isupper() for v in values]
    assert str_s.str.istitle().tolist() == [v.istitle() for v in values]


def test_isnumeric():
    # 0x00bc: ¼ VULGAR FRACTION ONE QUARTER
    # 0x2605: ★ not number
    # 0x1378: ፸ ETHIOPIC NUMBER SEVENTY
    # 0xFF13: ３ Em 3
    values = ["A", "3", "¼", "★", "፸", "３", "four"]
    s = Series(values)
    numeric_e = [False, True, True, False, True, True, False]
    decimal_e = [False, True, False, False, False, True, False]
    tm.assert_series_equal(s.str.isnumeric(), Series(numeric_e))
    tm.assert_series_equal(s.str.isdecimal(), Series(decimal_e))

    unicodes = ["A", "3", "¼", "★", "፸", "３", "four"]
    assert s.str.isnumeric().tolist() == [v.isnumeric() for v in unicodes]
    assert s.str.isdecimal().tolist() == [v.isdecimal() for v in unicodes]

    values = ["A", np.nan, "¼", "★", np.nan, "３", "four"]
    s = Series(values)
    numeric_e = [False, np.nan, True, False, np.nan, True, False]
    decimal_e = [False, np.nan, False, False, np.nan, True, False]
    tm.assert_series_equal(s.str.isnumeric(), Series(numeric_e))
    tm.assert_series_equal(s.str.isdecimal(), Series(decimal_e))


def test_get_dummies():
    s = Series(["a|b", "a|c", np.nan])
    result = s.str.get_dummies("|")
    expected = DataFrame([[1, 1, 0], [1, 0, 1], [0, 0, 0]], columns=list("abc"))
    tm.assert_frame_equal(result, expected)

    s = Series(["a;b", "a", 7])
    result = s.str.get_dummies(";")
    expected = DataFrame([[0, 1, 1], [0, 1, 0], [1, 0, 0]], columns=list("7ab"))
    tm.assert_frame_equal(result, expected)

    # GH9980, GH8028
    idx = Index(["a|b", "a|c", "b|c"])
    result = idx.str.get_dummies("|")

    expected = MultiIndex.from_tuples(
        [(1, 1, 0), (1, 0, 1), (0, 1, 1)], names=("a", "b", "c")
    )
    tm.assert_index_equal(result, expected)


def test_get_dummies_with_name_dummy():
    # GH 12180
    # Dummies named 'name' should work as expected
    s = Series(["a", "b,name", "b"])
    result = s.str.get_dummies(",")
    expected = DataFrame([[1, 0, 0], [0, 1, 1], [0, 1, 0]], columns=["a", "b", "name"])
    tm.assert_frame_equal(result, expected)

    idx = Index(["a|b", "name|c", "b|name"])
    result = idx.str.get_dummies("|")

    expected = MultiIndex.from_tuples(
        [(1, 1, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1)], names=("a", "b", "c", "name")
    )
    tm.assert_index_equal(result, expected)


def test_join():
    values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])
    result = values.str.split("_").str.join("_")
    tm.assert_series_equal(values, result)

    # mixed
    mixed = Series(
        [
            "a_b",
            np.nan,
            "asdf_cas_asdf",
            True,
            datetime.today(),
            "foo",
            None,
            1,
            2.0,
        ]
    )

    rs = Series(mixed).str.split("_").str.join("_")
    xp = Series(
        [
            "a_b",
            np.nan,
            "asdf_cas_asdf",
            np.nan,
            np.nan,
            "foo",
            np.nan,
            np.nan,
            np.nan,
        ]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_len():
    values = Series(["foo", "fooo", "fooooo", np.nan, "fooooooo"])

    result = values.str.len()
    exp = values.map(lambda x: len(x) if notna(x) else np.nan)
    tm.assert_series_equal(result, exp)

    # mixed
    mixed = Series(
        [
            "a_b",
            np.nan,
            "asdf_cas_asdf",
            True,
            datetime.today(),
            "foo",
            None,
            1,
            2.0,
        ]
    )

    rs = Series(mixed).str.len()
    xp = Series([3, np.nan, 13, np.nan, np.nan, 3, np.nan, np.nan, np.nan])

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_index():
    def _check(result, expected):
        if isinstance(result, Series):
            tm.assert_series_equal(result, expected)
        else:
            tm.assert_index_equal(result, expected)

    for klass in [Series, Index]:
        s = klass(["ABCDEFG", "BCDEFEF", "DEFGHIJEF", "EFGHEF"])

        result = s.str.index("EF")
        _check(result, klass([4, 3, 1, 0]))
        expected = np.array([v.index("EF") for v in s.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = s.str.rindex("EF")
        _check(result, klass([4, 5, 7, 4]))
        expected = np.array([v.rindex("EF") for v in s.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = s.str.index("EF", 3)
        _check(result, klass([4, 3, 7, 4]))
        expected = np.array([v.index("EF", 3) for v in s.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = s.str.rindex("EF", 3)
        _check(result, klass([4, 5, 7, 4]))
        expected = np.array([v.rindex("EF", 3) for v in s.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = s.str.index("E", 4, 8)
        _check(result, klass([4, 5, 7, 4]))
        expected = np.array([v.index("E", 4, 8) for v in s.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = s.str.rindex("E", 0, 5)
        _check(result, klass([4, 3, 1, 4]))
        expected = np.array([v.rindex("E", 0, 5) for v in s.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        with pytest.raises(ValueError, match="substring not found"):
            result = s.str.index("DE")

        msg = "expected a string object, not int"
        with pytest.raises(TypeError, match=msg):
            result = s.str.index(0)

        with pytest.raises(TypeError, match=msg):
            result = s.str.rindex(0)

    # test with nan
    s = Series(["abcb", "ab", "bcbe", np.nan])
    result = s.str.index("b")
    tm.assert_series_equal(result, Series([1, 1, 0, np.nan]))
    result = s.str.rindex("b")
    tm.assert_series_equal(result, Series([3, 1, 2, np.nan]))


def test_pad():
    values = Series(["a", "b", np.nan, "c", np.nan, "eeeeee"])

    result = values.str.pad(5, side="left")
    exp = Series(["    a", "    b", np.nan, "    c", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="right")
    exp = Series(["a    ", "b    ", np.nan, "c    ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="both")
    exp = Series(["  a  ", "  b  ", np.nan, "  c  ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "ee", None, 1, 2.0])

    rs = Series(mixed).str.pad(5, side="left")
    xp = Series(
        ["    a", np.nan, "    b", np.nan, np.nan, "   ee", np.nan, np.nan, np.nan]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    mixed = Series(["a", np.nan, "b", True, datetime.today(), "ee", None, 1, 2.0])

    rs = Series(mixed).str.pad(5, side="right")
    xp = Series(
        ["a    ", np.nan, "b    ", np.nan, np.nan, "ee   ", np.nan, np.nan, np.nan]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    mixed = Series(["a", np.nan, "b", True, datetime.today(), "ee", None, 1, 2.0])

    rs = Series(mixed).str.pad(5, side="both")
    xp = Series(
        ["  a  ", np.nan, "  b  ", np.nan, np.nan, "  ee ", np.nan, np.nan, np.nan]
    )

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_pad_fillchar():

    values = Series(["a", "b", np.nan, "c", np.nan, "eeeeee"])

    result = values.str.pad(5, side="left", fillchar="X")
    exp = Series(["XXXXa", "XXXXb", np.nan, "XXXXc", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="right", fillchar="X")
    exp = Series(["aXXXX", "bXXXX", np.nan, "cXXXX", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.pad(5, side="both", fillchar="X")
    exp = Series(["XXaXX", "XXbXX", np.nan, "XXcXX", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    msg = "fillchar must be a character, not str"
    with pytest.raises(TypeError, match=msg):
        result = values.str.pad(5, fillchar="XY")

    msg = "fillchar must be a character, not int"
    with pytest.raises(TypeError, match=msg):
        result = values.str.pad(5, fillchar=5)


@pytest.mark.parametrize("f", ["center", "ljust", "rjust", "zfill", "pad"])
def test_pad_width(f):
    # see gh-13598
    s = Series(["1", "22", "a", "bb"])
    msg = "width must be of integer type, not*"

    with pytest.raises(TypeError, match=msg):
        getattr(s.str, f)("f")


def test_center_ljust_rjust():
    values = Series(["a", "b", np.nan, "c", np.nan, "eeeeee"])

    result = values.str.center(5)
    exp = Series(["  a  ", "  b  ", np.nan, "  c  ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.ljust(5)
    exp = Series(["a    ", "b    ", np.nan, "c    ", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    result = values.str.rjust(5)
    exp = Series(["    a", "    b", np.nan, "    c", np.nan, "eeeeee"])
    tm.assert_almost_equal(result, exp)

    # mixed
    mixed = Series(["a", np.nan, "b", True, datetime.today(), "c", "eee", None, 1, 2.0])

    rs = Series(mixed).str.center(5)
    xp = Series(
        [
            "  a  ",
            np.nan,
            "  b  ",
            np.nan,
            np.nan,
            "  c  ",
            " eee ",
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.ljust(5)
    xp = Series(
        [
            "a    ",
            np.nan,
            "b    ",
            np.nan,
            np.nan,
            "c    ",
            "eee  ",
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.rjust(5)
    xp = Series(
        [
            "    a",
            np.nan,
            "    b",
            np.nan,
            np.nan,
            "    c",
            "  eee",
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_center_ljust_rjust_fillchar():
    values = Series(["a", "bb", "cccc", "ddddd", "eeeeee"])

    result = values.str.center(5, fillchar="X")
    expected = Series(["XXaXX", "XXbbX", "Xcccc", "ddddd", "eeeeee"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.center(5, "X") for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.ljust(5, fillchar="X")
    expected = Series(["aXXXX", "bbXXX", "ccccX", "ddddd", "eeeeee"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.ljust(5, "X") for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.rjust(5, fillchar="X")
    expected = Series(["XXXXa", "XXXbb", "Xcccc", "ddddd", "eeeeee"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.rjust(5, "X") for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    # If fillchar is not a charatter, normal str raises TypeError
    # 'aaa'.ljust(5, 'XY')
    # TypeError: must be char, not str
    template = "fillchar must be a character, not {dtype}"

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.center(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.ljust(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="str")):
        values.str.rjust(5, fillchar="XY")

    with pytest.raises(TypeError, match=template.format(dtype="int")):
        values.str.center(5, fillchar=1)

    with pytest.raises(TypeError, match=template.format(dtype="int")):
        values.str.ljust(5, fillchar=1)

    with pytest.raises(TypeError, match=template.format(dtype="int")):
        values.str.rjust(5, fillchar=1)


def test_zfill():
    values = Series(["1", "22", "aaa", "333", "45678"])

    result = values.str.zfill(5)
    expected = Series(["00001", "00022", "00aaa", "00333", "45678"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.zfill(5) for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    result = values.str.zfill(3)
    expected = Series(["001", "022", "aaa", "333", "45678"])
    tm.assert_series_equal(result, expected)
    expected = np.array([v.zfill(3) for v in values.values], dtype=np.object_)
    tm.assert_numpy_array_equal(result.values, expected)

    values = Series(["1", np.nan, "aaa", np.nan, "45678"])
    result = values.str.zfill(5)
    expected = Series(["00001", np.nan, "00aaa", np.nan, "45678"])
    tm.assert_series_equal(result, expected)


def test_pipe_failures():
    # #2119
    s = Series(["A|B|C"])

    result = s.str.split("|")
    exp = Series([["A", "B", "C"]])

    tm.assert_series_equal(result, exp)

    result = s.str.replace("|", " ", regex=False)
    exp = Series(["A B C"])

    tm.assert_series_equal(result, exp)


@pytest.mark.parametrize(
    "start, stop, step, expected",
    [
        (2, 5, None, Series(["foo", "bar", np.nan, "baz"])),
        (0, 3, -1, Series(["", "", np.nan, ""])),
        (None, None, -1, Series(["owtoofaa", "owtrabaa", np.nan, "xuqzabaa"])),
        (3, 10, 2, Series(["oto", "ato", np.nan, "aqx"])),
        (3, 0, -1, Series(["ofa", "aba", np.nan, "aba"])),
    ],
)
def test_slice(start, stop, step, expected):
    values = Series(["aafootwo", "aabartwo", np.nan, "aabazqux"])
    result = values.str.slice(start, stop, step)
    tm.assert_series_equal(result, expected)

    # mixed
    mixed = Series(
        ["aafootwo", np.nan, "aabartwo", True, datetime.today(), None, 1, 2.0]
    )

    rs = Series(mixed).str.slice(2, 5)
    xp = Series(["foo", np.nan, "bar", np.nan, np.nan, np.nan, np.nan, np.nan])

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.slice(2, 5, -1)
    xp = Series(["oof", np.nan, "rab", np.nan, np.nan, np.nan, np.nan, np.nan])


def test_slice_replace():
    values = Series(["short", "a bit longer", "evenlongerthanthat", "", np.nan])

    exp = Series(["shrt", "a it longer", "evnlongerthanthat", "", np.nan])
    result = values.str.slice_replace(2, 3)
    tm.assert_series_equal(result, exp)

    exp = Series(["shzrt", "a zit longer", "evznlongerthanthat", "z", np.nan])
    result = values.str.slice_replace(2, 3, "z")
    tm.assert_series_equal(result, exp)

    exp = Series(["shzort", "a zbit longer", "evzenlongerthanthat", "z", np.nan])
    result = values.str.slice_replace(2, 2, "z")
    tm.assert_series_equal(result, exp)

    exp = Series(["shzort", "a zbit longer", "evzenlongerthanthat", "z", np.nan])
    result = values.str.slice_replace(2, 1, "z")
    tm.assert_series_equal(result, exp)

    exp = Series(["shorz", "a bit longez", "evenlongerthanthaz", "z", np.nan])
    result = values.str.slice_replace(-1, None, "z")
    tm.assert_series_equal(result, exp)

    exp = Series(["zrt", "zer", "zat", "z", np.nan])
    result = values.str.slice_replace(None, -2, "z")
    tm.assert_series_equal(result, exp)

    exp = Series(["shortz", "a bit znger", "evenlozerthanthat", "z", np.nan])
    result = values.str.slice_replace(6, 8, "z")
    tm.assert_series_equal(result, exp)

    exp = Series(["zrt", "a zit longer", "evenlongzerthanthat", "z", np.nan])
    result = values.str.slice_replace(-10, 3, "z")
    tm.assert_series_equal(result, exp)


def test_strip_lstrip_rstrip():
    values = Series(["  aa   ", " bb \n", np.nan, "cc  "])

    result = values.str.strip()
    exp = Series(["aa", "bb", np.nan, "cc"])
    tm.assert_series_equal(result, exp)

    result = values.str.lstrip()
    exp = Series(["aa   ", "bb \n", np.nan, "cc  "])
    tm.assert_series_equal(result, exp)

    result = values.str.rstrip()
    exp = Series(["  aa", " bb", np.nan, "cc"])
    tm.assert_series_equal(result, exp)


def test_strip_lstrip_rstrip_mixed():
    # mixed
    mixed = Series(["  aa  ", np.nan, " bb \t\n", True, datetime.today(), None, 1, 2.0])

    rs = Series(mixed).str.strip()
    xp = Series(["aa", np.nan, "bb", np.nan, np.nan, np.nan, np.nan, np.nan])

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.lstrip()
    xp = Series(["aa  ", np.nan, "bb \t\n", np.nan, np.nan, np.nan, np.nan, np.nan])

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)

    rs = Series(mixed).str.rstrip()
    xp = Series(["  aa", np.nan, " bb", np.nan, np.nan, np.nan, np.nan, np.nan])

    assert isinstance(rs, Series)
    tm.assert_almost_equal(rs, xp)


def test_strip_lstrip_rstrip_args():
    values = Series(["xxABCxx", "xx BNSD", "LDFJH xx"])

    rs = values.str.strip("x")
    xp = Series(["ABC", " BNSD", "LDFJH "])
    tm.assert_series_equal(rs, xp)

    rs = values.str.lstrip("x")
    xp = Series(["ABCxx", " BNSD", "LDFJH xx"])
    tm.assert_series_equal(rs, xp)

    rs = values.str.rstrip("x")
    xp = Series(["xxABC", "xx BNSD", "LDFJH "])
    tm.assert_series_equal(rs, xp)


def test_wrap():
    # test values are: two words less than width, two words equal to width,
    # two words greater than width, one word less than width, one word
    # equal to width, one word greater than width, multiple tokens with
    # trailing whitespace equal to width
    values = Series(
        [
            "hello world",
            "hello world!",
            "hello world!!",
            "abcdefabcde",
            "abcdefabcdef",
            "abcdefabcdefa",
            "ab ab ab ab ",
            "ab ab ab ab a",
            "\t",
        ]
    )

    # expected values
    xp = Series(
        [
            "hello world",
            "hello world!",
            "hello\nworld!!",
            "abcdefabcde",
            "abcdefabcdef",
            "abcdefabcdef\na",
            "ab ab ab ab",
            "ab ab ab ab\na",
            "",
        ]
    )

    rs = values.str.wrap(12, break_long_words=True)
    tm.assert_series_equal(rs, xp)

    # test with pre and post whitespace (non-unicode), NaN, and non-ascii
    # Unicode
    values = Series(["  pre  ", np.nan, "\xac\u20ac\U00008000 abadcafe"])
    xp = Series(["  pre", np.nan, "\xac\u20ac\U00008000 ab\nadcafe"])
    rs = values.str.wrap(6)
    tm.assert_series_equal(rs, xp)


def test_string_slice_get_syntax():
    s = Series(
        [
            "YYY",
            "B",
            "C",
            "YYYYYYbYYY",
            "BYYYcYYY",
            np.nan,
            "CYYYBYYY",
            "dog",
            "cYYYt",
        ]
    )

    result = s.str[0]
    expected = s.str.get(0)
    tm.assert_series_equal(result, expected)

    result = s.str[:3]
    expected = s.str.slice(stop=3)
    tm.assert_series_equal(result, expected)

    result = s.str[2::-1]
    expected = s.str.slice(start=2, step=-1)
    tm.assert_series_equal(result, expected)


def test_string_slice_out_of_bounds():
    s = Series([(1, 2), (1,), (3, 4, 5)])

    result = s.str[1]
    expected = Series([2, np.nan, 4])

    tm.assert_series_equal(result, expected)

    s = Series(["foo", "b", "ba"])
    result = s.str[1]
    expected = Series(["o", np.nan, "a"])
    tm.assert_series_equal(result, expected)


def test_encode_decode():
    base = Series(["a", "b", "a\xe4"])
    series = base.str.encode("utf-8")

    f = lambda x: x.decode("utf-8")
    result = series.str.decode("utf-8")
    exp = series.map(f)

    tm.assert_series_equal(result, exp)


def test_encode_decode_errors():
    encodeBase = Series(["a", "b", "a\x9d"])

    msg = (
        r"'charmap' codec can't encode character '\\x9d' in position 1: "
        "character maps to <undefined>"
    )
    with pytest.raises(UnicodeEncodeError, match=msg):
        encodeBase.str.encode("cp1252")

    f = lambda x: x.encode("cp1252", "ignore")
    result = encodeBase.str.encode("cp1252", "ignore")
    exp = encodeBase.map(f)
    tm.assert_series_equal(result, exp)

    decodeBase = Series([b"a", b"b", b"a\x9d"])

    msg = (
        "'charmap' codec can't decode byte 0x9d in position 1: "
        "character maps to <undefined>"
    )
    with pytest.raises(UnicodeDecodeError, match=msg):
        decodeBase.str.decode("cp1252")

    f = lambda x: x.decode("cp1252", "ignore")
    result = decodeBase.str.decode("cp1252", "ignore")
    exp = decodeBase.map(f)

    tm.assert_series_equal(result, exp)


def test_normalize():
    values = ["ABC", "ＡＢＣ", "１２３", np.nan, "ｱｲｴ"]
    s = Series(values, index=["a", "b", "c", "d", "e"])

    normed = ["ABC", "ABC", "123", np.nan, "アイエ"]
    expected = Series(normed, index=["a", "b", "c", "d", "e"])

    result = s.str.normalize("NFKC")
    tm.assert_series_equal(result, expected)

    expected = Series(
        ["ABC", "ＡＢＣ", "１２３", np.nan, "ｱｲｴ"], index=["a", "b", "c", "d", "e"]
    )

    result = s.str.normalize("NFC")
    tm.assert_series_equal(result, expected)

    with pytest.raises(ValueError, match="invalid normalization form"):
        s.str.normalize("xxx")

    s = Index(["ＡＢＣ", "１２３", "ｱｲｴ"])
    expected = Index(["ABC", "123", "アイエ"])
    result = s.str.normalize("NFKC")
    tm.assert_index_equal(result, expected)


def test_index_str_accessor_visibility():
    from pandas.core.strings import StringMethods

    cases = [
        (["a", "b"], "string"),
        (["a", "b", 1], "mixed-integer"),
        (["a", "b", 1.3], "mixed"),
        (["a", "b", 1.3, 1], "mixed-integer"),
        (["aa", datetime(2011, 1, 1)], "mixed"),
    ]
    for values, tp in cases:
        idx = Index(values)
        assert isinstance(Series(values).str, StringMethods)
        assert isinstance(idx.str, StringMethods)
        assert idx.inferred_type == tp

    for values, tp in cases:
        idx = Index(values)
        assert isinstance(Series(values).str, StringMethods)
        assert isinstance(idx.str, StringMethods)
        assert idx.inferred_type == tp

    cases = [
        ([1, np.nan], "floating"),
        ([datetime(2011, 1, 1)], "datetime64"),
        ([timedelta(1)], "timedelta64"),
    ]
    for values, tp in cases:
        idx = Index(values)
        message = "Can only use .str accessor with string values"
        with pytest.raises(AttributeError, match=message):
            Series(values).str
        with pytest.raises(AttributeError, match=message):
            idx.str
        assert idx.inferred_type == tp

    # MultiIndex has mixed dtype, but not allow to use accessor
    idx = MultiIndex.from_tuples([("a", "b"), ("a", "b")])
    assert idx.inferred_type == "mixed"
    message = "Can only use .str accessor with Index, not MultiIndex"
    with pytest.raises(AttributeError, match=message):
        idx.str


def test_str_accessor_no_new_attributes():
    # https://github.com/pandas-dev/pandas/issues/10673
    s = Series(list("aabbcde"))
    with pytest.raises(AttributeError, match="You cannot add any new attribute"):
        s.str.xlabel = "a"


def test_method_on_bytes():
    lhs = Series(np.array(list("abc"), "S1").astype(object))
    rhs = Series(np.array(list("def"), "S1").astype(object))
    with pytest.raises(TypeError, match="Cannot use .str.cat with values of.*"):
        lhs.str.cat(rhs)


def test_casefold():
    # GH25405
    expected = Series(["ss", np.nan, "case", "ssd"])
    s = Series(["ß", np.nan, "case", "ßd"])
    result = s.str.casefold()

    tm.assert_series_equal(result, expected)


def test_string_array(any_string_method):
    method_name, args, kwargs = any_string_method
    if method_name == "decode":
        pytest.skip("decode requires bytes.")

    data = ["a", "bb", np.nan, "ccc"]
    a = Series(data, dtype=object)
    b = Series(data, dtype="string")

    expected = getattr(a.str, method_name)(*args, **kwargs)
    result = getattr(b.str, method_name)(*args, **kwargs)

    if isinstance(expected, Series):
        if expected.dtype == "object" and lib.is_string_array(
            expected.dropna().values,
        ):
            assert result.dtype == "string"
            result = result.astype(object)

        elif expected.dtype == "object" and lib.is_bool_array(
            expected.values, skipna=True
        ):
            assert result.dtype == "boolean"
            result = result.astype(object)

        elif expected.dtype == "bool":
            assert result.dtype == "boolean"
            result = result.astype("bool")

        elif expected.dtype == "float" and expected.isna().any():
            assert result.dtype == "Int64"
            result = result.astype("float")

    elif isinstance(expected, DataFrame):
        columns = expected.select_dtypes(include="object").columns
        assert all(result[columns].dtypes == "string")
        result[columns] = result[columns].astype(object)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "method,expected",
    [
        ("count", [2, None]),
        ("find", [0, None]),
        ("index", [0, None]),
        ("rindex", [2, None]),
    ],
)
def test_string_array_numeric_integer_array(method, expected):
    s = Series(["aba", None], dtype="string")
    result = getattr(s.str, method)("a")
    expected = Series(expected, dtype="Int64")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method,expected",
    [
        ("isdigit", [False, None, True]),
        ("isalpha", [True, None, False]),
        ("isalnum", [True, None, True]),
        ("isdigit", [False, None, True]),
    ],
)
def test_string_array_boolean_array(method, expected):
    s = Series(["a", None, "1"], dtype="string")
    result = getattr(s.str, method)()
    expected = Series(expected, dtype="boolean")
    tm.assert_series_equal(result, expected)


def test_string_array_extract():
    # https://github.com/pandas-dev/pandas/issues/30969
    # Only expand=False & multiple groups was failing
    a = Series(["a1", "b2", "cc"], dtype="string")
    b = Series(["a1", "b2", "cc"], dtype="object")
    pat = r"(\w)(\d)"

    result = a.str.extract(pat, expand=False)
    expected = b.str.extract(pat, expand=False)
    assert all(result.dtypes == "string")

    result = result.astype(object)
    tm.assert_equal(result, expected)


def test_str_get_stringarray_multiple_nans():
    s = Series(pd.array(["a", "ab", pd.NA, "abc"]))
    result = s.str.get(2)
    expected = Series(pd.array([pd.NA, pd.NA, pd.NA, "c"]))
    tm.assert_series_equal(result, expected)


def test_str_accessor_in_apply_func():
    # https://github.com/pandas-dev/pandas/issues/38979
    df = DataFrame(zip("abc", "def"))
    expected = Series(["A/D", "B/E", "C/F"])
    result = df.apply(lambda f: "/".join(f.str.upper()), axis=1)
    tm.assert_series_equal(result, expected)
