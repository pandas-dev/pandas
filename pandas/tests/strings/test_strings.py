from datetime import (
    datetime,
    timedelta,
)

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    isna,
)
import pandas._testing as tm


def assert_series_or_index_equal(left, right):
    if isinstance(left, Series):
        tm.assert_series_equal(left, right)
    else:  # Index
        tm.assert_index_equal(left, right)


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


def test_repeat_with_null(nullable_string_dtype):
    # GH: 31632
    ser = Series(["a", None], dtype=nullable_string_dtype)
    result = ser.str.repeat([3, 4])
    expected = Series(["aaa", None], dtype=nullable_string_dtype)
    tm.assert_series_equal(result, expected)

    ser = Series(["a", "b"], dtype=nullable_string_dtype)
    result = ser.str.repeat([3, None])
    expected = Series(["aaa", None], dtype=nullable_string_dtype)
    tm.assert_series_equal(result, expected)


def test_empty_str_methods(any_string_dtype):
    empty_str = empty = Series(dtype=any_string_dtype)
    if any_string_dtype == "object":
        empty_int = Series(dtype="int64")
        empty_bool = Series(dtype=bool)
    else:
        empty_int = Series(dtype="Int64")
        empty_bool = Series(dtype="boolean")
    empty_object = Series(dtype=object)
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
    tm.assert_series_equal(empty_object, empty_str.str.findall("a"))
    tm.assert_series_equal(empty_int, empty.str.find("a"))
    tm.assert_series_equal(empty_int, empty.str.rfind("a"))
    tm.assert_series_equal(empty_str, empty.str.pad(42))
    tm.assert_series_equal(empty_str, empty.str.center(42))
    tm.assert_series_equal(empty_object, empty.str.split("a"))
    tm.assert_series_equal(empty_object, empty.str.rsplit("a"))
    tm.assert_series_equal(empty_object, empty.str.partition("a", expand=False))
    tm.assert_series_equal(empty_object, empty.str.rpartition("a", expand=False))
    tm.assert_series_equal(empty_str, empty.str.slice(stop=1))
    tm.assert_series_equal(empty_str, empty.str.slice(step=1))
    tm.assert_series_equal(empty_str, empty.str.strip())
    tm.assert_series_equal(empty_str, empty.str.lstrip())
    tm.assert_series_equal(empty_str, empty.str.rstrip())
    tm.assert_series_equal(empty_str, empty.str.wrap(42))
    tm.assert_series_equal(empty_str, empty.str.get(0))
    tm.assert_series_equal(empty_object, empty_bytes.str.decode("ascii"))
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


def test_ismethods(any_string_dtype):
    values = ["A", "b", "Xy", "4", "3A", "", "TT", "55", "-", "  "]
    str_s = Series(values, dtype=any_string_dtype)
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

    dtype = "bool" if any_string_dtype == "object" else "boolean"
    tm.assert_series_equal(str_s.str.isalnum(), Series(alnum_e, dtype=dtype))
    tm.assert_series_equal(str_s.str.isalpha(), Series(alpha_e, dtype=dtype))
    tm.assert_series_equal(str_s.str.isdigit(), Series(digit_e, dtype=dtype))
    tm.assert_series_equal(str_s.str.isspace(), Series(space_e, dtype=dtype))
    tm.assert_series_equal(str_s.str.islower(), Series(lower_e, dtype=dtype))
    tm.assert_series_equal(str_s.str.isupper(), Series(upper_e, dtype=dtype))
    tm.assert_series_equal(str_s.str.istitle(), Series(title_e, dtype=dtype))

    assert str_s.str.isalnum().tolist() == [v.isalnum() for v in values]
    assert str_s.str.isalpha().tolist() == [v.isalpha() for v in values]
    assert str_s.str.isdigit().tolist() == [v.isdigit() for v in values]
    assert str_s.str.isspace().tolist() == [v.isspace() for v in values]
    assert str_s.str.islower().tolist() == [v.islower() for v in values]
    assert str_s.str.isupper().tolist() == [v.isupper() for v in values]
    assert str_s.str.istitle().tolist() == [v.istitle() for v in values]


def test_isnumeric(any_string_dtype):
    # 0x00bc: ¼ VULGAR FRACTION ONE QUARTER
    # 0x2605: ★ not number
    # 0x1378: ፸ ETHIOPIC NUMBER SEVENTY
    # 0xFF13: ３ Em 3
    values = ["A", "3", "¼", "★", "፸", "３", "four"]
    s = Series(values, dtype=any_string_dtype)
    numeric_e = [False, True, True, False, True, True, False]
    decimal_e = [False, True, False, False, False, True, False]
    dtype = "bool" if any_string_dtype == "object" else "boolean"
    tm.assert_series_equal(s.str.isnumeric(), Series(numeric_e, dtype=dtype))
    tm.assert_series_equal(s.str.isdecimal(), Series(decimal_e, dtype=dtype))

    unicodes = ["A", "3", "¼", "★", "፸", "３", "four"]
    assert s.str.isnumeric().tolist() == [v.isnumeric() for v in unicodes]
    assert s.str.isdecimal().tolist() == [v.isdecimal() for v in unicodes]

    values = ["A", np.nan, "¼", "★", np.nan, "３", "four"]
    s = Series(values, dtype=any_string_dtype)
    numeric_e = [False, np.nan, True, False, np.nan, True, False]
    decimal_e = [False, np.nan, False, False, np.nan, True, False]
    dtype = "object" if any_string_dtype == "object" else "boolean"
    tm.assert_series_equal(s.str.isnumeric(), Series(numeric_e, dtype=dtype))
    tm.assert_series_equal(s.str.isdecimal(), Series(decimal_e, dtype=dtype))


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


def test_len(any_string_dtype):
    values = Series(
        ["foo", "fooo", "fooooo", np.nan, "fooooooo", "foo\n", "あ"],
        dtype=any_string_dtype,
    )

    result = values.str.len()
    expected_dtype = "float64" if any_string_dtype == "object" else "Int64"
    expected = Series([3, 4, 6, np.nan, 8, 4, 1], dtype=expected_dtype)
    tm.assert_series_equal(result, expected)


def test_len_mixed():
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


def test_strip_lstrip_rstrip(any_string_dtype):
    values = Series(["  aa   ", " bb \n", np.nan, "cc  "], dtype=any_string_dtype)

    result = values.str.strip()
    exp = Series(["aa", "bb", np.nan, "cc"], dtype=any_string_dtype)
    tm.assert_series_equal(result, exp)

    result = values.str.lstrip()
    exp = Series(["aa   ", "bb \n", np.nan, "cc  "], dtype=any_string_dtype)
    tm.assert_series_equal(result, exp)

    result = values.str.rstrip()
    exp = Series(["  aa", " bb", np.nan, "cc"], dtype=any_string_dtype)
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


def test_strip_lstrip_rstrip_args(any_string_dtype):
    values = Series(["xxABCxx", "xx BNSD", "LDFJH xx"], dtype=any_string_dtype)

    rs = values.str.strip("x")
    xp = Series(["ABC", " BNSD", "LDFJH "], dtype=any_string_dtype)
    tm.assert_series_equal(rs, xp)

    rs = values.str.lstrip("x")
    xp = Series(["ABCxx", " BNSD", "LDFJH xx"], dtype=any_string_dtype)
    tm.assert_series_equal(rs, xp)

    rs = values.str.rstrip("x")
    xp = Series(["xxABC", "xx BNSD", "LDFJH "], dtype=any_string_dtype)
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


def test_str_accessor_in_apply_func():
    # https://github.com/pandas-dev/pandas/issues/38979
    df = DataFrame(zip("abc", "def"))
    expected = Series(["A/D", "B/E", "C/F"])
    result = df.apply(lambda f: "/".join(f.str.upper()), axis=1)
    tm.assert_series_equal(result, expected)
