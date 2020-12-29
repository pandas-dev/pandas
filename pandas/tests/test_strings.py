from datetime import datetime, timedelta
import re

import numpy as np
import pytest

from pandas._libs import lib

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series, concat, isna, notna
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


class TestStringMethods:
    def test_api(self):

        # GH 6106, GH 9322
        assert Series.str is strings.StringMethods
        assert isinstance(Series([""]).str, strings.StringMethods)

    def test_api_mi_raises(self):
        # GH 23679
        mi = MultiIndex.from_arrays([["a", "b", "c"]])
        msg = "Can only use .str accessor with Index, not MultiIndex"
        with pytest.raises(AttributeError, match=msg):
            mi.str
        assert not hasattr(mi, "str")

    @pytest.mark.parametrize("dtype", [object, "category"])
    def test_api_per_dtype(self, index_or_series, dtype, any_skipna_inferred_dtype):
        # one instance of parametrized fixture
        box = index_or_series
        inferred_dtype, values = any_skipna_inferred_dtype

        t = box(values, dtype=dtype)  # explicit dtype to avoid casting

        types_passing_constructor = [
            "string",
            "unicode",
            "empty",
            "bytes",
            "mixed",
            "mixed-integer",
        ]
        if inferred_dtype in types_passing_constructor:
            # GH 6106
            assert isinstance(t.str, strings.StringMethods)
        else:
            # GH 9184, GH 23011, GH 23163
            msg = "Can only use .str accessor with string values.*"
            with pytest.raises(AttributeError, match=msg):
                t.str
            assert not hasattr(t, "str")

    @pytest.mark.parametrize("dtype", [object, "category"])
    def test_api_per_method(
        self,
        index_or_series,
        dtype,
        any_allowed_skipna_inferred_dtype,
        any_string_method,
        request,
    ):
        # this test does not check correctness of the different methods,
        # just that the methods work on the specified (inferred) dtypes,
        # and raise on all others
        box = index_or_series

        # one instance of each parametrized fixture
        inferred_dtype, values = any_allowed_skipna_inferred_dtype
        method_name, args, kwargs = any_string_method

        # TODO: get rid of these xfails
        reason = None
        if box is Index and values.size == 0:
            if method_name in ["partition", "rpartition"] and kwargs.get(
                "expand", True
            ):
                reason = "Method cannot deal with empty Index"
            elif method_name == "split" and kwargs.get("expand", None):
                reason = "Split fails on empty Series when expand=True"
            elif method_name == "get_dummies":
                reason = "Need to fortify get_dummies corner cases"

        elif box is Index and inferred_dtype == "empty" and dtype == object:
            if method_name == "get_dummies":
                reason = "Need to fortify get_dummies corner cases"

        if reason is not None:
            mark = pytest.mark.xfail(reason=reason)
            request.node.add_marker(mark)

        t = box(values, dtype=dtype)  # explicit dtype to avoid casting
        method = getattr(t.str, method_name)

        bytes_allowed = method_name in ["decode", "get", "len", "slice"]
        # as of v0.23.4, all methods except 'cat' are very lenient with the
        # allowed data types, just returning NaN for entries that error.
        # This could be changed with an 'errors'-kwarg to the `str`-accessor,
        # see discussion in GH 13877
        mixed_allowed = method_name not in ["cat"]

        allowed_types = (
            ["string", "unicode", "empty"]
            + ["bytes"] * bytes_allowed
            + ["mixed", "mixed-integer"] * mixed_allowed
        )

        if inferred_dtype in allowed_types:
            # xref GH 23555, GH 23556
            method(*args, **kwargs)  # works!
        else:
            # GH 23011, GH 23163
            msg = (
                f"Cannot use .str.{method_name} with values of "
                f"inferred dtype {repr(inferred_dtype)}."
            )
            with pytest.raises(TypeError, match=msg):
                method(*args, **kwargs)

    def test_api_for_categorical(self, any_string_method):
        # https://github.com/pandas-dev/pandas/issues/10661
        s = Series(list("aabb"))
        s = s + " " + s
        c = s.astype("category")
        assert isinstance(c.str, strings.StringMethods)

        method_name, args, kwargs = any_string_method

        result = getattr(c.str, method_name)(*args, **kwargs)
        expected = getattr(s.str, method_name)(*args, **kwargs)

        if isinstance(result, DataFrame):
            tm.assert_frame_equal(result, expected)
        elif isinstance(result, Series):
            tm.assert_series_equal(result, expected)
        else:
            # str.cat(others=None) returns string, for example
            assert result == expected

    def test_iter(self):
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

    def test_iter_empty(self):
        ds = Series([], dtype=object)

        i, s = 100, 1

        with tm.assert_produces_warning(FutureWarning):
            for i, s in enumerate(ds.str):
                pass

        # nothing to iterate over so nothing defined values should remain
        # unchanged
        assert i == 100
        assert s == 1

    def test_iter_single_element(self):
        ds = Series(["a"])

        with tm.assert_produces_warning(FutureWarning):
            for i, s in enumerate(ds.str):
                pass

        assert not i
        tm.assert_series_equal(ds, s)

    def test_iter_object_try_string(self):
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

    @pytest.mark.parametrize("other", [None, Series, Index])
    def test_str_cat_name(self, index_or_series, other):
        # GH 21053
        box = index_or_series
        values = ["a", "b"]
        if other:
            other = other(values)
        else:
            other = values
        result = box(values, name="name").str.cat(other, sep=",")
        assert result.name == "name"

    def test_str_cat(self, index_or_series):
        box = index_or_series
        # test_cat above tests "str_cat" from ndarray;
        # here testing "str.cat" from Series/Indext to ndarray/list
        s = box(["a", "a", "b", "b", "c", np.nan])

        # single array
        result = s.str.cat()
        expected = "aabbc"
        assert result == expected

        result = s.str.cat(na_rep="-")
        expected = "aabbc-"
        assert result == expected

        result = s.str.cat(sep="_", na_rep="NA")
        expected = "a_a_b_b_c_NA"
        assert result == expected

        t = np.array(["a", np.nan, "b", "d", "foo", np.nan], dtype=object)
        expected = box(["aa", "a-", "bb", "bd", "cfoo", "--"])

        # Series/Index with array
        result = s.str.cat(t, na_rep="-")
        assert_series_or_index_equal(result, expected)

        # Series/Index with list
        result = s.str.cat(list(t), na_rep="-")
        assert_series_or_index_equal(result, expected)

        # errors for incorrect lengths
        rgx = r"If `others` contains arrays or lists \(or other list-likes.*"
        z = Series(["1", "2", "3"])

        with pytest.raises(ValueError, match=rgx):
            s.str.cat(z.values)

        with pytest.raises(ValueError, match=rgx):
            s.str.cat(list(z))

    def test_str_cat_raises_intuitive_error(self, index_or_series):
        # GH 11334
        box = index_or_series
        s = box(["a", "b", "c", "d"])
        message = "Did you mean to supply a `sep` keyword?"
        with pytest.raises(ValueError, match=message):
            s.str.cat("|")
        with pytest.raises(ValueError, match=message):
            s.str.cat("    ")

    @pytest.mark.parametrize("sep", ["", None])
    @pytest.mark.parametrize("dtype_target", ["object", "category"])
    @pytest.mark.parametrize("dtype_caller", ["object", "category"])
    def test_str_cat_categorical(
        self, index_or_series, dtype_caller, dtype_target, sep
    ):
        box = index_or_series

        s = Index(["a", "a", "b", "a"], dtype=dtype_caller)
        s = s if box == Index else Series(s, index=s)
        t = Index(["b", "a", "b", "c"], dtype=dtype_target)

        expected = Index(["ab", "aa", "bb", "ac"])
        expected = expected if box == Index else Series(expected, index=s)

        # Series/Index with unaligned Index -> t.values
        result = s.str.cat(t.values, sep=sep)
        assert_series_or_index_equal(result, expected)

        # Series/Index with Series having matching Index
        t = Series(t.values, index=s)
        result = s.str.cat(t, sep=sep)
        assert_series_or_index_equal(result, expected)

        # Series/Index with Series.values
        result = s.str.cat(t.values, sep=sep)
        assert_series_or_index_equal(result, expected)

        # Series/Index with Series having different Index
        t = Series(t.values, index=t.values)
        expected = Index(["aa", "aa", "aa", "bb", "bb"])
        expected = (
            expected if box == Index else Series(expected, index=expected.str[:1])
        )

        result = s.str.cat(t, sep=sep)
        assert_series_or_index_equal(result, expected)

    # test integer/float dtypes (inferred by constructor) and mixed
    @pytest.mark.parametrize(
        "data",
        [[1, 2, 3], [0.1, 0.2, 0.3], [1, 2, "b"]],
        ids=["integers", "floats", "mixed"],
    )
    # without dtype=object, np.array would cast [1, 2, 'b'] to ['1', '2', 'b']
    @pytest.mark.parametrize(
        "box",
        [Series, Index, list, lambda x: np.array(x, dtype=object)],
        ids=["Series", "Index", "list", "np.array"],
    )
    def test_str_cat_wrong_dtype_raises(self, box, data):
        # GH 22722
        s = Series(["a", "b", "c"])
        t = box(data)

        msg = "Concatenation requires list-likes containing only strings.*"
        with pytest.raises(TypeError, match=msg):
            # need to use outer and na_rep, as otherwise Index would not raise
            s.str.cat(t, join="outer", na_rep="-")

    def test_str_cat_mixed_inputs(self, index_or_series):
        box = index_or_series
        s = Index(["a", "b", "c", "d"])
        s = s if box == Index else Series(s, index=s)

        t = Series(["A", "B", "C", "D"], index=s.values)
        d = concat([t, Series(s, index=s)], axis=1)

        expected = Index(["aAa", "bBb", "cCc", "dDd"])
        expected = expected if box == Index else Series(expected.values, index=s.values)

        # Series/Index with DataFrame
        result = s.str.cat(d)
        assert_series_or_index_equal(result, expected)

        # Series/Index with two-dimensional ndarray
        result = s.str.cat(d.values)
        assert_series_or_index_equal(result, expected)

        # Series/Index with list of Series
        result = s.str.cat([t, s])
        assert_series_or_index_equal(result, expected)

        # Series/Index with mixed list of Series/array
        result = s.str.cat([t, s.values])
        assert_series_or_index_equal(result, expected)

        # Series/Index with list of Series; different indexes
        t.index = ["b", "c", "d", "a"]
        expected = box(["aDa", "bAb", "cBc", "dCd"])
        expected = expected if box == Index else Series(expected.values, index=s.values)
        result = s.str.cat([t, s])
        assert_series_or_index_equal(result, expected)

        # Series/Index with mixed list; different index
        result = s.str.cat([t, s.values])
        assert_series_or_index_equal(result, expected)

        # Series/Index with DataFrame; different indexes
        d.index = ["b", "c", "d", "a"]
        expected = box(["aDd", "bAa", "cBb", "dCc"])
        expected = expected if box == Index else Series(expected.values, index=s.values)
        result = s.str.cat(d)
        assert_series_or_index_equal(result, expected)

        # errors for incorrect lengths
        rgx = r"If `others` contains arrays or lists \(or other list-likes.*"
        z = Series(["1", "2", "3"])
        e = concat([z, z], axis=1)

        # two-dimensional ndarray
        with pytest.raises(ValueError, match=rgx):
            s.str.cat(e.values)

        # list of list-likes
        with pytest.raises(ValueError, match=rgx):
            s.str.cat([z.values, s.values])

        # mixed list of Series/list-like
        with pytest.raises(ValueError, match=rgx):
            s.str.cat([z.values, s])

        # errors for incorrect arguments in list-like
        rgx = "others must be Series, Index, DataFrame,.*"
        # make sure None/NaN do not crash checks in _get_series_list
        u = Series(["a", np.nan, "c", None])

        # mix of string and Series
        with pytest.raises(TypeError, match=rgx):
            s.str.cat([u, "u"])

        # DataFrame in list
        with pytest.raises(TypeError, match=rgx):
            s.str.cat([u, d])

        # 2-dim ndarray in list
        with pytest.raises(TypeError, match=rgx):
            s.str.cat([u, d.values])

        # nested lists
        with pytest.raises(TypeError, match=rgx):
            s.str.cat([u, [u, d]])

        # forbidden input type: set
        # GH 23009
        with pytest.raises(TypeError, match=rgx):
            s.str.cat(set(u))

        # forbidden input type: set in list
        # GH 23009
        with pytest.raises(TypeError, match=rgx):
            s.str.cat([u, set(u)])

        # other forbidden input type, e.g. int
        with pytest.raises(TypeError, match=rgx):
            s.str.cat(1)

        # nested list-likes
        with pytest.raises(TypeError, match=rgx):
            s.str.cat(iter([t.values, list(s)]))

    @pytest.mark.parametrize("join", ["left", "outer", "inner", "right"])
    def test_str_cat_align_indexed(self, index_or_series, join):
        # https://github.com/pandas-dev/pandas/issues/18657
        box = index_or_series

        s = Series(["a", "b", "c", "d"], index=["a", "b", "c", "d"])
        t = Series(["D", "A", "E", "B"], index=["d", "a", "e", "b"])
        sa, ta = s.align(t, join=join)
        # result after manual alignment of inputs
        expected = sa.str.cat(ta, na_rep="-")

        if box == Index:
            s = Index(s)
            sa = Index(sa)
            expected = Index(expected)

        result = s.str.cat(t, join=join, na_rep="-")
        assert_series_or_index_equal(result, expected)

    @pytest.mark.parametrize("join", ["left", "outer", "inner", "right"])
    def test_str_cat_align_mixed_inputs(self, join):
        s = Series(["a", "b", "c", "d"])
        t = Series(["d", "a", "e", "b"], index=[3, 0, 4, 1])
        d = concat([t, t], axis=1)

        expected_outer = Series(["aaa", "bbb", "c--", "ddd", "-ee"])
        expected = expected_outer.loc[s.index.join(t.index, how=join)]

        # list of Series
        result = s.str.cat([t, t], join=join, na_rep="-")
        tm.assert_series_equal(result, expected)

        # DataFrame
        result = s.str.cat(d, join=join, na_rep="-")
        tm.assert_series_equal(result, expected)

        # mixed list of indexed/unindexed
        u = np.array(["A", "B", "C", "D"])
        expected_outer = Series(["aaA", "bbB", "c-C", "ddD", "-e-"])
        # joint index of rhs [t, u]; u will be forced have index of s
        rhs_idx = (
            t.index.intersection(s.index) if join == "inner" else t.index.union(s.index)
        )

        expected = expected_outer.loc[s.index.join(rhs_idx, how=join)]
        result = s.str.cat([t, u], join=join, na_rep="-")
        tm.assert_series_equal(result, expected)

        with pytest.raises(TypeError, match="others must be Series,.*"):
            # nested lists are forbidden
            s.str.cat([t, list(u)], join=join)

        # errors for incorrect lengths
        rgx = r"If `others` contains arrays or lists \(or other list-likes.*"
        z = Series(["1", "2", "3"]).values

        # unindexed object of wrong length
        with pytest.raises(ValueError, match=rgx):
            s.str.cat(z, join=join)

        # unindexed object of wrong length in list
        with pytest.raises(ValueError, match=rgx):
            s.str.cat([t, z], join=join)

    def test_str_cat_all_na(self, index_or_series, index_or_series2):
        # GH 24044
        box = index_or_series
        other = index_or_series2

        # check that all NaNs in caller / target work
        s = Index(["a", "b", "c", "d"])
        s = s if box == Index else Series(s, index=s)
        t = other([np.nan] * 4, dtype=object)
        # add index of s for alignment
        t = t if other == Index else Series(t, index=s)

        # all-NA target
        if box == Series:
            expected = Series([np.nan] * 4, index=s.index, dtype=object)
        else:  # box == Index
            expected = Index([np.nan] * 4, dtype=object)
        result = s.str.cat(t, join="left")
        assert_series_or_index_equal(result, expected)

        # all-NA caller (only for Series)
        if other == Series:
            expected = Series([np.nan] * 4, dtype=object, index=t.index)
            result = t.str.cat(s, join="left")
            tm.assert_series_equal(result, expected)

    def test_str_cat_special_cases(self):
        s = Series(["a", "b", "c", "d"])
        t = Series(["d", "a", "e", "b"], index=[3, 0, 4, 1])

        # iterator of elements with different types
        expected = Series(["aaa", "bbb", "c-c", "ddd", "-e-"])
        result = s.str.cat(iter([t, s.values]), join="outer", na_rep="-")
        tm.assert_series_equal(result, expected)

        # right-align with different indexes in others
        expected = Series(["aa-", "d-d"], index=[0, 3])
        result = s.str.cat([t.loc[[0]], t.loc[[3]]], join="right", na_rep="-")
        tm.assert_series_equal(result, expected)

    def test_cat_on_filtered_index(self):
        df = DataFrame(
            index=MultiIndex.from_product(
                [[2011, 2012], [1, 2, 3]], names=["year", "month"]
            )
        )

        df = df.reset_index()
        df = df[df.month > 1]

        str_year = df.year.astype("str")
        str_month = df.month.astype("str")
        str_both = str_year.str.cat(str_month, sep=" ")

        assert str_both.loc[1] == "2011 2"

        str_multiple = str_year.str.cat([str_month, str_month], sep=" ")

        assert str_multiple.loc[1] == "2011 2 2"

    def test_count(self):
        values = np.array(
            ["foo", "foofoo", np.nan, "foooofooofommmfoo"], dtype=np.object_
        )

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

    def test_contains(self):
        values = np.array(
            ["foo", np.nan, "fooommm__foo", "mmm_", "foommm[_]+bar"], dtype=np.object_
        )
        values = Series(values)
        pat = "mmm[_]+"

        result = values.str.contains(pat)
        expected = Series(
            np.array([False, np.nan, True, True, False], dtype=np.object_)
        )
        tm.assert_series_equal(result, expected)

        result = values.str.contains(pat, regex=False)
        expected = Series(
            np.array([False, np.nan, False, False, True], dtype=np.object_)
        )
        tm.assert_series_equal(result, expected)

        values = Series(np.array(["foo", "xyz", "fooommm__foo", "mmm_"], dtype=object))
        result = values.str.contains(pat)
        expected = Series(np.array([False, False, True, True]))
        assert result.dtype == np.bool_
        tm.assert_series_equal(result, expected)

        # case insensitive using regex
        values = Series(np.array(["Foo", "xYz", "fOOomMm__fOo", "MMM_"], dtype=object))
        result = values.str.contains("FOO|mmm", case=False)
        expected = Series(np.array([True, False, True, True]))
        tm.assert_series_equal(result, expected)

        # case insensitive without regex
        result = Series(values).str.contains("foo", regex=False, case=False)
        expected = Series(np.array([True, False, True, False]))
        tm.assert_series_equal(result, expected)

        # mixed
        mixed = Series(
            np.array(
                ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
                dtype=object,
            )
        )
        rs = mixed.str.contains("o")
        xp = Series(
            np.array(
                [False, np.nan, False, np.nan, np.nan, True, np.nan, np.nan, np.nan],
                dtype=np.object_,
            )
        )
        tm.assert_series_equal(rs, xp)

        rs = mixed.str.contains("o")
        xp = Series(
            [False, np.nan, False, np.nan, np.nan, True, np.nan, np.nan, np.nan]
        )
        assert isinstance(rs, Series)
        tm.assert_series_equal(rs, xp)

        # unicode
        values = Series(
            np.array(["foo", np.nan, "fooommm__foo", "mmm_"], dtype=np.object_)
        )
        pat = "mmm[_]+"

        result = values.str.contains(pat)
        expected = Series(np.array([False, np.nan, True, True], dtype=np.object_))
        tm.assert_series_equal(result, expected)

        result = values.str.contains(pat, na=False)
        expected = Series(np.array([False, False, True, True]))
        tm.assert_series_equal(result, expected)

        values = Series(
            np.array(["foo", "xyz", "fooommm__foo", "mmm_"], dtype=np.object_)
        )
        result = values.str.contains(pat)
        expected = Series(np.array([False, False, True, True]))
        assert result.dtype == np.bool_
        tm.assert_series_equal(result, expected)

    def test_contains_for_object_category(self):
        # gh 22158

        # na for category
        values = Series(["a", "b", "c", "a", np.nan], dtype="category")
        result = values.str.contains("a", na=True)
        expected = Series([True, False, False, True, True])
        tm.assert_series_equal(result, expected)

        result = values.str.contains("a", na=False)
        expected = Series([True, False, False, True, False])
        tm.assert_series_equal(result, expected)

        # na for objects
        values = Series(["a", "b", "c", "a", np.nan])
        result = values.str.contains("a", na=True)
        expected = Series([True, False, False, True, True])
        tm.assert_series_equal(result, expected)

        result = values.str.contains("a", na=False)
        expected = Series([True, False, False, True, False])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [None, "category"])
    @pytest.mark.parametrize("null_value", [None, np.nan, pd.NA])
    @pytest.mark.parametrize("na", [True, False])
    def test_startswith(self, dtype, null_value, na):
        # add category dtype parametrizations for GH-36241
        values = Series(
            ["om", null_value, "foo_nom", "nom", "bar_foo", null_value, "foo"],
            dtype=dtype,
        )

        result = values.str.startswith("foo")
        exp = Series([False, np.nan, True, False, False, np.nan, True])
        tm.assert_series_equal(result, exp)

        result = values.str.startswith("foo", na=na)
        exp = Series([False, na, True, False, False, na, True])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = np.array(
            ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
            dtype=np.object_,
        )
        rs = Series(mixed).str.startswith("f")
        xp = Series(
            [False, np.nan, False, np.nan, np.nan, True, np.nan, np.nan, np.nan]
        )
        tm.assert_series_equal(rs, xp)

    @pytest.mark.parametrize("dtype", [None, "category"])
    @pytest.mark.parametrize("null_value", [None, np.nan, pd.NA])
    @pytest.mark.parametrize("na", [True, False])
    def test_endswith(self, dtype, null_value, na):
        # add category dtype parametrizations for GH-36241
        values = Series(
            ["om", null_value, "foo_nom", "nom", "bar_foo", null_value, "foo"],
            dtype=dtype,
        )

        result = values.str.endswith("foo")
        exp = Series([False, np.nan, False, False, True, np.nan, True])
        tm.assert_series_equal(result, exp)

        result = values.str.endswith("foo", na=na)
        exp = Series([False, na, False, False, True, na, True])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = np.array(
            ["a", np.nan, "b", True, datetime.today(), "foo", None, 1, 2.0],
            dtype=object,
        )
        rs = Series(mixed).str.endswith("f")
        xp = Series(
            [False, np.nan, False, np.nan, np.nan, False, np.nan, np.nan, np.nan]
        )
        tm.assert_series_equal(rs, xp)

    def test_title(self):
        values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])

        result = values.str.title()
        exp = Series(["Foo", "Bar", np.nan, "Blah", "Blurg"])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(
            ["FOO", np.nan, "bar", True, datetime.today(), "blah", None, 1, 2.0]
        )
        mixed = mixed.str.title()
        exp = Series(
            ["Foo", np.nan, "Bar", np.nan, np.nan, "Blah", np.nan, np.nan, np.nan]
        )
        tm.assert_almost_equal(mixed, exp)

    def test_lower_upper(self):
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

    def test_capitalize(self):
        values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])
        result = values.str.capitalize()
        exp = Series(["Foo", "Bar", np.nan, "Blah", "Blurg"])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(
            ["FOO", np.nan, "bar", True, datetime.today(), "blah", None, 1, 2.0]
        )
        mixed = mixed.str.capitalize()
        exp = Series(
            ["Foo", np.nan, "Bar", np.nan, np.nan, "Blah", np.nan, np.nan, np.nan]
        )
        tm.assert_almost_equal(mixed, exp)

    def test_swapcase(self):
        values = Series(["FOO", "BAR", np.nan, "Blah", "blurg"])
        result = values.str.swapcase()
        exp = Series(["foo", "bar", np.nan, "bLAH", "BLURG"])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(
            ["FOO", np.nan, "bar", True, datetime.today(), "Blah", None, 1, 2.0]
        )
        mixed = mixed.str.swapcase()
        exp = Series(
            ["foo", np.nan, "BAR", np.nan, np.nan, "bLAH", np.nan, np.nan, np.nan]
        )
        tm.assert_almost_equal(mixed, exp)

    def test_casemethods(self):
        values = ["aaa", "bbb", "CCC", "Dddd", "eEEE"]
        s = Series(values)
        assert s.str.lower().tolist() == [v.lower() for v in values]
        assert s.str.upper().tolist() == [v.upper() for v in values]
        assert s.str.title().tolist() == [v.title() for v in values]
        assert s.str.capitalize().tolist() == [v.capitalize() for v in values]
        assert s.str.swapcase().tolist() == [v.swapcase() for v in values]

    def test_replace(self):
        values = Series(["fooBAD__barBAD", np.nan])

        result = values.str.replace("BAD[_]*", "", regex=True)
        exp = Series(["foobar", np.nan])
        tm.assert_series_equal(result, exp)

        result = values.str.replace("BAD[_]*", "", n=1, regex=True)
        exp = Series(["foobarBAD", np.nan])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(
            ["aBAD", np.nan, "bBAD", True, datetime.today(), "fooBAD", None, 1, 2.0]
        )

        rs = Series(mixed).str.replace("BAD[_]*", "", regex=True)
        xp = Series(["a", np.nan, "b", np.nan, np.nan, "foo", np.nan, np.nan, np.nan])
        assert isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # flags + unicode
        values = Series([b"abcd,\xc3\xa0".decode("utf-8")])
        exp = Series([b"abcd, \xc3\xa0".decode("utf-8")])
        result = values.str.replace(
            r"(?<=\w),(?=\w)", ", ", flags=re.UNICODE, regex=True
        )
        tm.assert_series_equal(result, exp)

        # GH 13438
        msg = "repl must be a string or callable"
        for klass in (Series, Index):
            for repl in (None, 3, {"a": "b"}):
                for data in (["a", "b", None], ["a", "b", "c", "ad"]):
                    values = klass(data)
                    with pytest.raises(TypeError, match=msg):
                        values.str.replace("a", repl)

    def test_replace_callable(self):
        # GH 15055
        values = Series(["fooBAD__barBAD", np.nan])

        # test with callable
        repl = lambda m: m.group(0).swapcase()
        result = values.str.replace("[a-z][A-Z]{2}", repl, n=2, regex=True)
        exp = Series(["foObaD__baRbaD", np.nan])
        tm.assert_series_equal(result, exp)

        # test with wrong number of arguments, raising an error
        p_err = (
            r"((takes)|(missing)) (?(2)from \d+ to )?\d+ "
            r"(?(3)required )positional arguments?"
        )

        repl = lambda: None
        with pytest.raises(TypeError, match=p_err):
            values.str.replace("a", repl)

        repl = lambda m, x: None
        with pytest.raises(TypeError, match=p_err):
            values.str.replace("a", repl)

        repl = lambda m, x, y=None: None
        with pytest.raises(TypeError, match=p_err):
            values.str.replace("a", repl)

        # test regex named groups
        values = Series(["Foo Bar Baz", np.nan])
        pat = r"(?P<first>\w+) (?P<middle>\w+) (?P<last>\w+)"
        repl = lambda m: m.group("middle").swapcase()
        result = values.str.replace(pat, repl, regex=True)
        exp = Series(["bAR", np.nan])
        tm.assert_series_equal(result, exp)

    def test_replace_compiled_regex(self):
        # GH 15446
        values = Series(["fooBAD__barBAD", np.nan])

        # test with compiled regex
        pat = re.compile(r"BAD[_]*")
        result = values.str.replace(pat, "", regex=True)
        exp = Series(["foobar", np.nan])
        tm.assert_series_equal(result, exp)

        result = values.str.replace(pat, "", n=1, regex=True)
        exp = Series(["foobarBAD", np.nan])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(
            ["aBAD", np.nan, "bBAD", True, datetime.today(), "fooBAD", None, 1, 2.0]
        )

        rs = Series(mixed).str.replace(pat, "", regex=True)
        xp = Series(["a", np.nan, "b", np.nan, np.nan, "foo", np.nan, np.nan, np.nan])
        assert isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # flags + unicode
        values = Series([b"abcd,\xc3\xa0".decode("utf-8")])
        exp = Series([b"abcd, \xc3\xa0".decode("utf-8")])
        pat = re.compile(r"(?<=\w),(?=\w)", flags=re.UNICODE)
        result = values.str.replace(pat, ", ")
        tm.assert_series_equal(result, exp)

        # case and flags provided to str.replace will have no effect
        # and will produce warnings
        values = Series(["fooBAD__barBAD__bad", np.nan])
        pat = re.compile(r"BAD[_]*")

        with pytest.raises(ValueError, match="case and flags cannot be"):
            result = values.str.replace(pat, "", flags=re.IGNORECASE)

        with pytest.raises(ValueError, match="case and flags cannot be"):
            result = values.str.replace(pat, "", case=False)

        with pytest.raises(ValueError, match="case and flags cannot be"):
            result = values.str.replace(pat, "", case=True)

        # test with callable
        values = Series(["fooBAD__barBAD", np.nan])
        repl = lambda m: m.group(0).swapcase()
        pat = re.compile("[a-z][A-Z]{2}")
        result = values.str.replace(pat, repl, n=2)
        exp = Series(["foObaD__baRbaD", np.nan])
        tm.assert_series_equal(result, exp)

    def test_replace_literal(self):
        # GH16808 literal replace (regex=False vs regex=True)
        values = Series(["f.o", "foo", np.nan])
        exp = Series(["bao", "bao", np.nan])
        result = values.str.replace("f.", "ba", regex=True)
        tm.assert_series_equal(result, exp)

        exp = Series(["bao", "foo", np.nan])
        result = values.str.replace("f.", "ba", regex=False)
        tm.assert_series_equal(result, exp)

        # Cannot do a literal replace if given a callable repl or compiled
        # pattern
        callable_repl = lambda m: m.group(0).swapcase()
        compiled_pat = re.compile("[a-z][A-Z]{2}")

        msg = "Cannot use a callable replacement when regex=False"
        with pytest.raises(ValueError, match=msg):
            values.str.replace("abc", callable_repl, regex=False)

        msg = "Cannot use a compiled regex as replacement pattern with regex=False"
        with pytest.raises(ValueError, match=msg):
            values.str.replace(compiled_pat, "", regex=False)

    def test_repeat(self):
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

    def test_repeat_with_null(self):
        # GH: 31632
        values = Series(["a", None], dtype="string")
        result = values.str.repeat([3, 4])
        exp = Series(["aaa", None], dtype="string")
        tm.assert_series_equal(result, exp)

        values = Series(["a", "b"], dtype="string")
        result = values.str.repeat([3, None])
        exp = Series(["aaa", None], dtype="string")
        tm.assert_series_equal(result, exp)

    def test_match(self):
        # New match behavior introduced in 0.13
        values = Series(["fooBAD__barBAD", np.nan, "foo"])
        result = values.str.match(".*(BAD[_]+).*(BAD)")
        exp = Series([True, np.nan, False])
        tm.assert_series_equal(result, exp)

        values = Series(["fooBAD__barBAD", "BAD_BADleroybrown", np.nan, "foo"])
        result = values.str.match(".*BAD[_]+.*BAD")
        exp = Series([True, True, np.nan, False])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(
            [
                "aBAD_BAD",
                np.nan,
                "BAD_b_BAD",
                True,
                datetime.today(),
                "foo",
                None,
                1,
                2.0,
            ]
        )
        rs = Series(mixed).str.match(".*(BAD[_]+).*(BAD)")
        xp = Series([True, np.nan, True, np.nan, np.nan, False, np.nan, np.nan, np.nan])
        assert isinstance(rs, Series)
        tm.assert_series_equal(rs, xp)

        # na GH #6609
        res = Series(["a", 0, np.nan]).str.match("a", na=False)
        exp = Series([True, False, False])
        tm.assert_series_equal(exp, res)
        res = Series(["a", 0, np.nan]).str.match("a")
        exp = Series([True, np.nan, np.nan])
        tm.assert_series_equal(exp, res)

        values = Series(["ab", "AB", "abc", "ABC"])
        result = values.str.match("ab", case=False)
        expected = Series([True, True, True, True])
        tm.assert_series_equal(result, expected)

    def test_fullmatch(self):
        # GH 32806
        values = Series(["fooBAD__barBAD", "BAD_BADleroybrown", np.nan, "foo"])
        result = values.str.fullmatch(".*BAD[_]+.*BAD")
        exp = Series([True, False, np.nan, False])
        tm.assert_series_equal(result, exp)

        # Make sure that the new string arrays work
        string_values = Series(
            ["fooBAD__barBAD", "BAD_BADleroybrown", np.nan, "foo"], dtype="string"
        )
        result = string_values.str.fullmatch(".*BAD[_]+.*BAD")
        # Result is nullable boolean with StringDtype
        string_exp = Series([True, False, np.nan, False], dtype="boolean")
        tm.assert_series_equal(result, string_exp)

        values = Series(["ab", "AB", "abc", "ABC"])
        result = values.str.fullmatch("ab", case=False)
        expected = Series([True, True, False, False])
        tm.assert_series_equal(result, expected)

    def test_extract_expand_None(self):
        values = Series(["fooBAD__barBAD", np.nan, "foo"])
        with pytest.raises(ValueError, match="expand must be True or False"):
            values.str.extract(".*(BAD[_]+).*(BAD)", expand=None)

    def test_extract_expand_unspecified(self):
        values = Series(["fooBAD__barBAD", np.nan, "foo"])
        result_unspecified = values.str.extract(".*(BAD[_]+).*")
        assert isinstance(result_unspecified, DataFrame)
        result_true = values.str.extract(".*(BAD[_]+).*", expand=True)
        tm.assert_frame_equal(result_unspecified, result_true)

    def test_extract_expand_False(self):
        # Contains tests like those in test_match and some others.
        values = Series(["fooBAD__barBAD", np.nan, "foo"])
        er = [np.nan, np.nan]  # empty row

        result = values.str.extract(".*(BAD[_]+).*(BAD)", expand=False)
        exp = DataFrame([["BAD__", "BAD"], er, er])
        tm.assert_frame_equal(result, exp)

        # mixed
        mixed = Series(
            [
                "aBAD_BAD",
                np.nan,
                "BAD_b_BAD",
                True,
                datetime.today(),
                "foo",
                None,
                1,
                2.0,
            ]
        )

        rs = Series(mixed).str.extract(".*(BAD[_]+).*(BAD)", expand=False)
        exp = DataFrame([["BAD_", "BAD"], er, ["BAD_", "BAD"], er, er, er, er, er, er])
        tm.assert_frame_equal(rs, exp)

        # unicode
        values = Series(["fooBAD__barBAD", np.nan, "foo"])

        result = values.str.extract(".*(BAD[_]+).*(BAD)", expand=False)
        exp = DataFrame([["BAD__", "BAD"], er, er])
        tm.assert_frame_equal(result, exp)

        # GH9980
        # Index only works with one regex group since
        # multi-group would expand to a frame
        idx = Index(["A1", "A2", "A3", "A4", "B5"])
        with pytest.raises(ValueError, match="supported"):
            idx.str.extract("([AB])([123])", expand=False)

        # these should work for both Series and Index
        for klass in [Series, Index]:
            # no groups
            s_or_idx = klass(["A1", "B2", "C3"])
            msg = "pattern contains no capture groups"
            with pytest.raises(ValueError, match=msg):
                s_or_idx.str.extract("[ABC][123]", expand=False)

            # only non-capturing groups
            with pytest.raises(ValueError, match=msg):
                s_or_idx.str.extract("(?:[AB]).*", expand=False)

            # single group renames series/index properly
            s_or_idx = klass(["A1", "A2"])
            result = s_or_idx.str.extract(r"(?P<uno>A)\d", expand=False)
            assert result.name == "uno"

            exp = klass(["A", "A"], name="uno")
            if klass == Series:
                tm.assert_series_equal(result, exp)
            else:
                tm.assert_index_equal(result, exp)

        s = Series(["A1", "B2", "C3"])
        # one group, no matches
        result = s.str.extract("(_)", expand=False)
        exp = Series([np.nan, np.nan, np.nan], dtype=object)
        tm.assert_series_equal(result, exp)

        # two groups, no matches
        result = s.str.extract("(_)(_)", expand=False)
        exp = DataFrame(
            [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]], dtype=object
        )
        tm.assert_frame_equal(result, exp)

        # one group, some matches
        result = s.str.extract("([AB])[123]", expand=False)
        exp = Series(["A", "B", np.nan])
        tm.assert_series_equal(result, exp)

        # two groups, some matches
        result = s.str.extract("([AB])([123])", expand=False)
        exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
        tm.assert_frame_equal(result, exp)

        # one named group
        result = s.str.extract("(?P<letter>[AB])", expand=False)
        exp = Series(["A", "B", np.nan], name="letter")
        tm.assert_series_equal(result, exp)

        # two named groups
        result = s.str.extract("(?P<letter>[AB])(?P<number>[123])", expand=False)
        exp = DataFrame(
            [["A", "1"], ["B", "2"], [np.nan, np.nan]], columns=["letter", "number"]
        )
        tm.assert_frame_equal(result, exp)

        # mix named and unnamed groups
        result = s.str.extract("([AB])(?P<number>[123])", expand=False)
        exp = DataFrame(
            [["A", "1"], ["B", "2"], [np.nan, np.nan]], columns=[0, "number"]
        )
        tm.assert_frame_equal(result, exp)

        # one normal group, one non-capturing group
        result = s.str.extract("([AB])(?:[123])", expand=False)
        exp = Series(["A", "B", np.nan])
        tm.assert_series_equal(result, exp)

        # two normal groups, one non-capturing group
        result = Series(["A11", "B22", "C33"]).str.extract(
            "([AB])([123])(?:[123])", expand=False
        )
        exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
        tm.assert_frame_equal(result, exp)

        # one optional group followed by one normal group
        result = Series(["A1", "B2", "3"]).str.extract(
            "(?P<letter>[AB])?(?P<number>[123])", expand=False
        )
        exp = DataFrame(
            [["A", "1"], ["B", "2"], [np.nan, "3"]], columns=["letter", "number"]
        )
        tm.assert_frame_equal(result, exp)

        # one normal group followed by one optional group
        result = Series(["A1", "B2", "C"]).str.extract(
            "(?P<letter>[ABC])(?P<number>[123])?", expand=False
        )
        exp = DataFrame(
            [["A", "1"], ["B", "2"], ["C", np.nan]], columns=["letter", "number"]
        )
        tm.assert_frame_equal(result, exp)

        # GH6348
        # not passing index to the extractor
        def check_index(index):
            data = ["A1", "B2", "C"]
            index = index[: len(data)]
            s = Series(data, index=index)
            result = s.str.extract(r"(\d)", expand=False)
            exp = Series(["1", "2", np.nan], index=index)
            tm.assert_series_equal(result, exp)

            result = Series(data, index=index).str.extract(
                r"(?P<letter>\D)(?P<number>\d)?", expand=False
            )
            e_list = [["A", "1"], ["B", "2"], ["C", np.nan]]
            exp = DataFrame(e_list, columns=["letter", "number"], index=index)
            tm.assert_frame_equal(result, exp)

        i_funs = [
            tm.makeStringIndex,
            tm.makeUnicodeIndex,
            tm.makeIntIndex,
            tm.makeDateIndex,
            tm.makePeriodIndex,
            tm.makeRangeIndex,
        ]
        for index in i_funs:
            check_index(index())

        # single_series_name_is_preserved.
        s = Series(["a3", "b3", "c2"], name="bob")
        r = s.str.extract(r"(?P<sue>[a-z])", expand=False)
        e = Series(["a", "b", "c"], name="sue")
        tm.assert_series_equal(r, e)
        assert r.name == e.name

    def test_extract_expand_True(self):
        # Contains tests like those in test_match and some others.
        values = Series(["fooBAD__barBAD", np.nan, "foo"])
        er = [np.nan, np.nan]  # empty row

        result = values.str.extract(".*(BAD[_]+).*(BAD)", expand=True)
        exp = DataFrame([["BAD__", "BAD"], er, er])
        tm.assert_frame_equal(result, exp)

        # mixed
        mixed = Series(
            [
                "aBAD_BAD",
                np.nan,
                "BAD_b_BAD",
                True,
                datetime.today(),
                "foo",
                None,
                1,
                2.0,
            ]
        )

        rs = Series(mixed).str.extract(".*(BAD[_]+).*(BAD)", expand=True)
        exp = DataFrame([["BAD_", "BAD"], er, ["BAD_", "BAD"], er, er, er, er, er, er])
        tm.assert_frame_equal(rs, exp)

        # these should work for both Series and Index
        for klass in [Series, Index]:
            # no groups
            s_or_idx = klass(["A1", "B2", "C3"])
            msg = "pattern contains no capture groups"
            with pytest.raises(ValueError, match=msg):
                s_or_idx.str.extract("[ABC][123]", expand=True)

            # only non-capturing groups
            with pytest.raises(ValueError, match=msg):
                s_or_idx.str.extract("(?:[AB]).*", expand=True)

            # single group renames series/index properly
            s_or_idx = klass(["A1", "A2"])
            result_df = s_or_idx.str.extract(r"(?P<uno>A)\d", expand=True)
            assert isinstance(result_df, DataFrame)
            result_series = result_df["uno"]
            tm.assert_series_equal(result_series, Series(["A", "A"], name="uno"))

    def test_extract_series(self):
        # extract should give the same result whether or not the
        # series has a name.
        for series_name in None, "series_name":
            s = Series(["A1", "B2", "C3"], name=series_name)
            # one group, no matches
            result = s.str.extract("(_)", expand=True)
            exp = DataFrame([np.nan, np.nan, np.nan], dtype=object)
            tm.assert_frame_equal(result, exp)

            # two groups, no matches
            result = s.str.extract("(_)(_)", expand=True)
            exp = DataFrame(
                [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]], dtype=object
            )
            tm.assert_frame_equal(result, exp)

            # one group, some matches
            result = s.str.extract("([AB])[123]", expand=True)
            exp = DataFrame(["A", "B", np.nan])
            tm.assert_frame_equal(result, exp)

            # two groups, some matches
            result = s.str.extract("([AB])([123])", expand=True)
            exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
            tm.assert_frame_equal(result, exp)

            # one named group
            result = s.str.extract("(?P<letter>[AB])", expand=True)
            exp = DataFrame({"letter": ["A", "B", np.nan]})
            tm.assert_frame_equal(result, exp)

            # two named groups
            result = s.str.extract("(?P<letter>[AB])(?P<number>[123])", expand=True)
            e_list = [["A", "1"], ["B", "2"], [np.nan, np.nan]]
            exp = DataFrame(e_list, columns=["letter", "number"])
            tm.assert_frame_equal(result, exp)

            # mix named and unnamed groups
            result = s.str.extract("([AB])(?P<number>[123])", expand=True)
            exp = DataFrame(e_list, columns=[0, "number"])
            tm.assert_frame_equal(result, exp)

            # one normal group, one non-capturing group
            result = s.str.extract("([AB])(?:[123])", expand=True)
            exp = DataFrame(["A", "B", np.nan])
            tm.assert_frame_equal(result, exp)

    def test_extract_optional_groups(self):

        # two normal groups, one non-capturing group
        result = Series(["A11", "B22", "C33"]).str.extract(
            "([AB])([123])(?:[123])", expand=True
        )
        exp = DataFrame([["A", "1"], ["B", "2"], [np.nan, np.nan]])
        tm.assert_frame_equal(result, exp)

        # one optional group followed by one normal group
        result = Series(["A1", "B2", "3"]).str.extract(
            "(?P<letter>[AB])?(?P<number>[123])", expand=True
        )
        e_list = [["A", "1"], ["B", "2"], [np.nan, "3"]]
        exp = DataFrame(e_list, columns=["letter", "number"])
        tm.assert_frame_equal(result, exp)

        # one normal group followed by one optional group
        result = Series(["A1", "B2", "C"]).str.extract(
            "(?P<letter>[ABC])(?P<number>[123])?", expand=True
        )
        e_list = [["A", "1"], ["B", "2"], ["C", np.nan]]
        exp = DataFrame(e_list, columns=["letter", "number"])
        tm.assert_frame_equal(result, exp)

        # GH6348
        # not passing index to the extractor
        def check_index(index):
            data = ["A1", "B2", "C"]
            index = index[: len(data)]
            result = Series(data, index=index).str.extract(r"(\d)", expand=True)
            exp = DataFrame(["1", "2", np.nan], index=index)
            tm.assert_frame_equal(result, exp)

            result = Series(data, index=index).str.extract(
                r"(?P<letter>\D)(?P<number>\d)?", expand=True
            )
            e_list = [["A", "1"], ["B", "2"], ["C", np.nan]]
            exp = DataFrame(e_list, columns=["letter", "number"], index=index)
            tm.assert_frame_equal(result, exp)

        i_funs = [
            tm.makeStringIndex,
            tm.makeUnicodeIndex,
            tm.makeIntIndex,
            tm.makeDateIndex,
            tm.makePeriodIndex,
            tm.makeRangeIndex,
        ]
        for index in i_funs:
            check_index(index())

    def test_extract_single_group_returns_frame(self):
        # GH11386 extract should always return DataFrame, even when
        # there is only one group. Prior to v0.18.0, extract returned
        # Series when there was only one group in the regex.
        s = Series(["a3", "b3", "c2"], name="series_name")
        r = s.str.extract(r"(?P<letter>[a-z])", expand=True)
        e = DataFrame({"letter": ["a", "b", "c"]})
        tm.assert_frame_equal(r, e)

    def test_extractall(self):
        subject_list = [
            "dave@google.com",
            "tdhock5@gmail.com",
            "maudelaperriere@gmail.com",
            "rob@gmail.com some text steve@gmail.com",
            "a@b.com some text c@d.com and e@f.com",
            np.nan,
            "",
        ]
        expected_tuples = [
            ("dave", "google", "com"),
            ("tdhock5", "gmail", "com"),
            ("maudelaperriere", "gmail", "com"),
            ("rob", "gmail", "com"),
            ("steve", "gmail", "com"),
            ("a", "b", "com"),
            ("c", "d", "com"),
            ("e", "f", "com"),
        ]
        named_pattern = r"""
        (?P<user>[a-z0-9]+)
        @
        (?P<domain>[a-z]+)
        \.
        (?P<tld>[a-z]{2,4})
        """
        expected_columns = ["user", "domain", "tld"]
        S = Series(subject_list)
        # extractall should return a DataFrame with one row for each
        # match, indexed by the subject from which the match came.
        expected_index = MultiIndex.from_tuples(
            [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (4, 0), (4, 1), (4, 2)],
            names=(None, "match"),
        )
        expected_df = DataFrame(expected_tuples, expected_index, expected_columns)
        computed_df = S.str.extractall(named_pattern, re.VERBOSE)
        tm.assert_frame_equal(computed_df, expected_df)

        # The index of the input Series should be used to construct
        # the index of the output DataFrame:
        series_index = MultiIndex.from_tuples(
            [
                ("single", "Dave"),
                ("single", "Toby"),
                ("single", "Maude"),
                ("multiple", "robAndSteve"),
                ("multiple", "abcdef"),
                ("none", "missing"),
                ("none", "empty"),
            ]
        )
        Si = Series(subject_list, series_index)
        expected_index = MultiIndex.from_tuples(
            [
                ("single", "Dave", 0),
                ("single", "Toby", 0),
                ("single", "Maude", 0),
                ("multiple", "robAndSteve", 0),
                ("multiple", "robAndSteve", 1),
                ("multiple", "abcdef", 0),
                ("multiple", "abcdef", 1),
                ("multiple", "abcdef", 2),
            ],
            names=(None, None, "match"),
        )
        expected_df = DataFrame(expected_tuples, expected_index, expected_columns)
        computed_df = Si.str.extractall(named_pattern, re.VERBOSE)
        tm.assert_frame_equal(computed_df, expected_df)

        # MultiIndexed subject with names.
        Sn = Series(subject_list, series_index)
        Sn.index.names = ("matches", "description")
        expected_index.names = ("matches", "description", "match")
        expected_df = DataFrame(expected_tuples, expected_index, expected_columns)
        computed_df = Sn.str.extractall(named_pattern, re.VERBOSE)
        tm.assert_frame_equal(computed_df, expected_df)

        # optional groups.
        subject_list = ["", "A1", "32"]
        named_pattern = "(?P<letter>[AB])?(?P<number>[123])"
        computed_df = Series(subject_list).str.extractall(named_pattern)
        expected_index = MultiIndex.from_tuples(
            [(1, 0), (2, 0), (2, 1)], names=(None, "match")
        )
        expected_df = DataFrame(
            [("A", "1"), (np.nan, "3"), (np.nan, "2")],
            expected_index,
            columns=["letter", "number"],
        )
        tm.assert_frame_equal(computed_df, expected_df)

        # only one of two groups has a name.
        pattern = "([AB])?(?P<number>[123])"
        computed_df = Series(subject_list).str.extractall(pattern)
        expected_df = DataFrame(
            [("A", "1"), (np.nan, "3"), (np.nan, "2")],
            expected_index,
            columns=[0, "number"],
        )
        tm.assert_frame_equal(computed_df, expected_df)

    def test_extractall_single_group(self):
        # extractall(one named group) returns DataFrame with one named
        # column.
        s = Series(["a3", "b3", "d4c2"], name="series_name")
        r = s.str.extractall(r"(?P<letter>[a-z])")
        i = MultiIndex.from_tuples(
            [(0, 0), (1, 0), (2, 0), (2, 1)], names=(None, "match")
        )
        e = DataFrame({"letter": ["a", "b", "d", "c"]}, i)
        tm.assert_frame_equal(r, e)

        # extractall(one un-named group) returns DataFrame with one
        # un-named column.
        r = s.str.extractall(r"([a-z])")
        e = DataFrame(["a", "b", "d", "c"], i)
        tm.assert_frame_equal(r, e)

    def test_extractall_single_group_with_quantifier(self):
        # extractall(one un-named group with quantifier) returns
        # DataFrame with one un-named column (GH13382).
        s = Series(["ab3", "abc3", "d4cd2"], name="series_name")
        r = s.str.extractall(r"([a-z]+)")
        i = MultiIndex.from_tuples(
            [(0, 0), (1, 0), (2, 0), (2, 1)], names=(None, "match")
        )
        e = DataFrame(["ab", "abc", "d", "cd"], i)
        tm.assert_frame_equal(r, e)

    @pytest.mark.parametrize(
        "data, names",
        [
            ([], (None,)),
            ([], ("i1",)),
            ([], (None, "i2")),
            ([], ("i1", "i2")),
            (["a3", "b3", "d4c2"], (None,)),
            (["a3", "b3", "d4c2"], ("i1", "i2")),
            (["a3", "b3", "d4c2"], (None, "i2")),
            (["a3", "b3", "d4c2"], ("i1", "i2")),
        ],
    )
    def test_extractall_no_matches(self, data, names):
        # GH19075 extractall with no matches should return a valid MultiIndex
        n = len(data)
        if len(names) == 1:
            i = Index(range(n), name=names[0])
        else:
            a = (tuple([i] * (n - 1)) for i in range(n))
            i = MultiIndex.from_tuples(a, names=names)
        s = Series(data, name="series_name", index=i, dtype="object")
        ei = MultiIndex.from_tuples([], names=(names + ("match",)))

        # one un-named group.
        r = s.str.extractall("(z)")
        e = DataFrame(columns=[0], index=ei)
        tm.assert_frame_equal(r, e)

        # two un-named groups.
        r = s.str.extractall("(z)(z)")
        e = DataFrame(columns=[0, 1], index=ei)
        tm.assert_frame_equal(r, e)

        # one named group.
        r = s.str.extractall("(?P<first>z)")
        e = DataFrame(columns=["first"], index=ei)
        tm.assert_frame_equal(r, e)

        # two named groups.
        r = s.str.extractall("(?P<first>z)(?P<second>z)")
        e = DataFrame(columns=["first", "second"], index=ei)
        tm.assert_frame_equal(r, e)

        # one named, one un-named.
        r = s.str.extractall("(z)(?P<second>z)")
        e = DataFrame(columns=[0, "second"], index=ei)
        tm.assert_frame_equal(r, e)

    def test_extractall_stringindex(self):
        s = Series(["a1a2", "b1", "c1"], name="xxx")
        res = s.str.extractall(r"[ab](?P<digit>\d)")
        exp_idx = MultiIndex.from_tuples(
            [(0, 0), (0, 1), (1, 0)], names=[None, "match"]
        )
        exp = DataFrame({"digit": ["1", "2", "1"]}, index=exp_idx)
        tm.assert_frame_equal(res, exp)

        # index should return the same result as the default index without name
        # thus index.name doesn't affect to the result
        for idx in [
            Index(["a1a2", "b1", "c1"]),
            Index(["a1a2", "b1", "c1"], name="xxx"),
        ]:

            res = idx.str.extractall(r"[ab](?P<digit>\d)")
            tm.assert_frame_equal(res, exp)

        s = Series(
            ["a1a2", "b1", "c1"],
            name="s_name",
            index=Index(["XX", "yy", "zz"], name="idx_name"),
        )
        res = s.str.extractall(r"[ab](?P<digit>\d)")
        exp_idx = MultiIndex.from_tuples(
            [("XX", 0), ("XX", 1), ("yy", 0)], names=["idx_name", "match"]
        )
        exp = DataFrame({"digit": ["1", "2", "1"]}, index=exp_idx)
        tm.assert_frame_equal(res, exp)

    def test_extractall_errors(self):
        # Does not make sense to use extractall with a regex that has
        # no capture groups. (it returns DataFrame with one column for
        # each capture group)
        s = Series(["a3", "b3", "d4c2"], name="series_name")
        with pytest.raises(ValueError, match="no capture groups"):
            s.str.extractall(r"[a-z]")

    def test_extract_index_one_two_groups(self):
        s = Series(["a3", "b3", "d4c2"], index=["A3", "B3", "D4"], name="series_name")
        r = s.index.str.extract(r"([A-Z])", expand=True)
        e = DataFrame(["A", "B", "D"])
        tm.assert_frame_equal(r, e)

        # Prior to v0.18.0, index.str.extract(regex with one group)
        # returned Index. With more than one group, extract raised an
        # error (GH9980). Now extract always returns DataFrame.
        r = s.index.str.extract(r"(?P<letter>[A-Z])(?P<digit>[0-9])", expand=True)
        e_list = [("A", "3"), ("B", "3"), ("D", "4")]
        e = DataFrame(e_list, columns=["letter", "digit"])
        tm.assert_frame_equal(r, e)

    def test_extractall_same_as_extract(self):
        s = Series(["a3", "b3", "c2"], name="series_name")

        pattern_two_noname = r"([a-z])([0-9])"
        extract_two_noname = s.str.extract(pattern_two_noname, expand=True)
        has_multi_index = s.str.extractall(pattern_two_noname)
        no_multi_index = has_multi_index.xs(0, level="match")
        tm.assert_frame_equal(extract_two_noname, no_multi_index)

        pattern_two_named = r"(?P<letter>[a-z])(?P<digit>[0-9])"
        extract_two_named = s.str.extract(pattern_two_named, expand=True)
        has_multi_index = s.str.extractall(pattern_two_named)
        no_multi_index = has_multi_index.xs(0, level="match")
        tm.assert_frame_equal(extract_two_named, no_multi_index)

        pattern_one_named = r"(?P<group_name>[a-z])"
        extract_one_named = s.str.extract(pattern_one_named, expand=True)
        has_multi_index = s.str.extractall(pattern_one_named)
        no_multi_index = has_multi_index.xs(0, level="match")
        tm.assert_frame_equal(extract_one_named, no_multi_index)

        pattern_one_noname = r"([a-z])"
        extract_one_noname = s.str.extract(pattern_one_noname, expand=True)
        has_multi_index = s.str.extractall(pattern_one_noname)
        no_multi_index = has_multi_index.xs(0, level="match")
        tm.assert_frame_equal(extract_one_noname, no_multi_index)

    def test_extractall_same_as_extract_subject_index(self):
        # same as above tests, but s has an MultiIndex.
        i = MultiIndex.from_tuples(
            [("A", "first"), ("B", "second"), ("C", "third")],
            names=("capital", "ordinal"),
        )
        s = Series(["a3", "b3", "c2"], i, name="series_name")

        pattern_two_noname = r"([a-z])([0-9])"
        extract_two_noname = s.str.extract(pattern_two_noname, expand=True)
        has_match_index = s.str.extractall(pattern_two_noname)
        no_match_index = has_match_index.xs(0, level="match")
        tm.assert_frame_equal(extract_two_noname, no_match_index)

        pattern_two_named = r"(?P<letter>[a-z])(?P<digit>[0-9])"
        extract_two_named = s.str.extract(pattern_two_named, expand=True)
        has_match_index = s.str.extractall(pattern_two_named)
        no_match_index = has_match_index.xs(0, level="match")
        tm.assert_frame_equal(extract_two_named, no_match_index)

        pattern_one_named = r"(?P<group_name>[a-z])"
        extract_one_named = s.str.extract(pattern_one_named, expand=True)
        has_match_index = s.str.extractall(pattern_one_named)
        no_match_index = has_match_index.xs(0, level="match")
        tm.assert_frame_equal(extract_one_named, no_match_index)

        pattern_one_noname = r"([a-z])"
        extract_one_noname = s.str.extract(pattern_one_noname, expand=True)
        has_match_index = s.str.extractall(pattern_one_noname)
        no_match_index = has_match_index.xs(0, level="match")
        tm.assert_frame_equal(extract_one_noname, no_match_index)

    def test_empty_str_methods(self):
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

    def test_empty_str_methods_to_frame(self):
        empty = Series(dtype=str)
        empty_df = DataFrame()
        tm.assert_frame_equal(empty_df, empty.str.partition("a"))
        tm.assert_frame_equal(empty_df, empty.str.rpartition("a"))

    def test_ismethods(self):
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

    def test_isnumeric(self):
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

    def test_get_dummies(self):
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

    def test_get_dummies_with_name_dummy(self):
        # GH 12180
        # Dummies named 'name' should work as expected
        s = Series(["a", "b,name", "b"])
        result = s.str.get_dummies(",")
        expected = DataFrame(
            [[1, 0, 0], [0, 1, 1], [0, 1, 0]], columns=["a", "b", "name"]
        )
        tm.assert_frame_equal(result, expected)

        idx = Index(["a|b", "name|c", "b|name"])
        result = idx.str.get_dummies("|")

        expected = MultiIndex.from_tuples(
            [(1, 1, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1)], names=("a", "b", "c", "name")
        )
        tm.assert_index_equal(result, expected)

    def test_join(self):
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

    def test_len(self):
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

    def test_findall(self):
        values = Series(["fooBAD__barBAD", np.nan, "foo", "BAD"])

        result = values.str.findall("BAD[_]*")
        exp = Series([["BAD__", "BAD"], np.nan, [], ["BAD"]])
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = Series(
            [
                "fooBAD__barBAD",
                np.nan,
                "foo",
                True,
                datetime.today(),
                "BAD",
                None,
                1,
                2.0,
            ]
        )

        rs = Series(mixed).str.findall("BAD[_]*")
        xp = Series(
            [
                ["BAD__", "BAD"],
                np.nan,
                [],
                np.nan,
                np.nan,
                ["BAD"],
                np.nan,
                np.nan,
                np.nan,
            ]
        )

        assert isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

    def test_find(self):
        values = Series(["ABCDEFG", "BCDEFEF", "DEFGHIJEF", "EFGHEF", "XXXX"])
        result = values.str.find("EF")
        tm.assert_series_equal(result, Series([4, 3, 1, 0, -1]))
        expected = np.array([v.find("EF") for v in values.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rfind("EF")
        tm.assert_series_equal(result, Series([4, 5, 7, 4, -1]))
        expected = np.array([v.rfind("EF") for v in values.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.find("EF", 3)
        tm.assert_series_equal(result, Series([4, 3, 7, 4, -1]))
        expected = np.array([v.find("EF", 3) for v in values.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rfind("EF", 3)
        tm.assert_series_equal(result, Series([4, 5, 7, 4, -1]))
        expected = np.array([v.rfind("EF", 3) for v in values.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.find("EF", 3, 6)
        tm.assert_series_equal(result, Series([4, 3, -1, 4, -1]))
        expected = np.array([v.find("EF", 3, 6) for v in values.values], dtype=np.int64)
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rfind("EF", 3, 6)
        tm.assert_series_equal(result, Series([4, 3, -1, 4, -1]))
        expected = np.array(
            [v.rfind("EF", 3, 6) for v in values.values], dtype=np.int64
        )
        tm.assert_numpy_array_equal(result.values, expected)

        with pytest.raises(TypeError, match="expected a string object, not int"):
            result = values.str.find(0)

        with pytest.raises(TypeError, match="expected a string object, not int"):
            result = values.str.rfind(0)

    def test_find_nan(self):
        values = Series(["ABCDEFG", np.nan, "DEFGHIJEF", np.nan, "XXXX"])
        result = values.str.find("EF")
        tm.assert_series_equal(result, Series([4, np.nan, 1, np.nan, -1]))

        result = values.str.rfind("EF")
        tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

        result = values.str.find("EF", 3)
        tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

        result = values.str.rfind("EF", 3)
        tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

        result = values.str.find("EF", 3, 6)
        tm.assert_series_equal(result, Series([4, np.nan, -1, np.nan, -1]))

        result = values.str.rfind("EF", 3, 6)
        tm.assert_series_equal(result, Series([4, np.nan, -1, np.nan, -1]))

    def test_index(self):
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

    def test_pad(self):
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

    def test_pad_fillchar(self):

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
    def test_pad_width(self, f):
        # see gh-13598
        s = Series(["1", "22", "a", "bb"])
        msg = "width must be of integer type, not*"

        with pytest.raises(TypeError, match=msg):
            getattr(s.str, f)("f")

    def test_translate(self):
        def _check(result, expected):
            if isinstance(result, Series):
                tm.assert_series_equal(result, expected)
            else:
                tm.assert_index_equal(result, expected)

        for klass in [Series, Index]:
            s = klass(["abcdefg", "abcc", "cdddfg", "cdefggg"])
            table = str.maketrans("abc", "cde")
            result = s.str.translate(table)
            expected = klass(["cdedefg", "cdee", "edddfg", "edefggg"])
            _check(result, expected)

        # Series with non-string values
        s = Series(["a", "b", "c", 1.2])
        expected = Series(["c", "d", "e", np.nan])
        result = s.str.translate(table)
        tm.assert_series_equal(result, expected)

    def test_center_ljust_rjust(self):
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
        mixed = Series(
            ["a", np.nan, "b", True, datetime.today(), "c", "eee", None, 1, 2.0]
        )

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

    def test_center_ljust_rjust_fillchar(self):
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

    def test_zfill(self):
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

    def test_split(self):
        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

        result = values.str.split("_")
        exp = Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(["a__b__c", "c__d__e", np.nan, "f__g__h"])
        result = values.str.split("__")
        tm.assert_series_equal(result, exp)

        result = values.str.split("__", expand=False)
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(["a_b_c", np.nan, "d_e_f", True, datetime.today(), None, 1, 2.0])
        result = mixed.str.split("_")
        exp = Series(
            [
                ["a", "b", "c"],
                np.nan,
                ["d", "e", "f"],
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ]
        )
        assert isinstance(result, Series)
        tm.assert_almost_equal(result, exp)

        result = mixed.str.split("_", expand=False)
        assert isinstance(result, Series)
        tm.assert_almost_equal(result, exp)

        # regex split
        values = Series(["a,b_c", "c_d,e", np.nan, "f,g,h"])
        result = values.str.split("[,_]")
        exp = Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
        tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize("dtype", [object, "string"])
    @pytest.mark.parametrize("method", ["split", "rsplit"])
    def test_split_n(self, dtype, method):
        s = Series(["a b", pd.NA, "b c"], dtype=dtype)
        expected = Series([["a", "b"], pd.NA, ["b", "c"]])

        result = getattr(s.str, method)(" ", n=None)
        tm.assert_series_equal(result, expected)

        result = getattr(s.str, method)(" ", n=0)
        tm.assert_series_equal(result, expected)

    def test_rsplit(self):
        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])
        result = values.str.rsplit("_")
        exp = Series([["a", "b", "c"], ["c", "d", "e"], np.nan, ["f", "g", "h"]])
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(["a__b__c", "c__d__e", np.nan, "f__g__h"])
        result = values.str.rsplit("__")
        tm.assert_series_equal(result, exp)

        result = values.str.rsplit("__", expand=False)
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(["a_b_c", np.nan, "d_e_f", True, datetime.today(), None, 1, 2.0])
        result = mixed.str.rsplit("_")
        exp = Series(
            [
                ["a", "b", "c"],
                np.nan,
                ["d", "e", "f"],
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ]
        )
        assert isinstance(result, Series)
        tm.assert_almost_equal(result, exp)

        result = mixed.str.rsplit("_", expand=False)
        assert isinstance(result, Series)
        tm.assert_almost_equal(result, exp)

        # regex split is not supported by rsplit
        values = Series(["a,b_c", "c_d,e", np.nan, "f,g,h"])
        result = values.str.rsplit("[,_]")
        exp = Series([["a,b_c"], ["c_d,e"], np.nan, ["f,g,h"]])
        tm.assert_series_equal(result, exp)

        # setting max number of splits, make sure it's from reverse
        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])
        result = values.str.rsplit("_", n=1)
        exp = Series([["a_b", "c"], ["c_d", "e"], np.nan, ["f_g", "h"]])
        tm.assert_series_equal(result, exp)

    def test_split_blank_string(self):
        # expand blank split GH 20067
        values = Series([""], name="test")
        result = values.str.split(expand=True)
        exp = DataFrame([[]])  # NOTE: this is NOT an empty DataFrame
        tm.assert_frame_equal(result, exp)

        values = Series(["a b c", "a b", "", " "], name="test")
        result = values.str.split(expand=True)
        exp = DataFrame(
            [
                ["a", "b", "c"],
                ["a", "b", np.nan],
                [np.nan, np.nan, np.nan],
                [np.nan, np.nan, np.nan],
            ]
        )
        tm.assert_frame_equal(result, exp)

    def test_split_noargs(self):
        # #1859
        s = Series(["Wes McKinney", "Travis  Oliphant"])
        result = s.str.split()
        expected = ["Travis", "Oliphant"]
        assert result[1] == expected
        result = s.str.rsplit()
        assert result[1] == expected

    def test_split_maxsplit(self):
        # re.split 0, str.split -1
        s = Series(["bd asdf jfg", "kjasdflqw asdfnfk"])

        result = s.str.split(n=-1)
        xp = s.str.split()
        tm.assert_series_equal(result, xp)

        result = s.str.split(n=0)
        tm.assert_series_equal(result, xp)

        xp = s.str.split("asdf")
        result = s.str.split("asdf", n=0)
        tm.assert_series_equal(result, xp)

        result = s.str.split("asdf", n=-1)
        tm.assert_series_equal(result, xp)

    def test_split_no_pat_with_nonzero_n(self):
        s = Series(["split once", "split once too!"])
        result = s.str.split(n=1)
        expected = Series({0: ["split", "once"], 1: ["split", "once too!"]})
        tm.assert_series_equal(expected, result, check_index_type=False)

    def test_split_to_dataframe(self):
        s = Series(["nosplit", "alsonosplit"])
        result = s.str.split("_", expand=True)
        exp = DataFrame({0: Series(["nosplit", "alsonosplit"])})
        tm.assert_frame_equal(result, exp)

        s = Series(["some_equal_splits", "with_no_nans"])
        result = s.str.split("_", expand=True)
        exp = DataFrame(
            {0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]}
        )
        tm.assert_frame_equal(result, exp)

        s = Series(["some_unequal_splits", "one_of_these_things_is_not"])
        result = s.str.split("_", expand=True)
        exp = DataFrame(
            {
                0: ["some", "one"],
                1: ["unequal", "of"],
                2: ["splits", "these"],
                3: [np.nan, "things"],
                4: [np.nan, "is"],
                5: [np.nan, "not"],
            }
        )
        tm.assert_frame_equal(result, exp)

        s = Series(["some_splits", "with_index"], index=["preserve", "me"])
        result = s.str.split("_", expand=True)
        exp = DataFrame(
            {0: ["some", "with"], 1: ["splits", "index"]}, index=["preserve", "me"]
        )
        tm.assert_frame_equal(result, exp)

        with pytest.raises(ValueError, match="expand must be"):
            s.str.split("_", expand="not_a_boolean")

    def test_split_to_multiindex_expand(self):
        # https://github.com/pandas-dev/pandas/issues/23677

        idx = Index(["nosplit", "alsonosplit", np.nan])
        result = idx.str.split("_", expand=True)
        exp = idx
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 1

        idx = Index(["some_equal_splits", "with_no_nans", np.nan, None])
        result = idx.str.split("_", expand=True)
        exp = MultiIndex.from_tuples(
            [
                ("some", "equal", "splits"),
                ("with", "no", "nans"),
                [np.nan, np.nan, np.nan],
                [None, None, None],
            ]
        )
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 3

        idx = Index(["some_unequal_splits", "one_of_these_things_is_not", np.nan, None])
        result = idx.str.split("_", expand=True)
        exp = MultiIndex.from_tuples(
            [
                ("some", "unequal", "splits", np.nan, np.nan, np.nan),
                ("one", "of", "these", "things", "is", "not"),
                (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
                (None, None, None, None, None, None),
            ]
        )
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 6

        with pytest.raises(ValueError, match="expand must be"):
            idx.str.split("_", expand="not_a_boolean")

    def test_rsplit_to_dataframe_expand(self):
        s = Series(["nosplit", "alsonosplit"])
        result = s.str.rsplit("_", expand=True)
        exp = DataFrame({0: Series(["nosplit", "alsonosplit"])})
        tm.assert_frame_equal(result, exp)

        s = Series(["some_equal_splits", "with_no_nans"])
        result = s.str.rsplit("_", expand=True)
        exp = DataFrame(
            {0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]}
        )
        tm.assert_frame_equal(result, exp)

        result = s.str.rsplit("_", expand=True, n=2)
        exp = DataFrame(
            {0: ["some", "with"], 1: ["equal", "no"], 2: ["splits", "nans"]}
        )
        tm.assert_frame_equal(result, exp)

        result = s.str.rsplit("_", expand=True, n=1)
        exp = DataFrame({0: ["some_equal", "with_no"], 1: ["splits", "nans"]})
        tm.assert_frame_equal(result, exp)

        s = Series(["some_splits", "with_index"], index=["preserve", "me"])
        result = s.str.rsplit("_", expand=True)
        exp = DataFrame(
            {0: ["some", "with"], 1: ["splits", "index"]}, index=["preserve", "me"]
        )
        tm.assert_frame_equal(result, exp)

    def test_rsplit_to_multiindex_expand(self):
        idx = Index(["nosplit", "alsonosplit"])
        result = idx.str.rsplit("_", expand=True)
        exp = idx
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 1

        idx = Index(["some_equal_splits", "with_no_nans"])
        result = idx.str.rsplit("_", expand=True)
        exp = MultiIndex.from_tuples(
            [("some", "equal", "splits"), ("with", "no", "nans")]
        )
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 3

        idx = Index(["some_equal_splits", "with_no_nans"])
        result = idx.str.rsplit("_", expand=True, n=1)
        exp = MultiIndex.from_tuples([("some_equal", "splits"), ("with_no", "nans")])
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 2

    def test_split_nan_expand(self):
        # gh-18450
        s = Series(["foo,bar,baz", np.nan])
        result = s.str.split(",", expand=True)
        exp = DataFrame([["foo", "bar", "baz"], [np.nan, np.nan, np.nan]])
        tm.assert_frame_equal(result, exp)

        # check that these are actually np.nan and not None
        # TODO see GH 18463
        # tm.assert_frame_equal does not differentiate
        assert all(np.isnan(x) for x in result.iloc[1])

    def test_split_with_name(self):
        # GH 12617

        # should preserve name
        s = Series(["a,b", "c,d"], name="xxx")
        res = s.str.split(",")
        exp = Series([["a", "b"], ["c", "d"]], name="xxx")
        tm.assert_series_equal(res, exp)

        res = s.str.split(",", expand=True)
        exp = DataFrame([["a", "b"], ["c", "d"]])
        tm.assert_frame_equal(res, exp)

        idx = Index(["a,b", "c,d"], name="xxx")
        res = idx.str.split(",")
        exp = Index([["a", "b"], ["c", "d"]], name="xxx")
        assert res.nlevels == 1
        tm.assert_index_equal(res, exp)

        res = idx.str.split(",", expand=True)
        exp = MultiIndex.from_tuples([("a", "b"), ("c", "d")])
        assert res.nlevels == 2
        tm.assert_index_equal(res, exp)

    def test_partition_series(self):
        # https://github.com/pandas-dev/pandas/issues/23558

        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None])

        result = values.str.partition("_", expand=False)
        exp = Series(
            [("a", "_", "b_c"), ("c", "_", "d_e"), np.nan, ("f", "_", "g_h"), None]
        )
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition("_", expand=False)
        exp = Series(
            [("a_b", "_", "c"), ("c_d", "_", "e"), np.nan, ("f_g", "_", "h"), None]
        )
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(["a__b__c", "c__d__e", np.nan, "f__g__h", None])
        result = values.str.partition("__", expand=False)
        exp = Series(
            [
                ("a", "__", "b__c"),
                ("c", "__", "d__e"),
                np.nan,
                ("f", "__", "g__h"),
                None,
            ]
        )
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition("__", expand=False)
        exp = Series(
            [
                ("a__b", "__", "c"),
                ("c__d", "__", "e"),
                np.nan,
                ("f__g", "__", "h"),
                None,
            ]
        )
        tm.assert_series_equal(result, exp)

        # None
        values = Series(["a b c", "c d e", np.nan, "f g h", None])
        result = values.str.partition(expand=False)
        exp = Series(
            [("a", " ", "b c"), ("c", " ", "d e"), np.nan, ("f", " ", "g h"), None]
        )
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition(expand=False)
        exp = Series(
            [("a b", " ", "c"), ("c d", " ", "e"), np.nan, ("f g", " ", "h"), None]
        )
        tm.assert_series_equal(result, exp)

        # Not split
        values = Series(["abc", "cde", np.nan, "fgh", None])
        result = values.str.partition("_", expand=False)
        exp = Series([("abc", "", ""), ("cde", "", ""), np.nan, ("fgh", "", ""), None])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition("_", expand=False)
        exp = Series([("", "", "abc"), ("", "", "cde"), np.nan, ("", "", "fgh"), None])
        tm.assert_series_equal(result, exp)

        # unicode
        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

        result = values.str.partition("_", expand=False)
        exp = Series([("a", "_", "b_c"), ("c", "_", "d_e"), np.nan, ("f", "_", "g_h")])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition("_", expand=False)
        exp = Series([("a_b", "_", "c"), ("c_d", "_", "e"), np.nan, ("f_g", "_", "h")])
        tm.assert_series_equal(result, exp)

        # compare to standard lib
        values = Series(["A_B_C", "B_C_D", "E_F_G", "EFGHEF"])
        result = values.str.partition("_", expand=False).tolist()
        assert result == [v.partition("_") for v in values]
        result = values.str.rpartition("_", expand=False).tolist()
        assert result == [v.rpartition("_") for v in values]

    def test_partition_index(self):
        # https://github.com/pandas-dev/pandas/issues/23558

        values = Index(["a_b_c", "c_d_e", "f_g_h", np.nan, None])

        result = values.str.partition("_", expand=False)
        exp = Index(
            np.array(
                [("a", "_", "b_c"), ("c", "_", "d_e"), ("f", "_", "g_h"), np.nan, None],
                dtype=object,
            )
        )
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 1

        result = values.str.rpartition("_", expand=False)
        exp = Index(
            np.array(
                [("a_b", "_", "c"), ("c_d", "_", "e"), ("f_g", "_", "h"), np.nan, None],
                dtype=object,
            )
        )
        tm.assert_index_equal(result, exp)
        assert result.nlevels == 1

        result = values.str.partition("_")
        exp = Index(
            [
                ("a", "_", "b_c"),
                ("c", "_", "d_e"),
                ("f", "_", "g_h"),
                (np.nan, np.nan, np.nan),
                (None, None, None),
            ]
        )
        tm.assert_index_equal(result, exp)
        assert isinstance(result, MultiIndex)
        assert result.nlevels == 3

        result = values.str.rpartition("_")
        exp = Index(
            [
                ("a_b", "_", "c"),
                ("c_d", "_", "e"),
                ("f_g", "_", "h"),
                (np.nan, np.nan, np.nan),
                (None, None, None),
            ]
        )
        tm.assert_index_equal(result, exp)
        assert isinstance(result, MultiIndex)
        assert result.nlevels == 3

    def test_partition_to_dataframe(self):
        # https://github.com/pandas-dev/pandas/issues/23558

        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None])
        result = values.str.partition("_")
        exp = DataFrame(
            {
                0: ["a", "c", np.nan, "f", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["b_c", "d_e", np.nan, "g_h", None],
            }
        )
        tm.assert_frame_equal(result, exp)

        result = values.str.rpartition("_")
        exp = DataFrame(
            {
                0: ["a_b", "c_d", np.nan, "f_g", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["c", "e", np.nan, "h", None],
            }
        )
        tm.assert_frame_equal(result, exp)

        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h", None])
        result = values.str.partition("_", expand=True)
        exp = DataFrame(
            {
                0: ["a", "c", np.nan, "f", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["b_c", "d_e", np.nan, "g_h", None],
            }
        )
        tm.assert_frame_equal(result, exp)

        result = values.str.rpartition("_", expand=True)
        exp = DataFrame(
            {
                0: ["a_b", "c_d", np.nan, "f_g", None],
                1: ["_", "_", np.nan, "_", None],
                2: ["c", "e", np.nan, "h", None],
            }
        )
        tm.assert_frame_equal(result, exp)

    def test_partition_with_name(self):
        # GH 12617

        s = Series(["a,b", "c,d"], name="xxx")
        res = s.str.partition(",")
        exp = DataFrame({0: ["a", "c"], 1: [",", ","], 2: ["b", "d"]})
        tm.assert_frame_equal(res, exp)

        # should preserve name
        res = s.str.partition(",", expand=False)
        exp = Series([("a", ",", "b"), ("c", ",", "d")], name="xxx")
        tm.assert_series_equal(res, exp)

        idx = Index(["a,b", "c,d"], name="xxx")
        res = idx.str.partition(",")
        exp = MultiIndex.from_tuples([("a", ",", "b"), ("c", ",", "d")])
        assert res.nlevels == 3
        tm.assert_index_equal(res, exp)

        # should preserve name
        res = idx.str.partition(",", expand=False)
        exp = Index(np.array([("a", ",", "b"), ("c", ",", "d")]), name="xxx")
        assert res.nlevels == 1
        tm.assert_index_equal(res, exp)

    def test_partition_sep_kwarg(self):
        # GH 22676; depr kwarg "pat" in favor of "sep"
        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

        expected = values.str.partition(sep="_")
        result = values.str.partition("_")
        tm.assert_frame_equal(result, expected)

        expected = values.str.rpartition(sep="_")
        result = values.str.rpartition("_")
        tm.assert_frame_equal(result, expected)

    def test_pipe_failures(self):
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
    def test_slice(self, start, stop, step, expected):
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

    def test_slice_replace(self):
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

    def test_strip_lstrip_rstrip(self):
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

    def test_strip_lstrip_rstrip_mixed(self):
        # mixed
        mixed = Series(
            ["  aa  ", np.nan, " bb \t\n", True, datetime.today(), None, 1, 2.0]
        )

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

    def test_strip_lstrip_rstrip_args(self):
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

    def test_wrap(self):
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

    def test_get(self):
        values = Series(["a_b_c", "c_d_e", np.nan, "f_g_h"])

        result = values.str.split("_").str.get(1)
        expected = Series(["b", "d", np.nan, "g"])
        tm.assert_series_equal(result, expected)

        # mixed
        mixed = Series(["a_b_c", np.nan, "c_d_e", True, datetime.today(), None, 1, 2.0])

        rs = Series(mixed).str.split("_").str.get(1)
        xp = Series(["b", np.nan, "d", np.nan, np.nan, np.nan, np.nan, np.nan])

        assert isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # bounds testing
        values = Series(["1_2_3_4_5", "6_7_8_9_10", "11_12"])

        # positive index
        result = values.str.split("_").str.get(2)
        expected = Series(["3", "8", np.nan])
        tm.assert_series_equal(result, expected)

        # negative index
        result = values.str.split("_").str.get(-3)
        expected = Series(["3", "8", np.nan])
        tm.assert_series_equal(result, expected)

    def test_get_complex(self):
        # GH 20671, getting value not in dict raising `KeyError`
        values = Series([(1, 2, 3), [1, 2, 3], {1, 2, 3}, {1: "a", 2: "b", 3: "c"}])

        result = values.str.get(1)
        expected = Series([2, 2, np.nan, "a"])
        tm.assert_series_equal(result, expected)

        result = values.str.get(-1)
        expected = Series([3, 3, np.nan, np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("to_type", [tuple, list, np.array])
    def test_get_complex_nested(self, to_type):
        values = Series([to_type([to_type([1, 2])])])

        result = values.str.get(0)
        expected = Series([to_type([1, 2])])
        tm.assert_series_equal(result, expected)

        result = values.str.get(1)
        expected = Series([np.nan])
        tm.assert_series_equal(result, expected)

    def test_contains_moar(self):
        # PR #1179
        s = Series(["A", "B", "C", "Aaba", "Baca", "", np.nan, "CABA", "dog", "cat"])

        result = s.str.contains("a")
        expected = Series(
            [False, False, False, True, True, False, np.nan, False, False, True]
        )
        tm.assert_series_equal(result, expected)

        result = s.str.contains("a", case=False)
        expected = Series(
            [True, False, False, True, True, False, np.nan, True, False, True]
        )
        tm.assert_series_equal(result, expected)

        result = s.str.contains("Aa")
        expected = Series(
            [False, False, False, True, False, False, np.nan, False, False, False]
        )
        tm.assert_series_equal(result, expected)

        result = s.str.contains("ba")
        expected = Series(
            [False, False, False, True, False, False, np.nan, False, False, False]
        )
        tm.assert_series_equal(result, expected)

        result = s.str.contains("ba", case=False)
        expected = Series(
            [False, False, False, True, True, False, np.nan, True, False, False]
        )
        tm.assert_series_equal(result, expected)

    def test_contains_nan(self):
        # PR #14171
        s = Series([np.nan, np.nan, np.nan], dtype=np.object_)

        result = s.str.contains("foo", na=False)
        expected = Series([False, False, False], dtype=np.bool_)
        tm.assert_series_equal(result, expected)

        result = s.str.contains("foo", na=True)
        expected = Series([True, True, True], dtype=np.bool_)
        tm.assert_series_equal(result, expected)

        result = s.str.contains("foo", na="foo")
        expected = Series(["foo", "foo", "foo"], dtype=np.object_)
        tm.assert_series_equal(result, expected)

        result = s.str.contains("foo")
        expected = Series([np.nan, np.nan, np.nan], dtype=np.object_)
        tm.assert_series_equal(result, expected)

    def test_replace_moar(self):
        # PR #1179
        s = Series(["A", "B", "C", "Aaba", "Baca", "", np.nan, "CABA", "dog", "cat"])

        result = s.str.replace("A", "YYY")
        expected = Series(
            ["YYY", "B", "C", "YYYaba", "Baca", "", np.nan, "CYYYBYYY", "dog", "cat"]
        )
        tm.assert_series_equal(result, expected)

        result = s.str.replace("A", "YYY", case=False)
        expected = Series(
            [
                "YYY",
                "B",
                "C",
                "YYYYYYbYYY",
                "BYYYcYYY",
                "",
                np.nan,
                "CYYYBYYY",
                "dog",
                "cYYYt",
            ]
        )
        tm.assert_series_equal(result, expected)

        result = s.str.replace("^.a|dog", "XX-XX ", case=False, regex=True)
        expected = Series(
            [
                "A",
                "B",
                "C",
                "XX-XX ba",
                "XX-XX ca",
                "",
                np.nan,
                "XX-XX BA",
                "XX-XX ",
                "XX-XX t",
            ]
        )
        tm.assert_series_equal(result, expected)

    def test_string_slice_get_syntax(self):
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

    def test_string_slice_out_of_bounds(self):
        s = Series([(1, 2), (1,), (3, 4, 5)])

        result = s.str[1]
        expected = Series([2, np.nan, 4])

        tm.assert_series_equal(result, expected)

        s = Series(["foo", "b", "ba"])
        result = s.str[1]
        expected = Series(["o", np.nan, "a"])
        tm.assert_series_equal(result, expected)

    def test_match_findall_flags(self):
        data = {
            "Dave": "dave@google.com",
            "Steve": "steve@gmail.com",
            "Rob": "rob@gmail.com",
            "Wes": np.nan,
        }
        data = Series(data)

        pat = r"([A-Z0-9._%+-]+)@([A-Z0-9.-]+)\.([A-Z]{2,4})"

        result = data.str.extract(pat, flags=re.IGNORECASE, expand=True)
        assert result.iloc[0].tolist() == ["dave", "google", "com"]

        result = data.str.match(pat, flags=re.IGNORECASE)
        assert result[0]

        result = data.str.fullmatch(pat, flags=re.IGNORECASE)
        assert result[0]

        result = data.str.findall(pat, flags=re.IGNORECASE)
        assert result[0][0] == ("dave", "google", "com")

        result = data.str.count(pat, flags=re.IGNORECASE)
        assert result[0] == 1

        with tm.assert_produces_warning(UserWarning):
            result = data.str.contains(pat, flags=re.IGNORECASE)
        assert result[0]

    def test_encode_decode(self):
        base = Series(["a", "b", "a\xe4"])
        series = base.str.encode("utf-8")

        f = lambda x: x.decode("utf-8")
        result = series.str.decode("utf-8")
        exp = series.map(f)

        tm.assert_series_equal(result, exp)

    def test_encode_decode_errors(self):
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

    def test_normalize(self):
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

    def test_index_str_accessor_visibility(self):
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

    def test_str_accessor_no_new_attributes(self):
        # https://github.com/pandas-dev/pandas/issues/10673
        s = Series(list("aabbcde"))
        with pytest.raises(AttributeError, match="You cannot add any new attribute"):
            s.str.xlabel = "a"

    def test_method_on_bytes(self):
        lhs = Series(np.array(list("abc"), "S1").astype(object))
        rhs = Series(np.array(list("def"), "S1").astype(object))
        with pytest.raises(TypeError, match="Cannot use .str.cat with values of.*"):
            lhs.str.cat(rhs)

    def test_casefold(self):
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


@pytest.mark.parametrize("klass", [tuple, list, np.array, pd.Series, pd.Index])
def test_cat_different_classes(klass):
    # https://github.com/pandas-dev/pandas/issues/33425
    s = Series(["a", "b", "c"])
    result = s.str.cat(klass(["x", "y", "z"]))
    expected = Series(["ax", "by", "cz"])
    tm.assert_series_equal(result, expected)


def test_str_get_stringarray_multiple_nans():
    s = Series(pd.array(["a", "ab", pd.NA, "abc"]))
    result = s.str.get(2)
    expected = Series(pd.array([pd.NA, pd.NA, pd.NA, "c"]))
    tm.assert_series_equal(result, expected)
