import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series
from pandas.core.util.hashing import _hash_scalar, hash_tuple, hash_tuples
from pandas.util import hash_array, hash_pandas_object
import pandas.util.testing as tm


@pytest.fixture(
    params=[
        Series([1, 2, 3] * 3, dtype="int32"),
        Series([None, 2.5, 3.5] * 3, dtype="float32"),
        Series(["a", "b", "c"] * 3, dtype="category"),
        Series(["d", "e", "f"] * 3),
        Series([True, False, True] * 3),
        Series(pd.date_range("20130101", periods=9)),
        Series(pd.date_range("20130101", periods=9, tz="US/Eastern")),
        Series(pd.timedelta_range("2000", periods=9)),
    ]
)
def series(request):
    return request.param


@pytest.fixture(params=[True, False])
def index(request):
    return request.param


def _check_equal(obj, **kwargs):
    """
    Check that hashing an objects produces the same value each time.

    Parameters
    ----------
    obj : object
        The object to hash.
    kwargs : kwargs
        Keyword arguments to pass to the hashing function.
    """
    a = hash_pandas_object(obj, **kwargs)
    b = hash_pandas_object(obj, **kwargs)
    tm.assert_series_equal(a, b)


def _check_not_equal_with_index(obj):
    """
    Check the hash of an object with and without its index is not the same.

    Parameters
    ----------
    obj : object
        The object to hash.
    """
    if not isinstance(obj, Index):
        a = hash_pandas_object(obj, index=True)
        b = hash_pandas_object(obj, index=False)

        if len(obj):
            assert not (a == b).all()


def test_consistency():
    # Check that our hash doesn't change because of a mistake
    # in the actual code; this is the ground truth.
    result = hash_pandas_object(Index(["foo", "bar", "baz"]))
    expected = Series(
        np.array(
            [3600424527151052760, 1374399572096150070, 477881037637427054],
            dtype="uint64",
        ),
        index=["foo", "bar", "baz"],
    )
    tm.assert_series_equal(result, expected)


def test_hash_array(series):
    arr = series.values
    tm.assert_numpy_array_equal(hash_array(arr), hash_array(arr))


@pytest.mark.parametrize(
    "arr2", [np.array([3, 4, "All"]), np.array([3, 4, "All"], dtype=object)]
)
def test_hash_array_mixed(arr2):
    result1 = hash_array(np.array(["3", "4", "All"]))
    result2 = hash_array(arr2)

    tm.assert_numpy_array_equal(result1, result2)


@pytest.mark.parametrize("val", [5, "foo", pd.Timestamp("20130101")])
def test_hash_array_errors(val):
    msg = "must pass a ndarray-like"
    with pytest.raises(TypeError, match=msg):
        hash_array(val)


def test_hash_tuples():
    tuples = [(1, "one"), (1, "two"), (2, "one")]
    result = hash_tuples(tuples)

    expected = hash_pandas_object(MultiIndex.from_tuples(tuples)).values
    tm.assert_numpy_array_equal(result, expected)

    result = hash_tuples(tuples[0])
    assert result == expected[0]


@pytest.mark.parametrize(
    "tup",
    [(1, "one"), (1, np.nan), (1.0, pd.NaT, "A"), ("A", pd.Timestamp("2012-01-01"))],
)
def test_hash_tuple(tup):
    # Test equivalence between
    # hash_tuples and hash_tuple.
    result = hash_tuple(tup)
    expected = hash_tuples([tup])[0]

    assert result == expected


@pytest.mark.parametrize(
    "val",
    [
        1,
        1.4,
        "A",
        b"A",
        pd.Timestamp("2012-01-01"),
        pd.Timestamp("2012-01-01", tz="Europe/Brussels"),
        datetime.datetime(2012, 1, 1),
        pd.Timestamp("2012-01-01", tz="EST").to_pydatetime(),
        pd.Timedelta("1 days"),
        datetime.timedelta(1),
        pd.Period("2012-01-01", freq="D"),
        pd.Interval(0, 1),
        np.nan,
        pd.NaT,
        None,
    ],
)
def test_hash_scalar(val):
    result = _hash_scalar(val)
    expected = hash_array(np.array([val], dtype=object), categorize=True)

    assert result[0] == expected[0]


@pytest.mark.parametrize("val", [5, "foo", pd.Timestamp("20130101")])
def test_hash_tuples_err(val):
    msg = "must be convertible to a list-of-tuples"
    with pytest.raises(TypeError, match=msg):
        hash_tuples(val)


def test_multiindex_unique():
    mi = MultiIndex.from_tuples([(118, 472), (236, 118), (51, 204), (102, 51)])
    assert mi.is_unique is True

    result = hash_pandas_object(mi)
    assert result.is_unique is True


def test_multiindex_objects():
    mi = MultiIndex(
        levels=[["b", "d", "a"], [1, 2, 3]],
        codes=[[0, 1, 0, 2], [2, 0, 0, 1]],
        names=["col1", "col2"],
    )
    recons = mi._sort_levels_monotonic()

    # These are equal.
    assert mi.equals(recons)
    assert Index(mi.values).equals(Index(recons.values))

    # _hashed_values and hash_pandas_object(..., index=False) equivalency.
    expected = hash_pandas_object(mi, index=False).values
    result = mi._hashed_values

    tm.assert_numpy_array_equal(result, expected)

    expected = hash_pandas_object(recons, index=False).values
    result = recons._hashed_values

    tm.assert_numpy_array_equal(result, expected)

    expected = mi._hashed_values
    result = recons._hashed_values

    # Values should match, but in different order.
    tm.assert_numpy_array_equal(np.sort(result), np.sort(expected))


@pytest.mark.parametrize(
    "obj",
    [
        Series([1, 2, 3]),
        Series([1.0, 1.5, 3.2]),
        Series([1.0, 1.5, np.nan]),
        Series([1.0, 1.5, 3.2], index=[1.5, 1.1, 3.3]),
        Series(["a", "b", "c"]),
        Series(["a", np.nan, "c"]),
        Series(["a", None, "c"]),
        Series([True, False, True]),
        Series(dtype=object),
        Index([1, 2, 3]),
        Index([True, False, True]),
        DataFrame({"x": ["a", "b", "c"], "y": [1, 2, 3]}),
        DataFrame(),
        tm.makeMissingDataframe(),
        tm.makeMixedDataFrame(),
        tm.makeTimeDataFrame(),
        tm.makeTimeSeries(),
        tm.makeTimedeltaIndex(),
        tm.makePeriodIndex(),
        Series(tm.makePeriodIndex()),
        Series(pd.date_range("20130101", periods=3, tz="US/Eastern")),
        MultiIndex.from_product(
            [range(5), ["foo", "bar", "baz"], pd.date_range("20130101", periods=2)]
        ),
        MultiIndex.from_product([pd.CategoricalIndex(list("aabc")), range(3)]),
    ],
)
def test_hash_pandas_object(obj, index):
    _check_equal(obj, index=index)
    _check_not_equal_with_index(obj)


def test_hash_pandas_object2(series, index):
    _check_equal(series, index=index)
    _check_not_equal_with_index(series)


@pytest.mark.parametrize(
    "obj", [Series([], dtype="float64"), Series([], dtype="object"), Index([])]
)
def test_hash_pandas_empty_object(obj, index):
    # These are by-definition the same with
    # or without the index as the data is empty.
    _check_equal(obj, index=index)


@pytest.mark.parametrize(
    "s1",
    [
        Series(["a", "b", "c", "d"]),
        Series([1000, 2000, 3000, 4000]),
        Series(pd.date_range(0, periods=4)),
    ],
)
@pytest.mark.parametrize("categorize", [True, False])
def test_categorical_consistency(s1, categorize):
    # see gh-15143
    #
    # Check that categoricals hash consistent with their values,
    # not codes. This should work for categoricals of any dtype.
    s2 = s1.astype("category").cat.set_categories(s1)
    s3 = s2.cat.set_categories(list(reversed(s1)))

    # These should all hash identically.
    h1 = hash_pandas_object(s1, categorize=categorize)
    h2 = hash_pandas_object(s2, categorize=categorize)
    h3 = hash_pandas_object(s3, categorize=categorize)

    tm.assert_series_equal(h1, h2)
    tm.assert_series_equal(h1, h3)


def test_categorical_with_nan_consistency():
    c = pd.Categorical.from_codes(
        [-1, 0, 1, 2, 3, 4], categories=pd.date_range("2012-01-01", periods=5, name="B")
    )
    expected = hash_array(c, categorize=False)

    c = pd.Categorical.from_codes([-1, 0], categories=[pd.Timestamp("2012-01-01")])
    result = hash_array(c, categorize=False)

    assert result[0] in expected
    assert result[1] in expected


@pytest.mark.parametrize("obj", [pd.Timestamp("20130101")])
def test_pandas_errors(obj):
    msg = "Unexpected type for hashing"
    with pytest.raises(TypeError, match=msg):
        hash_pandas_object(obj)


def test_hash_keys():
    # Using different hash keys, should have
    # different hashes for the same data.
    #
    # This only matters for object dtypes.
    obj = Series(list("abc"))

    a = hash_pandas_object(obj, hash_key="9876543210123456")
    b = hash_pandas_object(obj, hash_key="9876543210123465")

    assert (a != b).all()


def test_invalid_key():
    # This only matters for object dtypes.
    msg = "key should be a 16-byte string encoded"

    with pytest.raises(ValueError, match=msg):
        hash_pandas_object(Series(list("abc")), hash_key="foo")


def test_already_encoded(index):
    # If already encoded, then ok.
    obj = Series(list("abc")).str.encode("utf8")
    _check_equal(obj, index=index)


def test_alternate_encoding(index):
    obj = Series(list("abc"))
    _check_equal(obj, index=index, encoding="ascii")


@pytest.mark.parametrize("l_exp", range(8))
@pytest.mark.parametrize("l_add", [0, 1])
def test_same_len_hash_collisions(l_exp, l_add):
    length = 2 ** (l_exp + 8) + l_add
    s = tm.rands_array(length, 2)

    result = hash_array(s, "utf8")
    assert not result[0] == result[1]


def test_hash_collisions():
    # Hash collisions are bad.
    #
    # https://github.com/pandas-dev/pandas/issues/14711#issuecomment-264885726
    hashes = [
        "Ingrid-9Z9fKIZmkO7i7Cn51Li34pJm44fgX6DYGBNj3VPlOH50m7HnBlPxfIwFMrcNJNMP6PSgLmwWnInciMWrCSAlLEvt7JkJl4IxiMrVbXSa8ZQoVaq5xoQPjltuJEfwdNlO6jo8qRRHvD8sBEBMQASrRa6TsdaPTPCBo3nwIBpE7YzzmyH0vMBhjQZLx1aCT7faSEx7PgFxQhHdKFWROcysamgy9iVj8DO2Fmwg1NNl93rIAqC3mdqfrCxrzfvIY8aJdzin2cHVzy3QUJxZgHvtUtOLxoqnUHsYbNTeq0xcLXpTZEZCxD4PGubIuCNf32c33M7HFsnjWSEjE2yVdWKhmSVodyF8hFYVmhYnMCztQnJrt3O8ZvVRXd5IKwlLexiSp4h888w7SzAIcKgc3g5XQJf6MlSMftDXm9lIsE1mJNiJEv6uY6pgvC3fUPhatlR5JPpVAHNSbSEE73MBzJrhCAbOLXQumyOXigZuPoME7QgJcBalliQol7YZ9",  # noqa: E501
        "Tim-b9MddTxOWW2AT1Py6vtVbZwGAmYCjbp89p8mxsiFoVX4FyDOF3wFiAkyQTUgwg9sVqVYOZo09Dh1AzhFHbgij52ylF0SEwgzjzHH8TGY8Lypart4p4onnDoDvVMBa0kdthVGKl6K0BDVGzyOXPXKpmnMF1H6rJzqHJ0HywfwS4XYpVwlAkoeNsiicHkJUFdUAhG229INzvIAiJuAHeJDUoyO4DCBqtoZ5TDend6TK7Y914yHlfH3g1WZu5LksKv68VQHJriWFYusW5e6ZZ6dKaMjTwEGuRgdT66iU5nqWTHRH8WSzpXoCFwGcTOwyuqPSe0fTe21DVtJn1FKj9F9nEnR9xOvJUO7E0piCIF4Ad9yAIDY4DBimpsTfKXCu1vdHpKYerzbndfuFe5AhfMduLYZJi5iAw8qKSwR5h86ttXV0Mc0QmXz8dsRvDgxjXSmupPxBggdlqUlC828hXiTPD7am0yETBV0F3bEtvPiNJfremszcV8NcqAoARMe",  # noqa: E501
    ]

    # These should be different.
    result1 = hash_array(np.asarray(hashes[0:1], dtype=object), "utf8")
    expected1 = np.array([14963968704024874985], dtype=np.uint64)
    tm.assert_numpy_array_equal(result1, expected1)

    result2 = hash_array(np.asarray(hashes[1:2], dtype=object), "utf8")
    expected2 = np.array([16428432627716348016], dtype=np.uint64)
    tm.assert_numpy_array_equal(result2, expected2)

    result = hash_array(np.asarray(hashes, dtype=object), "utf8")
    tm.assert_numpy_array_equal(result, np.concatenate([expected1, expected2], axis=0))
