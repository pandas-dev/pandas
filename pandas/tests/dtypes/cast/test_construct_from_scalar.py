import numpy as np

from pandas.core.dtypes.cast import construct_1d_arraylike_from_scalar
from pandas.core.dtypes.dtypes import CategoricalDtype

from pandas import Categorical, Timedelta, Timestamp
import pandas._testing as tm


def test_cast_1d_array_like_from_scalar_categorical():
    # see gh-19565
    #
    # Categorical result from scalar did not maintain
    # categories and ordering of the passed dtype.
    cats = ["a", "b", "c"]
    cat_type = CategoricalDtype(categories=cats, ordered=False)
    expected = Categorical(["a", "a"], categories=cats)

    result = construct_1d_arraylike_from_scalar("a", len(expected), cat_type)
    tm.assert_categorical_equal(result, expected)


def test_cast_1d_array_like_from_timestamp():
    # check we dont lose nanoseconds
    ts = Timestamp.now() + Timedelta(1)
    res = construct_1d_arraylike_from_scalar(ts, 2, np.dtype("M8[ns]"))
    assert res[0] == ts


def test_cast_1d_array_like_from_timedelta():
    # check we dont lose nanoseconds
    td = Timedelta(1)
    res = construct_1d_arraylike_from_scalar(td, 2, np.dtype("m8[ns]"))
    assert res[0] == td
