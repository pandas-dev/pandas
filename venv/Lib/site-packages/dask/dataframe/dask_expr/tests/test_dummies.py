from __future__ import annotations

import pandas
import pytest

from dask.dataframe.dask_expr import from_pandas, get_dummies
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

pd = _backend_library()


@pytest.mark.parametrize(
    "data",
    [
        pd.Series([1, 1, 1, 2, 2, 1, 3, 4], dtype="category"),
        pd.Series(
            pandas.Categorical([1, 1, 1, 2, 2, 1, 3, 4], categories=[4, 3, 2, 1])
        ),
        pd.DataFrame(
            {"a": [1, 2, 3, 4, 4, 3, 2, 1], "b": pandas.Categorical(list("abcdabcd"))}
        ),
    ],
)
def test_get_dummies(data):
    exp = pd.get_dummies(data)

    ddata = from_pandas(data, 2)
    res = get_dummies(ddata)
    assert_eq(res, exp)
    pandas.testing.assert_index_equal(res.columns, exp.columns)
