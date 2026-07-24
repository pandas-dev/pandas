from __future__ import annotations

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


def test_repartition_quantiles():
    pdf = pd.DataFrame({"a": [1, 2, 3, 4, 5, 15, 7, 8, 9, 10, 11], "d": 3})
    df = from_pandas(pdf, npartitions=5)
    result = df.a._repartition_quantiles(npartitions=5)
    expected = pd.Series(
        [1, 1, 3, 7, 9, 15], index=[0, 0.2, 0.4, 0.6, 0.8, 1], name="a"
    )
    assert_eq(result, expected, check_exact=False)

    result = df.a._repartition_quantiles(npartitions=4)
    expected = pd.Series([1, 2, 5, 8, 15], index=[0, 0.25, 0.5, 0.75, 1], name="a")
    assert_eq(result, expected, check_exact=False)
