from __future__ import annotations

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


def test_cumulative_empty_partitions():
    pdf = pd.DataFrame(
        {"x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]},
        index=pd.date_range("1995-02-26", periods=8, freq="5min"),
        dtype=float,
    )
    pdf2 = pdf.drop(pdf.between_time("00:10", "00:20").index)

    df = from_pandas(pdf, npartitions=8)
    df2 = from_pandas(pdf2, npartitions=1).repartition(df.divisions)

    assert_eq(df2.cumprod(), pdf2.cumprod())
    assert_eq(df2.cumsum(), pdf2.cumsum())
