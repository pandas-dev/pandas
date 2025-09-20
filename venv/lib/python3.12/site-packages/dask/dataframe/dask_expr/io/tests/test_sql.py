from __future__ import annotations

import pytest

from dask.dataframe import read_sql_table
from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq
from dask.utils import tmpfile

pd = _backend_library()

pytest.importorskip("sqlalchemy")


def test_shuffle_after_read_sql():
    with tmpfile() as f:
        uri = "sqlite:///%s" % f

        df = pd.DataFrame(
            {
                "id": [1, 2, 3, 4, 5, 6, 7, 8],
                "value": [
                    "value1",
                    "value2",
                    "value3",
                    "value3",
                    "value4",
                    "value4",
                    "value4",
                    "value5",
                ],
            }
        ).set_index("id")
        ddf = from_pandas(df, npartitions=1)

        ddf.to_sql("test_table", uri, if_exists="append")
        result = read_sql_table("test_table", con=uri, index_col="id")
        assert_eq(
            result["value"].unique(), pd.Series(df["value"].unique(), name="value")
        )
        assert_eq(result.shuffle(on_index=True), df)
