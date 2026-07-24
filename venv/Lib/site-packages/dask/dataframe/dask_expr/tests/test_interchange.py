from __future__ import annotations

from pandas.core.interchange.dataframe_protocol import DtypeKind

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library

# Set DataFrame backend for this module
pd = _backend_library()


def test_interchange_protocol():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    df = from_pandas(pdf, npartitions=2)
    df_int = df.__dataframe__()
    pd.testing.assert_index_equal(pdf.columns, df_int.column_names())
    assert df_int.num_columns() == 2
    assert df_int.num_rows() == 3
    column = df_int.get_columns()[0]
    assert column.dtype()[0] == DtypeKind.INT
