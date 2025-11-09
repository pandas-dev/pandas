import numpy as np
import pytest

import pandas as pd

pa = pytest.importorskip("pyarrow")

xfail_sparse = pytest.mark.xfail(
    reason="pending implementation of preserve_sparse for Parquet",
    strict=False,
)


@xfail_sparse
@pytest.mark.parametrize(
    "subtype, fill_value, data",
    [
        ("int64", 0, [0, 0, 3, 0, 5]),
        ("float64", 0.0, [0.0, 0.0, 1.5, 0.0, 2.5]),
        ("boolean", False, [False, False, True, False, True]),
    ],
)
def test_parquet_sparse_roundtrip(tmp_path, subtype, fill_value, data):
    path = tmp_path / "out.parquet"
    s = pd.Series(pd.arrays.SparseArray(data, fill_value=fill_value))
    df = pd.DataFrame({"s": s, "x": np.arange(len(s))})

    df.to_parquet(path, preserve_sparse=True)
    df2 = pd.read_parquet(path, preserve_sparse=True)

    assert isinstance(df2["s"].dtype, pd.SparseDtype)
    assert df2["s"].dtype.fill_value == fill_value
    pd.testing.assert_series_equal(
        df2["s"].sparse.to_dense(), s.sparse.to_dense(), check_dtype=False
    )
