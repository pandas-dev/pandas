import pytest

import pandas as pd

schema = {
    "id": "int64[pyarrow]",
    "time": "timestamp[s][pyarrow]",
    "value": "float[pyarrow]",
}


@pytest.mark.parametrize("dtype", ["timestamp[s][pyarrow]"])
def test_concat_preserves_pyarrow_timestamp(dtype):
    dfA = (
        pd.DataFrame(
            [
                (0, "2021-01-01 00:00:00", 5.3),
                (1, "2021-01-01 00:01:00", 5.4),
                (2, "2021-01-01 00:01:00", 5.4),
                (3, "2021-01-01 00:02:00", 5.5),
            ],
            columns=schema,
        )
        .astype(schema)
        .set_index(["id", "time"])
    )

    dfB = (
        pd.DataFrame(
            [
                (1, "2022-01-01 08:00:00", 6.3),
                (2, "2022-01-01 08:01:00", 6.4),
                (3, "2022-01-01 08:02:00", 6.5),
            ],
            columns=schema,
        )
        .astype(schema)
        .set_index(["id", "time"])
    )

    df = pd.concat([dfA, dfB], keys=[0, 1], names=["run"])

    # check whether df.index is multiIndex
    assert isinstance(df.index, pd.MultiIndex), (
        f"Expected MultiIndex, but received {type(df.index)}"
    )

    # Verifying special dtype timestamp[s][pyarrow] stays intact after concat
    assert df.index.levels[2].dtype == dtype, (
        f"Expected {dtype}, but received {df.index.levels[2].dtype}"
    )
