import io

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.xfail(
    reason=(
        "fastparquet bug: index is corrupted when reading multiple "
        "BytesIO parquets (GH#64007)"
    ),
    strict=True,
)
@pytest.mark.filterwarnings("ignore:.*fastparquet.*")
def test_fastparquet_bytesio_multiple_files_preserve_index():
    # GH#64007

    def make_df(start):
        return pd.DataFrame(
            {"A": [1, 2, 3], "B": [1, 2, 3]},
            index=[start, start + 1, start + 2],
        )

    df1 = make_df(1)
    df2 = make_df(2)

    buf1 = io.BytesIO()
    buf2 = io.BytesIO()

    df1.to_parquet(buf1, engine="fastparquet")
    df2.to_parquet(buf2, engine="fastparquet")

    buf1.seek(0)
    buf2.seek(0)

    result1 = pd.read_parquet(buf1, engine="fastparquet")
    result2 = pd.read_parquet(buf2, engine="fastparquet")

    tm.assert_index_equal(result1.index, df1.index)
    tm.assert_index_equal(result2.index, df2.index)

