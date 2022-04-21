from typing import Union

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def dataframe(attr_param: Union[list, None, int, float, str]) -> pd.DataFrame:
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df.attrs["test_attr"] = attr_param
    return df


@pytest.mark.parametrize(
    "attr_param", [[1, 2, 3], ["a", "b", "c"], None, "test_attr_str", 1, 1.0]
)
def test_read_write_pyarrow(dataframe):
    with tm.ensure_clean("__tmp_to_parquet_from_parquet__") as path:
        dataframe.to_parquet(path, engine="pyarrow")
        dataframe_read = pd.read_parquet(path, engine="pyarrow")
    assert dataframe.attrs == dataframe_read.attrs


@pytest.mark.parametrize(
    "attr_param", [[1, 2, 3], ["a", "b", "c"], None, "test_attr_str", 1, 1.0]
)
def test_read_write_fastparquet(dataframe):
    with tm.ensure_clean("__tmp_to_parquet_from_parquet__") as path:
        dataframe.to_parquet(path, engine="fastparquet")
        dataframe_read = pd.read_parquet(path, engine="fastparquet")
    assert dataframe.attrs == dataframe_read.attrs
