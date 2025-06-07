import pytest

import pandas as pd


@pytest.mark.parametrize("tz", ["US/Eastern", "Europe/Berlin"])
def test_hdf_where_query_on_tzindex_raises(tmp_path, tz):
    df = pd.DataFrame({"x": range(5)})
    df["dt"] = pd.date_range("2020-01-01", periods=5, tz=tz)
    df = df.set_index("dt")

    file_path = tmp_path / "test.h5"
    df.to_hdf(file_path, key="df", format="table")

    with pytest.raises(ValueError, match="invalid variable reference"):
        pd.read_hdf(file_path, key="df", where='dt=="2020-01-03 00:00:00-05:00"')
