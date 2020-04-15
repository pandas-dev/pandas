from datetime import datetime, timezone
import pandas as pd


def test_at_timezone():
    # https://github.com/pandas-dev/pandas/issues/33544
    df = pd.DataFrame({"foo": [datetime(2000, 1, 1)]})
    df.at[0, "foo"] = datetime(2000, 1, 2, tzinfo=timezone.utc)
    assert df.at[0, "foo"].tzinfo is not None
