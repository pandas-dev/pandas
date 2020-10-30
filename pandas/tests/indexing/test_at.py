from datetime import datetime, timezone

import pandas as pd
import pandas._testing as tm


def test_at_timezone():
    # https://github.com/pandas-dev/pandas/issues/33544
    result = pd.DataFrame({"foo": [datetime(2000, 1, 1)]})
    result.at[0, "foo"] = datetime(2000, 1, 2, tzinfo=timezone.utc)
    expected = pd.DataFrame(
        {"foo": [datetime(2000, 1, 2, tzinfo=timezone.utc)]}, dtype=object
    )
    tm.assert_frame_equal(result, expected)
