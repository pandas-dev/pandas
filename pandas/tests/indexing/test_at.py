from datetime import datetime, timezone

import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm


def test_at_timezone():
    # https://github.com/pandas-dev/pandas/issues/33544
    result = DataFrame({"foo": [datetime(2000, 1, 1)]})
    result.at[0, "foo"] = datetime(2000, 1, 2, tzinfo=timezone.utc)
    expected = DataFrame(
        {"foo": [datetime(2000, 1, 2, tzinfo=timezone.utc)]}, dtype=object
    )
    tm.assert_frame_equal(result, expected)


class TestAtWithDuplicates:
    def test_at_with_duplicate_axes_requires_scalar_lookup(self):
        # GH#33041 check that falling back to loc doesn't allow non-scalar
        #  args to slip in

        arr = np.random.randn(6).reshape(3, 2)
        df = DataFrame(arr, columns=["A", "A"])

        msg = "Invalid call for scalar access"
        with pytest.raises(ValueError, match=msg):
            df.at[[1, 2]]
        with pytest.raises(ValueError, match=msg):
            df.at[1, ["A"]]
        with pytest.raises(ValueError, match=msg):
            df.at[:, "A"]

        with pytest.raises(ValueError, match=msg):
            df.at[[1, 2]] = 1
        with pytest.raises(ValueError, match=msg):
            df.at[1, ["A"]] = 1
        with pytest.raises(ValueError, match=msg):
            df.at[:, "A"] = 1
