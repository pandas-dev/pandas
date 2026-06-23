import numpy as np
import pytest

import pandas as pd
from pandas import MultiIndex
import pandas._testing as tm


def test_pickle_compat_construction():
    # this is testing for pickle compat
    # need an object to create with
    with pytest.raises(TypeError, match="Must pass both levels and codes"):
        MultiIndex()


def test_multiindex_datetime64ns_pickle_roundtrip():
    # GH#63078 - pickling a MultiIndex with datetime64[ns] level raised
    # NotImplementedError in NDArrayBacked.__setstate__
    df = pd.DataFrame(
        {
            "date": [np.datetime64("20110101", "ns")],
            "id": [1],
            "val": [1],
        }
    ).set_index(["date", "id"])

    result = tm.round_trip_pickle(df)
    tm.assert_frame_equal(result, df)
