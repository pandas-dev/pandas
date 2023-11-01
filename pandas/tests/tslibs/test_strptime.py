from datetime import (
    datetime,
    timezone,
)

import numpy as np
import pytest

from pandas._libs.tslibs.dtypes import NpyDatetimeUnit
from pandas._libs.tslibs.strptime import array_strptime

from pandas import Timestamp
import pandas._testing as tm

creso_infer = NpyDatetimeUnit.NPY_FR_GENERIC.value


class TestArrayStrptimeResolutionInference:
    @pytest.mark.parametrize("tz", [None, timezone.utc])
    def test_array_strptime_resolution_inference_homogeneous_strings(self, tz):
        dt = datetime(2016, 1, 2, 3, 4, 5, 678900, tzinfo=tz)

        fmt = "%Y-%m-%d %H:%M:%S"
        dtstr = dt.strftime(fmt)
        arr = np.array([dtstr] * 3, dtype=object)
        expected = np.array([dt.replace(tzinfo=None)] * 3, dtype="M8[s]")

        res, _ = array_strptime(arr, fmt=fmt, utc=False, creso=creso_infer)
        tm.assert_numpy_array_equal(res, expected)

        fmt = "%Y-%m-%d %H:%M:%S.%f"
        dtstr = dt.strftime(fmt)
        arr = np.array([dtstr] * 3, dtype=object)
        expected = np.array([dt.replace(tzinfo=None)] * 3, dtype="M8[us]")

        res, _ = array_strptime(arr, fmt=fmt, utc=False, creso=creso_infer)
        tm.assert_numpy_array_equal(res, expected)

        fmt = "ISO8601"
        res, _ = array_strptime(arr, fmt=fmt, utc=False, creso=creso_infer)
        tm.assert_numpy_array_equal(res, expected)

    @pytest.mark.parametrize("tz", [None, timezone.utc])
    def test_array_strptime_resolution_mixed(self, tz):
        dt = datetime(2016, 1, 2, 3, 4, 5, 678900, tzinfo=tz)

        ts = Timestamp(dt).as_unit("ns")

        arr = np.array([dt, ts], dtype=object)
        expected = np.array(
            [Timestamp(dt).as_unit("ns").asm8, ts.asm8],
            dtype="M8[ns]",
        )

        fmt = "%Y-%m-%d %H:%M:%S"
        res, _ = array_strptime(arr, fmt=fmt, utc=False, creso=creso_infer)
        tm.assert_numpy_array_equal(res, expected)

        fmt = "ISO8601"
        res, _ = array_strptime(arr, fmt=fmt, utc=False, creso=creso_infer)
        tm.assert_numpy_array_equal(res, expected)
