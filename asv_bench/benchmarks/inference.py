"""
The functions benchmarked in this file depend _almost_ exclusively on
_libs, but not in a way that is easy to formalize.

If a PR does not change anything in pandas/_libs/ or pandas/core/tools/, then
it is likely that these benchmarks will be unaffected.
"""

import numpy as np

from pandas import (
    NaT,
    Series,
    date_range,
    to_datetime,
    to_numeric,
    to_timedelta,
)

from .pandas_vb_common import (
    lib,
    tm,
)


class ToNumeric:

    params = ["ignore", "coerce"]
    param_names = ["errors"]

    def setup(self, errors):
        N = 10000
        self.float = Series(np.random.randn(N))
        self.numstr = self.float.astype("str")
        self.str = Series(tm.makeStringIndex(N))

    def time_from_float(self, errors):
        to_numeric(self.float, errors=errors)

    def time_from_numeric_str(self, errors):
        to_numeric(self.numstr, errors=errors)

    def time_from_str(self, errors):
        to_numeric(self.str, errors=errors)


class ToNumericDowncast:

    param_names = ["dtype", "downcast"]
    params = [
        [
            "string-float",
            "string-int",
            "string-nint",
            "datetime64",
            "int-list",
            "int32",
        ],
        [None, "integer", "signed", "unsigned", "float"],
    ]

    N = 500000
    N2 = N // 2

    data_dict = {
        "string-int": ["1"] * N2 + [2] * N2,
        "string-nint": ["-1"] * N2 + [2] * N2,
        "datetime64": np.repeat(
            np.array(["1970-01-01", "1970-01-02"], dtype="datetime64[D]"), N
        ),
        "string-float": ["1.1"] * N2 + [2] * N2,
        "int-list": [1] * N2 + [2] * N2,
        "int32": np.repeat(np.int32(1), N),
    }

    def setup(self, dtype, downcast):
        self.data = self.data_dict[dtype]

    def time_downcast(self, dtype, downcast):
        to_numeric(self.data, downcast=downcast)


class MaybeConvertNumeric:
    # maybe_convert_numeric depends _exclusively_ on _libs, could
    #  go in benchmarks/libs.py

    def setup_cache(self):
        N = 10 ** 6
        arr = np.repeat([2 ** 63], N) + np.arange(N).astype("uint64")
        data = arr.astype(object)
        data[1::2] = arr[1::2].astype(str)
        data[-1] = -1
        return data

    def time_convert(self, data):
        lib.maybe_convert_numeric(data, set(), coerce_numeric=False)


class MaybeConvertObjects:
    # maybe_convert_objects depends _almost_ exclusively on _libs, but
    #  does have some run-time imports from outside of _libs

    def setup(self):
        N = 10 ** 5

        data = list(range(N))
        data[0] = NaT
        data = np.array(data)
        self.data = data

    def time_maybe_convert_objects(self):
        lib.maybe_convert_objects(self.data)


class ToDatetimeFromIntsFloats:
    def setup(self):
        self.ts_sec = Series(range(1521080307, 1521685107), dtype="int64")
        self.ts_sec_float = self.ts_sec.astype("float64")

        self.ts_nanosec = 1_000_000 * self.ts_sec
        self.ts_nanosec_float = self.ts_nanosec.astype("float64")

    # speed of int64 and float64 paths should be comparable

    def time_nanosec_int64(self):
        to_datetime(self.ts_nanosec, unit="ns")

    def time_nanosec_float64(self):
        to_datetime(self.ts_nanosec_float, unit="ns")

    def time_sec_int64(self):
        to_datetime(self.ts_sec, unit="s")

    def time_sec_float64(self):
        to_datetime(self.ts_sec_float, unit="s")


class ToDatetimeYYYYMMDD:
    def setup(self):
        rng = date_range(start="1/1/2000", periods=10000, freq="D")
        self.stringsD = Series(rng.strftime("%Y%m%d"))

    def time_format_YYYYMMDD(self):
        to_datetime(self.stringsD, format="%Y%m%d")


class ToDatetimeCacheSmallCount:

    params = ([True, False], [50, 500, 5000, 100000])
    param_names = ["cache", "count"]

    def setup(self, cache, count):
        rng = date_range(start="1/1/1971", periods=count)
        self.unique_date_strings = rng.strftime("%Y-%m-%d").tolist()

    def time_unique_date_strings(self, cache, count):
        to_datetime(self.unique_date_strings, cache=cache)


class ToDatetimeISO8601:
    def setup(self):
        rng = date_range(start="1/1/2000", periods=20000, freq="H")
        self.strings = rng.strftime("%Y-%m-%d %H:%M:%S").tolist()
        self.strings_nosep = rng.strftime("%Y%m%d %H:%M:%S").tolist()
        self.strings_tz_space = [
            x.strftime("%Y-%m-%d %H:%M:%S") + " -0800" for x in rng
        ]

    def time_iso8601(self):
        to_datetime(self.strings)

    def time_iso8601_nosep(self):
        to_datetime(self.strings_nosep)

    def time_iso8601_format(self):
        to_datetime(self.strings, format="%Y-%m-%d %H:%M:%S")

    def time_iso8601_format_no_sep(self):
        to_datetime(self.strings_nosep, format="%Y%m%d %H:%M:%S")

    def time_iso8601_tz_spaceformat(self):
        to_datetime(self.strings_tz_space)


class ToDatetimeNONISO8601:
    def setup(self):
        N = 10000
        half = N // 2
        ts_string_1 = "March 1, 2018 12:00:00+0400"
        ts_string_2 = "March 1, 2018 12:00:00+0500"
        self.same_offset = [ts_string_1] * N
        self.diff_offset = [ts_string_1] * half + [ts_string_2] * half

    def time_same_offset(self):
        to_datetime(self.same_offset)

    def time_different_offset(self):
        to_datetime(self.diff_offset)


class ToDatetimeFormatQuarters:
    def setup(self):
        self.s = Series(["2Q2005", "2Q05", "2005Q1", "05Q1"] * 10000)

    def time_infer_quarter(self):
        to_datetime(self.s)


class ToDatetimeFormat:
    def setup(self):
        N = 100000
        self.s = Series(["19MAY11", "19MAY11:00:00:00"] * N)
        self.s2 = self.s.str.replace(":\\S+$", "")

        self.same_offset = ["10/11/2018 00:00:00.045-07:00"] * N
        self.diff_offset = [
            f"10/11/2018 00:00:00.045-0{offset}:00" for offset in range(10)
        ] * (N // 10)

    def time_exact(self):
        to_datetime(self.s2, format="%d%b%y")

    def time_no_exact(self):
        to_datetime(self.s, format="%d%b%y", exact=False)

    def time_same_offset(self):
        to_datetime(self.same_offset, format="%m/%d/%Y %H:%M:%S.%f%z")

    def time_different_offset(self):
        to_datetime(self.diff_offset, format="%m/%d/%Y %H:%M:%S.%f%z")

    def time_same_offset_to_utc(self):
        to_datetime(self.same_offset, format="%m/%d/%Y %H:%M:%S.%f%z", utc=True)

    def time_different_offset_to_utc(self):
        to_datetime(self.diff_offset, format="%m/%d/%Y %H:%M:%S.%f%z", utc=True)


class ToDatetimeCache:

    params = [True, False]
    param_names = ["cache"]

    def setup(self, cache):
        N = 10000
        self.unique_numeric_seconds = list(range(N))
        self.dup_numeric_seconds = [1000] * N
        self.dup_string_dates = ["2000-02-11"] * N
        self.dup_string_with_tz = ["2000-02-11 15:00:00-0800"] * N

    def time_unique_seconds_and_unit(self, cache):
        to_datetime(self.unique_numeric_seconds, unit="s", cache=cache)

    def time_dup_seconds_and_unit(self, cache):
        to_datetime(self.dup_numeric_seconds, unit="s", cache=cache)

    def time_dup_string_dates(self, cache):
        to_datetime(self.dup_string_dates, cache=cache)

    def time_dup_string_dates_and_format(self, cache):
        to_datetime(self.dup_string_dates, format="%Y-%m-%d", cache=cache)

    def time_dup_string_tzoffset_dates(self, cache):
        to_datetime(self.dup_string_with_tz, cache=cache)


class ToTimedelta:
    def setup(self):
        self.ints = np.random.randint(0, 60, size=10000)
        self.str_days = []
        self.str_seconds = []
        for i in self.ints:
            self.str_days.append(f"{i} days")
            self.str_seconds.append(f"00:00:{i:02d}")

    def time_convert_int(self):
        to_timedelta(self.ints, unit="s")

    def time_convert_string_days(self):
        to_timedelta(self.str_days)

    def time_convert_string_seconds(self):
        to_timedelta(self.str_seconds)


class ToTimedeltaErrors:

    params = ["coerce", "ignore"]
    param_names = ["errors"]

    def setup(self, errors):
        ints = np.random.randint(0, 60, size=10000)
        self.arr = [f"{i} days" for i in ints]
        self.arr[-1] = "apple"

    def time_convert(self, errors):
        to_timedelta(self.arr, errors=errors)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
