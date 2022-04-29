"""Tests that the tslibs API is locked down"""

from pandas._libs import tslibs


def test_namespace():

    submodules = [
        "base",
        "ccalendar",
        "conversion",
        "dtypes",
        "fields",
        "nattype",
        "np_datetime",
        "offsets",
        "parsing",
        "period",
        "strptime",
        "vectorized",
        "timedeltas",
        "timestamps",
        "timezones",
        "tzconversion",
    ]

    api = [
        "BaseOffset",
        "NaT",
        "NaTType",
        "iNaT",
        "nat_strings",
        "OutOfBoundsDatetime",
        "OutOfBoundsTimedelta",
        "Period",
        "IncompatibleFrequency",
        "Resolution",
        "Tick",
        "Timedelta",
        "dt64arr_to_periodarr",
        "Timestamp",
        "is_date_array_normalized",
        "ints_to_pydatetime",
        "normalize_i8_timestamps",
        "get_resolution",
        "delta_to_nanoseconds",
        "ints_to_pytimedelta",
        "localize_pydatetime",
        "tz_convert_from_utc",
        "tz_convert_from_utc_single",
        "to_offset",
        "tz_compare",
    ]

    expected = set(submodules + api)
    # exclude "ctime" bc it is not (yet) imported outside of tests
    names = [x for x in dir(tslibs) if not x.startswith("__") and x != "ctime"]
    assert set(names) == expected
