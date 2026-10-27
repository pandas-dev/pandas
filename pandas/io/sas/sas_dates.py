"""
Conversion of SAS date and datetime values, shared by the SAS7BDAT and XPORT
readers.

SAS stores a date as a float64 count of days, and a datetime as a float64 count
of seconds, since 1960-01-01. Nothing in the value itself says which; the
column's display format does (``DATE9.`` makes it a date, ``DATETIME20.`` a
datetime), so the readers look the format name up here to decide whether and
how to convert a numeric column.
"""

from __future__ import annotations

import numpy as np

from pandas._libs.tslibs.conversion import cast_from_unit_vectorized
from pandas.errors import OutOfBoundsDatetime

from pandas import Timestamp

import pandas.io.sas.sas_constants as const

_unix_origin = Timestamp("1970-01-01")
_sas_origin = Timestamp("1960-01-01")

# SAS uses a modified Gregorian calendar where years divisible by 4000 are
# not leap years (unlike proleptic Gregorian). These are the SAS day counts
# for the first day affected by each 4000-year boundary.
# See https://communities.sas.com/t5/SAS-Programming/Leap-Years-divisible-by-4000/td-p/663467
_SAS_MARCH1_4000 = 745154  # SAS day count for March 1, 4000
_SAS_MARCH1_8000 = 2206123  # SAS day count for March 1, 8000


def convert_type_for_format(fmt: str) -> str | None:
    """
    Which conversion a numeric column with SAS display format `fmt` needs.

    Parameters
    ----------
    fmt : str
        The format name without width or decimals, as the file records it
        (``"DATE"`` for ``DATE9.``).

    Returns
    -------
    {"d", "s"} or None
        "d" if the column holds day counts (a date format), "s" if it holds
        second counts (a datetime format), None if it is not a date at all.
    """
    if fmt in const.sas_date_formats:
        return "d"
    if fmt in const.sas_datetime_formats:
        return "s"
    return None


def _sas_to_gregorian_correction(values: np.ndarray, unit: str) -> np.ndarray:
    """
    Compute the additive correction (in `unit`) to convert SAS day/second counts
    to proleptic Gregorian day/second counts.

    SAS omits Feb 29 for years divisible by 4000 (unlike proleptic Gregorian);
    this adds back the missing days. `unit` must be "d" (days) or "s" (seconds).
    """
    scale = 86400 if unit == "s" else 1
    thresholds = np.array([_SAS_MARCH1_4000, _SAS_MARCH1_8000], dtype=np.int64) * scale
    correction = np.zeros(len(values), dtype=np.float64)
    valid = ~np.isnan(values)
    for threshold in thresholds:
        correction[valid] += (values[valid] >= threshold).astype(np.float64) * scale
    return correction


def convert_datetimes(sas_datetimes: np.ndarray, unit: str) -> np.ndarray:
    """
    Convert SAS day or second counts to a datetime64 array.

    Parameters
    ----------
    sas_datetimes : ndarray of float64
       Dates or datetimes in SAS; NaN marks a missing value.
    unit : {'d', 's'}
       "d" if the floats represent dates, "s" for datetimes

    Returns
    -------
    ndarray
       datetime64[s] for unit="d", datetime64[ms] for unit="s". Missing
       values become NaT.
    """
    td = (_sas_origin - _unix_origin).as_unit("s")
    # SAS's own date range tops out near 6e6 days, so a count this size is not
    # a date the file could legitimately hold -- it is corrupt bytes, or a
    # numeric column carrying a date format. Casting it does not overflow, it
    # saturates: a negative one lands on the NaT sentinel and is read as
    # missing, and a positive one does raise below, but naming the date the
    # saturated cast landed on rather than anything the file holds.
    # A day count is scaled by 86400 below, so its own limit is that much lower.
    limit = 2.0**63 / 86400 if unit == "d" else 2.0**63
    too_large = np.abs(sas_datetimes) >= limit
    if too_large.any():
        value = sas_datetimes[too_large][0]
        what = "date" if unit == "d" else "datetime"
        raise OutOfBoundsDatetime(
            f"Out of bounds SAS {what} value: {value}; no SAS {what} can be this "
            f"large, so the file is corrupt or the column is not a {what}"
        )
    if unit == "s":
        corrected = sas_datetimes + _sas_to_gregorian_correction(
            sas_datetimes, unit="s"
        )
        millis = cast_from_unit_vectorized(corrected, unit="s", out_unit="ms")
        return millis.view("M8[ms]") + td
    else:
        corrected = sas_datetimes + _sas_to_gregorian_correction(
            sas_datetimes, unit="d"
        )
        # A date-formatted column is a float64 day count that SAS does not force
        # whole, so scale the fraction in rather than truncating it with an M8[D]
        # cast. Round to seconds here instead of scaling from "D" inside
        # cast_from_unit_vectorized, whose rounding precision is a power of ten:
        # 1e-4 of a day, coarser than the seconds this returns.
        secs = cast_from_unit_vectorized(
            np.round(corrected * 86400.0), unit="s", out_unit="s"
        )
        return secs.view("M8[s]") + td
