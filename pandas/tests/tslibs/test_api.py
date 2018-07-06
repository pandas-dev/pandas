# -*- coding: utf-8 -*-
"""Tests that the tslibs API is locked down"""

from pandas._libs import tslibs


def test_namespace():
    # just check that none of these raise NameErrors
    tslibs.normalize_date
    tslibs.localize_pydatetime
    tslibs.tz_convert_single
    tslibs.NaT
    tslibs.iNaT
    tslibs.OutOfBoundsDatetime
    tslibs.Timestamp
    tslibs.Timedelta
    tslibs.delta_to_nanoseconds
