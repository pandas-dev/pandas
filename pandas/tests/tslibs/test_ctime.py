"""
These tests are adapted from cython's tests
https://github.com/cython/cython/blob/master/tests/run/time_pxd.pyx
"""
# TODO(cython3): use cython's cpython.time implementation


import time

# error: Module "pandas._libs.tslibs" has no attribute "ctime"
from pandas._libs.tslibs import ctime  # type: ignore[attr-defined]


def test_time():
    # check that ctime.time() matches time.time() to within call-time tolerance
    tic1 = time.time()
    tic2 = ctime.pytime()
    tic3 = time.time()

    assert tic1 <= tic3  # sanity check
    assert tic1 <= tic2
    assert tic2 <= tic3


def test_localtime():
    ltp = time.localtime()
    ltc = ctime.pylocaltime()

    if ltp.tm_sec != ltc["tm_sec"]:
        # If the time.localtime call is just before the end of a second and the
        #  ctime.localtime call is just after the beginning of the next second,
        #  re-call.  This should not occur twice in a row.
        ltp = time.localtime()
        ltc = ctime.pylocaltime()

    assert ltp.tm_year == ltc["tm_year"] or (ltp.tm_year, ltc["tm_year"])
    assert ltp.tm_mon == ltc["tm_mon"] or (ltp.tm_mon, ltc["tm_mon"])
    assert ltp.tm_mday == ltc["tm_mday"] or (ltp.tm_mday, ltc["tm_mday"])
    assert ltp.tm_hour == ltc["tm_hour"] or (ltp.tm_hour, ltc["tm_hour"])
    assert ltp.tm_min == ltc["tm_min"] or (ltp.tm_min, ltc["tm_min"])
    assert ltp.tm_sec == ltc["tm_sec"] or (ltp.tm_sec, ltc["tm_sec"])
    assert ltp.tm_wday == ltc["tm_wday"] or (ltp.tm_wday, ltc["tm_wday"])
    assert ltp.tm_yday == ltc["tm_yday"] or (ltp.tm_yday, ltc["tm_yday"])
    assert ltp.tm_isdst == ltc["tm_isdst"] or (ltp.tm_isdst, ltc["tm_isdst"])
