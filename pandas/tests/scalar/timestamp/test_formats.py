import pytest

from pandas import Timestamp


def test_isoformat():
    ts = Timestamp(
        year=2019, month=5, day=18, hour=15, minute=17, second=8, microsecond=132263
    )
    assert ts.isoformat() == "2019-05-18T15:17:08.132263"
    assert ts.isoformat(timespec="seconds") == "2019-05-18T15:17:08"
