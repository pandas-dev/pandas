import numpy as np
import pytest

from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds

from pandas import (
    Timedelta,
    offsets,
)

from pandas.errors import OutOfBoundsTimedelta


@pytest.mark.parametrize(
    "obj,expected",
    [
        (np.timedelta64(14, "D"), 14 * 24 * 3600 * 1e9),
        (Timedelta(minutes=-7), -7 * 60 * 1e9),
        (Timedelta(minutes=-7).to_pytimedelta(), -7 * 60 * 1e9),
        (offsets.Nano(125), 125),
        (1, 1),
        (np.int64(2), 2),
        (np.int32(3), 3),
    ],
)
def test_delta_to_nanoseconds(obj, expected):
    result = delta_to_nanoseconds(obj)
    assert result == expected


def test_delta_to_nanoseconds_error():
    obj = np.array([123456789], dtype="m8[ns]")

    with pytest.raises(TypeError, match="<class 'numpy.ndarray'>"):
        delta_to_nanoseconds(obj)


def test_huge_nanoseconds_overflow():
    # GH 32402
    assert delta_to_nanoseconds(Timedelta(1e10)) == 1e10
    assert delta_to_nanoseconds(Timedelta(nanoseconds=1e10)) == 1e10


def test_out_of_bounds():
    # GH 36615
    with pytest.raises(OutOfBoundsTimedelta, match="200000 days"):
        Timedelta(np.timedelta64(200000, "D"))
