import operator

import numpy as np
import pytest

from pandas._libs import ops


@pytest.fixture(name="int_max")
def fixture_int_max() -> int:
    return np.iinfo(np.int64).max


@pytest.fixture(name="int_min")
def fixture_int_min() -> int:
    return np.iinfo(np.int64).min


@pytest.fixture(name="overflow_msg")
def fixture_overflow_msg() -> str:
    return "Python int too large to convert to C long"


def test_raises_for_too_large_arg(int_max: int, overflow_msg: str):
    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.add, int_max + 1, 1)

    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.add, 1, int_max + 1)


def test_raises_for_too_small_arg(int_min: int, overflow_msg: str):
    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.add, int_min - 1, 1)

    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.add, 1, int_min - 1)


def test_raises_for_too_large_result(int_max: int, overflow_msg: str):
    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.add, int_max, 1)

    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.add, 1, int_max)


def test_raises_for_too_small_result(int_min: int, overflow_msg: str):
    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.sub, int_min, 1)

    with pytest.raises(OverflowError, match=overflow_msg):
        ops.calculate(operator.sub, 1, int_min)
