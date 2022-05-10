import operator

import numpy as np
import pytest

from pandas._libs import ops

from pandas.core.ops import (
    radd,
    rsub,
)


@pytest.fixture(name="int_max", scope="module")
def fixture_int_max() -> int:
    return np.iinfo(np.int64).max


@pytest.fixture(name="int_min", scope="module")
def fixture_int_min() -> int:
    return np.iinfo(np.int64).min


@pytest.fixture(name="overflow_msg", scope="module")
def fixture_overflow_msg() -> str:
    return "|".join(
        (
            "Python int too large to convert to C long",
            "int too big to convert",
        )
    )


@pytest.fixture(name="add_op", params=(operator.add, radd))
def fixture_add_op(request):
    return request.param


@pytest.fixture(name="sub_op", params=(operator.sub, rsub))
def fixture_sub_op(request):
    return request.param


class TestCalcIntInt:
    def test_raises_for_too_large_arg(self, int_max: int, add_op, overflow_msg: str):
        add_op(int_max + 1, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(add_op, int_max + 1, 1)

    def test_raises_for_too_small_arg(self, int_min: int, sub_op, overflow_msg: str):
        sub_op(int_min - 1, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(sub_op, int_min - 1, 1)

    def test_raises_for_too_large_result(self, int_max: int, add_op, overflow_msg: str):
        assert add_op(int_max, 1) == int_max + 1

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(add_op, int_max, 1)

    def test_raises_for_too_small_result(self, int_min: int, sub_op, overflow_msg: str):
        assert abs(sub_op(int_min, 1)) == abs(int_min - 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(sub_op, int_min, 1)
