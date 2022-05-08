import operator
from platform import architecture

import numpy as np
import pytest

from pandas._libs import ops


@pytest.fixture(name="int_max", scope="module")
def fixture_int_max() -> int:
    return np.iinfo(np.int64).max


@pytest.fixture(name="int_min", scope="module")
def fixture_int_min() -> int:
    return np.iinfo(np.int64).min


@pytest.fixture(name="float_max", scope="module")
def fixture_float_max() -> np.float64:
    return np.finfo(np.float64).max


@pytest.fixture(name="float_min", scope="module")
def fixture_float_min() -> np.float64:
    return np.finfo(np.float64).min


@pytest.fixture(name="overflow_msg", scope="module")
def fixture_overflow_msg() -> str:
    return "|".join(
        (
            "Python int too large to convert to C long",
            "int too big to convert",
        )
    )


class TestCalcIntInt:
    def test_raises_for_too_large_arg(self, int_max: int, overflow_msg: str):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.add, int_max + 1, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.add, 1, int_max + 1)

    def test_raises_for_too_small_arg(self, int_min: int, overflow_msg: str):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.add, int_min - 1, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.add, 1, int_min - 1)

    def test_raises_for_too_large_result(self, int_max: int, overflow_msg: str):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.add, int_max, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.add, 1, int_max)

    def test_raises_for_too_small_result(self, int_min: int, overflow_msg: str):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.sub, int_min, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_int(operator.sub, 1, int_min)


class TestCalcIntFloat:
    @pytest.mark.parametrize(
        "op,lval,rval,expected",
        (
            (operator.add, 1, 1.0, 2),
            (operator.sub, 2, 1.0, 1),
            (operator.mul, 1, 2.0, 2),
            (operator.truediv, 1, 0.5, 2),
        ),
        ids=("+", "-", "*", "/"),
    )
    def test_arithmetic_ops(self, op, lval: int, rval: float, expected: int):
        result = ops.calc_int_float(op, lval, rval)

        assert result == expected
        assert isinstance(result, int)

    def test_raises_for_too_large_arg(
        self,
        int_max: int,
        float_max: float,
        overflow_msg: str,
    ):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.add, int_max + 1, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.add, 1, float_max + 1)

    def test_raises_for_too_small_arg(
        self,
        int_min: int,
        float_min: float,
        overflow_msg: str,
    ):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.add, int_min - 1, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.add, 1, float_min - 1)

    def test_raises_for_too_large_result(
        self,
        int_max: int,
        float_max: float,
        overflow_msg: str,
    ):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.add, int_max, 1)

        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.add, 1, float_max)

    @pytest.mark.parametrize(
        "value",
        (
            pytest.param(
                1024,
                marks=pytest.mark.xfail(
                    reason="TBD",
                    raises=pytest.fail.Exception,
                    strict=True,
                ),
            ),
            pytest.param(
                1024.1,
                marks=pytest.mark.xfail(
                    condition=architecture()[0] == "32bit",
                    reason="overflows earlier",
                    raises=pytest.fail.Exception,
                    strict=True,
                ),
            ),
        ),
    )
    def test_raises_for_most_too_small_results(
        self,
        value: float,
        int_min: int,
        overflow_msg: str,
    ):
        with pytest.raises(OverflowError, match=overflow_msg):
            ops.calc_int_float(operator.sub, int_min, value)
