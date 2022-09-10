import struct
import sys

from hypothesis import (
    assume,
    example,
    given,
    strategies as st,
)
import numpy as np
import pytest

from pandas.io.sas._byteswap import (
    read_double_with_byteswap,
    read_float_with_byteswap,
    read_uint16_with_byteswap,
    read_uint32_with_byteswap,
    read_uint64_with_byteswap,
)

_swapped_byte_order = {"big": "<", "little": ">"}[sys.byteorder]


@given(read_offset=st.integers(0, 11), number=st.integers(min_value=0))
@example(number=2**16, read_offset=0)
@example(number=2**32, read_offset=0)
@example(number=2**64, read_offset=0)
@pytest.mark.parametrize("int_type", ["H", "I", "Q"])
@pytest.mark.parametrize("should_byteswap", [True, False])
def test_int_byteswap(read_offset, number, int_type, should_byteswap):
    int_type_nbytes = struct.calcsize(int_type)
    assume(number < 2 ** (8 * int_type_nbytes))
    number_bytes = struct.pack(int_type, number)
    data = bytearray(np.random.default_rng().bytes(20))
    data[read_offset : read_offset + int_type_nbytes] = number_bytes
    read_uintxx_with_byteswap = {
        "H": read_uint16_with_byteswap,
        "I": read_uint32_with_byteswap,
        "Q": read_uint64_with_byteswap,
    }[int_type]
    output_number = read_uintxx_with_byteswap(bytes(data), read_offset, should_byteswap)
    if should_byteswap:
        (number_bytes_swapped,) = struct.unpack(
            _swapped_byte_order + int_type, number_bytes
        )
        assert output_number == number_bytes_swapped
    else:
        assert output_number == number
