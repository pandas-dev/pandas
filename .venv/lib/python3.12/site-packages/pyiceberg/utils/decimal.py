#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.

"""Helper methods for working with Python Decimals."""

import math
from decimal import Decimal


def decimal_to_unscaled(value: Decimal) -> int:
    """Get an unscaled value given a Decimal value.

    Args:
        value (Decimal): A Decimal instance.

    Returns:
        int: The unscaled value.
    """
    sign, digits, _ = value.as_tuple()
    return int(Decimal((sign, digits, 0)).to_integral_value())


def unscaled_to_decimal(unscaled: int, scale: int) -> Decimal:
    """Get a scaled Decimal value given an unscaled value and a scale.

    Args:
        unscaled (int): An unscaled value.
        scale (int): A scale to set for the returned Decimal instance.

    Returns:
        Decimal: A scaled Decimal instance.
    """
    sign, digits, _ = Decimal(unscaled).as_tuple()
    return Decimal((sign, digits, -scale))


def bytes_required(value: int | Decimal) -> int:
    """Return the minimum number of bytes needed to serialize a decimal or unscaled value.

    Args:
        value (int | Decimal): a Decimal value or unscaled int value.

    Returns:
        int: the minimum number of bytes needed to serialize the value.
    """
    if isinstance(value, int):
        return (value.bit_length() + 8) // 8
    elif isinstance(value, Decimal):
        return (decimal_to_unscaled(value).bit_length() + 8) // 8

    raise ValueError(f"Unsupported value: {value}")


def decimal_to_bytes(value: Decimal, byte_length: int | None = None) -> bytes:
    """Return a byte representation of a decimal.

    Args:
        value (Decimal): a decimal value.
        byte_length (int): The number of bytes.
    Returns:
        bytes: the unscaled value of the Decimal as bytes.
    """
    unscaled_value = decimal_to_unscaled(value)
    if byte_length is None:
        byte_length = bytes_required(unscaled_value)
    return unscaled_value.to_bytes(byte_length, byteorder="big", signed=True)


def bytes_to_decimal(value: bytes, scale: int) -> Decimal:
    """Return a decimal from the bytes.

    Args:
        value (bytes): the bytes to be converted into a decimal.
        scale (int): the scale of the decimal.

    Returns:
        Decimal: the scaled decimal.
    """
    unscaled_datum = int.from_bytes(value, byteorder="big", signed=True)
    return unscaled_to_decimal(unscaled_datum, scale)


def truncate_decimal(value: Decimal, width: int) -> Decimal:
    """Get a truncated Decimal value given a decimal value and a width.

    Args:
        value (Decimal): a decimal value.
        width (int): A width for the returned Decimal instance.
    Returns:
        Decimal: A truncated Decimal instance.
    """
    unscaled_value = decimal_to_unscaled(value)
    applied_value = unscaled_value - (((unscaled_value % width) + width) % width)
    return unscaled_to_decimal(applied_value, abs(int(value.as_tuple().exponent)))


MAX_PRECISION = tuple(math.floor(math.log10(math.fabs(math.pow(2, 8 * pos - 1) - 1))) for pos in range(24))
REQUIRED_LENGTH = tuple(next(pos for pos in range(24) if p <= MAX_PRECISION[pos]) for p in range(40))


def decimal_required_bytes(precision: int) -> int:
    """Compute the number of bytes required to store a precision.

    Args:
        precision: The number of digits to store.

    Returns:
        The number of bytes required to store a decimal with a certain precision.
    """
    if precision <= 0 or precision >= 40:
        raise ValueError(f"Unsupported precision, outside of (0, 40]: {precision}")

    return REQUIRED_LENGTH[precision]
