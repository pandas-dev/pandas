# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
"""Utility module for various conversions around PrimitiveType implementations.

This module enables:
    - Converting partition strings to built-in python objects.
    - Converting a value to a byte buffer.
    - Converting a byte buffer to a value.

Note:
    Conversion logic varies based on the PrimitiveType implementation. Therefore conversion functions
    are defined here as generic functions using the @singledispatch decorator. For each PrimitiveType
    implementation, a concrete function is registered for each generic conversion function. For PrimitiveType
    implementations that share the same conversion logic, registrations can be stacked.
"""

import uuid
from datetime import date, datetime, time
from decimal import Decimal
from functools import singledispatch
from struct import Struct
from typing import (
    Any,
    Callable,
    Optional,
    Union,
)

from pyiceberg.typedef import UTF8, L
from pyiceberg.types import (
    BinaryType,
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    FixedType,
    FloatType,
    IntegerType,
    LongType,
    PrimitiveType,
    StringType,
    TimestampType,
    TimestamptzType,
    TimeType,
    UUIDType,
    strtobool,
)
from pyiceberg.utils.datetime import date_to_days, datetime_to_micros, time_to_micros
from pyiceberg.utils.decimal import decimal_to_bytes, unscaled_to_decimal

_BOOL_STRUCT = Struct("<?")
_INT_STRUCT = Struct("<i")
_LONG_STRUCT = Struct("<q")
_FLOAT_STRUCT = Struct("<f")
_DOUBLE_STRUCT = Struct("<d")


def handle_none(func: Callable) -> Callable:  # type: ignore
    """Handle cases where partition values are `None` or "__HIVE_DEFAULT_PARTITION__".

    Args:
        func (Callable): A function registered to the singledispatch function `partition_to_py`.
    """

    def wrapper(primitive_type: PrimitiveType, value_str: Optional[str]) -> Any:
        if value_str is None:
            return None
        elif value_str == "__HIVE_DEFAULT_PARTITION__":
            return None
        return func(primitive_type, value_str)

    return wrapper


@singledispatch
def partition_to_py(primitive_type: PrimitiveType, value_str: str) -> Union[int, float, str, uuid.UUID, bytes, Decimal]:
    """Convert a partition string to a python built-in.

    Args:
        primitive_type (PrimitiveType): An implementation of the PrimitiveType base class.
        value_str (str): A string representation of a partition value.
    """
    raise TypeError(f"Cannot convert '{value_str}' to unsupported type: {primitive_type}")


@partition_to_py.register(BooleanType)
@handle_none
def _(primitive_type: BooleanType, value_str: str) -> Union[int, float, str, uuid.UUID]:
    return strtobool(value_str)


@partition_to_py.register(IntegerType)
@partition_to_py.register(LongType)
@partition_to_py.register(DateType)
@partition_to_py.register(TimeType)
@partition_to_py.register(TimestampType)
@partition_to_py.register(TimestamptzType)
@handle_none
def _(primitive_type: PrimitiveType, value_str: str) -> int:
    """Convert a string to an integer value.

    Raises:
        ValueError: If the scale/exponent is not 0.
    """
    _, _, exponent = Decimal(value_str).as_tuple()
    if exponent != 0:  # Raise if there are digits to the right of the decimal
        raise ValueError(f"Cannot convert partition value, value cannot have fractional digits for {primitive_type} partition")
    return int(float(value_str))


@partition_to_py.register(FloatType)
@partition_to_py.register(DoubleType)
@handle_none
def _(_: PrimitiveType, value_str: str) -> float:
    return float(value_str)


@partition_to_py.register(StringType)
@handle_none
def _(_: StringType, value_str: str) -> str:
    return value_str


@partition_to_py.register(UUIDType)
@handle_none
def _(_: UUIDType, value_str: str) -> uuid.UUID:
    return uuid.UUID(value_str)


@partition_to_py.register(FixedType)
@partition_to_py.register(BinaryType)
@handle_none
def _(_: PrimitiveType, value_str: str) -> bytes:
    return bytes(value_str, UTF8)


@partition_to_py.register(DecimalType)
@handle_none
def _(_: DecimalType, value_str: str) -> Decimal:
    return Decimal(value_str)


@singledispatch
def to_bytes(
    primitive_type: PrimitiveType, _: Union[bool, bytes, Decimal, date, datetime, float, int, str, time, uuid.UUID]
) -> bytes:
    """Convert a built-in python value to bytes.

    This conversion follows the serialization scheme for storing single values as individual binary values defined in the Iceberg specification that
    can be found at https://iceberg.apache.org/spec/#appendix-d-single-value-serialization

    Args:
        primitive_type (PrimitiveType): An implementation of the PrimitiveType base class.
        _: The value to convert to bytes (The type of this value depends on which dispatched function is
            used--check dispatchable functions for type hints).
    """
    raise TypeError(f"scale does not match {primitive_type}")


@to_bytes.register(BooleanType)
def _(_: BooleanType, value: bool) -> bytes:
    return _BOOL_STRUCT.pack(1 if value else 0)


@to_bytes.register(IntegerType)
def _(_: PrimitiveType, value: int) -> bytes:
    return _INT_STRUCT.pack(value)


@to_bytes.register(LongType)
def _(_: PrimitiveType, value: int) -> bytes:
    return _LONG_STRUCT.pack(value)


@to_bytes.register(TimestampType)
@to_bytes.register(TimestamptzType)
def _(_: TimestampType, value: Union[datetime, int]) -> bytes:
    if isinstance(value, datetime):
        value = datetime_to_micros(value)
    return _LONG_STRUCT.pack(value)


@to_bytes.register(DateType)
def _(_: DateType, value: Union[date, int]) -> bytes:
    if isinstance(value, date):
        value = date_to_days(value)
    return _INT_STRUCT.pack(value)


@to_bytes.register(TimeType)
def _(_: TimeType, value: Union[time, int]) -> bytes:
    if isinstance(value, time):
        value = time_to_micros(value)
    return _LONG_STRUCT.pack(value)


@to_bytes.register(FloatType)
def _(_: FloatType, value: float) -> bytes:
    """Convert a float value into bytes.

    Note: float in python is implemented using a double in C. Therefore this involves a conversion of a 32-bit (single precision)
    float to a 64-bit (double precision) float which introduces some imprecision.
    """
    return _FLOAT_STRUCT.pack(value)


@to_bytes.register(DoubleType)
def _(_: DoubleType, value: float) -> bytes:
    return _DOUBLE_STRUCT.pack(value)


@to_bytes.register(StringType)
def _(_: StringType, value: str) -> bytes:
    return value.encode(UTF8)


@to_bytes.register(UUIDType)
def _(_: UUIDType, value: Union[uuid.UUID, bytes]) -> bytes:
    if isinstance(value, bytes):
        return value
    return value.bytes


@to_bytes.register(BinaryType)
@to_bytes.register(FixedType)
def _(_: PrimitiveType, value: bytes) -> bytes:
    return value


@to_bytes.register(DecimalType)
def _(primitive_type: DecimalType, value: Decimal) -> bytes:
    """Convert a Decimal value to bytes given a DecimalType instance with defined precision and scale.

    Args:
        primitive_type (DecimalType): A DecimalType instance with precision and scale.
        value (Decimal): A Decimal instance.

    Raises:
        ValueError: If either the precision or scale of `value` does not match that defined in the DecimalType instance.


    Returns:
        bytes: The byte representation of `value`.
    """
    _, digits, exponent = value.as_tuple()
    exponent = abs(int(exponent))
    if exponent != primitive_type.scale:
        raise ValueError(f"Cannot serialize value, scale of value does not match type {primitive_type}: {exponent}")
    elif len(digits) > primitive_type.precision:
        raise ValueError(
            f"Cannot serialize value, precision of value is greater than precision of type {primitive_type}: {len(digits)}"
        )

    return decimal_to_bytes(value)


@singledispatch  # type: ignore
def from_bytes(primitive_type: PrimitiveType, b: bytes) -> L:  # type: ignore
    """Convert bytes to a built-in python value.

    Args:
        primitive_type (PrimitiveType): An implementation of the PrimitiveType base class.
        b (bytes): The bytes to convert.
    """
    raise TypeError(f"Cannot deserialize bytes, type {primitive_type} not supported: {str(b)}")


@from_bytes.register(BooleanType)
def _(_: BooleanType, b: bytes) -> bool:
    return _BOOL_STRUCT.unpack(b)[0] != 0


@from_bytes.register(IntegerType)
@from_bytes.register(DateType)
def _(_: PrimitiveType, b: bytes) -> int:
    return _INT_STRUCT.unpack(b)[0]


@from_bytes.register(LongType)
@from_bytes.register(TimeType)
@from_bytes.register(TimestampType)
@from_bytes.register(TimestamptzType)
def _(_: PrimitiveType, b: bytes) -> int:
    return _LONG_STRUCT.unpack(b)[0]


@from_bytes.register(FloatType)
def _(_: FloatType, b: bytes) -> float:
    return _FLOAT_STRUCT.unpack(b)[0]


@from_bytes.register(DoubleType)
def _(_: DoubleType, b: bytes) -> float:
    return _DOUBLE_STRUCT.unpack(b)[0]


@from_bytes.register(StringType)
def _(_: StringType, b: bytes) -> str:
    return bytes(b).decode(UTF8)


@from_bytes.register(BinaryType)
@from_bytes.register(FixedType)
@from_bytes.register(UUIDType)
def _(_: PrimitiveType, b: bytes) -> bytes:
    return b


@from_bytes.register(DecimalType)
def _(primitive_type: DecimalType, buf: bytes) -> Decimal:
    unscaled = int.from_bytes(buf, "big", signed=True)
    return unscaled_to_decimal(unscaled, primitive_type.scale)
