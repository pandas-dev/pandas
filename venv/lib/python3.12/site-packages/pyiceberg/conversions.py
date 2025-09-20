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
    - Converting a json-single field serialized field

Note:
    Conversion logic varies based on the PrimitiveType implementation. Therefore conversion functions
    are defined here as generic functions using the @singledispatch decorator. For each PrimitiveType
    implementation, a concrete function is registered for each generic conversion function. For PrimitiveType
    implementations that share the same conversion logic, registrations can be stacked.
"""

import codecs
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
    TimestampNanoType,
    TimestampType,
    TimestamptzNanoType,
    TimestamptzType,
    TimeType,
    UnknownType,
    UUIDType,
    strtobool,
)
from pyiceberg.utils.datetime import (
    date_str_to_days,
    date_to_days,
    datetime_to_micros,
    datetime_to_nanos,
    days_to_date,
    micros_to_time,
    micros_to_timestamp,
    micros_to_timestamptz,
    time_str_to_micros,
    time_to_micros,
    timestamp_to_micros,
    timestamptz_to_micros,
    to_human_day,
    to_human_time,
    to_human_timestamp,
    to_human_timestamptz,
)
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
@partition_to_py.register(TimestampNanoType)
@partition_to_py.register(TimestamptzType)
@partition_to_py.register(TimestamptzNanoType)
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


@partition_to_py.register(UnknownType)
@handle_none
def _(type_: UnknownType, _: str) -> None:
    return None


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
def _(_: PrimitiveType, value: Union[datetime, int]) -> bytes:
    if isinstance(value, datetime):
        value = datetime_to_micros(value)
    return _LONG_STRUCT.pack(value)


@to_bytes.register(TimestampNanoType)
@to_bytes.register(TimestamptzNanoType)
def _(_: PrimitiveType, value: Union[datetime, int]) -> bytes:
    if isinstance(value, datetime):
        value = datetime_to_nanos(value)
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
    raise TypeError(f"Cannot deserialize bytes, type {primitive_type} not supported: {b!r}")


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
@from_bytes.register(TimestampNanoType)
@from_bytes.register(TimestamptzNanoType)
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


@from_bytes.register(UnknownType)
def _(type_: UnknownType, buf: bytes) -> None:
    return None


@singledispatch  # type: ignore
def to_json(primitive_type: PrimitiveType, val: Any) -> L:  # type: ignore
    """Convert built-in python values into JSON value types.

    https://iceberg.apache.org/spec/#json-single-value-serialization

    Args:
        primitive_type (PrimitiveType): An implementation of the PrimitiveType base class.
        val (Any): The arbitrary built-in value to convert into the right form
    """
    raise TypeError(f"Cannot deserialize bytes, type {primitive_type} not supported: {val}")


@to_json.register(BooleanType)
def _(_: BooleanType, val: bool) -> bool:
    """Python bool automatically converts into a JSON bool."""
    return val


@to_json.register(IntegerType)
@to_json.register(LongType)
def _(_: Union[IntegerType, LongType], val: int) -> int:
    """Python int automatically converts to a JSON int."""
    return val


@to_json.register(DateType)
def _(_: DateType, val: Union[date, int]) -> str:
    """JSON date is string encoded."""
    if isinstance(val, date):
        val = date_to_days(val)
    return to_human_day(val)


@to_json.register(TimeType)
def _(_: TimeType, val: Union[int, time]) -> str:
    """Python time or microseconds since epoch serializes into an ISO8601 time."""
    if isinstance(val, time):
        val = time_to_micros(val)
    return to_human_time(val)


@to_json.register(TimestampType)
def _(_: PrimitiveType, val: Union[int, datetime]) -> str:
    """Python datetime (without timezone) or microseconds since epoch serializes into an ISO8601 timestamp."""
    if isinstance(val, datetime):
        val = datetime_to_micros(val)

    return to_human_timestamp(val)


@to_json.register(TimestamptzType)
def _(_: TimestamptzType, val: Union[int, datetime]) -> str:
    """Python datetime (with timezone) or microseconds since epoch serializes into an ISO8601 timestamp."""
    if isinstance(val, datetime):
        val = datetime_to_micros(val)
    return to_human_timestamptz(val)


@to_json.register(FloatType)
@to_json.register(DoubleType)
def _(_: Union[FloatType, DoubleType], val: float) -> float:
    """Float serializes into JSON float."""
    return val


@to_json.register(StringType)
def _(_: StringType, val: str) -> str:
    """Python string serializes into JSON string."""
    return val


@to_json.register(FixedType)
def _(t: FixedType, b: bytes) -> str:
    """Python bytes serializes into hexadecimal encoded string."""
    if len(t) != len(b):
        raise ValueError(f"FixedType has length {len(t)}, which is different from the value: {len(b)}")

    return codecs.encode(b, "hex").decode(UTF8)


@to_json.register(BinaryType)
def _(_: BinaryType, b: bytes) -> str:
    """Python bytes serializes into hexadecimal encoded string."""
    return codecs.encode(b, "hex").decode(UTF8)


@to_json.register(DecimalType)
def _(_: DecimalType, val: Decimal) -> str:
    """Python decimal serializes into string.

    Stores the string representation of the decimal value, specifically, for
    values with a positive scale, the number of digits to the right of the
    decimal point is used to indicate scale, for values with a negative scale,
    the scientific notation is used and the exponent must equal the negated scale.
    """
    return str(val)


@to_json.register(UUIDType)
def _(_: UUIDType, val: uuid.UUID) -> str:
    """Serialize into a JSON string."""
    return str(val)


@singledispatch  # type: ignore
def from_json(primitive_type: PrimitiveType, val: Any) -> L:  # type: ignore
    """Convert JSON value types into built-in python values.

    https://iceberg.apache.org/spec/#json-single-value-serialization

    Args:
        primitive_type (PrimitiveType): An implementation of the PrimitiveType base class.
        val (Any): The arbitrary JSON value to convert into the right form
    """
    raise TypeError(f"Cannot deserialize bytes, type {primitive_type} not supported: {str(val)}")


@from_json.register(BooleanType)
def _(_: BooleanType, val: bool) -> bool:
    """JSON bool automatically converts into a Python bool."""
    return val


@from_json.register(IntegerType)
@from_json.register(LongType)
def _(_: Union[IntegerType, LongType], val: int) -> int:
    """JSON int automatically converts to a Python int."""
    return val


@from_json.register(DateType)
def _(_: DateType, val: Union[str, int, date]) -> date:
    """JSON date is string encoded."""
    if isinstance(val, str):
        val = date_str_to_days(val)
    if isinstance(val, int):
        return days_to_date(val)
    else:
        return val


@from_json.register(TimeType)
def _(_: TimeType, val: Union[str, int, time]) -> time:
    """JSON ISO8601 string into Python time."""
    if isinstance(val, str):
        val = time_str_to_micros(val)
    if isinstance(val, int):
        return micros_to_time(val)
    else:
        return val


@from_json.register(TimestampType)
def _(_: PrimitiveType, val: Union[str, int, datetime]) -> datetime:
    """JSON ISO8601 string into Python datetime."""
    if isinstance(val, str):
        val = timestamp_to_micros(val)
    if isinstance(val, int):
        return micros_to_timestamp(val)
    else:
        return val


@from_json.register(TimestamptzType)
def _(_: TimestamptzType, val: Union[str, int, datetime]) -> datetime:
    """JSON ISO8601 string into Python datetime."""
    if isinstance(val, str):
        val = timestamptz_to_micros(val)
    if isinstance(val, int):
        return micros_to_timestamptz(val)
    else:
        return val


@from_json.register(FloatType)
@from_json.register(DoubleType)
def _(_: Union[FloatType, DoubleType], val: float) -> float:
    """JSON float deserializes into a Python float."""
    return val


@from_json.register(StringType)
def _(_: StringType, val: str) -> str:
    """JSON string serializes into a Python string."""
    return val


@from_json.register(FixedType)
def _(t: FixedType, val: Union[str, bytes]) -> bytes:
    """JSON hexadecimal encoded string into bytes."""
    if isinstance(val, str):
        val = codecs.decode(val.encode(UTF8), "hex")

    if len(t) != len(val):
        raise ValueError(f"FixedType has length {len(t)}, which is different from the value: {len(val)}")

    return val


@from_json.register(BinaryType)
def _(_: BinaryType, val: Union[bytes, str]) -> bytes:
    """JSON hexadecimal encoded string into bytes."""
    if isinstance(val, str):
        return codecs.decode(val.encode(UTF8), "hex")
    else:
        return val


@from_json.register(DecimalType)
def _(_: DecimalType, val: str) -> Decimal:
    """Convert JSON string into a Python Decimal."""
    return Decimal(val)


@from_json.register(UUIDType)
def _(_: UUIDType, val: Union[str, bytes, uuid.UUID]) -> uuid.UUID:
    """Convert JSON string into Python UUID."""
    if isinstance(val, str):
        return uuid.UUID(val)
    elif isinstance(val, bytes):
        return uuid.UUID(bytes=val)
    else:
        return val
