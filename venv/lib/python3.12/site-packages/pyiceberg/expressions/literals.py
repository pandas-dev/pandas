#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
# pylint: disable=W0613
from __future__ import annotations

import struct
from abc import ABC, abstractmethod
from datetime import date, datetime, time
from decimal import ROUND_HALF_UP, Decimal
from functools import singledispatchmethod
from math import isnan
from typing import Any, Generic, Type
from uuid import UUID

from pyiceberg.typedef import L
from pyiceberg.types import (
    BinaryType,
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    FixedType,
    FloatType,
    IcebergType,
    IntegerType,
    LongType,
    StringType,
    TimestampType,
    TimestamptzType,
    TimeType,
    UUIDType,
)
from pyiceberg.utils.datetime import (
    date_str_to_days,
    date_to_days,
    datetime_to_micros,
    micros_to_days,
    time_str_to_micros,
    time_to_micros,
    timestamp_to_micros,
    timestamptz_to_micros,
)
from pyiceberg.utils.decimal import decimal_to_unscaled, unscaled_to_decimal
from pyiceberg.utils.singleton import Singleton

UUID_BYTES_LENGTH = 16


class Literal(Generic[L], ABC):
    """Literal which has a value and can be converted between types."""

    _value: L

    def __init__(self, value: L, value_type: Type[L]):
        if value is None or not isinstance(value, value_type):
            raise TypeError(f"Invalid literal value: {value!r} (not a {value_type})")
        if isinstance(value, float) and isnan(value):
            raise ValueError("Cannot create expression literal from NaN.")
        self._value = value

    @property
    def value(self) -> L:
        return self._value

    @singledispatchmethod
    @abstractmethod
    def to(self, type_var: IcebergType) -> Literal[L]: ...  # pragma: no cover

    def __repr__(self) -> str:
        """Return the string representation of the Literal class."""
        return f"{type(self).__name__}({self.value!r})"

    def __str__(self) -> str:
        """Return the string representation of the Literal class."""
        return str(self.value)

    def __hash__(self) -> int:
        """Return a hashed representation of the Literal class."""
        return hash(self.value)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Literal class."""
        if not isinstance(other, Literal):
            return False
        return self.value == other.value

    def __ne__(self, other: Any) -> bool:
        """Return the inequality of two instances of the Literal class."""
        return not self.__eq__(other)

    def __lt__(self, other: Any) -> bool:
        """Return if one instance of the Literal class is less than another instance."""
        return self.value < other.value

    def __gt__(self, other: Any) -> bool:
        """Return if one instance of the Literal class is greater than another instance."""
        return self.value > other.value

    def __le__(self, other: Any) -> bool:
        """Return if one instance of the Literal class is less than or equal to another instance."""
        return self.value <= other.value

    def __ge__(self, other: Any) -> bool:
        """Return if one instance of the Literal class is greater than or equal to another instance."""
        return self.value >= other.value


def literal(value: L) -> Literal[L]:
    """
    Construct an Iceberg Literal based on Python primitive data type.

    Args:
        value (Python primitive type): the value to be associated with literal.

    Example:
        from pyiceberg.expressions.literals import literal.
        >>> literal(123)
        LongLiteral(123)
    """
    if isinstance(value, float):
        return DoubleLiteral(value)  # type: ignore
    elif isinstance(value, bool):
        return BooleanLiteral(value)
    elif isinstance(value, int):
        return LongLiteral(value)
    elif isinstance(value, str):
        return StringLiteral(value)
    elif isinstance(value, UUID):
        return UUIDLiteral(value.bytes)  # type: ignore
    elif isinstance(value, bytes):
        return BinaryLiteral(value)
    elif isinstance(value, Decimal):
        return DecimalLiteral(value)
    elif isinstance(value, datetime):
        return TimestampLiteral(datetime_to_micros(value))  # type: ignore
    elif isinstance(value, date):
        return DateLiteral(date_to_days(value))  # type: ignore
    elif isinstance(value, time):
        return TimeLiteral(time_to_micros(value))  # type: ignore
    else:
        raise TypeError(f"Invalid literal value: {repr(value)}")


class AboveMax(Literal[L]):
    def __repr__(self) -> str:
        """Return the string representation of the AboveMax class."""
        return f"{self.__class__.__name__}()"

    def __str__(self) -> str:
        """Return the string representation of the AboveMax class."""
        return self.__class__.__name__


class BelowMin(Literal[L]):
    def __repr__(self) -> str:
        """Return the string representation of the BelowMin class."""
        return f"{self.__class__.__name__}()"

    def __str__(self) -> str:
        """Return the string representation of the BelowMin class."""
        return self.__class__.__name__


class FloatAboveMax(AboveMax[float], Singleton):
    def __init__(self) -> None:
        super().__init__(FloatType.max, float)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError("Cannot change the type of FloatAboveMax")

    @to.register(FloatType)
    def _(self, _: FloatType) -> Literal[float]:
        return self


class FloatBelowMin(BelowMin[float], Singleton):
    def __init__(self) -> None:
        super().__init__(FloatType.min, float)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError("Cannot change the type of FloatBelowMin")

    @to.register(FloatType)
    def _(self, _: FloatType) -> Literal[float]:
        return self


class IntAboveMax(AboveMax[int], Singleton):
    def __init__(self) -> None:
        super().__init__(IntegerType.max, int)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError("Cannot change the type of IntAboveMax")

    @to.register(IntegerType)
    def _(self, _: IntegerType) -> Literal[int]:
        return self


class IntBelowMin(BelowMin[int], Singleton):
    def __init__(self) -> None:
        super().__init__(IntegerType.min, int)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError("Cannot change the type of IntBelowMin")

    @to.register(IntegerType)
    def _(self, _: IntegerType) -> Literal[int]:
        return self


class LongAboveMax(AboveMax[int], Singleton):
    def __init__(self) -> None:
        super().__init__(LongType.max, int)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError("Cannot change the type of IntAboveMax")

    @to.register(LongType)
    def _(self, _: LongType) -> Literal[int]:
        return self


class LongBelowMin(BelowMin[int], Singleton):
    def __init__(self) -> None:
        super().__init__(LongType.min, int)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError("Cannot change the type of IntBelowMin")

    @to.register(LongType)
    def _(self, _: LongType) -> Literal[int]:
        return self


class BooleanLiteral(Literal[bool]):
    def __init__(self, value: bool) -> None:
        super().__init__(value, bool)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal[bool]:
        raise TypeError(f"Cannot convert BooleanLiteral into {type_var}")

    @to.register(BooleanType)
    def _(self, _: BooleanType) -> Literal[bool]:
        return self


class LongLiteral(Literal[int]):
    def __init__(self, value: int) -> None:
        super().__init__(value, int)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert LongLiteral into {type_var}")

    def increment(self) -> Literal[int]:
        return LongLiteral(self.value + 1)

    def decrement(self) -> Literal[int]:
        return LongLiteral(self.value - 1)

    @to.register(LongType)
    def _(self, _: LongType) -> Literal[int]:
        if LongType.max < self.value:
            return LongAboveMax()
        elif LongType.min > self.value:
            return LongBelowMin()
        else:
            return self

    @to.register(IntegerType)
    def _(self, _: IntegerType) -> Literal[int]:
        if IntegerType.max < self.value:
            return IntAboveMax()
        elif IntegerType.min > self.value:
            return IntBelowMin()
        return self

    @to.register(FloatType)
    def _(self, _: FloatType) -> Literal[float]:
        return FloatLiteral(float(self.value))

    @to.register(DoubleType)
    def _(self, _: DoubleType) -> Literal[float]:
        return DoubleLiteral(float(self.value))

    @to.register(DateType)
    def _(self, _: DateType) -> Literal[int]:
        return DateLiteral(self.value)

    @to.register(TimeType)
    def _(self, _: TimeType) -> Literal[int]:
        return TimeLiteral(self.value)

    @to.register(TimestampType)
    def _(self, _: TimestampType) -> Literal[int]:
        return TimestampLiteral(self.value)

    @to.register(TimestamptzType)
    def _(self, _: TimestamptzType) -> Literal[int]:
        return TimestampLiteral(self.value)

    @to.register(DecimalType)
    def _(self, type_var: DecimalType) -> Literal[Decimal]:
        unscaled = Decimal(self.value)
        if type_var.scale == 0:
            return DecimalLiteral(unscaled)
        else:
            sign, digits, _ = unscaled.as_tuple()
            zeros = (0,) * type_var.scale
            return DecimalLiteral(Decimal((sign, digits + zeros, -type_var.scale)))


class FloatLiteral(Literal[float]):
    def __init__(self, value: float) -> None:
        super().__init__(value, float)
        self._value32 = struct.unpack("<f", struct.pack("<f", value))[0]

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the FloatLiteral class."""
        return self._value32 == other

    def __lt__(self, other: Any) -> bool:
        """Return if one instance of the FloatLiteral class is less than another instance."""
        return self._value32 < other

    def __gt__(self, other: Any) -> bool:
        """Return if one instance of the FloatLiteral class is greater than another instance."""
        return self._value32 > other

    def __le__(self, other: Any) -> bool:
        """Return if one instance of the FloatLiteral class is less than or equal to another instance."""
        return self._value32 <= other

    def __ge__(self, other: Any) -> bool:
        """Return if one instance of the FloatLiteral class is greater than or equal to another instance."""
        return self._value32 >= other

    def __hash__(self) -> int:
        """Return a hashed representation of the FloatLiteral class."""
        return hash(self._value32)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert FloatLiteral into {type_var}")

    @to.register(FloatType)
    def _(self, _: FloatType) -> Literal[float]:
        return self

    @to.register(DoubleType)
    def _(self, _: DoubleType) -> Literal[float]:
        return DoubleLiteral(self.value)

    @to.register(DecimalType)
    def _(self, type_var: DecimalType) -> Literal[Decimal]:
        return DecimalLiteral(Decimal(self.value).quantize(Decimal((0, (1,), -type_var.scale)), rounding=ROUND_HALF_UP))


class DoubleLiteral(Literal[float]):
    def __init__(self, value: float) -> None:
        super().__init__(value, float)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert DoubleLiteral into {type_var}")

    @to.register(DoubleType)
    def _(self, _: DoubleType) -> Literal[float]:
        return self

    @to.register(FloatType)
    def _(self, _: FloatType) -> Literal[float]:
        if FloatType.max < self.value:
            return FloatAboveMax()
        elif FloatType.min > self.value:
            return FloatBelowMin()
        return FloatLiteral(self.value)

    @to.register(DecimalType)
    def _(self, type_var: DecimalType) -> Literal[Decimal]:
        return DecimalLiteral(Decimal(self.value).quantize(Decimal((0, (1,), -type_var.scale)), rounding=ROUND_HALF_UP))


class DateLiteral(Literal[int]):
    def __init__(self, value: int) -> None:
        super().__init__(value, int)

    def increment(self) -> Literal[int]:
        return DateLiteral(self.value + 1)

    def decrement(self) -> Literal[int]:
        return DateLiteral(self.value - 1)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert DateLiteral into {type_var}")

    @to.register(DateType)
    def _(self, _: DateType) -> Literal[int]:
        return self


class TimeLiteral(Literal[int]):
    def __init__(self, value: int) -> None:
        super().__init__(value, int)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert TimeLiteral into {type_var}")

    @to.register(TimeType)
    def _(self, _: TimeType) -> Literal[int]:
        return self


class TimestampLiteral(Literal[int]):
    def __init__(self, value: int) -> None:
        super().__init__(value, int)

    def increment(self) -> Literal[int]:
        return TimestampLiteral(self.value + 1)

    def decrement(self) -> Literal[int]:
        return TimestampLiteral(self.value - 1)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert TimestampLiteral into {type_var}")

    @to.register(TimestampType)
    def _(self, _: TimestampType) -> Literal[int]:
        return self

    @to.register(TimestamptzType)
    def _(self, _: TimestamptzType) -> Literal[int]:
        return self

    @to.register(DateType)
    def _(self, _: DateType) -> Literal[int]:
        return DateLiteral(micros_to_days(self.value))


class DecimalLiteral(Literal[Decimal]):
    def __init__(self, value: Decimal) -> None:
        super().__init__(value, Decimal)

    def increment(self) -> Literal[Decimal]:
        original_scale = abs(int(self.value.as_tuple().exponent))
        unscaled = decimal_to_unscaled(self.value)
        return DecimalLiteral(unscaled_to_decimal(unscaled + 1, original_scale))

    def decrement(self) -> Literal[Decimal]:
        original_scale = abs(int(self.value.as_tuple().exponent))
        unscaled = decimal_to_unscaled(self.value)
        return DecimalLiteral(unscaled_to_decimal(unscaled - 1, original_scale))

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert DecimalLiteral into {type_var}")

    @to.register(DecimalType)
    def _(self, type_var: DecimalType) -> Literal[Decimal]:
        if type_var.scale == abs(int(self.value.as_tuple().exponent)):
            return self
        raise ValueError(f"Could not convert {self.value} into a {type_var}")

    @to.register(IntegerType)
    def _(self, _: IntegerType) -> Literal[int]:
        value_int = int(self.value.to_integral_value())
        if value_int > IntegerType.max:
            return IntAboveMax()
        elif value_int < IntegerType.min:
            return IntBelowMin()
        else:
            return LongLiteral(value_int)

    @to.register(LongType)
    def _(self, _: LongType) -> Literal[int]:
        value_int = int(self.value.to_integral_value())
        if value_int > LongType.max:
            return IntAboveMax()
        elif value_int < LongType.min:
            return IntBelowMin()
        else:
            return LongLiteral(value_int)

    @to.register(FloatType)
    def _(self, _: FloatType) -> Literal[float]:
        value_float = float(self.value)
        if value_float > FloatType.max:
            return FloatAboveMax()
        elif value_float < FloatType.min:
            return FloatBelowMin()
        else:
            return FloatLiteral(value_float)

    @to.register(DoubleType)
    def _(self, _: DoubleLiteral) -> Literal[float]:
        return DoubleLiteral(float(self.value))


class StringLiteral(Literal[str]):
    def __init__(self, value: str) -> None:
        super().__init__(value, str)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert StringLiteral into {type_var}")

    @to.register(StringType)
    def _(self, _: StringType) -> Literal[str]:
        return self

    @to.register(IntegerType)
    def _(self, type_var: IntegerType) -> Literal[int]:
        try:
            number = int(float(self.value))

            if IntegerType.max < number:
                return IntAboveMax()
            elif IntegerType.min > number:
                return IntBelowMin()
            return LongLiteral(number)
        except ValueError as e:
            raise ValueError(f"Could not convert {self.value} into a {type_var}") from e

    @to.register(LongType)
    def _(self, type_var: LongType) -> Literal[int]:
        try:
            long_value = int(float(self.value))
            if LongType.max < long_value:
                return LongAboveMax()
            elif LongType.min > long_value:
                return LongBelowMin()
            else:
                return LongLiteral(long_value)
        except (TypeError, ValueError) as e:
            raise ValueError(f"Could not convert {self.value} into a {type_var}") from e

    @to.register(DateType)
    def _(self, type_var: DateType) -> Literal[int]:
        try:
            return DateLiteral(date_str_to_days(self.value))
        except (TypeError, ValueError) as e:
            raise ValueError(f"Could not convert {self.value} into a {type_var}") from e

    @to.register(TimeType)
    def _(self, type_var: TimeType) -> Literal[int]:
        try:
            return TimeLiteral(time_str_to_micros(self.value))
        except (TypeError, ValueError) as e:
            raise ValueError(f"Could not convert {self.value} into a {type_var}") from e

    @to.register(TimestampType)
    def _(self, _: TimestampType) -> Literal[int]:
        return TimestampLiteral(timestamp_to_micros(self.value))

    @to.register(TimestamptzType)
    def _(self, _: TimestamptzType) -> Literal[int]:
        return TimestampLiteral(timestamptz_to_micros(self.value))

    @to.register(UUIDType)
    def _(self, _: UUIDType) -> Literal[bytes]:
        return UUIDLiteral(UUID(self.value).bytes)

    @to.register(DecimalType)
    def _(self, type_var: DecimalType) -> Literal[Decimal]:
        dec = Decimal(self.value)
        scale = abs(int(dec.as_tuple().exponent))
        if type_var.scale == scale:
            return DecimalLiteral(dec)
        else:
            raise ValueError(f"Could not convert {self.value} into a {type_var}, scales differ {type_var.scale} <> {scale}")

    @to.register(BooleanType)
    def _(self, type_var: BooleanType) -> Literal[bool]:
        value_upper = self.value.upper()
        if value_upper in ["TRUE", "FALSE"]:
            return BooleanLiteral(value_upper == "TRUE")
        else:
            raise ValueError(f"Could not convert {self.value} into a {type_var}")

    @to.register(FloatType)
    def _(self, type_var: FloatType) -> Literal[float]:
        try:
            number = float(self.value)
            if FloatType.max < number:
                return FloatAboveMax()
            elif FloatType.min > number:
                return FloatBelowMin()
            return FloatLiteral(number)
        except ValueError as e:
            raise ValueError(f"Could not convert {self.value} into a {type_var}") from e

    @to.register(DoubleType)
    def _(self, type_var: DoubleType) -> Literal[float]:
        try:
            number = float(self.value)
            return DoubleLiteral(number)
        except ValueError as e:
            raise ValueError(f"Could not convert {self.value} into a {type_var}") from e

    def __repr__(self) -> str:
        """Return the string representation of the StringLiteral class."""
        return f"literal({repr(self.value)})"


class UUIDLiteral(Literal[bytes]):
    def __init__(self, value: bytes) -> None:
        super().__init__(value, bytes)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert UUIDLiteral into {type_var}")

    @to.register(UUIDType)
    def _(self, _: UUIDType) -> Literal[bytes]:
        return self

    @to.register(FixedType)
    def _(self, type_var: FixedType) -> Literal[bytes]:
        if len(type_var) == UUID_BYTES_LENGTH:
            return FixedLiteral(self.value)
        else:
            raise TypeError(
                f"Cannot convert UUIDLiteral into {type_var}, different length: {len(type_var)} <> {UUID_BYTES_LENGTH}"
            )

    @to.register(BinaryType)
    def _(self, _: BinaryType) -> Literal[bytes]:
        return BinaryLiteral(self.value)


class FixedLiteral(Literal[bytes]):
    def __init__(self, value: bytes) -> None:
        super().__init__(value, bytes)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert FixedLiteral into {type_var}")

    @to.register(FixedType)
    def _(self, type_var: FixedType) -> Literal[bytes]:
        if len(self.value) == len(type_var):
            return self
        else:
            raise ValueError(
                f"Could not convert {self.value!r} into a {type_var}, lengths differ {len(self.value)} <> {len(type_var)}"
            )

    @to.register(BinaryType)
    def _(self, _: BinaryType) -> Literal[bytes]:
        return BinaryLiteral(self.value)

    @to.register(UUIDType)
    def _(self, type_var: UUIDType) -> Literal[bytes]:
        if len(self.value) == UUID_BYTES_LENGTH:
            return UUIDLiteral(self.value)
        else:
            raise TypeError(
                f"Could not convert {self.value!r} into a {type_var}, lengths differ {len(self.value)} <> {UUID_BYTES_LENGTH}"
            )


class BinaryLiteral(Literal[bytes]):
    def __init__(self, value: bytes) -> None:
        super().__init__(value, bytes)

    @singledispatchmethod
    def to(self, type_var: IcebergType) -> Literal:  # type: ignore
        raise TypeError(f"Cannot convert BinaryLiteral into {type_var}")

    @to.register(BinaryType)
    def _(self, _: BinaryType) -> Literal[bytes]:
        return self

    @to.register(FixedType)
    def _(self, type_var: FixedType) -> Literal[bytes]:
        if len(type_var) == len(self.value):
            return FixedLiteral(self.value)
        else:
            raise TypeError(
                f"Cannot convert BinaryLiteral into {type_var}, different length: {len(type_var)} <> {len(self.value)}"
            )

    @to.register(UUIDType)
    def _(self, type_var: UUIDType) -> Literal[bytes]:
        if len(self.value) == UUID_BYTES_LENGTH:
            return UUIDLiteral(self.value)
        else:
            raise TypeError(
                f"Cannot convert BinaryLiteral into {type_var}, different length: {UUID_BYTES_LENGTH} <> {len(self.value)}"
            )
