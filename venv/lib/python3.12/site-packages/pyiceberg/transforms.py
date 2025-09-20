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

import base64
import datetime as py_datetime
import importlib
import struct
import types
from abc import ABC, abstractmethod
from enum import IntEnum
from functools import singledispatch
from typing import TYPE_CHECKING, Any, Callable, Generic, Optional, TypeVar
from typing import Literal as LiteralType
from uuid import UUID

import mmh3
from pydantic import Field, PositiveInt, PrivateAttr

from pyiceberg.exceptions import NotInstalledError
from pyiceberg.expressions import (
    BoundEqualTo,
    BoundGreaterThan,
    BoundGreaterThanOrEqual,
    BoundIn,
    BoundLessThan,
    BoundLessThanOrEqual,
    BoundLiteralPredicate,
    BoundNotEqualTo,
    BoundNotIn,
    BoundNotStartsWith,
    BoundPredicate,
    BoundSetPredicate,
    BoundStartsWith,
    BoundTerm,
    BoundUnaryPredicate,
    EqualTo,
    GreaterThan,
    GreaterThanOrEqual,
    LessThan,
    LessThanOrEqual,
    NotEqualTo,
    NotStartsWith,
    Reference,
    StartsWith,
    UnboundPredicate,
)
from pyiceberg.expressions.literals import (
    DateLiteral,
    DecimalLiteral,
    Literal,
    LongLiteral,
    TimestampLiteral,
    literal,
)
from pyiceberg.typedef import IcebergRootModel, L
from pyiceberg.types import (
    BinaryType,
    DateType,
    DecimalType,
    FixedType,
    IcebergType,
    IntegerType,
    LongType,
    StringType,
    TimestampNanoType,
    TimestampType,
    TimestamptzNanoType,
    TimestamptzType,
    TimeType,
    UUIDType,
)
from pyiceberg.utils import datetime
from pyiceberg.utils.decimal import decimal_to_bytes, truncate_decimal
from pyiceberg.utils.parsing import ParseNumberFromBrackets
from pyiceberg.utils.singleton import Singleton

if TYPE_CHECKING:
    import pyarrow as pa

    ArrayLike = TypeVar("ArrayLike", pa.Array, pa.ChunkedArray)

S = TypeVar("S")
T = TypeVar("T")

IDENTITY = "identity"
VOID = "void"
BUCKET = "bucket"
TRUNCATE = "truncate"
YEAR = "year"
MONTH = "month"
DAY = "day"
HOUR = "hour"

BUCKET_PARSER = ParseNumberFromBrackets(BUCKET)
TRUNCATE_PARSER = ParseNumberFromBrackets(TRUNCATE)


def _try_import(module_name: str, extras_name: Optional[str] = None) -> types.ModuleType:
    try:
        return importlib.import_module(module_name)
    except ImportError:
        if extras_name:
            msg = f'{module_name} needs to be installed. pip install "pyiceberg[{extras_name}]"'
        else:
            msg = f"{module_name} needs to be installed."
        raise NotInstalledError(msg) from None


def _transform_literal(func: Callable[[L], L], lit: Literal[L]) -> Literal[L]:
    """Small helper to upwrap the value from the literal, and wrap it again."""
    return literal(func(lit.value))


def _pyiceberg_transform_wrapper(
    transform_func: Callable[["ArrayLike", Any], "ArrayLike"],
    *args: Any,
    expected_type: Optional["pa.DataType"] = None,
) -> Callable[["ArrayLike"], "ArrayLike"]:
    try:
        import pyarrow as pa
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError("For partition transforms, PyArrow needs to be installed") from e

    def _transform(array: "ArrayLike") -> "ArrayLike":
        def _cast_if_needed(arr: "ArrayLike") -> "ArrayLike":
            if expected_type is not None:
                return arr.cast(expected_type)
            else:
                return arr

        if isinstance(array, pa.Array):
            return _cast_if_needed(transform_func(array, *args))
        elif isinstance(array, pa.ChunkedArray):
            result_chunks = []
            for arr in array.iterchunks():
                result_chunks.append(_cast_if_needed(transform_func(arr, *args)))
            return pa.chunked_array(result_chunks)
        else:
            raise ValueError(f"PyArrow array can only be of type pa.Array or pa.ChunkedArray, but found {type(array)}")

    return _transform


class Transform(IcebergRootModel[str], ABC, Generic[S, T]):
    """Transform base class for concrete transforms.

    A base class to transform values and project predicates on partition values.
    This class is not used directly. Instead, use one of module method to create the child classes.
    """

    root: str = Field()

    @abstractmethod
    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[T]]: ...

    @abstractmethod
    def can_transform(self, source: IcebergType) -> bool:
        return False

    @abstractmethod
    def result_type(self, source: IcebergType) -> IcebergType:
        """Return the `IcebergType` produced by this transform given a source type.

        This method defines both the physical and display representation of the partition field.

        The physical representation must conform to the Iceberg spec. The display representation
        can deviate from the spec, such as by transforming the value into a more human-readable format.
        """
        ...

    @abstractmethod
    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]: ...

    @abstractmethod
    def strict_project(self, name: str, pred: BoundPredicate[Any]) -> Optional[UnboundPredicate[Any]]: ...

    @property
    def preserves_order(self) -> bool:
        return False

    def satisfies_order_of(self, other: Any) -> bool:
        return self == other

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        return str(value) if value is not None else "null"

    @property
    def dedup_name(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        """Return the string representation of the Transform class."""
        return self.root

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Transform class."""
        if isinstance(other, Transform):
            return self.root == other.root
        return False

    @abstractmethod
    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]": ...


def parse_transform(v: Any) -> Transform[Any, Any]:
    if isinstance(v, str):
        if v == IDENTITY:
            return IdentityTransform()
        elif v == VOID:
            return VoidTransform()
        elif v.startswith(BUCKET):
            return BucketTransform(num_buckets=BUCKET_PARSER.match(v))
        elif v.startswith(TRUNCATE):
            return TruncateTransform(width=TRUNCATE_PARSER.match(v))
        elif v == YEAR:
            return YearTransform()
        elif v == MONTH:
            return MonthTransform()
        elif v == DAY:
            return DayTransform()
        elif v == HOUR:
            return HourTransform()
        else:
            return UnknownTransform(transform=v)
    return v


class BucketTransform(Transform[S, int]):
    """Base Transform class to transform a value into a bucket partition value.

    Transforms are parameterized by a number of buckets. Bucket partition transforms use a 32-bit
    hash of the source value to produce a positive value by mod the bucket number.

    Args:
      num_buckets (int): The number of buckets.
    """

    root: str = Field()
    _num_buckets: PositiveInt = PrivateAttr()

    def __init__(self, num_buckets: int, **data: Any) -> None:
        super().__init__(f"bucket[{num_buckets}]", **data)
        self._num_buckets = num_buckets

    @property
    def num_buckets(self) -> int:
        return self._num_buckets

    def hash(self, value: S) -> int:
        raise NotImplementedError()

    def apply(self, value: Optional[S]) -> Optional[int]:
        return (self.hash(value) & IntegerType.max) % self._num_buckets if value else None

    def result_type(self, source: IcebergType) -> IcebergType:
        return IntegerType()

    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        transformer = self.transform(pred.term.ref().field.field_type)

        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)
        elif isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundEqualTo):
            return pred.as_unbound(Reference(name), _transform_literal(transformer, pred.literal))
        elif isinstance(pred, BoundIn):  # NotIn can't be projected
            return pred.as_unbound(Reference(name), {_transform_literal(transformer, literal) for literal in pred.literals})
        else:
            # - Comparison predicates can't be projected, notEq can't be projected
            # - Small ranges can be projected:
            #   For example, (x > 0) and (x < 3) can be turned into in({1, 2}) and projected.
            return None

    def strict_project(self, name: str, pred: BoundPredicate[Any]) -> Optional[UnboundPredicate[Any]]:
        transformer = self.transform(pred.term.ref().field.field_type)

        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)
        elif isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundNotEqualTo):
            return pred.as_unbound(Reference(name), _transform_literal(transformer, pred.literal))
        elif isinstance(pred, BoundNotIn):
            return pred.as_unbound(Reference(name), {_transform_literal(transformer, literal) for literal in pred.literals})
        else:
            # no strict projection for comparison or equality
            return None

    def can_transform(self, source: IcebergType) -> bool:
        return isinstance(
            source,
            (
                IntegerType,
                DateType,
                LongType,
                TimeType,
                TimestampType,
                TimestamptzType,
                TimestampNanoType,
                TimestamptzNanoType,
                DecimalType,
                StringType,
                FixedType,
                BinaryType,
                UUIDType,
            ),
        )

    def transform(self, source: IcebergType, bucket: bool = True) -> Callable[[Optional[Any]], Optional[int]]:
        if isinstance(source, TimeType):

            def hash_func(v: Any) -> int:
                if isinstance(v, py_datetime.time):
                    v = datetime.time_to_micros(v)

                return mmh3.hash(struct.pack("<q", v))

        elif isinstance(source, DateType):

            def hash_func(v: Any) -> int:
                if isinstance(v, py_datetime.date):
                    v = datetime.date_to_days(v)

                return mmh3.hash(struct.pack("<q", v))

        elif isinstance(source, (TimestampType, TimestamptzType)):

            def hash_func(v: Any) -> int:
                if isinstance(v, py_datetime.datetime):
                    v = datetime.datetime_to_micros(v)

                return mmh3.hash(struct.pack("<q", v))

        elif isinstance(source, (TimestampNanoType, TimestamptzNanoType)):

            def hash_func(v: Any) -> int:
                # In order to bucket TimestampNano the same as Timestamp
                # convert to micros before hashing.
                if isinstance(v, py_datetime.datetime):
                    v = datetime.datetime_to_micros(v)
                else:
                    v = datetime.nanos_to_micros(v)

                return mmh3.hash(struct.pack("<q", v))

        elif isinstance(source, (IntegerType, LongType)):

            def hash_func(v: Any) -> int:
                return mmh3.hash(struct.pack("<q", v))

        elif isinstance(source, DecimalType):

            def hash_func(v: Any) -> int:
                return mmh3.hash(decimal_to_bytes(v))

        elif isinstance(source, (StringType, FixedType, BinaryType)):

            def hash_func(v: Any) -> int:
                return mmh3.hash(v)

        elif isinstance(source, UUIDType):

            def hash_func(v: Any) -> int:
                if isinstance(v, UUID):
                    return mmh3.hash(v.bytes)
                return mmh3.hash(v)

        else:
            raise ValueError(f"Unknown type {source}")

        if bucket:
            return lambda v: (hash_func(v) & IntegerType.max) % self._num_buckets if v is not None else None
        return hash_func

    def __repr__(self) -> str:
        """Return the string representation of the BucketTransform class."""
        return f"BucketTransform(num_buckets={self._num_buckets})"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        pyiceberg_core_transform = _try_import("pyiceberg_core", extras_name="pyiceberg-core").transform
        return _pyiceberg_transform_wrapper(pyiceberg_core_transform.bucket, self._num_buckets)


class TimeResolution(IntEnum):
    YEAR = 6
    MONTH = 5
    WEEK = 4
    DAY = 3
    HOUR = 2
    MINUTE = 1
    SECOND = 0


class TimeTransform(Transform[S, int], Generic[S], Singleton):
    @property
    @abstractmethod
    def granularity(self) -> TimeResolution: ...

    def satisfies_order_of(self, other: Transform[S, T]) -> bool:
        return self.granularity <= other.granularity if hasattr(other, "granularity") else False

    def result_type(self, source: IcebergType) -> IntegerType:
        return IntegerType()

    @abstractmethod
    def transform(self, source: IcebergType) -> Callable[[Optional[Any]], Optional[int]]: ...

    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        transformer = self.transform(pred.term.ref().field.field_type)
        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)
        elif isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundLiteralPredicate):
            return _truncate_number(name, pred, transformer)
        elif isinstance(pred, BoundIn):  # NotIn can't be projected
            return _set_apply_transform(name, pred, transformer)
        else:
            return None

    def strict_project(self, name: str, pred: BoundPredicate[Any]) -> Optional[UnboundPredicate[Any]]:
        transformer = self.transform(pred.term.ref().field.field_type)
        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)
        elif isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundLiteralPredicate):
            return _truncate_number_strict(name, pred, transformer)
        elif isinstance(pred, BoundNotIn):
            return _set_apply_transform(name, pred, transformer)
        else:
            return None

    @property
    def dedup_name(self) -> str:
        return "time"

    @property
    def preserves_order(self) -> bool:
        return True


class YearTransform(TimeTransform[S]):
    """Transforms a datetime value into a year value.

    Example:
        >>> transform = YearTransform()
        >>> transform.transform(TimestampType())(1512151975038194)
        47
    """

    root: LiteralType["year"] = Field(default="year")  # noqa: F821

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[int]]:
        if isinstance(source, DateType):

            def year_func(v: Any) -> int:
                if isinstance(v, py_datetime.date):
                    v = datetime.date_to_days(v)

                return datetime.days_to_years(v)

        elif isinstance(source, (TimestampType, TimestamptzType)):

            def year_func(v: Any) -> int:
                if isinstance(v, py_datetime.datetime):
                    v = datetime.datetime_to_micros(v)

                return datetime.micros_to_years(v)

        elif isinstance(source, (TimestampNanoType, TimestamptzNanoType)):

            def year_func(v: Any) -> int:
                # python datetime has no nanoseconds support.
                # nanosecond datetimes will be expressed as int as a workaround
                return datetime.nanos_to_years(v)

        else:
            raise ValueError(f"Cannot apply year transform for type: {source}")

        return lambda v: year_func(v) if v is not None else None

    def can_transform(self, source: IcebergType) -> bool:
        return isinstance(source, (DateType, TimestampType, TimestamptzType, TimestampNanoType, TimestamptzNanoType))

    @property
    def granularity(self) -> TimeResolution:
        return TimeResolution.YEAR

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        return datetime.to_human_year(value) if isinstance(value, int) else "null"

    def __repr__(self) -> str:
        """Return the string representation of the YearTransform class."""
        return "YearTransform()"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        pa = _try_import("pyarrow")
        pyiceberg_core_transform = _try_import("pyiceberg_core", extras_name="pyiceberg-core").transform
        return _pyiceberg_transform_wrapper(pyiceberg_core_transform.year, expected_type=pa.int32())


class MonthTransform(TimeTransform[S]):
    """Transforms a datetime value into a month value.

    Example:
        >>> transform = MonthTransform()
        >>> transform.transform(DateType())(17501)
        575
    """

    root: LiteralType["month"] = Field(default="month")  # noqa: F821

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[int]]:
        if isinstance(source, DateType):

            def month_func(v: Any) -> int:
                if isinstance(v, py_datetime.date):
                    v = datetime.date_to_days(v)

                return datetime.days_to_months(v)

        elif isinstance(source, (TimestampType, TimestamptzType)):

            def month_func(v: Any) -> int:
                if isinstance(v, py_datetime.datetime):
                    v = datetime.datetime_to_micros(v)

                return datetime.micros_to_months(v)

        elif isinstance(source, (TimestampNanoType, TimestamptzNanoType)):

            def month_func(v: Any) -> int:
                # python datetime has no nanoseconds support.
                # nanosecond datetimes will be expressed as int as a workaround
                return datetime.nanos_to_months(v)

        else:
            raise ValueError(f"Cannot apply month transform for type: {source}")

        return lambda v: month_func(v) if v is not None else None

    def can_transform(self, source: IcebergType) -> bool:
        return isinstance(source, (DateType, TimestampType, TimestamptzType, TimestampNanoType, TimestamptzNanoType))

    @property
    def granularity(self) -> TimeResolution:
        return TimeResolution.MONTH

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        return datetime.to_human_month(value) if isinstance(value, int) else "null"

    def __repr__(self) -> str:
        """Return the string representation of the MonthTransform class."""
        return "MonthTransform()"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        pa = _try_import("pyarrow")
        pyiceberg_core_transform = _try_import("pyiceberg_core", extras_name="pyiceberg-core").transform

        return _pyiceberg_transform_wrapper(pyiceberg_core_transform.month, expected_type=pa.int32())


class DayTransform(TimeTransform[S]):
    """Transforms a datetime value into a day value.

    Example:
        >>> transform = DayTransform()
        >>> transform.transform(DateType())(17501)
        17501
    """

    root: LiteralType["day"] = Field(default="day")  # noqa: F821

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[int]]:
        if isinstance(source, DateType):

            def day_func(v: Any) -> int:
                if isinstance(v, py_datetime.date):
                    v = datetime.date_to_days(v)

                return v

        elif isinstance(source, (TimestampType, TimestamptzType)):

            def day_func(v: Any) -> int:
                if isinstance(v, py_datetime.datetime):
                    v = datetime.datetime_to_micros(v)

                return datetime.micros_to_days(v)

        elif isinstance(source, (TimestampNanoType, TimestamptzNanoType)):

            def day_func(v: Any) -> int:
                # python datetime has no nanoseconds support.
                # nanosecond datetimes will be expressed as int as a workaround
                return datetime.nanos_to_days(v)

        else:
            raise ValueError(f"Cannot apply day transform for type: {source}")

        return lambda v: day_func(v) if v is not None else None

    def can_transform(self, source: IcebergType) -> bool:
        return isinstance(source, (DateType, TimestampType, TimestamptzType, TimestampNanoType, TimestamptzNanoType))

    def result_type(self, source: IcebergType) -> IcebergType:
        """Return the result type of a day transform.

        The physical representation conforms to the Iceberg spec as DateType is internally converted to int.
        The DateType returned here provides a more human-readable way to display the partition field.
        """
        return DateType()

    @property
    def granularity(self) -> TimeResolution:
        return TimeResolution.DAY

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        return datetime.to_human_day(value) if isinstance(value, int) else "null"

    def __repr__(self) -> str:
        """Return the string representation of the DayTransform class."""
        return "DayTransform()"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        pa = _try_import("pyarrow", extras_name="pyarrow")
        pyiceberg_core_transform = _try_import("pyiceberg_core", extras_name="pyiceberg-core").transform

        return _pyiceberg_transform_wrapper(pyiceberg_core_transform.day, expected_type=pa.int32())


class HourTransform(TimeTransform[S]):
    """Transforms a datetime value into a hour value.

    Example:
        >>> transform = HourTransform()
        >>> transform.transform(TimestampType())(1512151975038194)
        420042
    """

    root: LiteralType["hour"] = Field(default="hour")  # noqa: F821

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[int]]:
        if isinstance(source, (TimestampType, TimestamptzType)):

            def hour_func(v: Any) -> int:
                if isinstance(v, py_datetime.datetime):
                    v = datetime.datetime_to_micros(v)

                return datetime.micros_to_hours(v)

        elif isinstance(source, (TimestampNanoType, TimestamptzNanoType)):

            def hour_func(v: Any) -> int:
                # python datetime has no nanoseconds support.
                # nanosecond datetimes will be expressed as int as a workaround
                return datetime.nanos_to_hours(v)

        else:
            raise ValueError(f"Cannot apply hour transform for type: {source}")

        return lambda v: hour_func(v) if v is not None else None

    def can_transform(self, source: IcebergType) -> bool:
        return isinstance(source, (TimestampType, TimestamptzType, TimestampNanoType, TimestamptzNanoType))

    @property
    def granularity(self) -> TimeResolution:
        return TimeResolution.HOUR

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        return datetime.to_human_hour(value) if isinstance(value, int) else "null"

    def __repr__(self) -> str:
        """Return the string representation of the HourTransform class."""
        return "HourTransform()"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        pyiceberg_core_transform = _try_import("pyiceberg_core", extras_name="pyiceberg-core").transform

        return _pyiceberg_transform_wrapper(pyiceberg_core_transform.hour)


def _base64encode(buffer: bytes) -> str:
    """Convert bytes to base64 string."""
    return base64.b64encode(buffer).decode("ISO-8859-1")


class IdentityTransform(Transform[S, S]):
    """Transforms a value into itself.

    Example:
        >>> transform = IdentityTransform()
        >>> transform.transform(StringType())('hello-world')
        'hello-world'
    """

    root: LiteralType["identity"] = Field(default="identity")  # noqa: F821

    def __init__(self) -> None:
        super().__init__("identity")

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[S]]:
        return lambda v: v

    def can_transform(self, source: IcebergType) -> bool:
        return source.is_primitive

    def result_type(self, source: IcebergType) -> IcebergType:
        return source

    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)
        elif isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundLiteralPredicate):
            return pred.as_unbound(Reference(name), pred.literal)
        elif isinstance(pred, BoundSetPredicate):
            return pred.as_unbound(Reference(name), pred.literals)
        else:
            return None

    def strict_project(self, name: str, pred: BoundPredicate[Any]) -> Optional[UnboundPredicate[Any]]:
        if isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundLiteralPredicate):
            return pred.as_unbound(Reference(name), pred.literal)
        elif isinstance(pred, BoundSetPredicate):
            return pred.as_unbound(Reference(name), pred.literals)
        else:
            return None

    @property
    def preserves_order(self) -> bool:
        return True

    def satisfies_order_of(self, other: Transform[S, T]) -> bool:
        """Ordering by value is the same as long as the other preserves order."""
        return other.preserves_order

    def to_human_string(self, source_type: IcebergType, value: Optional[S]) -> str:
        return _human_string(value, source_type) if value is not None else "null"

    def __str__(self) -> str:
        """Return the string representation of the IdentityTransform class."""
        return "identity"

    def __repr__(self) -> str:
        """Return the string representation of the IdentityTransform class."""
        return "IdentityTransform()"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        return lambda v: v


class TruncateTransform(Transform[S, S]):
    """A transform for truncating a value to a specified width.

    Args:
      width (int): The truncate width, should be positive.
    Raises:
      ValueError: If a type is provided that is incompatible with a Truncate transform.
    """

    root: str = Field()
    _source_type: IcebergType = PrivateAttr()
    _width: PositiveInt = PrivateAttr()

    def __init__(self, width: int, **data: Any):
        super().__init__(root=f"truncate[{width}]", **data)
        self._width = width

    def can_transform(self, source: IcebergType) -> bool:
        return isinstance(source, (IntegerType, LongType, StringType, BinaryType, DecimalType))

    def result_type(self, source: IcebergType) -> IcebergType:
        return source

    @property
    def preserves_order(self) -> bool:
        return True

    @property
    def source_type(self) -> IcebergType:
        return self._source_type

    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        field_type = pred.term.ref().field.field_type

        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)

        if isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))
        elif isinstance(pred, BoundIn):
            return _set_apply_transform(name, pred, self.transform(field_type))
        elif isinstance(field_type, (IntegerType, LongType, DecimalType)):
            if isinstance(pred, BoundLiteralPredicate):
                return _truncate_number(name, pred, self.transform(field_type))
        elif isinstance(field_type, (BinaryType, StringType)):
            if isinstance(pred, BoundLiteralPredicate):
                return _truncate_array(name, pred, self.transform(field_type))
        return None

    def strict_project(self, name: str, pred: BoundPredicate[Any]) -> Optional[UnboundPredicate[Any]]:
        field_type = pred.term.ref().field.field_type

        if isinstance(pred.term, BoundTransform):
            return _project_transform_predicate(self, name, pred)

        if isinstance(pred, BoundUnaryPredicate):
            return pred.as_unbound(Reference(name))

        if isinstance(field_type, (IntegerType, LongType, DecimalType)):
            if isinstance(pred, BoundLiteralPredicate):
                return _truncate_number_strict(name, pred, self.transform(field_type))
            elif isinstance(pred, BoundNotIn):
                return _set_apply_transform(name, pred, self.transform(field_type))
            else:
                return None

        if isinstance(pred, BoundLiteralPredicate):
            if isinstance(pred, BoundStartsWith):
                literal_width = len(pred.literal.value)
                if literal_width < self.width:
                    return pred.as_unbound(name, pred.literal.value)
                elif literal_width == self.width:
                    return EqualTo(name, pred.literal.value)
                else:
                    return None
            elif isinstance(pred, BoundNotStartsWith):
                literal_width = len(pred.literal.value)
                if literal_width < self.width:
                    return pred.as_unbound(name, pred.literal.value)
                elif literal_width == self.width:
                    return NotEqualTo(name, pred.literal.value)
                else:
                    return pred.as_unbound(name, self.transform(field_type)(pred.literal.value))
            else:
                # ProjectionUtil.truncateArrayStrict(name, pred, this);
                return _truncate_array_strict(name, pred, self.transform(field_type))
        elif isinstance(pred, BoundNotIn):
            return _set_apply_transform(name, pred, self.transform(field_type))
        else:
            return None

    @property
    def width(self) -> int:
        return self._width

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[S]]:
        if isinstance(source, (IntegerType, LongType)):

            def truncate_func(v: Any) -> Any:
                return v - v % self._width

        elif isinstance(source, (StringType, BinaryType)):

            def truncate_func(v: Any) -> Any:
                return v[0 : min(self._width, len(v))]

        elif isinstance(source, DecimalType):

            def truncate_func(v: Any) -> Any:
                return truncate_decimal(v, self._width)

        else:
            raise ValueError(f"Cannot truncate for type: {source}")

        return lambda v: truncate_func(v) if v is not None else None

    def satisfies_order_of(self, other: Transform[S, T]) -> bool:
        if self == other:
            return True
        elif (
            isinstance(self.source_type, StringType)
            and isinstance(other, TruncateTransform)
            and isinstance(other.source_type, StringType)
        ):
            return self.width >= other.width

        return False

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        if value is None:
            return "null"
        elif isinstance(value, bytes):
            return _base64encode(value)
        else:
            return str(value)

    def __repr__(self) -> str:
        """Return the string representation of the TruncateTransform class."""
        return f"TruncateTransform(width={self._width})"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        pyiceberg_core_transform = _try_import("pyiceberg_core", extras_name="pyiceberg-core").transform

        return _pyiceberg_transform_wrapper(pyiceberg_core_transform.truncate, self._width)


@singledispatch
def _human_string(value: Any, _type: IcebergType) -> str:
    return str(value)


@_human_string.register(bytes)
def _(value: bytes, _type: IcebergType) -> str:
    return _base64encode(value)


@_human_string.register(int)
def _(value: int, _type: IcebergType) -> str:
    return _int_to_human_string(_type, value)


@_human_string.register(bool)
def _(value: bool, _type: IcebergType) -> str:
    return str(value).lower()


@singledispatch
def _int_to_human_string(_type: IcebergType, value: int) -> str:
    return str(value)


@_int_to_human_string.register(DateType)
def _(_type: IcebergType, value: int) -> str:
    return datetime.to_human_day(value)


@_int_to_human_string.register(TimeType)
def _(_type: IcebergType, value: int) -> str:
    return datetime.to_human_time(value)


@_int_to_human_string.register(TimestampType)
def _(_type: IcebergType, value: int) -> str:
    return datetime.to_human_timestamp(value)


@_int_to_human_string.register(TimestamptzType)
def _(_type: IcebergType, value: int) -> str:
    return datetime.to_human_timestamptz(value)


class UnknownTransform(Transform[S, T]):
    """A transform that represents when an unknown transform is provided.

    Args:
      transform (str): A string name of a transform.

    Keyword Args:
      source_type (IcebergType): An Iceberg `Type`.
    """

    root: LiteralType["unknown"] = Field(default="unknown")  # noqa: F821
    _transform: str = PrivateAttr()

    def __init__(self, transform: str, **data: Any):
        super().__init__(**data)
        self._transform = transform

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[T]]:
        raise AttributeError(f"Cannot apply unsupported transform: {self}")

    def can_transform(self, source: IcebergType) -> bool:
        return False

    def result_type(self, source: IcebergType) -> StringType:
        return StringType()

    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        return None

    def strict_project(self, name: str, pred: BoundPredicate[Any]) -> Optional[UnboundPredicate[Any]]:
        return None

    def __repr__(self) -> str:
        """Return the string representation of the UnknownTransform class."""
        return f"UnknownTransform(transform={repr(self._transform)})"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        raise NotImplementedError()


class VoidTransform(Transform[S, None], Singleton):
    """A transform that always returns None."""

    root: str = "void"

    def transform(self, source: IcebergType) -> Callable[[Optional[S]], Optional[T]]:
        return lambda v: None

    def can_transform(self, _: IcebergType) -> bool:
        return True

    def result_type(self, source: IcebergType) -> IcebergType:
        return source

    def project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        return None

    def strict_project(self, name: str, pred: BoundPredicate[L]) -> Optional[UnboundPredicate[Any]]:
        return None

    def to_human_string(self, _: IcebergType, value: Optional[S]) -> str:
        return "null"

    def __repr__(self) -> str:
        """Return the string representation of the VoidTransform class."""
        return "VoidTransform()"

    def pyarrow_transform(self, source: IcebergType) -> "Callable[[pa.Array], pa.Array]":
        try:
            import pyarrow as pa
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("For partition transforms, PyArrow needs to be installed") from e

        return lambda arr: pa.nulls(len(arr), type=arr.type)


def _truncate_number(
    name: str, pred: BoundLiteralPredicate[L], transform: Callable[[Optional[L]], Optional[L]]
) -> Optional[UnboundPredicate[Any]]:
    boundary = pred.literal

    if not isinstance(boundary, (LongLiteral, DecimalLiteral, DateLiteral, TimestampLiteral)):
        raise ValueError(f"Expected a numeric literal, got: {type(boundary)}")

    if isinstance(pred, BoundLessThan):
        return LessThanOrEqual(Reference(name), _transform_literal(transform, boundary.decrement()))
    elif isinstance(pred, BoundLessThanOrEqual):
        return LessThanOrEqual(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundGreaterThan):
        return GreaterThanOrEqual(Reference(name), _transform_literal(transform, boundary.increment()))
    elif isinstance(pred, BoundGreaterThanOrEqual):
        return GreaterThanOrEqual(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundEqualTo):
        return EqualTo(Reference(name), _transform_literal(transform, boundary))
    else:
        return None


def _truncate_number_strict(
    name: str, pred: BoundLiteralPredicate[L], transform: Callable[[Optional[L]], Optional[L]]
) -> Optional[UnboundPredicate[Any]]:
    boundary = pred.literal

    if not isinstance(boundary, (LongLiteral, DecimalLiteral, DateLiteral, TimestampLiteral)):
        raise ValueError(f"Expected a numeric literal, got: {type(boundary)}")

    if isinstance(pred, BoundLessThan):
        return LessThan(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundLessThanOrEqual):
        return LessThan(Reference(name), _transform_literal(transform, boundary.increment()))
    elif isinstance(pred, BoundGreaterThan):
        return GreaterThan(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundGreaterThanOrEqual):
        return GreaterThan(Reference(name), _transform_literal(transform, boundary.decrement()))
    elif isinstance(pred, BoundNotEqualTo):
        return NotEqualTo(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundEqualTo):
        # there is no predicate that guarantees equality because adjacent longs transform to the
        # same value
        return None
    else:
        return None


def _truncate_array_strict(
    name: str, pred: BoundLiteralPredicate[L], transform: Callable[[Optional[L]], Optional[L]]
) -> Optional[UnboundPredicate[Any]]:
    boundary = pred.literal

    if isinstance(pred, (BoundLessThan, BoundLessThanOrEqual)):
        return LessThan(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, (BoundGreaterThan, BoundGreaterThanOrEqual)):
        return GreaterThan(Reference(name), _transform_literal(transform, boundary))
    if isinstance(pred, BoundNotEqualTo):
        return NotEqualTo(Reference(name), _transform_literal(transform, boundary))
    else:
        return None


def _truncate_array(
    name: str, pred: BoundLiteralPredicate[L], transform: Callable[[Optional[L]], Optional[L]]
) -> Optional[UnboundPredicate[Any]]:
    boundary = pred.literal

    if isinstance(pred, (BoundLessThan, BoundLessThanOrEqual)):
        return LessThanOrEqual(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, (BoundGreaterThan, BoundGreaterThanOrEqual)):
        return GreaterThanOrEqual(Reference(name), _transform_literal(transform, boundary))
    if isinstance(pred, BoundEqualTo):
        return EqualTo(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundStartsWith):
        return StartsWith(Reference(name), _transform_literal(transform, boundary))
    elif isinstance(pred, BoundNotStartsWith):
        return NotStartsWith(Reference(name), _transform_literal(transform, boundary))
    else:
        return None


def _project_transform_predicate(
    transform: Transform[Any, Any], partition_name: str, pred: BoundPredicate[L]
) -> Optional[UnboundPredicate[Any]]:
    term = pred.term
    if isinstance(term, BoundTransform) and transform == term.transform:
        return _remove_transform(partition_name, pred)
    return None


def _remove_transform(partition_name: str, pred: BoundPredicate[L]) -> UnboundPredicate[Any]:
    if isinstance(pred, BoundUnaryPredicate):
        return pred.as_unbound(Reference(partition_name))
    elif isinstance(pred, BoundLiteralPredicate):
        return pred.as_unbound(Reference(partition_name), pred.literal)
    elif isinstance(pred, (BoundIn, BoundNotIn)):
        return pred.as_unbound(Reference(partition_name), pred.literals)
    else:
        raise ValueError(f"Cannot replace transform in unknown predicate: {pred}")


def _set_apply_transform(name: str, pred: BoundSetPredicate[L], transform: Callable[[L], L]) -> UnboundPredicate[Any]:
    literals = pred.literals
    if isinstance(pred, BoundSetPredicate):
        transformed_literals = {_transform_literal(transform, literal) for literal in literals}
        return pred.as_unbound(Reference(name=name), literals=transformed_literals)
    else:
        raise ValueError(f"Unknown BoundSetPredicate: {pred}")


class BoundTransform(BoundTerm[L]):
    """A transform expression."""

    transform: Transform[L, Any]

    def __init__(self, term: BoundTerm[L], transform: Transform[L, Any]):
        self.term: BoundTerm[L] = term
        self.transform = transform
