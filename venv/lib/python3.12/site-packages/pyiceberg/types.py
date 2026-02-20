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
"""Data types used in describing Iceberg schemas.

This module implements the data types described in the Iceberg specification for Iceberg schemas. To
describe an Iceberg table schema, these classes can be used in the construction of a StructType instance.

Example:
    >>> str(StructType(
    ...     NestedField(1, "required_field", StringType(), True),
    ...     NestedField(2, "optional_field", IntegerType())
    ... ))
    'struct<1: required_field: required string, 2: optional_field: optional int>'

Notes:
  - https://iceberg.apache.org/spec/#primitive-types
"""

from __future__ import annotations

import re
from functools import cached_property
from typing import (
    Annotated,
    Any,
    ClassVar,
    Literal,
)

from pydantic import (
    BeforeValidator,
    Field,
    PrivateAttr,
    SerializeAsAny,
    field_validator,
    model_serializer,
    model_validator,
)
from pydantic_core.core_schema import ValidationInfo, ValidatorFunctionWrapHandler

from pyiceberg.exceptions import ValidationError
from pyiceberg.typedef import IcebergBaseModel, IcebergRootModel, L, TableVersion
from pyiceberg.utils.parsing import ParseNumberFromBrackets
from pyiceberg.utils.singleton import Singleton

DECIMAL_REGEX = re.compile(r"decimal\((\d+),\s*(\d+)\)")
FIXED = "fixed"
FIXED_PARSER = ParseNumberFromBrackets(FIXED)


def transform_dict_value_to_str(dict: dict[str, Any]) -> dict[str, str]:
    """Transform all values in the dictionary to string. Raise an error if any value is None."""
    for key, value in dict.items():
        if value is None:
            raise ValueError(f"None type is not a supported value in properties: {key}")
    return {k: str(v).lower() if isinstance(v, bool) else str(v) for k, v in dict.items()}


def _parse_decimal_type(decimal: Any) -> tuple[int, int]:
    if isinstance(decimal, str):
        matches = DECIMAL_REGEX.search(decimal)
        if matches:
            return int(matches.group(1)), int(matches.group(2))
        else:
            raise ValidationError(f"Could not parse {decimal} into a DecimalType")
    elif isinstance(decimal, dict):
        return decimal["precision"], decimal["scale"]
    else:
        return decimal


def _parse_fixed_type(fixed: Any) -> int:
    if isinstance(fixed, str):
        return FIXED_PARSER.match(fixed)
    elif isinstance(fixed, dict):
        return fixed["length"]
    else:
        return fixed


def strtobool(val: str) -> bool:
    """Convert a string representation of truth to true (1) or false (0).

    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"Invalid truth value: {val!r}")


class IcebergType(IcebergBaseModel):
    """Base type for all Iceberg Types.

    Example:
        >>> str(IcebergType())
        'IcebergType()'
        >>> repr(IcebergType())
        'IcebergType()'
    """

    @model_validator(mode="wrap")
    @classmethod
    def handle_primitive_type(cls, v: Any, handler: ValidatorFunctionWrapHandler) -> IcebergType:
        # Pydantic works mostly around dicts, and there seems to be something
        # by not serializing into a RootModel, might revisit this.
        if isinstance(v, str):
            if v == "boolean":
                return BooleanType()
            elif v == "string":
                return StringType()
            elif v == "int":
                return IntegerType()
            elif v == "long":
                return LongType()
            if v == "float":
                return FloatType()
            if v == "double":
                return DoubleType()
            if v == "timestamp":
                return TimestampType()
            if v == "timestamptz":
                return TimestamptzType()
            if v == "timestamp_ns":
                return TimestampNanoType()
            if v == "timestamptz_ns":
                return TimestamptzNanoType()
            if v == "date":
                return DateType()
            if v == "time":
                return TimeType()
            if v == "uuid":
                return UUIDType()
            if v == "binary":
                return BinaryType()
            if v == "unknown":
                return UnknownType()
            if v.startswith("fixed"):
                return FixedType(_parse_fixed_type(v))
            if v.startswith("decimal"):
                precision, scale = _parse_decimal_type(v)
                return DecimalType(precision, scale)
            else:
                raise ValueError(f"Type not recognized: {v}")
        if isinstance(v, dict) and cls == IcebergType:
            complex_type = v.get("type")
            if complex_type == "list":
                return ListType(**v)
            elif complex_type == "map":
                return MapType(**v)
            elif complex_type == "struct":
                return StructType(**v)
            else:
                return NestedField(**v)
        return handler(v)

    @property
    def is_primitive(self) -> bool:
        return isinstance(self, PrimitiveType)

    @property
    def is_struct(self) -> bool:
        return isinstance(self, StructType)

    def minimum_format_version(self) -> TableVersion:
        """Minimum Iceberg format version after which this type is supported."""
        return 1


class PrimitiveType(Singleton, IcebergRootModel[str], IcebergType):
    """Base class for all Iceberg Primitive Types."""

    root: Any = Field()

    def __repr__(self) -> str:
        """Return the string representation of the PrimitiveType class."""
        return f"{type(self).__name__}()"

    def __str__(self) -> str:
        """Return the string representation of the PrimitiveType class."""
        return self.root


class FixedType(PrimitiveType):
    """A fixed data type in Iceberg.

    Example:
        >>> FixedType(8)
        FixedType(length=8)
        >>> FixedType(8) == FixedType(8)
        True
        >>> FixedType(19) == FixedType(25)
        False
    """

    root: int = Field()

    def __init__(self, length: int) -> None:
        super().__init__(root=length)

    @model_serializer
    def ser_model(self) -> str:
        return f"fixed[{self.root}]"

    def __len__(self) -> int:
        """Return the length of an instance of the FixedType class."""
        return self.root

    def __str__(self) -> str:
        """Return the string representation."""
        return f"fixed[{self.root}]"

    def __repr__(self) -> str:
        """Return the string representation of the FixedType class."""
        return f"FixedType(length={self.root})"

    def __getnewargs__(self) -> tuple[int]:
        """Pickle the FixedType class."""
        return (self.root,)


class DecimalType(PrimitiveType):
    """A decimal data type in Iceberg.

    Example:
        >>> DecimalType(32, 3)
        DecimalType(precision=32, scale=3)
        >>> DecimalType(8, 3) == DecimalType(8, 3)
        True
    """

    root: tuple[int, int]

    def __init__(self, precision: int, scale: int) -> None:
        super().__init__(root=(precision, scale))

    @model_serializer
    def ser_model(self) -> str:
        """Serialize the model to a string."""
        return f"decimal({self.precision}, {self.scale})"

    @property
    def precision(self) -> int:
        """Return the precision of the decimal."""
        return self.root[0]

    @property
    def scale(self) -> int:
        """Return the scale of the decimal."""
        return self.root[1]

    def __repr__(self) -> str:
        """Return the string representation of the DecimalType class."""
        return f"DecimalType(precision={self.precision}, scale={self.scale})"

    def __str__(self) -> str:
        """Return the string representation."""
        return f"decimal({self.precision}, {self.scale})"

    def __hash__(self) -> int:
        """Return the hash of the tuple."""
        return hash(self.root)

    def __getnewargs__(self) -> tuple[int, int]:
        """Pickle the DecimalType class."""
        return self.precision, self.scale

    def __eq__(self, other: Any) -> bool:
        """Compare to root to another object."""
        return self.root == other.root if isinstance(other, DecimalType) else False


def _deserialize_default_value(v: Any, context: ValidationInfo) -> Any:
    if v is not None:
        from pyiceberg.conversions import from_json

        return from_json(context.data.get("field_type"), v)
    else:
        return None


DefaultValue = Annotated[L, BeforeValidator(_deserialize_default_value)]


class NestedField(IcebergType):
    """Represents a field of a struct, a map key, a map value, or a list element.

    This is where field IDs, names, docs, and nullability are tracked.

    Example:
        >>> str(NestedField(
        ...     field_id=1,
        ...     name='foo',
        ...     field_type=FixedType(22),
        ...     required=False,
        ... ))
        '1: foo: optional fixed[22]'
        >>> str(NestedField(
        ...     field_id=2,
        ...     name='bar',
        ...     field_type=LongType(),
        ...     is_optional=False,
        ...     doc="Just a long"
        ... ))
        '2: bar: required long (Just a long)'
        >>> str(NestedField(
        ...     field_id=3,
        ...     name='baz',
        ...     field_type="string",
        ...     required=True,
        ...     doc="A string field"
        ... ))
        '3: baz: required string (A string field)'
    """

    field_id: int = Field(alias="id")
    name: str = Field()
    field_type: SerializeAsAny[IcebergType] = Field(alias="type")
    required: bool = Field(default=False)
    doc: str | None = Field(default=None, repr=False)
    initial_default: DefaultValue | None = Field(alias="initial-default", default=None, repr=True)  # type: ignore
    write_default: DefaultValue | None = Field(alias="write-default", default=None, repr=True)  # type: ignore

    @field_validator("field_type", mode="before")
    def convert_field_type(cls, v: Any) -> IcebergType:
        """Convert string values into IcebergType instances."""
        if isinstance(v, str):
            try:
                return IcebergType.handle_primitive_type(v, None)
            except ValueError as e:
                raise ValueError(f"Unsupported field type: '{v}'") from e
        return v

    def __init__(
        self,
        field_id: int | None = None,
        name: str | None = None,
        field_type: IcebergType | str | None = None,
        required: bool = False,
        doc: str | None = None,
        initial_default: Any | None = None,
        write_default: L | None = None,
        **data: Any,
    ):
        # We need an init when we want to use positional arguments, but
        # need also to support the aliases.
        data["id"] = data["id"] if "id" in data else field_id
        data["name"] = name
        data["type"] = data["type"] if "type" in data else field_type
        data["required"] = required
        data["doc"] = doc
        data["initial-default"] = data["initial-default"] if "initial-default" in data else initial_default
        data["write-default"] = data["write-default"] if "write-default" in data else write_default
        super().__init__(**data)

    @model_serializer()
    def serialize_model(self) -> dict[str, Any]:
        from pyiceberg.conversions import to_json

        fields = {
            "id": self.field_id,
            "name": self.name,
            "type": self.field_type,
            "required": self.required,
        }

        if self.doc is not None:
            fields["doc"] = self.doc
        if self.initial_default is not None:
            fields["initial-default"] = to_json(self.field_type, self.initial_default)
        if self.write_default is not None:
            fields["write-default"] = to_json(self.field_type, self.write_default)

        return fields

    def __str__(self) -> str:
        """Return the string representation of the NestedField class."""
        doc = "" if not self.doc else f" ({self.doc})"
        req = "required" if self.required else "optional"
        return f"{self.field_id}: {self.name}: {req} {self.field_type}{doc}"

    def __repr__(self) -> str:
        """Return the string representation of the NestedField class."""
        parts = [
            f"field_id={self.field_id}",
            f"name={self.name!r}",
            f"field_type={self.field_type!r}",
            f"required={self.required}",
        ]
        if self.initial_default is not None:
            parts.append(f"initial_default={self.initial_default!r}")
        if self.write_default is not None:
            parts.append(f"write_default={self.write_default!r}")

        return f"NestedField({', '.join(parts)})"

    def __getnewargs__(self) -> tuple[int, str, IcebergType, bool, str | None]:
        """Pickle the NestedField class."""
        return (self.field_id, self.name, self.field_type, self.required, self.doc)

    @property
    def optional(self) -> bool:
        return not self.required


class StructType(IcebergType):
    """A struct type in Iceberg.

    Example:
        >>> str(StructType(
        ...     NestedField(1, "required_field", StringType(), True),
        ...     NestedField(2, "optional_field", IntegerType())
        ... ))
        'struct<1: required_field: optional string, 2: optional_field: optional int>'
    """

    type: Literal["struct"] = Field(default="struct")
    fields: tuple[NestedField, ...] = Field(default_factory=tuple)
    _hash: int = PrivateAttr()

    def __init__(self, *fields: NestedField, **data: Any):
        # In case we use positional arguments, instead of keyword args
        if fields:
            data["fields"] = fields
        super().__init__(**data)
        self._hash = hash(self.fields)

    def field(self, field_id: int) -> NestedField | None:
        for field in self.fields:
            if field.field_id == field_id:
                return field
        return None

    def field_by_name(self, name: str, case_sensitive: bool = True) -> NestedField | None:
        if case_sensitive:
            for field in self.fields:
                if field.name == name:
                    return field
        else:
            name_lower = name.lower()
            for field in self.fields:
                if field.name.lower() == name_lower:
                    return field
        return None

    def __str__(self) -> str:
        """Return the string representation of the StructType class."""
        return f"struct<{', '.join(map(str, self.fields))}>"

    def __repr__(self) -> str:
        """Return the string representation of the StructType class."""
        return f"StructType(fields=({', '.join(map(repr, self.fields))},))"

    def __len__(self) -> int:
        """Return the length of an instance of the StructType class."""
        return len(self.fields)

    def __getnewargs__(self) -> tuple[NestedField, ...]:
        """Pickle the StructType class."""
        return self.fields

    def __hash__(self) -> int:
        """Use the cache hash value of the StructType class."""
        return self._hash

    def __eq__(self, other: Any) -> bool:
        """Compare the object if it is equal to another object."""
        return self.fields == other.fields if isinstance(other, StructType) else False


class ListType(IcebergType):
    """A list type in Iceberg.

    Example:
        >>> ListType(element_id=3, element_type=StringType(), element_required=True)
        ListType(element_id=3, element_type=StringType(), element_required=True)
    """

    type: Literal["list"] = Field(default="list")
    element_id: int = Field(alias="element-id")
    element_type: SerializeAsAny[IcebergType] = Field(alias="element")
    element_required: bool = Field(alias="element-required", default=True)
    _element_field: NestedField = PrivateAttr()
    _hash: int = PrivateAttr()

    def __init__(
        self, element_id: int | None = None, element: IcebergType | None = None, element_required: bool = True, **data: Any
    ):
        data["element-id"] = data["element-id"] if "element-id" in data else element_id
        data["element"] = element or data["element_type"]
        data["element-required"] = data["element-required"] if "element-required" in data else element_required
        super().__init__(**data)
        self._hash = hash(data.values())

    @cached_property
    def element_field(self) -> NestedField:
        return NestedField(
            name="element",
            field_id=self.element_id,
            field_type=self.element_type,
            required=self.element_required,
        )

    def __str__(self) -> str:
        """Return the string representation of the ListType class."""
        return f"list<{self.element_type}>"

    def __getnewargs__(self) -> tuple[int, IcebergType, bool]:
        """Pickle the ListType class."""
        return (self.element_id, self.element_type, self.element_required)

    def __hash__(self) -> int:
        """Use the cache hash value of the StructType class."""
        return self._hash

    def __eq__(self, other: Any) -> bool:
        """Compare the list type to another list type."""
        return self.element_field == other.element_field if isinstance(other, ListType) else False


class MapType(IcebergType):
    """A map type in Iceberg.

    Example:
        >>> MapType(key_id=1, key_type=StringType(), value_id=2, value_type=IntegerType(), value_required=True)
        MapType(key_id=1, key_type=StringType(), value_id=2, value_type=IntegerType(), value_required=True)
    """

    type: Literal["map"] = Field(default="map")
    key_id: int = Field(alias="key-id")
    key_type: SerializeAsAny[IcebergType] = Field(alias="key")
    value_id: int = Field(alias="value-id")
    value_type: SerializeAsAny[IcebergType] = Field(alias="value")
    value_required: bool = Field(alias="value-required", default=True)
    _hash: int = PrivateAttr()

    def __init__(
        self,
        key_id: int | None = None,
        key_type: IcebergType | None = None,
        value_id: int | None = None,
        value_type: IcebergType | None = None,
        value_required: bool = True,
        **data: Any,
    ):
        data["key-id"] = data["key-id"] if "key-id" in data else key_id
        data["key"] = data["key"] if "key" in data else key_type
        data["value-id"] = data["value-id"] if "value-id" in data else value_id
        data["value"] = data["value"] if "value" in data else value_type
        data["value-required"] = data["value-required"] if "value-required" in data else value_required
        super().__init__(**data)
        self._hash = hash(self.__getnewargs__())

    @cached_property
    def key_field(self) -> NestedField:
        return NestedField(
            name="key",
            field_id=self.key_id,
            field_type=self.key_type,
            required=True,
        )

    @cached_property
    def value_field(self) -> NestedField:
        return NestedField(
            name="value",
            field_id=self.value_id,
            field_type=self.value_type,
            required=self.value_required,
        )

    def __str__(self) -> str:
        """Return the string representation of the MapType class."""
        return f"map<{self.key_type}, {self.value_type}>"

    def __getnewargs__(self) -> tuple[int, IcebergType, int, IcebergType, bool]:
        """Pickle the MapType class."""
        return (self.key_id, self.key_type, self.value_id, self.value_type, self.value_required)

    def __hash__(self) -> int:
        """Return the hash of the MapType."""
        return self._hash

    def __eq__(self, other: Any) -> bool:
        """Compare the MapType to another object."""
        return (
            self.key_field == other.key_field and self.value_field == other.value_field if isinstance(other, MapType) else False
        )


class BooleanType(PrimitiveType):
    """A boolean data type in Iceberg can be represented using an instance of this class.

    Example:
        >>> column_foo = BooleanType()
        >>> isinstance(column_foo, BooleanType)
        True
        >>> column_foo
        BooleanType()
    """

    root: Literal["boolean"] = Field(default="boolean")


class IntegerType(PrimitiveType):
    """An Integer data type in Iceberg can be represented using an instance of this class.

    Integers in Iceberg are 32-bit signed and can be promoted to Longs.

    Example:
        >>> column_foo = IntegerType()
        >>> isinstance(column_foo, IntegerType)
        True

    Attributes:
        max (int): The maximum allowed value for Integers, inherited from the canonical Iceberg implementation
            in Java (returns `2147483647`)
        min (int): The minimum allowed value for Integers, inherited from the canonical Iceberg implementation
            in Java (returns `-2147483648`)
    """

    root: Literal["int"] = Field(default="int")

    max: ClassVar[int] = 2147483647
    min: ClassVar[int] = -2147483648


class LongType(PrimitiveType):
    """A Long data type in Iceberg can be represented using an instance of this class.

    Longs in Iceberg are 64-bit signed integers.

    Example:
        >>> column_foo = LongType()
        >>> isinstance(column_foo, LongType)
        True
        >>> column_foo
        LongType()
        >>> str(column_foo)
        'long'

    Attributes:
        max (int): The maximum allowed value for Longs, inherited from the canonical Iceberg implementation
            in Java. (returns `9223372036854775807`)
        min (int): The minimum allowed value for Longs, inherited from the canonical Iceberg implementation
            in Java (returns `-9223372036854775808`)
    """

    root: Literal["long"] = Field(default="long")

    max: ClassVar[int] = 9223372036854775807
    min: ClassVar[int] = -9223372036854775808


class FloatType(PrimitiveType):
    """A Float data type in Iceberg can be represented using an instance of this class.

    Floats in Iceberg are 32-bit IEEE 754 floating points and can be promoted to Doubles.

    Example:
        >>> column_foo = FloatType()
        >>> isinstance(column_foo, FloatType)
        True
        >>> column_foo
        FloatType()

    Attributes:
        max (float): The maximum allowed value for Floats, inherited from the canonical Iceberg implementation
            in Java. (returns `3.4028235e38`)
        min (float): The minimum allowed value for Floats, inherited from the canonical Iceberg implementation
            in Java (returns `-3.4028235e38`)
    """

    max: ClassVar[float] = 3.4028235e38
    min: ClassVar[float] = -3.4028235e38

    root: Literal["float"] = Field(default="float")


class DoubleType(PrimitiveType):
    """A Double data type in Iceberg can be represented using an instance of this class.

    Doubles in Iceberg are 64-bit IEEE 754 floating points.

    Example:
        >>> column_foo = DoubleType()
        >>> isinstance(column_foo, DoubleType)
        True
        >>> column_foo
        DoubleType()
    """

    root: Literal["double"] = Field(default="double")


class DateType(PrimitiveType):
    """A Date data type in Iceberg can be represented using an instance of this class.

    Dates in Iceberg are calendar dates without a timezone or time.

    Example:
        >>> column_foo = DateType()
        >>> isinstance(column_foo, DateType)
        True
        >>> column_foo
        DateType()
    """

    root: Literal["date"] = Field(default="date")


class TimeType(PrimitiveType):
    """A Time data type in Iceberg can be represented using an instance of this class.

    Times in Iceberg have microsecond precision and are a time of day without a date or timezone.

    Example:
        >>> column_foo = TimeType()
        >>> isinstance(column_foo, TimeType)
        True
        >>> column_foo
        TimeType()
    """

    root: Literal["time"] = Field(default="time")


class TimestampType(PrimitiveType):
    """A Timestamp data type in Iceberg can be represented using an instance of this class.

    Timestamps in Iceberg have microsecond precision and include a date and a time of day without a timezone.

    Example:
        >>> column_foo = TimestampType()
        >>> isinstance(column_foo, TimestampType)
        True
        >>> column_foo
        TimestampType()
    """

    root: Literal["timestamp"] = Field(default="timestamp")


class TimestamptzType(PrimitiveType):
    """A Timestamptz data type in Iceberg can be represented using an instance of this class.

    Timestamptzs in Iceberg are stored as UTC and include a date and a time of day with a timezone.

    Example:
        >>> column_foo = TimestamptzType()
        >>> isinstance(column_foo, TimestamptzType)
        True
        >>> column_foo
        TimestamptzType()
    """

    root: Literal["timestamptz"] = Field(default="timestamptz")


class TimestampNanoType(PrimitiveType):
    """A TimestampNano data type in Iceberg can be represented using an instance of this class.

    TimestampNanos in Iceberg have nanosecond precision and include a date and a time of day without a timezone.

    Example:
        >>> column_foo = TimestampNanoType()
        >>> isinstance(column_foo, TimestampNanoType)
        True
        >>> column_foo
        TimestampNanoType()
    """

    root: Literal["timestamp_ns"] = Field(default="timestamp_ns")

    def minimum_format_version(self) -> TableVersion:
        return 3


class TimestamptzNanoType(PrimitiveType):
    """A TimestamptzNano data type in Iceberg can be represented using an instance of this class.

    TimestamptzNanos in Iceberg are stored as UTC and include a date and a time of day with a timezone.

    Example:
        >>> column_foo = TimestamptzNanoType()
        >>> isinstance(column_foo, TimestamptzNanoType)
        True
        >>> column_foo
        TimestamptzNanoType()
    """

    root: Literal["timestamptz_ns"] = Field(default="timestamptz_ns")

    def minimum_format_version(self) -> TableVersion:
        return 3


class StringType(PrimitiveType):
    """A String data type in Iceberg can be represented using an instance of this class.

    Strings in Iceberg are arbitrary-length character sequences and are encoded with UTF-8.

    Example:
        >>> column_foo = StringType()
        >>> isinstance(column_foo, StringType)
        True
        >>> column_foo
        StringType()
    """

    root: Literal["string"] = Field(default="string")


class UUIDType(PrimitiveType):
    """A UUID data type in Iceberg can be represented using an instance of this class.

    UUIDs in Iceberg are universally unique identifiers.

    Example:
        >>> column_foo = UUIDType()
        >>> isinstance(column_foo, UUIDType)
        True
        >>> column_foo
        UUIDType()
    """

    root: Literal["uuid"] = Field(default="uuid")


class BinaryType(PrimitiveType):
    """A Binary data type in Iceberg can be represented using an instance of this class.

    Binaries in Iceberg are arbitrary-length byte arrays.

    Example:
        >>> column_foo = BinaryType()
        >>> isinstance(column_foo, BinaryType)
        True
        >>> column_foo
        BinaryType()
    """

    root: Literal["binary"] = Field(default="binary")


class UnknownType(PrimitiveType):
    """An unknown data type in Iceberg can be represented using an instance of this class.

    Unknowns in Iceberg are used to represent data types that are not known at the time of writing.

    Example:
        >>> column_foo = UnknownType()
        >>> isinstance(column_foo, UnknownType)
        True
        >>> column_foo
        UnknownType()
    """

    root: Literal["unknown"] = Field(default="unknown")

    def minimum_format_version(self) -> TableVersion:
        return 3
