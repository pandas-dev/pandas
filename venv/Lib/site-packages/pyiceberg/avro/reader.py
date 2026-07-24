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
"""
Classes for building the Reader tree.

Constructing a reader tree from the schema makes it easy
to decouple the reader implementation from the schema.

The reader tree can be changed in such a way that the
read schema is different, while respecting the read schema.
"""

from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from dataclasses import field as dataclassfield
from decimal import Decimal
from typing import (
    Any,
)
from uuid import UUID

from pyiceberg.avro.decoder import BinaryDecoder
from pyiceberg.typedef import StructProtocol
from pyiceberg.types import StructType
from pyiceberg.utils.decimal import bytes_to_decimal, decimal_required_bytes
from pyiceberg.utils.lazydict import LazyDict
from pyiceberg.utils.singleton import Singleton


def _skip_map_array(decoder: BinaryDecoder, skip_entry: Callable[[], None]) -> None:
    """Skips over an array or map.

    Both the array and map are encoded similar, and we can reuse
    the logic of skipping in an efficient way.

    From the Avro spec:

    Maps (and arrays) are encoded as a series of blocks.
    Each block consists of a long count value, followed by that many key/value pairs in the case of a map,
    and followed by that many array items in the case of an array. A block with count zero indicates the
    end of the map. Each item is encoded per the map's value schema.

    If a block's count is negative, its absolute value is used, and the count is followed immediately by a
    long block size indicating the number of bytes in the block. This block size permits fast skipping
    through data, e.g., when projecting a record to a subset of its fields.

    Args:
        decoder:
            The decoder that reads the types from the underlying data.
        skip_entry:
            Function to skip over the underlying data, element in case of an array, and the
            key/value in the case of a map.
    """
    block_count = decoder.read_int()
    while block_count != 0:
        if block_count < 0:
            # The length in bytes in encoded, so we can skip over it right away
            block_size = decoder.read_int()
            decoder.skip(block_size)
        else:
            for _ in range(block_count):
                skip_entry()
        block_count = decoder.read_int()


class Reader(Singleton):
    @abstractmethod
    def read(self, decoder: BinaryDecoder) -> Any: ...

    @abstractmethod
    def skip(self, decoder: BinaryDecoder) -> None: ...

    def __repr__(self) -> str:
        """Return the string representation of the Reader class."""
        return f"{self.__class__.__name__}()"


class NoneReader(Reader):
    def read(self, _: BinaryDecoder) -> None:
        return None

    def skip(self, decoder: BinaryDecoder) -> None:
        return None


class DefaultReader(Reader):
    __slots__ = ("default_value",)
    default_value: Any

    def __init__(self, default_value: Any) -> None:
        self.default_value = default_value

    def read(self, _: BinaryDecoder) -> Any:
        return self.default_value

    def skip(self, decoder: BinaryDecoder) -> None:
        pass


class BooleanReader(Reader):
    def read(self, decoder: BinaryDecoder) -> bool:
        return decoder.read_boolean()

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_boolean()


class IntegerReader(Reader):
    """Longs and ints are encoded the same way, and there is no long in Python."""

    def read(self, decoder: BinaryDecoder) -> int:
        return decoder.read_int()

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_int()


class FloatReader(Reader):
    def read(self, decoder: BinaryDecoder) -> float:
        return decoder.read_float()

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_float()


class DoubleReader(Reader):
    def read(self, decoder: BinaryDecoder) -> float:
        return decoder.read_double()

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_double()


class DateReader(IntegerReader):
    """Reads a day granularity date from the stream.

    The number of days from 1 January 1970.
    """


class TimeReader(IntegerReader):
    """Reads a microsecond granularity timestamp from the stream.

    Long is decoded as an integer which represents
    the number of microseconds from the unix epoch, 1 January 1970.
    """


class TimestampReader(IntegerReader):
    """Reads a microsecond granularity timestamp from the stream.

    Long is decoded as python integer which represents
    the number of microseconds from the unix epoch, 1 January 1970.
    """


class TimestampNanoReader(IntegerReader):
    """Reads a nanosecond granularity timestamp from the stream.

    Long is decoded as python integer which represents
    the number of nanoseconds from the unix epoch, 1 January 1970.
    """


class TimestamptzReader(IntegerReader):
    """Reads a microsecond granularity timestamptz from the stream.

    Long is decoded as python integer which represents
    the number of microseconds from the unix epoch, 1 January 1970.

    Adjusted to UTC.
    """


class TimestamptzNanoReader(IntegerReader):
    """Reads a microsecond granularity timestamptz from the stream.

    Long is decoded as python integer which represents
    the number of nanoseconds from the unix epoch, 1 January 1970.

    Adjusted to UTC.
    """


class StringReader(Reader):
    def read(self, decoder: BinaryDecoder) -> str:
        return decoder.read_utf8()

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_utf8()


class UUIDReader(Reader):
    def read(self, decoder: BinaryDecoder) -> UUID:
        return UUID(bytes=decoder.read(16))

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip(16)


class UnknownReader(Reader):
    def read(self, decoder: BinaryDecoder) -> None:
        return None

    def skip(self, decoder: BinaryDecoder) -> None:
        pass


@dataclass(frozen=True)
class FixedReader(Reader):
    _len: int = dataclassfield()

    def read(self, decoder: BinaryDecoder) -> bytes:
        return decoder.read(len(self))

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip(len(self))

    def __len__(self) -> int:
        """Return the length of an instance of the FixedReader class."""
        return self._len

    def __repr__(self) -> str:
        """Return the string representation of the FixedReader class."""
        return f"FixedReader({self._len})"


class BinaryReader(Reader):
    """Read a binary value.

    First reads an integer, to get the length of the binary value,
    then reads the binary field itself.
    """

    def read(self, decoder: BinaryDecoder) -> bytes:
        return decoder.read_bytes()

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_bytes()


@dataclass(frozen=True, init=False)
class DecimalReader(Reader):
    """Reads a value as a decimal.

    Decimal bytes are decoded as signed short, int or long depending on the
    size of bytes.
    """

    precision: int = dataclassfield()
    scale: int = dataclassfield()
    _length: int

    def __init__(self, precision: int, scale: int):
        object.__setattr__(self, "precision", precision)
        object.__setattr__(self, "scale", scale)
        object.__setattr__(self, "_length", decimal_required_bytes(precision))

    def read(self, decoder: BinaryDecoder) -> Decimal:
        return bytes_to_decimal(decoder.read(self._length), self.scale)

    def skip(self, decoder: BinaryDecoder) -> None:
        decoder.skip_bytes()

    def __repr__(self) -> str:
        """Return the string representation of the DecimalReader class."""
        return f"DecimalReader({self.precision}, {self.scale})"


@dataclass(frozen=True)
class OptionReader(Reader):
    option: Reader = dataclassfield()

    def read(self, decoder: BinaryDecoder) -> Any | None:
        # For the Iceberg spec it is required to set the default value to null
        # From https://iceberg.apache.org/spec/#avro
        # Optional fields must always set the Avro field default value to null.
        #
        # This means that null has to come first:
        # https://avro.apache.org/docs/current/spec.html
        # type of the default value must match the first element of the union.
        # This is enforced in the schema conversion, which happens prior
        # to building the reader tree
        if decoder.read_int() > 0:
            return self.option.read(decoder)
        return None

    def skip(self, decoder: BinaryDecoder) -> None:
        if decoder.read_int() > 0:
            return self.option.skip(decoder)


class StructReader(Reader):
    __slots__ = (
        "field_readers",
        "create_struct",
        "struct",
        "_field_reader_functions",
        "_hash",
        "_max_pos",
    )
    field_readers: tuple[tuple[int | None, Reader], ...]
    create_struct: Callable[..., StructProtocol]
    struct: StructType
    field_reader_functions = tuple[tuple[str | None, int, Callable[[BinaryDecoder], Any] | None], ...]

    def __init__(
        self,
        field_readers: tuple[tuple[int | None, Reader], ...],
        create_struct: Callable[..., StructProtocol],
        struct: StructType,
    ) -> None:
        self.field_readers = field_readers
        self.create_struct = create_struct
        # TODO: Implement struct-reuse
        self.struct = struct

        if not isinstance(self.create_struct(), StructProtocol):
            raise ValueError(f"Incompatible with StructProtocol: {self.create_struct}")

        reading_callbacks: list[tuple[int | None, Callable[[BinaryDecoder], Any]]] = []
        max_pos = -1
        for pos, field in field_readers:
            if pos is not None:
                reading_callbacks.append((pos, field.read))
                max_pos = max(max_pos, pos)
            else:
                reading_callbacks.append((None, field.skip))

        self._field_reader_functions = tuple(reading_callbacks)
        self._hash = hash(self._field_reader_functions)
        self._max_pos = 1 + max_pos

    def read(self, decoder: BinaryDecoder) -> StructProtocol:
        # TODO: Implement struct-reuse
        struct = self.create_struct(*[None] * self._max_pos)
        for pos, field_reader in self._field_reader_functions:
            if pos is not None:
                struct[pos] = field_reader(decoder)  # later: pass reuse in here
            else:
                field_reader(decoder)

        return struct

    def skip(self, decoder: BinaryDecoder) -> None:
        for _, field in self.field_readers:
            field.skip(decoder)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the StructReader class."""
        return (
            self.field_readers == other.field_readers and self.create_struct == other.create_struct
            if isinstance(other, StructReader)
            else False
        )

    def __repr__(self) -> str:
        """Return the string representation of the StructReader class."""
        return f"StructReader(({','.join(repr(field) for field in self.field_readers)}), {repr(self.create_struct)})"

    def __hash__(self) -> int:
        """Return a hashed representation of the StructReader class."""
        return self._hash


@dataclass(frozen=False, init=False)
class ListReader(Reader):
    __slots__ = ("element", "_is_int_list", "_hash")
    element: Reader

    def __init__(self, element: Reader) -> None:
        super().__init__()
        self.element = element
        self._hash = hash(self.element)
        self._is_int_list = isinstance(self.element, IntegerReader)

    def read(self, decoder: BinaryDecoder) -> list[Any]:
        read_items: list[Any] = []
        block_count = decoder.read_int()
        while block_count != 0:
            if block_count < 0:
                block_count = -block_count
                _ = decoder.read_int()
            if self._is_int_list:
                read_items.extend(decoder.read_ints(block_count))
            else:
                for _ in range(block_count):
                    read_items.append(self.element.read(decoder))
            block_count = decoder.read_int()
        return read_items

    def skip(self, decoder: BinaryDecoder) -> None:
        _skip_map_array(decoder, lambda: self.element.skip(decoder))

    def __hash__(self) -> int:
        """Return a hashed representation of the ListReader class."""
        return self._hash


# Represent an empty dict as a singleton
EMPTY_DICT: dict[Any, Any] = {}


@dataclass(frozen=False, init=False)
class MapReader(Reader):
    __slots__ = ("key", "value", "_is_int_int", "_is_int_bytes", "_key_reader", "_value_reader", "_hash")
    key: Reader
    value: Reader

    def __init__(self, key: Reader, value: Reader) -> None:
        super().__init__()
        self.key = key
        self.value = value
        if isinstance(self.key, IntegerReader):
            self._is_int_int = isinstance(self.value, IntegerReader)
            self._is_int_bytes = isinstance(self.value, BinaryReader)
        else:
            self._is_int_int = False
            self._is_int_bytes = False
            self._key_reader = self.key.read
            self._value_reader = self.value.read
        self._hash = hash((self.key, self.value))

    def _read_int_int(self, decoder: BinaryDecoder) -> Mapping[int, int]:
        """Read a mapping from int to int from the decoder.

        Read a map of ints to ints from the decoder, since this is such a common
        data type, it is optimized to be faster than the generic map reader, by
        using a lazy dict.

        The time it takes to create the python dictionary is much larger than
        the time it takes to read the data from the decoder as an array, so the
        lazy dict defers creating the python dictionary until it is actually
        accessed.

        """
        block_count = decoder.read_int()

        # Often times the map is empty, so we can just return an empty dict without
        # instancing the LazyDict
        if block_count == 0:
            return EMPTY_DICT

        contents_array: list[tuple[int, ...]] = []

        while block_count != 0:
            if block_count < 0:
                block_count = -block_count
                # We ignore the block size for now
                decoder.skip_int()

            # Since the integers are encoding right next to each other
            # just read them all at once.
            contents_array.append(decoder.read_ints(block_count * 2))
            block_count = decoder.read_int()

        return LazyDict(contents_array)

    def read(self, decoder: BinaryDecoder) -> Mapping[Any, Any]:
        read_items: dict[Any, Any] = {}

        if self._is_int_int or self._is_int_bytes:
            if self._is_int_int:
                return self._read_int_int(decoder)

            block_count = decoder.read_int()
            while block_count != 0:
                if block_count < 0:
                    block_count = -block_count
                    # We ignore the block size for now
                    _ = decoder.read_int()
                decoder.read_int_bytes_dict(block_count, read_items)
                block_count = decoder.read_int()
        else:
            block_count = decoder.read_int()
            while block_count != 0:
                if block_count < 0:
                    block_count = -block_count
                    # We ignore the block size for now
                    _ = decoder.read_int()
                for _ in range(block_count):
                    key = self._key_reader(decoder)
                    read_items[key] = self._value_reader(decoder)
                block_count = decoder.read_int()

        return read_items

    def skip(self, decoder: BinaryDecoder) -> None:
        def skip() -> None:
            self.key.skip(decoder)
            self.value.skip(decoder)

        _skip_map_array(decoder, skip)

    def __hash__(self) -> int:
        """Return a hashed representation of the MapReader class."""
        return self._hash
