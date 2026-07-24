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
# pylint: disable=W0621
"""Avro reader for reading Avro files."""

from __future__ import annotations

import io
import json
import os
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from types import TracebackType
from typing import (
    Generic,
    TypeVar,
)

from pyiceberg.avro.codecs import AVRO_CODEC_KEY, CODEC_MAPPING_ICEBERG_TO_AVRO, KNOWN_CODECS
from pyiceberg.avro.codecs.codec import Codec
from pyiceberg.avro.decoder import BinaryDecoder, new_decoder
from pyiceberg.avro.encoder import BinaryEncoder
from pyiceberg.avro.reader import Reader
from pyiceberg.avro.resolver import construct_reader, construct_writer, resolve_reader, resolve_writer
from pyiceberg.avro.writer import Writer
from pyiceberg.io import InputFile, OutputFile, OutputStream
from pyiceberg.schema import Schema
from pyiceberg.typedef import EMPTY_DICT, Record, StructProtocol
from pyiceberg.types import (
    FixedType,
    MapType,
    NestedField,
    StringType,
    StructType,
)
from pyiceberg.utils.schema_conversion import AvroSchemaConversion

VERSION = 1
MAGIC = bytes(b"Obj" + bytearray([VERSION]))
MAGIC_SIZE = len(MAGIC)
SYNC_SIZE = 16
META_SCHEMA = StructType(
    NestedField(name="magic", field_id=100, field_type=FixedType(length=MAGIC_SIZE), required=True),
    NestedField(
        field_id=200,
        name="meta",
        field_type=MapType(key_id=201, key_type=StringType(), value_id=202, value_type=StringType(), value_required=True),
        required=True,
    ),
    NestedField(field_id=300, name="sync", field_type=FixedType(length=SYNC_SIZE), required=True),
)

_SCHEMA_KEY = "avro.schema"


class AvroFileHeader(Record):
    @property
    def magic(self) -> bytes:
        return self._data[0]

    @property
    def meta(self) -> dict[str, str]:
        return self._data[1]

    @property
    def sync(self) -> bytes:
        return self._data[2]

    def compression_codec(self) -> type[Codec] | None:
        """Get the file's compression codec algorithm from the file's metadata.

        In the case of a null codec, we return a None indicating that we
        don't need to compress/decompress.
        """
        from pyiceberg.table import TableProperties

        codec_name = self.meta.get(AVRO_CODEC_KEY, TableProperties.WRITE_AVRO_COMPRESSION_DEFAULT)
        if codec_name not in KNOWN_CODECS:
            raise ValueError(f"Unsupported codec: {codec_name}")

        return KNOWN_CODECS[codec_name]  # type: ignore

    def get_schema(self) -> Schema:
        if _SCHEMA_KEY in self.meta:
            avro_schema_string = self.meta[_SCHEMA_KEY]
            avro_schema = json.loads(avro_schema_string)
            return AvroSchemaConversion().avro_to_iceberg(avro_schema)
        else:
            raise ValueError("No schema found in Avro file headers")


D = TypeVar("D", bound=StructProtocol)


@dataclass
class Block(Generic[D]):
    reader: Reader
    block_records: int
    block_decoder: BinaryDecoder
    position: int = 0

    def __iter__(self) -> Block[D]:
        """Return an iterator for the Block class."""
        return self

    def has_next(self) -> bool:
        return self.position < self.block_records

    def __next__(self) -> D:
        """Return the next item when iterating over the Block class."""
        if self.has_next():
            self.position += 1
            return self.reader.read(self.block_decoder)
        raise StopIteration


class AvroFile(Generic[D]):
    __slots__ = (
        "input_file",
        "read_schema",
        "read_types",
        "read_enums",
        "header",
        "schema",
        "reader",
        "decoder",
        "block",
    )
    input_file: InputFile
    read_schema: Schema | None
    read_types: dict[int, Callable[..., StructProtocol]]
    read_enums: dict[int, Callable[..., Enum]]
    header: AvroFileHeader
    schema: Schema
    reader: Reader

    decoder: BinaryDecoder
    block: Block[D] | None

    def __init__(
        self,
        input_file: InputFile,
        read_schema: Schema | None = None,
        read_types: dict[int, Callable[..., StructProtocol]] = EMPTY_DICT,
        read_enums: dict[int, Callable[..., Enum]] = EMPTY_DICT,
    ) -> None:
        self.input_file = input_file
        self.read_schema = read_schema
        self.read_types = read_types
        self.read_enums = read_enums
        self.block = None

    def __enter__(self) -> AvroFile[D]:
        """Generate a reader tree for the payload within an avro file.

        Return:
            A generator returning the AvroStructs.
        """
        with self.input_file.open() as f:
            self.decoder = new_decoder(f.read())
        self.header = self._read_header()
        self.schema = self.header.get_schema()
        if not self.read_schema:
            self.read_schema = self.schema

        self.reader = resolve_reader(self.schema, self.read_schema, self.read_types, self.read_enums)

        return self

    def __exit__(self, exctype: type[BaseException] | None, excinst: BaseException | None, exctb: TracebackType | None) -> None:
        """Perform cleanup when exiting the scope of a 'with' statement."""

    def __iter__(self) -> AvroFile[D]:
        """Return an iterator for the AvroFile class."""
        return self

    def _read_block(self) -> int:
        # If there is already a block, we'll have the sync bytes
        if self.block:
            sync_marker = self.decoder.read(SYNC_SIZE)
            if sync_marker != self.header.sync:
                raise ValueError(f"Expected sync bytes {self.header.sync!r}, but got {sync_marker!r}")
        block_records = self.decoder.read_int()

        block_bytes = self.decoder.read_bytes()
        if codec := self.header.compression_codec():
            block_bytes = codec.decompress(block_bytes)

        self.block = Block(reader=self.reader, block_records=block_records, block_decoder=new_decoder(block_bytes))
        return block_records

    def __next__(self) -> D:
        """Return the next item when iterating over the AvroFile class."""
        if self.block and self.block.has_next():
            return next(self.block)

        try:
            new_block = self._read_block()
        except EOFError as exc:
            raise StopIteration from exc

        if new_block > 0:
            return self.__next__()
        raise StopIteration

    def _read_header(self) -> AvroFileHeader:
        return construct_reader(META_SCHEMA, {-1: AvroFileHeader}).read(self.decoder)


class AvroOutputFile(Generic[D]):
    output_file: OutputFile
    output_stream: OutputStream
    file_schema: Schema
    schema_name: str
    encoder: BinaryEncoder
    sync_bytes: bytes
    writer: Writer

    def __init__(
        self,
        output_file: OutputFile,
        file_schema: Schema,
        schema_name: str,
        record_schema: Schema | None = None,
        metadata: dict[str, str] = EMPTY_DICT,
    ) -> None:
        self.output_file = output_file
        self.file_schema = file_schema
        self.schema_name = schema_name
        self.sync_bytes = os.urandom(SYNC_SIZE)
        self.writer = (
            construct_writer(file_schema=self.file_schema)
            if record_schema is None
            else resolve_writer(record_schema=record_schema, file_schema=self.file_schema)
        )
        self.metadata = metadata

    def __enter__(self) -> AvroOutputFile[D]:
        """
        Open the file and writes the header.

        Returns:
            The file object to write records to
        """
        self.output_stream = self.output_file.create(overwrite=True)
        self.encoder = BinaryEncoder(self.output_stream)

        self._write_header()

        return self

    def __exit__(self, exctype: type[BaseException] | None, excinst: BaseException | None, exctb: TracebackType | None) -> None:
        """Perform cleanup when exiting the scope of a 'with' statement."""
        self.output_stream.close()

    def _write_header(self) -> None:
        from pyiceberg.table import TableProperties

        codec_name = self.metadata.get(AVRO_CODEC_KEY, TableProperties.WRITE_AVRO_COMPRESSION_DEFAULT)
        if avro_codec_name := CODEC_MAPPING_ICEBERG_TO_AVRO.get(codec_name):
            codec_name = avro_codec_name

        json_schema = json.dumps(AvroSchemaConversion().iceberg_to_avro(self.file_schema, schema_name=self.schema_name))

        meta = {**self.metadata, _SCHEMA_KEY: json_schema, AVRO_CODEC_KEY: codec_name}
        header = AvroFileHeader(MAGIC, meta, self.sync_bytes)
        construct_writer(META_SCHEMA).write(self.encoder, header)

    def compression_codec(self) -> type[Codec] | None:
        """Get the file's compression codec algorithm from the file's metadata.

        In the case of a null codec, we return a None indicating that we
        don't need to compress/decompress.
        """
        from pyiceberg.table import TableProperties

        codec_name = self.metadata.get(AVRO_CODEC_KEY, TableProperties.WRITE_AVRO_COMPRESSION_DEFAULT)

        if avro_codec_name := CODEC_MAPPING_ICEBERG_TO_AVRO.get(codec_name):
            codec_name = avro_codec_name

        if codec_name not in KNOWN_CODECS:
            raise ValueError(f"Unsupported codec: {codec_name}")

        return KNOWN_CODECS[codec_name]  # type: ignore

    def write_block(self, objects: list[D]) -> None:
        in_memory = io.BytesIO()
        block_content_encoder = BinaryEncoder(output_stream=in_memory)
        for obj in objects:
            self.writer.write(block_content_encoder, obj)
        block_content = in_memory.getvalue()

        self.encoder.write_int(len(objects))

        if codec := self.compression_codec():
            content, content_length = codec.compress(block_content)
            self.encoder.write_int(content_length)
            self.encoder.write(content)
        else:
            self.encoder.write_int(len(block_content))
            self.encoder.write(block_content)

        self.encoder.write(self.sync_bytes)
