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
from __future__ import annotations

import codecs
import gzip
from abc import ABC, abstractmethod
from collections.abc import Callable

from pyiceberg.io import InputFile, InputStream, OutputFile
from pyiceberg.table.metadata import TableMetadata, TableMetadataUtil
from pyiceberg.typedef import UTF8
from pyiceberg.utils.config import Config

GZIP = "gzip"


class Compressor(ABC):
    @staticmethod
    def get_compressor(location: str) -> Compressor:
        return GzipCompressor() if location.endswith(".gz.metadata.json") else NOOP_COMPRESSOR

    @abstractmethod
    def stream_decompressor(self, inp: InputStream) -> InputStream:
        """Return a stream decompressor.

        Args:
            inp: The input stream that needs decompressing.

        Returns:
            The wrapped stream
        """

    @abstractmethod
    def bytes_compressor(self) -> Callable[[bytes], bytes]:
        """Return a function to compress bytes.

        Returns:
            A function that can be used to compress bytes.
        """


class NoopCompressor(Compressor):
    def stream_decompressor(self, inp: InputStream) -> InputStream:
        return inp

    def bytes_compressor(self) -> Callable[[bytes], bytes]:
        return lambda b: b


NOOP_COMPRESSOR = NoopCompressor()


class GzipCompressor(Compressor):
    def stream_decompressor(self, inp: InputStream) -> InputStream:
        return gzip.open(inp)

    def bytes_compressor(self) -> Callable[[bytes], bytes]:
        return gzip.compress


class FromByteStream:
    """A collection of methods that deserialize dictionaries into Iceberg objects."""

    @staticmethod
    def table_metadata(
        byte_stream: InputStream, encoding: str = UTF8, compression: Compressor = NOOP_COMPRESSOR
    ) -> TableMetadata:
        """Instantiate a TableMetadata object from a byte stream.

        Args:
            byte_stream: A file-like byte stream object.
            encoding (default "utf-8"): The byte encoder to use for the reader.
            compression: Optional compression method
        """
        with compression.stream_decompressor(byte_stream) as byte_stream:
            reader = codecs.getreader(encoding)
            json_bytes = reader(byte_stream)
            metadata = json_bytes.read()

        return TableMetadataUtil.parse_raw(metadata)


class FromInputFile:
    """A collection of methods that deserialize InputFiles into Iceberg objects."""

    @staticmethod
    def table_metadata(input_file: InputFile, encoding: str = UTF8) -> TableMetadata:
        """Create a TableMetadata instance from an input file.

        Args:
            input_file (InputFile): A custom implementation of the iceberg.io.file.InputFile abstract base class.
            encoding (str): Encoding to use when loading bytestream.

        Returns:
            TableMetadata: A table metadata instance.

        """
        with input_file.open() as input_stream:
            return FromByteStream.table_metadata(
                byte_stream=input_stream, encoding=encoding, compression=Compressor.get_compressor(location=input_file.location)
            )


class ToOutputFile:
    """A collection of methods that serialize Iceberg objects into files given an OutputFile instance."""

    @staticmethod
    def table_metadata(metadata: TableMetadata, output_file: OutputFile, overwrite: bool = False) -> None:
        """Write a TableMetadata instance to an output file.

        Args:
            output_file (OutputFile): A custom implementation of the iceberg.io.file.OutputFile abstract base class.
            overwrite (bool): Where to overwrite the file if it already exists. Defaults to `False`.
        """
        with output_file.create(overwrite=overwrite) as output_stream:
            # We need to serialize None values, in order to dump `None` current-snapshot-id as `-1`
            exclude_none = False if Config().get_bool("legacy-current-snapshot-id") else True

            json_bytes = metadata.model_dump_json(exclude_none=exclude_none).encode(UTF8)
            json_bytes = Compressor.get_compressor(output_file.location).bytes_compressor()(json_bytes)
            output_stream.write(json_bytes)
