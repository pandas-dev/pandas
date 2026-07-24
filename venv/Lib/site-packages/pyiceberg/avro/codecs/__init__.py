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
Contains Codecs for Python Avro.

Note that the word "codecs" means "compression/decompression algorithms" in the
Avro world (https://avro.apache.org/docs/current/spec.html#Object+Container+Files),
so don't confuse it with the Python's "codecs", which is a package mainly for
converting character sets (https://docs.python.org/3/library/codecs.html).
"""

from __future__ import annotations

from typing import Literal, TypeAlias

from pyiceberg.avro.codecs.bzip2 import BZip2Codec
from pyiceberg.avro.codecs.codec import Codec
from pyiceberg.avro.codecs.deflate import DeflateCodec
from pyiceberg.avro.codecs.snappy_codec import SnappyCodec
from pyiceberg.avro.codecs.zstandard_codec import ZStandardCodec

AvroCompressionCodec: TypeAlias = Literal["null", "bzip2", "snappy", "zstandard", "deflate"]

AVRO_CODEC_KEY = "avro.codec"

KNOWN_CODECS: dict[AvroCompressionCodec, type[Codec] | None] = {
    "null": None,
    "bzip2": BZip2Codec,
    "snappy": SnappyCodec,
    "zstandard": ZStandardCodec,
    "deflate": DeflateCodec,
}

# Map to convert the naming from Iceberg to Avro
CODEC_MAPPING_ICEBERG_TO_AVRO: dict[str, str] = {"gzip": "deflate", "zstd": "zstandard"}
