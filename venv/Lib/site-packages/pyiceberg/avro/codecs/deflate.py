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

import zlib

from pyiceberg.avro.codecs.codec import Codec


class DeflateCodec(Codec):
    @staticmethod
    def compress(data: bytes) -> tuple[bytes, int]:
        # The first two characters and last character are zlib
        # wrappers around deflate data.
        compressed_data = zlib.compress(data)[2:-1]
        return compressed_data, len(compressed_data)

    @staticmethod
    def decompress(data: bytes) -> bytes:
        # -15 is the log of the window size; negative indicates
        # "raw" (no zlib headers) decompression.  See zlib.h.
        return zlib.decompress(data, -15)
