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

from pyiceberg.avro.decoder import BinaryDecoder

class CythonBinaryDecoder(BinaryDecoder):
    def __init__(self, input_contents: bytes) -> None:
        pass

    def tell(self) -> int:
        pass

    def read(self, n: int) -> bytes:
        pass

    def read_boolean(self) -> bool:
        pass

    def read_int(self) -> int:
        pass

    def read_ints(self, count: int) -> tuple[int, ...]:
        pass

    def read_int_bytes_dict(self, count: int, dest: dict[int, bytes]) -> None:
        pass

    def read_bytes(self) -> bytes:
        pass

    def read_float(self) -> float:
        pass

    def read_double(self) -> float:
        pass

    def read_utf8(self) -> str:
        pass

    def skip(self, n: int) -> None:
        pass

    def skip_int(self) -> None:
        pass

    def skip_boolean(self) -> None:
        pass

    def skip_float(self) -> None:
        pass

    def skip_double(self) -> None:
        pass

    def skip_bytes(self) -> None:
        pass

    def skip_utf8(self) -> None:
        pass
