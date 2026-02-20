#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.


def truncate_upper_bound_text_string(value: str, trunc_length: int | None) -> str | None:
    result = value[:trunc_length]
    if result != value:
        chars = [*result]

        for i in range(-1, -len(result) - 1, -1):
            try:
                to_inc = ord(chars[i])
                # will raise exception if the highest unicode code is reached
                _next = chr(to_inc + 1)
                chars[i] = _next
                return "".join(chars)
            except ValueError:
                pass
        return None  # didn't find a valid upper bound
    return result


def truncate_upper_bound_binary_string(value: bytes, trunc_length: int | None) -> bytes | None:
    result = value[:trunc_length]
    if result != value:
        _bytes = [*result]
        for i in range(-1, -len(result) - 1, -1):
            if _bytes[i] < 255:
                _bytes[i] += 1
                return b"".join([i.to_bytes(1, byteorder="little") for i in _bytes])
        return None

    return result
