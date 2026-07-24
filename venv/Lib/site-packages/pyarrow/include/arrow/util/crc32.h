// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#include <cstddef>
#include <cstdint>

#include "arrow/util/visibility.h"

namespace arrow {
namespace internal {

/// \brief Compute the CRC32 checksum of the given data
///
/// This function computes CRC32 with the polynomial 0x04C11DB7,
/// as used in zlib and others (note this is different from CRC32C).
/// To compute a running CRC32, pass the previous value in `prev`,
/// otherwise `prev` should be 0.
ARROW_EXPORT
uint32_t crc32(uint32_t prev, const void* data, size_t length);

}  // namespace internal
}  // namespace arrow
