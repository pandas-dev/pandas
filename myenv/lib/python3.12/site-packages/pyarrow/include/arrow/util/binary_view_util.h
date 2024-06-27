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

#pragma once

#include <string_view>
#include <utility>

#include "arrow/type.h"
#include "arrow/util/span.h"

namespace arrow::util {

inline BinaryViewType::c_type ToInlineBinaryView(const void* data, int32_t size) {
  // Small string: inlined. Bytes beyond size are zeroed
  BinaryViewType::c_type out;
  out.inlined = {size, {}};
  memcpy(&out.inlined.data, data, size);
  return out;
}

inline BinaryViewType::c_type ToInlineBinaryView(std::string_view v) {
  return ToInlineBinaryView(v.data(), static_cast<int32_t>(v.size()));
}

inline BinaryViewType::c_type ToBinaryView(const void* data, int32_t size,
                                           int32_t buffer_index, int32_t offset) {
  if (size <= BinaryViewType::kInlineSize) {
    return ToInlineBinaryView(data, size);
  }

  // Large string: store index/offset.
  BinaryViewType::c_type out;
  out.ref = {size, {}, buffer_index, offset};
  memcpy(&out.ref.prefix, data, sizeof(out.ref.prefix));
  return out;
}

inline BinaryViewType::c_type ToBinaryView(std::string_view v, int32_t buffer_index,
                                           int32_t offset) {
  return ToBinaryView(v.data(), static_cast<int32_t>(v.size()), buffer_index, offset);
}

template <typename BufferPtr>
std::string_view FromBinaryView(const BinaryViewType::c_type& v,
                                const BufferPtr* data_buffers) {
  auto* data = v.is_inline() ? v.inlined.data.data()
                             : data_buffers[v.ref.buffer_index]->data() + v.ref.offset;
  return {reinterpret_cast<const char*>(data), static_cast<size_t>(v.size())};
}
template <typename BufferPtr>
std::string_view FromBinaryView(BinaryViewType::c_type&&, const BufferPtr*) = delete;

template <typename BufferPtr>
bool EqualBinaryView(BinaryViewType::c_type l, BinaryViewType::c_type r,
                     const BufferPtr* l_buffers, const BufferPtr* r_buffers) {
  int64_t l_size_and_prefix, r_size_and_prefix;
  memcpy(&l_size_and_prefix, &l, sizeof(l_size_and_prefix));
  memcpy(&r_size_and_prefix, &r, sizeof(r_size_and_prefix));

  if (l_size_and_prefix != r_size_and_prefix) return false;

  if (l.is_inline()) {
    // The columnar spec mandates that the inlined part be zero-padded, so we can compare
    // a word at a time regardless of the exact size.
    int64_t l_inlined, r_inlined;
    memcpy(&l_inlined, l.inline_data() + BinaryViewType::kPrefixSize, sizeof(l_inlined));
    memcpy(&r_inlined, r.inline_data() + BinaryViewType::kPrefixSize, sizeof(r_inlined));
    return l_inlined == r_inlined;
  }

  // Sizes are equal and this is not inline, therefore both are out
  // of line and have kPrefixSize first in common.
  const uint8_t* l_data = l_buffers[l.ref.buffer_index]->data() + l.ref.offset;
  const uint8_t* r_data = r_buffers[r.ref.buffer_index]->data() + r.ref.offset;
  return memcmp(l_data + BinaryViewType::kPrefixSize,
                r_data + BinaryViewType::kPrefixSize,
                l.size() - BinaryViewType::kPrefixSize) == 0;
}

}  // namespace arrow::util
