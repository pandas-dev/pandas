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

// Contains utilities for making UBSan happy.

#pragma once

#include <cstring>
#include <memory>
#include <type_traits>

#include "arrow/util/aligned_storage.h"
#include "arrow/util/macros.h"

namespace arrow {
namespace util {

namespace internal {

constexpr uint8_t kNonNullFiller = 0;

}  // namespace internal

/// \brief Returns maybe_null if not null or a non-null pointer to an arbitrary memory
/// that shouldn't be dereferenced.
///
/// Memset/Memcpy are undefined when a nullptr is passed as an argument use this utility
/// method to wrap locations where this could happen.
///
/// Note: Flatbuffers has UBSan warnings if a zero length vector is passed.
/// https://github.com/google/flatbuffers/pull/5355 is trying to resolve
/// them.
template <typename T>
inline T* MakeNonNull(T* maybe_null = NULLPTR) {
  if (ARROW_PREDICT_TRUE(maybe_null != NULLPTR)) {
    return maybe_null;
  }

  return const_cast<T*>(reinterpret_cast<const T*>(&internal::kNonNullFiller));
}

template <typename T>
inline std::enable_if_t<std::is_trivially_copyable_v<T>, T> SafeLoadAs(
    const uint8_t* unaligned) {
  using Type = std::remove_const_t<T>;
  arrow::internal::AlignedStorage<Type> raw_data;
  std::memcpy(raw_data.get(), unaligned, sizeof(T));
  auto data = *raw_data.get();
  raw_data.destroy();
  return data;
}

template <typename T>
inline std::enable_if_t<std::is_trivially_copyable_v<T>, T> SafeLoad(const T* unaligned) {
  using Type = std::remove_const_t<T>;
  arrow::internal::AlignedStorage<Type> raw_data;
  std::memcpy(raw_data.get(), static_cast<const void*>(unaligned), sizeof(T));
  auto data = *raw_data.get();
  raw_data.destroy();
  return data;
}

template <typename U, typename T>
inline std::enable_if_t<std::is_trivially_copyable_v<T> &&
                            std::is_trivially_copyable_v<U> && sizeof(T) == sizeof(U),
                        U>
SafeCopy(T value) {
  using TypeU = std::remove_const_t<U>;
  arrow::internal::AlignedStorage<TypeU> raw_data;
  std::memcpy(raw_data.get(), static_cast<const void*>(&value), sizeof(T));
  auto data = *raw_data.get();
  raw_data.destroy();
  return data;
}

template <typename T>
inline std::enable_if_t<std::is_trivially_copyable_v<T>, void> SafeStore(void* unaligned,
                                                                         T value) {
  std::memcpy(unaligned, &value, sizeof(T));
}

}  // namespace util
}  // namespace arrow
