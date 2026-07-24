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

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <type_traits>

namespace arrow::util {

template <class T>
class span;

/// std::span polyfill.
///
/// Does not support static extents.
template <typename T>
class span {
  static_assert(sizeof(T),
                R"(
std::span allows contiguous_iterators instead of just pointers, the enforcement
of which requires T to be a complete type. arrow::util::span does not support
contiguous_iterators, but T is still required to be a complete type to prevent
writing code which would break when it is replaced by std::span.)");

 public:
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using iterator = T*;
  using const_iterator = T const*;

  span() = default;
  span(const span&) = default;
  span& operator=(const span&) = default;

  template <typename M, typename = std::enable_if_t<std::is_same_v<T, M const>>>
  // NOLINTNEXTLINE runtime/explicit
  constexpr span(span<M> mut) : span{mut.data(), mut.size()} {}

  constexpr span(T* data, size_t count) : data_{data}, size_{count} {}

  constexpr span(T* begin, T* end)
      : data_{begin}, size_{static_cast<size_t>(end - begin)} {}

  template <typename R, typename RD = decltype(std::data(std::declval<R>())),
            typename RS = decltype(std::size(std::declval<R>())),
            typename E = std::enable_if_t<std::is_constructible_v<T*, RD> &&
                                          std::is_constructible_v<size_t, RS>>>
  // NOLINTNEXTLINE runtime/explicit, non-const reference
  constexpr span(R&& range) : data_{std::data(range)}, size_{std::size(range)} {}

  constexpr T* begin() const { return data_; }
  constexpr T* end() const { return data_ + size_; }
  constexpr T* data() const { return data_; }

  constexpr size_t size() const { return size_; }
  constexpr size_t size_bytes() const { return size_ * sizeof(T); }
  constexpr bool empty() const { return size_ == 0; }

  constexpr T& operator[](size_t i) { return data_[i]; }
  constexpr const T& operator[](size_t i) const { return data_[i]; }

  constexpr span subspan(size_t offset) const {
    if (offset > size_) return {data_, data_};
    return {data_ + offset, size_ - offset};
  }

  constexpr span subspan(size_t offset, size_t count) const {
    auto out = subspan(offset);
    if (count < out.size_) {
      out.size_ = count;
    }
    return out;
  }

  constexpr bool operator==(const span& other) const {
    if (size_ != other.size_) return false;

    if constexpr (std::is_integral_v<T>) {
      if (size_ == 0) {
        return true;  // memcmp does not handle null pointers, even if size_ == 0
      }
      return std::memcmp(data_, other.data_, size_bytes()) == 0;
    } else {
      T* ptr = data_;
      for (T const& e : other) {
        if (*ptr++ != e) return false;
      }
      return true;
    }
  }
  constexpr bool operator!=(const span& other) const { return !(*this == other); }

 private:
  T* data_{};
  size_t size_{};
};

template <typename R>
span(R& range) -> span<std::remove_pointer_t<decltype(std::data(range))>>;

template <typename T>
span(T*, size_t) -> span<T>;

template <typename T>
constexpr span<std::byte const> as_bytes(span<T> s) {
  return {reinterpret_cast<const std::byte*>(s.data()), s.size_bytes()};
}

template <typename T>
constexpr span<std::byte> as_writable_bytes(span<T> s) {
  return {reinterpret_cast<std::byte*>(s.data()), s.size_bytes()};
}

}  // namespace arrow::util
