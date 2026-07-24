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

#include <cstdint>
#include <string>

#include "arrow/util/span.h"
#include "arrow/util/visibility.h"

namespace arrow::util {
/**
 * A secure string that ensures the wrapped string is cleared from memory on
 * deconstruction. This class can only be created from std::string that are securely
 * erased after creation.
 *
 * Note: This class does not provide a constructor / assignment operator that copies a
 * std::string because that would allow code to create a SecureString while accidentally
 * not noticing the need to securely erasing the argument after invoking the constructor /
 * calling the assignment operator.
 */
class ARROW_EXPORT SecureString {
 public:
  SecureString() = default;
  SecureString(SecureString&&) noexcept;
  SecureString(const SecureString&) = default;
  explicit SecureString(std::string&&) noexcept;
  explicit SecureString(size_t, char) noexcept;

  SecureString& operator=(SecureString&&) noexcept;
  SecureString& operator=(const SecureString&);
  SecureString& operator=(std::string&&) noexcept;

  bool operator==(const SecureString&) const;
  bool operator!=(const SecureString&) const;

  ~SecureString() { Dispose(); }

  [[nodiscard]] bool empty() const;
  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] std::size_t length() const;
  [[nodiscard]] std::size_t capacity() const;

  [[nodiscard]] span<uint8_t> as_span();
  [[nodiscard]] span<const uint8_t> as_span() const;
  [[nodiscard]] std::string_view as_view() const;

  void Dispose();

  static void SecureClear(std::string*);
  static void SecureClear(uint8_t* data, size_t size);

 private:
  std::string secret_;
};

}  // namespace arrow::util
