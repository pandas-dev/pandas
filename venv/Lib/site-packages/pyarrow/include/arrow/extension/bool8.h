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

#include "arrow/extension_type.h"

namespace arrow::extension {

/// \brief Bool8 is an alternate representation for boolean
/// arrays using 8 bits instead of 1 bit per value. The underlying
/// storage type is int8.
class ARROW_EXPORT Bool8Array : public ExtensionArray {
 public:
  using ExtensionArray::ExtensionArray;
};

/// \brief Bool8 is an alternate representation for boolean
/// arrays using 8 bits instead of 1 bit per value. The underlying
/// storage type is int8.
class ARROW_EXPORT Bool8Type : public ExtensionType {
 public:
  /// \brief Construct a Bool8Type.
  Bool8Type() : ExtensionType(int8()) {}

  std::string extension_name() const override { return "arrow.bool8"; }
  std::string ToString(bool show_metadata = false) const override;

  bool ExtensionEquals(const ExtensionType& other) const override;

  std::string Serialize() const override;

  Result<std::shared_ptr<DataType>> Deserialize(
      std::shared_ptr<DataType> storage_type,
      const std::string& serialized_data) const override;

  /// Create a Bool8Array from ArrayData
  std::shared_ptr<Array> MakeArray(std::shared_ptr<ArrayData> data) const override;

  static Result<std::shared_ptr<DataType>> Make();
};

/// \brief Return a Bool8Type instance.
ARROW_EXPORT std::shared_ptr<DataType> bool8();

}  // namespace arrow::extension
