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
#include "arrow/type.h"

namespace arrow::extension {

/// \brief Opaque is a placeholder for a type from an external (usually
///   non-Arrow) system that could not be interpreted.
class ARROW_EXPORT OpaqueType : public ExtensionType {
 public:
  /// \brief Construct an OpaqueType.
  ///
  /// \param[in] storage_type The underlying storage type.  Should be
  ///   arrow::null if there is no data.
  /// \param[in] type_name The name of the type in the external system.
  /// \param[in] vendor_name The name of the external system.
  explicit OpaqueType(std::shared_ptr<DataType> storage_type, std::string type_name,
                      std::string vendor_name)
      : ExtensionType(std::move(storage_type)),
        type_name_(std::move(type_name)),
        vendor_name_(std::move(vendor_name)) {}

  std::string extension_name() const override { return "arrow.opaque"; }
  std::string ToString(bool show_metadata) const override;
  bool ExtensionEquals(const ExtensionType& other) const override;
  std::string Serialize() const override;
  Result<std::shared_ptr<DataType>> Deserialize(
      std::shared_ptr<DataType> storage_type,
      const std::string& serialized_data) const override;
  /// Create an OpaqueArray from ArrayData
  std::shared_ptr<Array> MakeArray(std::shared_ptr<ArrayData> data) const override;

  std::string_view type_name() const { return type_name_; }
  std::string_view vendor_name() const { return vendor_name_; }

 private:
  std::string type_name_;
  std::string vendor_name_;
};

/// \brief Opaque is a wrapper for (usually binary) data from an external
///   (often non-Arrow) system that could not be interpreted.
class ARROW_EXPORT OpaqueArray : public ExtensionArray {
 public:
  using ExtensionArray::ExtensionArray;
};

/// \brief Return an OpaqueType instance.
ARROW_EXPORT std::shared_ptr<DataType> opaque(std::shared_ptr<DataType> storage_type,
                                              std::string type_name,
                                              std::string vendor_name);

}  // namespace arrow::extension
