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
#include <optional>
#include <string>
#include <variant>

#include "arrow/compare.h"
#include "arrow/type.h"
#include "arrow/util/visibility.h"

namespace arrow {

/// \class ArrayStatistics
/// \brief Statistics for an Array
///
/// Apache Arrow format doesn't have statistics but data source such
/// as Apache Parquet may have statistics. Statistics associated with
/// data source can be read unified API via this class.
struct ARROW_EXPORT ArrayStatistics {
  /// \brief The type for maximum and minimum values. If the target
  /// value exists, one of them is used. `std::nullopt` is used
  /// otherwise.
  using ValueType = std::variant<bool, int64_t, uint64_t, double, std::string>;

  static const std::shared_ptr<DataType>& ValueToArrowType(
      const std::optional<ValueType>& value,
      const std::shared_ptr<DataType>& array_type) {
    if (!value.has_value()) {
      return null();
    }

    struct Visitor {
      const std::shared_ptr<DataType>& array_type;

      const std::shared_ptr<DataType>& operator()(const bool&) { return boolean(); }
      const std::shared_ptr<DataType>& operator()(const int64_t&) { return int64(); }
      const std::shared_ptr<DataType>& operator()(const uint64_t&) { return uint64(); }
      const std::shared_ptr<DataType>& operator()(const double&) { return float64(); }
      const std::shared_ptr<DataType>& operator()(const std::string&) {
        switch (array_type->id()) {
          case Type::STRING:
          case Type::BINARY:
          case Type::FIXED_SIZE_BINARY:
          case Type::LARGE_STRING:
          case Type::LARGE_BINARY:
          case Type::BINARY_VIEW:
          case Type::STRING_VIEW:
            return array_type;
          default:
            return utf8();
        }
      }
    } visitor{array_type};
    return std::visit(visitor, value.value());
  }

  /// \brief The number of null values, may not be set
  std::optional<int64_t> null_count = std::nullopt;

  /// \brief The number of distinct values, may not be set
  std::optional<int64_t> distinct_count = std::nullopt;

  /// \brief The minimum value, may not be set
  std::optional<ValueType> min = std::nullopt;

  /// \brief Compute Arrow type of the minimum value.
  ///
  /// If \ref ValueType is `std::string`, `array_type` may be
  /// used. If `array_type` is a binary-like type such as \ref
  /// arrow::binary and \ref arrow::large_utf8, `array_type` is
  /// returned. \ref arrow::utf8 is returned otherwise.
  ///
  /// If \ref ValueType isn't `std::string`, `array_type` isn't used.
  ///
  /// \param array_type The Arrow type of the associated array.
  ///
  /// \return \ref arrow::null if the minimum value is `std::nullopt`,
  ///         Arrow type based on \ref ValueType of the \ref min
  ///         otherwise.
  const std::shared_ptr<DataType>& MinArrowType(
      const std::shared_ptr<DataType>& array_type) {
    return ValueToArrowType(min, array_type);
  }

  /// \brief Whether the minimum value is exact or not
  bool is_min_exact = false;

  /// \brief The maximum value, may not be set
  std::optional<ValueType> max = std::nullopt;

  /// \brief Compute Arrow type of the maximum value.
  ///
  /// If \ref ValueType is `std::string`, `array_type` may be
  /// used. If `array_type` is a binary-like type such as \ref
  /// arrow::binary and \ref arrow::large_utf8, `array_type` is
  /// returned. \ref arrow::utf8 is returned otherwise.
  ///
  /// If \ref ValueType isn't `std::string`, `array_type` isn't used.
  ///
  /// \param array_type The Arrow type of the associated array.
  ///
  /// \return \ref arrow::null if the maximum value is `std::nullopt`,
  ///         Arrow type based on \ref ValueType of the \ref max
  ///         otherwise.
  const std::shared_ptr<DataType>& MaxArrowType(
      const std::shared_ptr<DataType>& array_type) {
    return ValueToArrowType(max, array_type);
  }

  /// \brief Whether the maximum value is exact or not
  bool is_max_exact = false;

  /// \brief Check two \ref arrow::ArrayStatistics for equality
  ///
  /// \param other The \ref arrow::ArrayStatistics instance to compare against.
  ///
  /// \param equal_options Options used to compare double values for equality.
  ///
  /// \return True if the two \ref arrow::ArrayStatistics instances are equal; otherwise,
  /// false.
  bool Equals(const ArrayStatistics& other,
              const EqualOptions& equal_options = EqualOptions::Defaults()) const {
    return ArrayStatisticsEquals(*this, other, equal_options);
  }

  /// \brief Check two statistics for equality
  bool operator==(const ArrayStatistics& other) const { return Equals(other); }

  /// \brief Check two statistics for not equality
  bool operator!=(const ArrayStatistics& other) const { return !Equals(other); }
};

}  // namespace arrow
