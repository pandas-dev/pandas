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
#include <iosfwd>
#include <optional>
#include <vector>

#include "arrow/util/span.h"
#include "parquet/platform.h"
#include "parquet/type_fwd.h"

namespace parquet {

/// A structure for capturing metadata for estimating the unencoded,
/// uncompressed size of data written. This is useful for readers to estimate
/// how much memory is needed to reconstruct data in their memory model and for
/// fine-grained filter push down on nested structures (the histograms contained
/// in this structure can help determine the number of nulls at a particular
/// nesting level and maximum length of lists).
struct PARQUET_EXPORT SizeStatistics {
  /// When present, there is expected to be one element corresponding to each
  /// definition (i.e. size=max definition+1) where each element
  /// represents the number of times the definition level was observed in the
  /// data.
  ///
  /// This field may be omitted (a.k.a. zero-length vector) if max_definition_level
  /// is 0 without loss of information.
  std::vector<int64_t> definition_level_histogram;

  /// Same as definition_level_histogram except for repetition levels.
  ///
  /// This field may be omitted (a.k.a. zero-length vector) if max_repetition_level
  /// is 0 without loss of information.
  std::vector<int64_t> repetition_level_histogram;

  /// The number of physical bytes stored for BYTE_ARRAY data values assuming
  /// no encoding. This is exclusive of the bytes needed to store the length of
  /// each byte array. In other words, this field is equivalent to the `(size
  /// of PLAIN-ENCODING the byte array values) - (4 bytes * number of values
  /// written)`. To determine unencoded sizes of other types readers can use
  /// schema information multiplied by the number of non-null and null values.
  /// The number of null/non-null values can be inferred from the histograms
  /// below.
  ///
  /// For example, if a column chunk is dictionary-encoded with dictionary
  /// ["a", "bc", "cde"], and a data page contains the indices [0, 0, 1, 2],
  /// then this value for that data page should be 7 (1 + 1 + 2 + 3).
  ///
  /// This field should only be set for types that use BYTE_ARRAY as their
  /// physical type.
  std::optional<int64_t> unencoded_byte_array_data_bytes;

  /// \brief Check if the SizeStatistics is set.
  bool is_set() const {
    return !repetition_level_histogram.empty() || !definition_level_histogram.empty() ||
           unencoded_byte_array_data_bytes.has_value();
  }

  /// \brief Increment the unencoded byte array data bytes.
  void IncrementUnencodedByteArrayDataBytes(int64_t value);

  /// \brief Merge two SizeStatistics.
  /// \throws ParquetException if SizeStatistics to merge is not compatible.
  void Merge(const SizeStatistics& other);

  /// \brief Validate the SizeStatistics
  /// \throws ParquetException if the histograms don't have the right length,
  /// or if unencoded_byte_array_data_bytes is present for a non-BYTE_ARRAY column.
  void Validate(const ColumnDescriptor* descr) const;

  /// \brief Reset the SizeStatistics to be empty.
  void Reset();

  /// \brief Make an empty SizeStatistics object for specific type.
  static std::unique_ptr<SizeStatistics> Make(const ColumnDescriptor* descr);
};

PARQUET_EXPORT
std::ostream& operator<<(std::ostream&, const SizeStatistics&);

PARQUET_EXPORT
void UpdateLevelHistogram(::arrow::util::span<const int16_t> levels,
                          ::arrow::util::span<int64_t> histogram);

}  // namespace parquet
