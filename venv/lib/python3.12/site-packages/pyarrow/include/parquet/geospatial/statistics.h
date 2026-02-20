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
#include <memory>
#include <optional>

#include "parquet/platform.h"
#include "parquet/types.h"

namespace parquet::geospatial {

/// \brief The maximum number of dimensions represented by a geospatial type
/// (i.e., X, Y, Z, and M)
inline constexpr int kMaxDimensions = 4;

/// \brief NaN, used to represent bounds for which predicate pushdown cannnot
/// be applied (e.g., because a writer did not provide bounds for a given dimension)
inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// \brief Structure represented encoded statistics to be written to and read from Parquet
/// serialized metadata.
///
/// See the Parquet Thrift definition and GeoStatistics for the specific definition
/// of field values.
struct PARQUET_EXPORT EncodedGeoStatistics {
  bool xy_bounds_present{false};
  double xmin{kNaN};
  double xmax{kNaN};
  double ymin{kNaN};
  double ymax{kNaN};

  bool z_bounds_present{false};
  double zmin{kNaN};
  double zmax{kNaN};

  bool m_bounds_present{false};
  double mmin{kNaN};
  double mmax{kNaN};

  bool geospatial_types_present() const { return !geospatial_types.empty(); }
  std::vector<int32_t> geospatial_types;
};

class GeoStatisticsImpl;

/// \brief Base type for computing geospatial column statistics while writing a file
/// or representing them when reading a file
///
/// These statistics track the minimum and maximum value (omitting NaN values) of the
/// four possible dimensions (X, Y, Z, and M) and the distinct set of geometry
/// type/dimension combinations (e.g., point XY, linestring XYZM) present in the data.
/// Any of these individual components may be "invalid": for example, when reading a
/// Parquet file, information about individual components obtained from the column
/// chunk metadata may have been missing or deemed unusable. Orthogonally,
/// any of these individual components may be "empty": for example, when using
/// GeoStatistics to accumulate bounds whilst writing, if all geometries in a column chunk
/// are null, all ranges (X, Y, Z, and M) will be empty. If all geometries in a column
/// chunk contain only XY coordinates (the most common case), the Z and M ranges will
/// be empty but the X and Y ranges will contain finite bounds. Empty ranges are
/// considered "valid" because they are known to represent exactly zero values (in
/// contrast to an invalid range, whose contents is completely unknown). These concepts
/// are all necessary for this object to accurately represent (1) accumulated or partially
/// accumulated statistics during the writing process and (2) deserialized statistics read
/// from the column chunk metadata during the reading process.
///
/// EXPERIMENTAL
class PARQUET_EXPORT GeoStatistics {
 public:
  GeoStatistics();
  explicit GeoStatistics(const EncodedGeoStatistics& encoded);

  ~GeoStatistics();

  /// \brief Return true if bounds, geometry types, and validity are identical
  bool Equals(const GeoStatistics& other) const;

  /// \brief Update these statistics based on previously calculated or decoded statistics
  ///
  /// Merging statistics with wraparound X values is not currently supported. Merging
  /// two GeoStatistics where one or both has a wraparound X range will result in these
  /// statistics having an X dimension marked as invalid.
  void Merge(const GeoStatistics& other);

  /// \brief Update these statistics based on values
  void Update(const ByteArray* values, int64_t num_values);

  /// \brief Update these statistics based on the non-null elements of values
  void UpdateSpaced(const ByteArray* values, const uint8_t* valid_bits,
                    int64_t valid_bits_offset, int64_t num_spaced_values,
                    int64_t num_values);

  /// \brief Update these statistics based on the non-null elements of values
  ///
  /// Currently, BinaryArray and LargeBinaryArray input is supported.
  void Update(const ::arrow::Array& values);

  /// \brief Return these statistics to an empty state
  void Reset();

  /// \brief Encode the statistics for serializing to Thrift
  ///
  /// If invalid WKB was encountered or if the statistics contain NaN
  /// for any reason, Encode() will return nullopt to indicate that
  /// statistics should not be written to thrift.
  std::optional<EncodedGeoStatistics> Encode() const;

  /// \brief Returns false if invalid WKB was encountered
  bool is_valid() const;

  /// \brief Reset existing statistics and populate them from previously-encoded ones
  void Decode(const EncodedGeoStatistics& encoded);

  /// \brief Minimum values in XYZM order
  ///
  /// For dimensions where dimension_valid() is false, the value will be NaN. For
  /// dimensions where dimension_empty() is true, the value will be +Inf.
  ///
  /// For the first dimension (X) only, wraparound bounds apply where xmin > xmax. In this
  /// case, these bounds represent the union of the intervals [xmax, Inf] and [-Inf,
  /// xmin]. This implementation does not yet generate these types of bounds but they may
  /// be encountered in statistics when reading a Parquet file.
  std::array<double, kMaxDimensions> lower_bound() const;

  /// \brief Maximum values in XYZM order
  ///
  /// For dimensions where dimension_valid() is false, the value will be NaN. For
  /// dimensions where dimension_empty() is true, the value will be -Inf.
  ///
  /// For the first dimension (X) only, wraparound bounds apply where xmin > xmax. In this
  /// case, these bounds represent the union of the intervals [xmax, Inf] and [-Inf,
  /// xmin]. This implementation does not yet generate these types of bounds but they may
  /// be encountered in statistics when reading a Parquet file.
  std::array<double, kMaxDimensions> upper_bound() const;

  /// \brief Dimension emptiness in XYZM order
  ///
  /// True for a given dimension if and only if zero non-NaN values were encountered
  /// in that dimension and dimension_valid() is true for that dimension.
  ///
  /// When calculating statistics, zero or more of these values may be true because
  /// this implementation calculates bounds for all dimensions; however, it may be
  /// true that zero coordinates were encountered in a given dimension. For example,
  /// dimension_empty() will return four true values if Update() was not called
  /// or if Update() was called with only null values. If Update() was provided
  /// one or more geometries with X and Y dimensions but not Z or M dimensions,
  /// dimension_empty() will return true, true, false, false.
  ///
  /// For statistics read from a Parquet file, dimension_empty() will always contain
  /// false values because there is no mechanism to communicate an empty interval
  /// in the Thrift metadata.
  std::array<bool, kMaxDimensions> dimension_empty() const;

  /// \brief Dimension validity (i.e. presence) in XYZM order
  ///
  /// When calculating statistics, this will always be true because this implementation
  /// calculates statistics for all dimensions. When reading a Parquet file, one or more
  /// of these values may be false because the file may not have provided bounds for all
  /// dimensions.
  ///
  /// See documentation for dimension_empty(), lower_bound(), and/or upper_bound() for the
  /// canonical values of those outputs for the dimensions where dimension_valid() is
  /// false.
  std::array<bool, kMaxDimensions> dimension_valid() const;

  /// \brief Return the geometry type codes
  ///
  /// This implementation always returns sorted output with no duplicates. When
  /// calculating statistics, a value will always be returned (although the returned
  /// vector may be empty if Update() was never called or was only called with null
  /// values). When reading a Parquet file, std::nullopt may be returned because
  /// the file may not have provided this information.
  std::optional<std::vector<int32_t>> geometry_types() const;

  /// \brief Return a string representation of these statistics
  std::string ToString() const;

 private:
  std::unique_ptr<GeoStatisticsImpl> impl_;
};

}  // namespace parquet::geospatial
