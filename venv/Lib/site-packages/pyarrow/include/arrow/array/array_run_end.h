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

// Array accessor classes run-end encoded arrays

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "arrow/array/array_base.h"
#include "arrow/array/data.h"
#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/type.h"
#include "arrow/type_fwd.h"
#include "arrow/util/checked_cast.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

/// \addtogroup run-end-encoded-arrays
///
/// @{

// ----------------------------------------------------------------------
// RunEndEncoded

/// \brief Array type for run-end encoded data
class ARROW_EXPORT RunEndEncodedArray : public Array {
 private:
  std::shared_ptr<Array> run_ends_array_;
  std::shared_ptr<Array> values_array_;

 public:
  using TypeClass = RunEndEncodedType;

  explicit RunEndEncodedArray(const std::shared_ptr<ArrayData>& data);

  /// \brief Construct a RunEndEncodedArray from all parameters
  ///
  /// The length and offset parameters refer to the dimensions of the logical
  /// array which is the array we would get after expanding all the runs into
  /// repeated values. As such, length can be much greater than the length of
  /// the child run_ends and values arrays.
  RunEndEncodedArray(const std::shared_ptr<DataType>& type, int64_t length,
                     const std::shared_ptr<Array>& run_ends,
                     const std::shared_ptr<Array>& values, int64_t offset = 0);

  /// \brief Construct a RunEndEncodedArray from all parameters
  ///
  /// The length and offset parameters refer to the dimensions of the logical
  /// array which is the array we would get after expanding all the runs into
  /// repeated values. As such, length can be much greater than the length of
  /// the child run_ends and values arrays.
  static Result<std::shared_ptr<RunEndEncodedArray>> Make(
      const std::shared_ptr<DataType>& type, int64_t logical_length,
      const std::shared_ptr<Array>& run_ends, const std::shared_ptr<Array>& values,
      int64_t logical_offset = 0);

  /// \brief Construct a RunEndEncodedArray from values and run ends arrays
  ///
  /// The data type is automatically inferred from the arguments.
  /// The run_ends and values arrays must have the same length.
  static Result<std::shared_ptr<RunEndEncodedArray>> Make(
      int64_t logical_length, const std::shared_ptr<Array>& run_ends,
      const std::shared_ptr<Array>& values, int64_t logical_offset = 0);

 protected:
  void SetData(const std::shared_ptr<ArrayData>& data);

 public:
  /// \brief Returns an array holding the logical indexes of each run-end
  ///
  /// The physical offset to the array is applied.
  const std::shared_ptr<Array>& run_ends() const { return run_ends_array_; }

  /// \brief Returns an array holding the values of each run
  ///
  /// The physical offset to the array is applied.
  const std::shared_ptr<Array>& values() const { return values_array_; }

  /// \brief Returns an array holding the logical indexes of each run end
  ///
  /// If a non-zero logical offset is set, this function allocates a new
  /// array and rewrites all the run end values to be relative to the logical
  /// offset and cuts the end of the array to the logical length.
  Result<std::shared_ptr<Array>> LogicalRunEnds(MemoryPool* pool) const;

  /// \brief Returns an array holding the values of each run
  ///
  /// If a non-zero logical offset is set, this function allocates a new
  /// array containing only the values within the logical range.
  std::shared_ptr<Array> LogicalValues() const;

  /// \brief Find the physical offset of this REE array
  ///
  /// This function uses binary-search, so it has a O(log N) cost.
  int64_t FindPhysicalOffset() const;

  /// \brief Find the physical length of this REE array
  ///
  /// The physical length of an REE is the number of physical values (and
  /// run-ends) necessary to represent the logical range of values from offset
  /// to length.
  ///
  /// Avoid calling this function if the physical length can be established in
  /// some other way (e.g. when iterating over the runs sequentially until the
  /// end). This function uses binary-search, so it has a O(log N) cost.
  int64_t FindPhysicalLength() const;
};

/// @}

}  // namespace arrow
