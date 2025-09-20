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

#include "arrow/python/common.h"
#include "arrow/python/visibility.h"

namespace arrow::py {

/// \brief Create an array of evenly spaced values within a given interval.
/// This function is similar to Python's `range` function.
/// The resulting array will contain values starting from `start` up to but not
/// including `stop`, with a step size of `step`. If `step` is zero, the function
/// will return an error.
/// The resulting array will have a data type of `int64`.
/// \param[in] start initial value of the sequence.
/// \param[in] stop final value of the sequence (exclusive).
/// \param[in] step step size between consecutive values.
/// \param[in] pool Memory pool for any memory allocations.
/// \return Result Array
ARROW_PYTHON_EXPORT
Result<std::shared_ptr<Array>> Arange(int64_t start, int64_t stop, int64_t step,
                                      MemoryPool* pool);

}  // namespace arrow::py
