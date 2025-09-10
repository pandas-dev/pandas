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

#include <cstdint>
#include "arrow/array/data.h"

namespace arrow {
namespace union_util {

/// \brief Compute the number of of logical nulls in a sparse union array
int64_t LogicalSparseUnionNullCount(const ArraySpan& span);

/// \brief Compute the number of of logical nulls in a dense union array
int64_t LogicalDenseUnionNullCount(const ArraySpan& span);

}  // namespace union_util
}  // namespace arrow
