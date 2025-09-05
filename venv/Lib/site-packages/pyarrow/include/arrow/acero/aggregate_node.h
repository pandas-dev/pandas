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

// This API is EXPERIMENTAL.

#pragma once

#include <memory>
#include <vector>

#include "arrow/acero/visibility.h"
#include "arrow/compute/api_aggregate.h"
#include "arrow/compute/test_util_internal.h"
#include "arrow/compute/type_fwd.h"
#include "arrow/result.h"
#include "arrow/type_fwd.h"

namespace arrow {
namespace acero {
namespace aggregate {

using compute::Aggregate;
using compute::default_exec_context;
using compute::ExecContext;

/// \brief Make the output schema of an aggregate node
///
/// The output schema is determined by the aggregation kernels, which may depend on the
/// ExecContext argument. To guarantee correct results, the same ExecContext argument
/// should be used in execution.
///
/// \param[in] input_schema the schema of the input to the node
/// \param[in] keys the grouping keys for the aggregation
/// \param[in] segment_keys the segmenting keys for the aggregation
/// \param[in] aggregates the aggregates for the aggregation
/// \param[in] exec_ctx the execution context for the aggregation
ARROW_ACERO_EXPORT Result<std::shared_ptr<Schema>> MakeOutputSchema(
    const std::shared_ptr<Schema>& input_schema, const std::vector<FieldRef>& keys,
    const std::vector<FieldRef>& segment_keys, const std::vector<Aggregate>& aggregates,
    ExecContext* exec_ctx = default_exec_context());

}  // namespace aggregate
}  // namespace acero
}  // namespace arrow
