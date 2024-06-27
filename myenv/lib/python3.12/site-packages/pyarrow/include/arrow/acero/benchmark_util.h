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
#include <vector>

#include "benchmark/benchmark.h"

#include "arrow/acero/exec_plan.h"
#include "arrow/acero/test_util_internal.h"
#include "arrow/compute/exec.h"

namespace arrow {

namespace acero {

Status BenchmarkNodeOverhead(benchmark::State& state, int32_t num_batches,
                             int32_t batch_size, arrow::acero::BatchesWithSchema data,
                             std::vector<arrow::acero::Declaration>& node_declarations,
                             arrow::MemoryPool* pool = default_memory_pool());

Status BenchmarkIsolatedNodeOverhead(benchmark::State& state,
                                     arrow::compute::Expression expr, int32_t num_batches,
                                     int32_t batch_size,
                                     arrow::acero::BatchesWithSchema data,
                                     std::string factory_name,
                                     arrow::acero::ExecNodeOptions& options,
                                     arrow::MemoryPool* pool = default_memory_pool());

}  // namespace acero
}  // namespace arrow
