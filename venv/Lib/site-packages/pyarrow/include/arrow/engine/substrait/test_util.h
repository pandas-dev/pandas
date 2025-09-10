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

#include "arrow/testing/gtest_util.h"
#include "arrow/util/vector.h"

#include <functional>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include "arrow/acero/exec_plan.h"
#include "arrow/compute/exec.h"
#include "arrow/compute/kernel.h"
#include "arrow/testing/visibility.h"
#include "arrow/util/async_generator.h"
#include "arrow/util/pcg_random.h"

namespace arrow {
namespace engine {

Result<std::shared_ptr<Table>> SortTableOnAllFields(const std::shared_ptr<Table>& tab);

void AssertTablesEqualIgnoringOrder(const std::shared_ptr<Table>& exp,
                                    const std::shared_ptr<Table>& act);

}  // namespace engine
}  // namespace arrow
