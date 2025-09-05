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

#include <vector>

#include "arrow/acero/options.h"
#include "arrow/acero/visibility.h"
#include "arrow/compute/exec.h"
#include "arrow/type.h"

namespace arrow {
namespace acero {
namespace asofjoin {

using AsofJoinKeys = AsofJoinNodeOptions::Keys;

/// \brief Make the output schema of an as-of-join node
///
/// \param[in] input_schema the schema of each input to the node
/// \param[in] input_keys the key of each input to the node
ARROW_ACERO_EXPORT Result<std::shared_ptr<Schema>> MakeOutputSchema(
    const std::vector<std::shared_ptr<Schema>>& input_schema,
    const std::vector<AsofJoinKeys>& input_keys);

}  // namespace asofjoin
}  // namespace acero
}  // namespace arrow
