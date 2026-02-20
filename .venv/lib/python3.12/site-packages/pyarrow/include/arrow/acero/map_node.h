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
#include <functional>
#include <memory>
#include <vector>

#include "arrow/acero/exec_plan.h"
#include "arrow/acero/util.h"
#include "arrow/acero/visibility.h"
#include "arrow/compute/type_fwd.h"
#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/cancel.h"
#include "arrow/util/type_fwd.h"

namespace arrow {
namespace acero {

/// A utility base class for simple exec nodes with one input
///
/// Pause/Resume Producing are forwarded appropriately
/// There is nothing to do in StopProducingImpl
///
/// An AtomicCounter is used to keep track of when all data has arrived.  When it
/// has the Finish() method will be invoked
class ARROW_ACERO_EXPORT MapNode : public ExecNode, public TracedNode {
 public:
  MapNode(ExecPlan* plan, std::vector<ExecNode*> inputs,
          std::shared_ptr<Schema> output_schema);

  Status InputFinished(ExecNode* input, int total_batches) override;

  Status StartProducing() override;

  void PauseProducing(ExecNode* output, int32_t counter) override;

  void ResumeProducing(ExecNode* output, int32_t counter) override;

  Status InputReceived(ExecNode* input, ExecBatch batch) override;

  const Ordering& ordering() const override;

 protected:
  Status StopProducingImpl() override;

  /// Transform a batch
  ///
  /// The output batch will have the same guarantee as the input batch
  /// If this was the last batch this call may trigger Finish()
  virtual Result<ExecBatch> ProcessBatch(ExecBatch batch) = 0;

  /// Function called after all data has been received
  ///
  /// By default this does nothing.  Override this to provide a custom implementation.
  virtual void Finish();

 protected:
  // Counter for the number of batches received
  AtomicCounter input_counter_;
};

}  // namespace acero
}  // namespace arrow
