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
#include "arrow/acero/exec_plan.h"
#include "arrow/acero/options.h"

#include <memory>

namespace arrow::acero {

class BackpressureHandler {
 private:
  BackpressureHandler(ExecNode* input, size_t low_threshold, size_t high_threshold,
                      std::unique_ptr<BackpressureControl> backpressure_control)
      : input_(input),
        low_threshold_(low_threshold),
        high_threshold_(high_threshold),
        backpressure_control_(std::move(backpressure_control)) {}

 public:
  static Result<BackpressureHandler> Make(
      ExecNode* input, size_t low_threshold, size_t high_threshold,
      std::unique_ptr<BackpressureControl> backpressure_control) {
    if (low_threshold >= high_threshold) {
      return Status::Invalid("low threshold (", low_threshold,
                             ") must be less than high threshold (", high_threshold, ")");
    }
    if (backpressure_control == NULLPTR) {
      return Status::Invalid("null backpressure control parameter");
    }
    BackpressureHandler backpressure_handler(input, low_threshold, high_threshold,
                                             std::move(backpressure_control));
    return backpressure_handler;
  }

  void Handle(size_t start_level, size_t end_level) {
    if (start_level < high_threshold_ && end_level >= high_threshold_) {
      backpressure_control_->Pause();
    } else if (start_level > low_threshold_ && end_level <= low_threshold_) {
      backpressure_control_->Resume();
    }
  }

  Status ForceShutdown() {
    // It may be unintuitive to call Resume() here, but this is to avoid a deadlock.
    // Since acero's executor won't terminate if any one node is paused, we need to
    // force resume the node before stopping production.
    backpressure_control_->Resume();
    return input_->StopProducing();
  }

 private:
  ExecNode* input_;
  size_t low_threshold_;
  size_t high_threshold_;
  std::unique_ptr<BackpressureControl> backpressure_control_;
};

}  // namespace arrow::acero
