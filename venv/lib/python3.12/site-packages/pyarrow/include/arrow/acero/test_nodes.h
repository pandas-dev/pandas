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

#include <string>

#include "arrow/acero/options.h"
#include "arrow/acero/test_util_internal.h"
#include "arrow/testing/random.h"

namespace arrow {
namespace acero {

// \brief Make a delaying source that is optionally noisy (prints when it emits)
AsyncGenerator<std::optional<ExecBatch>> MakeDelayedGen(
    Iterator<std::optional<ExecBatch>> src, std::string label, double delay_sec,
    bool noisy = false);

// \brief Make a delaying source that is optionally noisy (prints when it emits)
AsyncGenerator<std::optional<ExecBatch>> MakeDelayedGen(
    AsyncGenerator<std::optional<ExecBatch>> src, std::string label, double delay_sec,
    bool noisy = false);

// \brief Make a delaying source that is optionally noisy (prints when it emits)
AsyncGenerator<std::optional<ExecBatch>> MakeDelayedGen(BatchesWithSchema src,
                                                        std::string label,
                                                        double delay_sec,
                                                        bool noisy = false);

/// A node that slightly resequences the input at random
struct JitterNodeOptions : public ExecNodeOptions {
  random::SeedType seed;
  /// The max amount to add to a node's "cost".
  int max_jitter_modifier;

  explicit JitterNodeOptions(random::SeedType seed, int max_jitter_modifier = 5)
      : seed(seed), max_jitter_modifier(max_jitter_modifier) {}
  static constexpr std::string_view kName = "jitter";
};

class GateImpl;

class Gate {
 public:
  static std::shared_ptr<Gate> Make();

  Gate();
  virtual ~Gate();

  void ReleaseAllBatches();
  void ReleaseOneBatch();
  Future<> WaitForNextReleasedBatch();

 private:
  ARROW_DISALLOW_COPY_AND_ASSIGN(Gate);

  GateImpl* impl_;
};

// A node that holds all input batches until a given gate is released
struct GatedNodeOptions : public ExecNodeOptions {
  explicit GatedNodeOptions(Gate* gate) : gate(gate) {}
  Gate* gate;

  static constexpr std::string_view kName = "gated";
};

void RegisterTestNodes();

}  // namespace acero
}  // namespace arrow
