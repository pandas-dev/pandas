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

#include "arrow/testing/visibility.h"
#include "arrow/type_fwd.h"

namespace arrow {

ARROW_TESTING_EXPORT
bool WithinUlp(util::Float16 left, util::Float16 right, int n_ulps);
ARROW_TESTING_EXPORT
bool WithinUlp(float left, float right, int n_ulps);
ARROW_TESTING_EXPORT
bool WithinUlp(double left, double right, int n_ulps);

ARROW_TESTING_EXPORT
void AssertWithinUlp(util::Float16 left, util::Float16 right, int n_ulps);
ARROW_TESTING_EXPORT
void AssertWithinUlp(float left, float right, int n_ulps);
ARROW_TESTING_EXPORT
void AssertWithinUlp(double left, double right, int n_ulps);

}  // namespace arrow
