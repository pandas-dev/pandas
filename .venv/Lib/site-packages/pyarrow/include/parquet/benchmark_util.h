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

#include <random>
#include <string>
#include <vector>

#include "parquet/types.h"

namespace parquet::benchmark {

template <typename T>
void GenerateBenchmarkData(uint32_t size, uint32_t seed, T* data,
                           std::vector<uint8_t>* heap, uint32_t data_string_length);

#define _GENERATE_BENCHMARK_DATA_DECL(KLASS)                            \
  template <>                                                           \
  void GenerateBenchmarkData(uint32_t size, uint32_t seed, KLASS* data, \
                             std::vector<uint8_t>* heap, uint32_t data_string_length);

_GENERATE_BENCHMARK_DATA_DECL(int32_t)
_GENERATE_BENCHMARK_DATA_DECL(int64_t)
_GENERATE_BENCHMARK_DATA_DECL(float)
_GENERATE_BENCHMARK_DATA_DECL(double)
_GENERATE_BENCHMARK_DATA_DECL(ByteArray)
_GENERATE_BENCHMARK_DATA_DECL(FLBA)
_GENERATE_BENCHMARK_DATA_DECL(Int96)

#undef _GENERATE_BENCHMARK_DATA_DECL

}  // namespace parquet::benchmark
