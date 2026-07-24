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

#include "arrow/python/util.h"

#include "arrow/array.h"
#include "arrow/python/common.h"

namespace arrow ::py {

Result<std::shared_ptr<Array>> Arange(int64_t start, int64_t stop, int64_t step,
                                      MemoryPool* pool) {
  int64_t size;
  if (step == 0) {
    return Status::Invalid("Step must not be zero");
  }
  if (step > 0 && stop > start) {
    // Ceiling division for positive step
    size = (stop - start + step - 1) / step;
  } else if (step < 0 && stop < start) {
    // Ceiling division for negative step
    size = (start - stop - step - 1) / (-step);
  } else {
    return MakeEmptyArray(int64());
  }
  std::shared_ptr<Buffer> data_buffer;
  ARROW_ASSIGN_OR_RAISE(data_buffer, AllocateBuffer(size * sizeof(int64_t), pool));
  auto values = reinterpret_cast<int64_t*>(data_buffer->mutable_data());
  for (int64_t i = 0; i < size; ++i) {
    values[i] = start + i * step;
  }
  auto data = ArrayData::Make(int64(), size, {nullptr, data_buffer}, 0);
  return MakeArray(data);
}

}  // namespace arrow::py
