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

#include <memory>

#include "arrow/python/common.h"
#include "arrow/python/visibility.h"
#include "arrow/record_batch.h"
#include "arrow/result.h"
#include "arrow/util/macros.h"

namespace arrow {
namespace py {

class ARROW_PYTHON_EXPORT PyRecordBatchReader : public RecordBatchReader {
 public:
  std::shared_ptr<Schema> schema() const override;

  Status ReadNext(std::shared_ptr<RecordBatch>* batch) override;

  // For use from Cython
  // Assumes that `iterable` is borrowed
  static Result<std::shared_ptr<RecordBatchReader>> Make(std::shared_ptr<Schema>,
                                                         PyObject* iterable);

 protected:
  PyRecordBatchReader();

  Status Init(std::shared_ptr<Schema>, PyObject* iterable);

  std::shared_ptr<Schema> schema_;
  OwnedRefNoGIL iterator_;
};

class ARROW_PYTHON_EXPORT CastingRecordBatchReader : public RecordBatchReader {
 public:
  std::shared_ptr<Schema> schema() const override;

  Status ReadNext(std::shared_ptr<RecordBatch>* batch) override;

  static Result<std::shared_ptr<RecordBatchReader>> Make(
      std::shared_ptr<RecordBatchReader> parent, std::shared_ptr<Schema> schema);

  Status Close() override;

 protected:
  CastingRecordBatchReader();

  Status Init(std::shared_ptr<RecordBatchReader> parent, std::shared_ptr<Schema> schema);

  std::shared_ptr<RecordBatchReader> parent_;
  std::shared_ptr<Schema> schema_;
};

}  // namespace py
}  // namespace arrow
