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

#include "ipc.h"

#include <memory>

#include "arrow/compute/cast.h"
#include "arrow/python/pyarrow.h"

namespace arrow {
namespace py {

PyRecordBatchReader::PyRecordBatchReader() {}

Status PyRecordBatchReader::Init(std::shared_ptr<Schema> schema, PyObject* iterable) {
  schema_ = std::move(schema);

  iterator_.reset(PyObject_GetIter(iterable));
  return CheckPyError();
}

std::shared_ptr<Schema> PyRecordBatchReader::schema() const { return schema_; }

Status PyRecordBatchReader::ReadNext(std::shared_ptr<RecordBatch>* batch) {
  PyAcquireGIL lock;

  if (!iterator_) {
    // End of stream
    batch->reset();
    return Status::OK();
  }

  OwnedRef py_batch(PyIter_Next(iterator_.obj()));
  if (!py_batch) {
    RETURN_IF_PYERROR();
    // End of stream
    batch->reset();
    iterator_.reset();
    return Status::OK();
  }

  return unwrap_batch(py_batch.obj()).Value(batch);
}

Result<std::shared_ptr<RecordBatchReader>> PyRecordBatchReader::Make(
    std::shared_ptr<Schema> schema, PyObject* iterable) {
  auto reader = std::shared_ptr<PyRecordBatchReader>(new PyRecordBatchReader());
  RETURN_NOT_OK(reader->Init(std::move(schema), iterable));
  return reader;
}

CastingRecordBatchReader::CastingRecordBatchReader() = default;

Status CastingRecordBatchReader::Init(std::shared_ptr<RecordBatchReader> parent,
                                      std::shared_ptr<Schema> schema) {
  std::shared_ptr<Schema> src = parent->schema();

  // The check for names has already been done in Python where it's easier to
  // generate a nice error message.
  int num_fields = schema->num_fields();
  if (src->num_fields() != num_fields) {
    return Status::Invalid("Number of fields not equal");
  }

  // Ensure all columns can be cast before succeeding
  for (int i = 0; i < num_fields; i++) {
    if (!compute::CanCast(*src->field(i)->type(), *schema->field(i)->type())) {
      return Status::TypeError("Field ", i, " cannot be cast from ",
                               src->field(i)->type()->ToString(), " to ",
                               schema->field(i)->type()->ToString());
    }
  }

  parent_ = std::move(parent);
  schema_ = std::move(schema);

  return Status::OK();
}

std::shared_ptr<Schema> CastingRecordBatchReader::schema() const { return schema_; }

Status CastingRecordBatchReader::ReadNext(std::shared_ptr<RecordBatch>* batch) {
  std::shared_ptr<RecordBatch> out;
  ARROW_RETURN_NOT_OK(parent_->ReadNext(&out));
  if (!out) {
    batch->reset();
    return Status::OK();
  }

  auto num_columns = out->num_columns();
  auto options = compute::CastOptions::Safe();
  ArrayVector columns(num_columns);
  for (int i = 0; i < num_columns; i++) {
    const Array& src = *out->column(i);
    if (!schema_->field(i)->nullable() && src.null_count() > 0) {
      return Status::Invalid(
          "Can't cast array that contains nulls to non-nullable field at index ", i);
    }

    ARROW_ASSIGN_OR_RAISE(columns[i],
                          compute::Cast(src, schema_->field(i)->type(), options));
  }

  *batch = RecordBatch::Make(schema_, out->num_rows(), std::move(columns));
  return Status::OK();
}

Result<std::shared_ptr<RecordBatchReader>> CastingRecordBatchReader::Make(
    std::shared_ptr<RecordBatchReader> parent, std::shared_ptr<Schema> schema) {
  auto reader = std::shared_ptr<CastingRecordBatchReader>(new CastingRecordBatchReader());
  ARROW_RETURN_NOT_OK(reader->Init(parent, schema));
  return reader;
}

Status CastingRecordBatchReader::Close() { return parent_->Close(); }

}  // namespace py
}  // namespace arrow
