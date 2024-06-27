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
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "arrow/array.h"
#include "arrow/array/builder_base.h"

namespace arrow {

/// \addtogroup run-end-encoded-builders
///
/// @{

namespace internal {

/// \brief An ArrayBuilder that deduplicates repeated values as they are
/// appended to the inner-ArrayBuilder and reports the length of the current run
/// of identical values.
///
/// The following sequence of calls
///
///     Append(2)
///     Append(2)
///     Append(2)
///     Append(7)
///     Append(7)
///     Append(2)
///     FinishInternal()
///
/// will cause the inner-builder to receive only 3 Append calls
///
///     Append(2)
///     Append(7)
///     Append(2)
///     FinishInternal()
///
/// Note that values returned by length(), null_count() and capacity() are
/// related to the compressed array built by the inner-ArrayBuilder.
class RunCompressorBuilder : public ArrayBuilder {
 public:
  RunCompressorBuilder(MemoryPool* pool, std::shared_ptr<ArrayBuilder> inner_builder,
                       std::shared_ptr<DataType> type);

  ~RunCompressorBuilder() override;

  ARROW_DISALLOW_COPY_AND_ASSIGN(RunCompressorBuilder);

  /// \brief Called right before a run is being closed
  ///
  /// Subclasses can override this function to perform an additional action when
  /// a run is closed (i.e. run-length is known and value is appended to the
  /// inner builder).
  ///
  /// \param value can be NULLPTR if closing a run of NULLs
  /// \param length the greater than 0 length of the value run being closed
  virtual Status WillCloseRun(const std::shared_ptr<const Scalar>& value,
                              int64_t length) {
    return Status::OK();
  }

  /// \brief Called right before a run of empty values is being closed
  ///
  /// Subclasses can override this function to perform an additional action when
  /// a run of empty values is appended (i.e. run-length is known and a single
  /// empty value is appended to the inner builder).
  ///
  /// \param length the greater than 0 length of the value run being closed
  virtual Status WillCloseRunOfEmptyValues(int64_t length) { return Status::OK(); }

  /// \brief Allocate enough memory for a given number of array elements.
  ///
  /// NOTE: Conservatively resizing a run-length compressed array for a given
  /// number of logical elements is not possible, since the physical length will
  /// vary depending on the values to be appended in the future. But we can
  /// pessimistically assume that each run will contain a single value and
  /// allocate that number of runs.
  Status Resize(int64_t capacity) override { return ResizePhysical(capacity); }

  /// \brief Allocate enough memory for a given number of runs.
  ///
  /// Like Resize on non-encoded builders, it does not account for variable size
  /// data.
  Status ResizePhysical(int64_t capacity);

  Status ReservePhysical(int64_t additional_capacity) {
    return Reserve(additional_capacity);
  }

  void Reset() override;

  Status AppendNull() final { return AppendNulls(1); }
  Status AppendNulls(int64_t length) override;

  Status AppendEmptyValue() final { return AppendEmptyValues(1); }
  Status AppendEmptyValues(int64_t length) override;

  Status AppendScalar(const Scalar& scalar, int64_t n_repeats) override;
  Status AppendScalars(const ScalarVector& scalars) override;

  // AppendArraySlice() is not implemented.

  /// \brief Append a slice of an array containing values from already
  /// compressed runs.
  ///
  /// NOTE: WillCloseRun() is not called as the length of each run cannot be
  /// determined at this point. Caller should ensure that !has_open_run() by
  /// calling FinishCurrentRun() before calling this.
  ///
  /// Pre-condition: !has_open_run()
  Status AppendRunCompressedArraySlice(const ArraySpan& array, int64_t offset,
                                       int64_t length);

  /// \brief Forces the closing of the current run if one is currently open.
  ///
  /// This can be called when one wants to ensure the current run will not be
  /// extended. This may cause identical values to appear close to each other in
  /// the underlying array (i.e. two runs that could be a single run) if more
  /// values are appended after this is called.
  ///
  /// Finish() and FinishInternal() call this automatically.
  virtual Status FinishCurrentRun();

  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  ArrayBuilder& inner_builder() const { return *inner_builder_; }

  std::shared_ptr<DataType> type() const override { return inner_builder_->type(); }

  bool has_open_run() const { return current_run_length_ > 0; }
  int64_t open_run_length() const { return current_run_length_; }

 private:
  inline void UpdateDimensions() {
    capacity_ = inner_builder_->capacity();
    length_ = inner_builder_->length();
    null_count_ = inner_builder_->null_count();
  }

 private:
  std::shared_ptr<ArrayBuilder> inner_builder_;
  std::shared_ptr<const Scalar> current_value_ = NULLPTR;
  int64_t current_run_length_ = 0;
};

}  // namespace internal

// ----------------------------------------------------------------------
// RunEndEncoded builder

/// \brief Run-end encoded array builder.
///
/// NOTE: the value returned by and capacity() is related to the
/// compressed array (physical) and not the decoded array (logical) that is
/// run-end encoded. null_count() always returns 0. length(), on the other hand,
/// returns the logical length of the run-end encoded array.
class ARROW_EXPORT RunEndEncodedBuilder : public ArrayBuilder {
 private:
  // An internal::RunCompressorBuilder that produces a run-end in the
  // RunEndEncodedBuilder every time a value-run is closed.
  class ValueRunBuilder : public internal::RunCompressorBuilder {
   public:
    ValueRunBuilder(MemoryPool* pool, const std::shared_ptr<ArrayBuilder>& value_builder,
                    const std::shared_ptr<DataType>& value_type,
                    RunEndEncodedBuilder& ree_builder);

    ~ValueRunBuilder() override = default;

    Status WillCloseRun(const std::shared_ptr<const Scalar>&, int64_t length) override {
      return ree_builder_.CloseRun(length);
    }

    Status WillCloseRunOfEmptyValues(int64_t length) override {
      return ree_builder_.CloseRun(length);
    }

   private:
    RunEndEncodedBuilder& ree_builder_;
  };

 public:
  RunEndEncodedBuilder(MemoryPool* pool,
                       const std::shared_ptr<ArrayBuilder>& run_end_builder,
                       const std::shared_ptr<ArrayBuilder>& value_builder,
                       std::shared_ptr<DataType> type);

  /// \brief Allocate enough memory for a given number of array elements.
  ///
  /// NOTE: Conservatively resizing an REE for a given number of logical
  /// elements is not possible, since the physical length will vary depending on
  /// the values to be appended in the future. But we can pessimistically assume
  /// that each run will contain a single value and allocate that number of
  /// runs.
  Status Resize(int64_t capacity) override { return ResizePhysical(capacity); }

  /// \brief Allocate enough memory for a given number of runs.
  Status ResizePhysical(int64_t capacity);

  /// \brief Ensure that there is enough space allocated to append the indicated
  /// number of run without any further reallocation. Overallocation is
  /// used in order to minimize the impact of incremental ReservePhysical() calls.
  /// Note that additional_capacity is relative to the current number of elements
  /// rather than to the current capacity, so calls to Reserve() which are not
  /// interspersed with addition of new elements may not increase the capacity.
  ///
  /// \param[in] additional_capacity the number of additional runs
  /// \return Status
  Status ReservePhysical(int64_t additional_capacity) {
    return Reserve(additional_capacity);
  }

  void Reset() override;

  Status AppendNull() final { return AppendNulls(1); }
  Status AppendNulls(int64_t length) override;

  Status AppendEmptyValue() final { return AppendEmptyValues(1); }
  Status AppendEmptyValues(int64_t length) override;
  Status AppendScalar(const Scalar& scalar, int64_t n_repeats) override;
  Status AppendScalars(const ScalarVector& scalars) override;
  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override;
  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<RunEndEncodedArray>* out) { return FinishTyped(out); }

  /// \brief Forces the closing of the current run if one is currently open.
  ///
  /// This can be called when one wants to ensure the current run will not be
  /// extended. This may cause identical values to appear close to each other in
  /// the values array (i.e. two runs that could be a single run) if more
  /// values are appended after this is called.
  Status FinishCurrentRun();

  std::shared_ptr<DataType> type() const override;

 private:
  /// \brief Update physical capacity and logical length
  ///
  /// \param committed_logical_length number of logical values that have been
  ///                                 committed to the values array
  /// \param open_run_length number of logical values in the currently open run if any
  inline void UpdateDimensions(int64_t committed_logical_length,
                               int64_t open_run_length) {
    capacity_ = run_end_builder().capacity();
    length_ = committed_logical_length + open_run_length;
    committed_logical_length_ = committed_logical_length;
  }

  // Pre-condition: !value_run_builder_.has_open_run()
  template <typename RunEndCType>
  Status DoAppendArraySlice(const ArraySpan& array, int64_t offset, int64_t length);

  template <typename RunEndCType>
  Status DoAppendRunEnd(int64_t run_end);

  /// \brief Cast run_end to the appropriate type and appends it to the run_ends
  /// array.
  Status AppendRunEnd(int64_t run_end);

  /// \brief Close a run by appending a value to the run_ends array and updating
  /// length_ to reflect the new run.
  ///
  /// Pre-condition: run_length > 0.
  [[nodiscard]] Status CloseRun(int64_t run_length);

  ArrayBuilder& run_end_builder();
  ArrayBuilder& value_builder();

 private:
  std::shared_ptr<RunEndEncodedType> type_;
  ValueRunBuilder* value_run_builder_;
  // The length not counting the current open run in the value_run_builder_
  int64_t committed_logical_length_ = 0;
};

/// @}

}  // namespace arrow
