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

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include "arrow/array/array_base.h"
#include "arrow/array/array_binary.h"
#include "arrow/array/builder_base.h"
#include "arrow/array/data.h"
#include "arrow/buffer.h"
#include "arrow/buffer_builder.h"
#include "arrow/status.h"
#include "arrow/type.h"
#include "arrow/util/binary_view_util.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

/// \addtogroup binary-builders
///
/// @{

// ----------------------------------------------------------------------
// Binary and String

template <typename TYPE>
class BaseBinaryBuilder
    : public ArrayBuilder,
      public internal::ArrayBuilderExtraOps<BaseBinaryBuilder<TYPE>, std::string_view> {
 public:
  using TypeClass = TYPE;
  using offset_type = typename TypeClass::offset_type;

  explicit BaseBinaryBuilder(MemoryPool* pool = default_memory_pool(),
                             int64_t alignment = kDefaultBufferAlignment)
      : ArrayBuilder(pool, alignment),
        offsets_builder_(pool, alignment),
        value_data_builder_(pool, alignment) {}

  BaseBinaryBuilder(const std::shared_ptr<DataType>& type, MemoryPool* pool)
      : BaseBinaryBuilder(pool) {}

  Status Append(const uint8_t* value, offset_type length) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppendNextOffset();
    // Safety check for UBSAN.
    if (ARROW_PREDICT_TRUE(length > 0)) {
      ARROW_RETURN_NOT_OK(ValidateOverflow(length));
      ARROW_RETURN_NOT_OK(value_data_builder_.Append(value, length));
    }

    UnsafeAppendToBitmap(true);
    return Status::OK();
  }

  Status Append(const char* value, offset_type length) {
    return Append(reinterpret_cast<const uint8_t*>(value), length);
  }

  Status Append(std::string_view value) {
    return Append(value.data(), static_cast<offset_type>(value.size()));
  }

  /// Extend the last appended value by appending more data at the end
  ///
  /// Unlike Append, this does not create a new offset.
  Status ExtendCurrent(const uint8_t* value, offset_type length) {
    // Safety check for UBSAN.
    if (ARROW_PREDICT_TRUE(length > 0)) {
      ARROW_RETURN_NOT_OK(ValidateOverflow(length));
      ARROW_RETURN_NOT_OK(value_data_builder_.Append(value, length));
    }
    return Status::OK();
  }

  Status ExtendCurrent(std::string_view value) {
    return ExtendCurrent(reinterpret_cast<const uint8_t*>(value.data()),
                         static_cast<offset_type>(value.size()));
  }

  Status AppendNulls(int64_t length) final {
    const int64_t num_bytes = value_data_builder_.length();
    ARROW_RETURN_NOT_OK(Reserve(length));
    for (int64_t i = 0; i < length; ++i) {
      offsets_builder_.UnsafeAppend(static_cast<offset_type>(num_bytes));
    }
    UnsafeAppendToBitmap(length, false);
    return Status::OK();
  }

  Status AppendNull() final {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppendNextOffset();
    UnsafeAppendToBitmap(false);
    return Status::OK();
  }

  Status AppendEmptyValue() final {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppendNextOffset();
    UnsafeAppendToBitmap(true);
    return Status::OK();
  }

  Status AppendEmptyValues(int64_t length) final {
    const int64_t num_bytes = value_data_builder_.length();
    ARROW_RETURN_NOT_OK(Reserve(length));
    for (int64_t i = 0; i < length; ++i) {
      offsets_builder_.UnsafeAppend(static_cast<offset_type>(num_bytes));
    }
    UnsafeAppendToBitmap(length, true);
    return Status::OK();
  }

  /// \brief Append without checking capacity
  ///
  /// Offsets and data should have been presized using Reserve() and
  /// ReserveData(), respectively.
  void UnsafeAppend(const uint8_t* value, offset_type length) {
    UnsafeAppendNextOffset();
    value_data_builder_.UnsafeAppend(value, length);
    UnsafeAppendToBitmap(true);
  }

  void UnsafeAppend(const char* value, offset_type length) {
    UnsafeAppend(reinterpret_cast<const uint8_t*>(value), length);
  }

  void UnsafeAppend(const std::string& value) {
    UnsafeAppend(value.c_str(), static_cast<offset_type>(value.size()));
  }

  void UnsafeAppend(std::string_view value) {
    UnsafeAppend(value.data(), static_cast<offset_type>(value.size()));
  }

  /// Like ExtendCurrent, but do not check capacity
  void UnsafeExtendCurrent(const uint8_t* value, offset_type length) {
    value_data_builder_.UnsafeAppend(value, length);
  }

  void UnsafeExtendCurrent(std::string_view value) {
    UnsafeExtendCurrent(reinterpret_cast<const uint8_t*>(value.data()),
                        static_cast<offset_type>(value.size()));
  }

  void UnsafeAppendNull() {
    const int64_t num_bytes = value_data_builder_.length();
    offsets_builder_.UnsafeAppend(static_cast<offset_type>(num_bytes));
    UnsafeAppendToBitmap(false);
  }

  void UnsafeAppendEmptyValue() {
    const int64_t num_bytes = value_data_builder_.length();
    offsets_builder_.UnsafeAppend(static_cast<offset_type>(num_bytes));
    UnsafeAppendToBitmap(true);
  }

  /// \brief Append a sequence of strings in one shot.
  ///
  /// \param[in] values a vector of strings
  /// \param[in] valid_bytes an optional sequence of bytes where non-zero
  /// indicates a valid (non-null) value
  /// \return Status
  Status AppendValues(const std::vector<std::string>& values,
                      const uint8_t* valid_bytes = NULLPTR) {
    std::size_t total_length = std::accumulate(
        values.begin(), values.end(), 0ULL,
        [](uint64_t sum, const std::string& str) { return sum + str.size(); });
    ARROW_RETURN_NOT_OK(Reserve(values.size()));
    ARROW_RETURN_NOT_OK(ReserveData(total_length));

    if (valid_bytes != NULLPTR) {
      for (std::size_t i = 0; i < values.size(); ++i) {
        UnsafeAppendNextOffset();
        if (valid_bytes[i]) {
          value_data_builder_.UnsafeAppend(
              reinterpret_cast<const uint8_t*>(values[i].data()), values[i].size());
        }
      }
    } else {
      for (const auto& value : values) {
        UnsafeAppendNextOffset();
        value_data_builder_.UnsafeAppend(reinterpret_cast<const uint8_t*>(value.data()),
                                         value.size());
      }
    }

    UnsafeAppendToBitmap(valid_bytes, values.size());
    return Status::OK();
  }

  /// \brief Append a sequence of nul-terminated strings in one shot.
  ///        If one of the values is NULL, it is processed as a null
  ///        value even if the corresponding valid_bytes entry is 1.
  ///
  /// \param[in] values a contiguous C array of nul-terminated char *
  /// \param[in] length the number of values to append
  /// \param[in] valid_bytes an optional sequence of bytes where non-zero
  /// indicates a valid (non-null) value
  /// \return Status
  Status AppendValues(const char** values, int64_t length,
                      const uint8_t* valid_bytes = NULLPTR) {
    std::size_t total_length = 0;
    std::vector<std::size_t> value_lengths(length);
    bool have_null_value = false;
    for (int64_t i = 0; i < length; ++i) {
      if (values[i] != NULLPTR) {
        auto value_length = strlen(values[i]);
        value_lengths[i] = value_length;
        total_length += value_length;
      } else {
        have_null_value = true;
      }
    }
    ARROW_RETURN_NOT_OK(Reserve(length));
    ARROW_RETURN_NOT_OK(ReserveData(total_length));

    if (valid_bytes) {
      int64_t valid_bytes_offset = 0;
      for (int64_t i = 0; i < length; ++i) {
        UnsafeAppendNextOffset();
        if (valid_bytes[i]) {
          if (values[i]) {
            value_data_builder_.UnsafeAppend(reinterpret_cast<const uint8_t*>(values[i]),
                                             value_lengths[i]);
          } else {
            UnsafeAppendToBitmap(valid_bytes + valid_bytes_offset,
                                 i - valid_bytes_offset);
            UnsafeAppendToBitmap(false);
            valid_bytes_offset = i + 1;
          }
        }
      }
      UnsafeAppendToBitmap(valid_bytes + valid_bytes_offset, length - valid_bytes_offset);
    } else {
      if (have_null_value) {
        std::vector<uint8_t> valid_vector(length, 0);
        for (int64_t i = 0; i < length; ++i) {
          UnsafeAppendNextOffset();
          if (values[i]) {
            value_data_builder_.UnsafeAppend(reinterpret_cast<const uint8_t*>(values[i]),
                                             value_lengths[i]);
            valid_vector[i] = 1;
          }
        }
        UnsafeAppendToBitmap(valid_vector.data(), length);
      } else {
        for (int64_t i = 0; i < length; ++i) {
          UnsafeAppendNextOffset();
          value_data_builder_.UnsafeAppend(reinterpret_cast<const uint8_t*>(values[i]),
                                           value_lengths[i]);
        }
        UnsafeAppendToBitmap(NULLPTR, length);
      }
    }
    return Status::OK();
  }

  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override {
    auto bitmap = array.GetValues<uint8_t>(0, 0);
    auto offsets = array.GetValues<offset_type>(1);
    auto data = array.GetValues<uint8_t>(2, 0);
    auto total_length = offsets[offset + length] - offsets[offset];
    ARROW_RETURN_NOT_OK(Reserve(length));
    ARROW_RETURN_NOT_OK(ReserveData(total_length));
    for (int64_t i = 0; i < length; i++) {
      if (!bitmap || bit_util::GetBit(bitmap, array.offset + offset + i)) {
        const offset_type start = offsets[offset + i];
        const offset_type end = offsets[offset + i + 1];
        UnsafeAppend(data + start, end - start);
      } else {
        UnsafeAppendNull();
      }
    }
    return Status::OK();
  }

  void Reset() override {
    ArrayBuilder::Reset();
    offsets_builder_.Reset();
    value_data_builder_.Reset();
  }

  Status ValidateOverflow(int64_t new_bytes) {
    auto new_size = value_data_builder_.length() + new_bytes;
    if (ARROW_PREDICT_FALSE(new_size > memory_limit())) {
      return Status::CapacityError("array cannot contain more than ", memory_limit(),
                                   " bytes, have ", new_size);
    } else {
      return Status::OK();
    }
  }

  Status Resize(int64_t capacity) override {
    ARROW_RETURN_NOT_OK(CheckCapacity(capacity));
    // One more than requested for offsets
    ARROW_RETURN_NOT_OK(offsets_builder_.Resize(capacity + 1));
    return ArrayBuilder::Resize(capacity);
  }

  /// \brief Ensures there is enough allocated capacity to append the indicated
  /// number of bytes to the value data buffer without additional allocations
  Status ReserveData(int64_t elements) {
    ARROW_RETURN_NOT_OK(ValidateOverflow(elements));
    return value_data_builder_.Reserve(elements);
  }

  Status FinishInternal(std::shared_ptr<ArrayData>* out) override {
    // Write final offset (values length)
    ARROW_RETURN_NOT_OK(AppendNextOffset());

    // These buffers' padding zeroed by BufferBuilder
    std::shared_ptr<Buffer> offsets, value_data, null_bitmap;
    ARROW_RETURN_NOT_OK(offsets_builder_.Finish(&offsets));
    ARROW_RETURN_NOT_OK(value_data_builder_.Finish(&value_data));
    ARROW_RETURN_NOT_OK(null_bitmap_builder_.Finish(&null_bitmap));

    *out = ArrayData::Make(type(), length_, {null_bitmap, offsets, value_data},
                           null_count_, 0);
    Reset();
    return Status::OK();
  }

  /// \return data pointer of the value date builder
  const uint8_t* value_data() const { return value_data_builder_.data(); }
  /// \return size of values buffer so far
  int64_t value_data_length() const { return value_data_builder_.length(); }
  /// \return capacity of values buffer
  int64_t value_data_capacity() const { return value_data_builder_.capacity(); }

  /// \return data pointer of the value date builder
  const offset_type* offsets_data() const { return offsets_builder_.data(); }

  /// Temporary access to a value.
  ///
  /// This pointer becomes invalid on the next modifying operation.
  const uint8_t* GetValue(int64_t i, offset_type* out_length) const {
    const offset_type* offsets = offsets_builder_.data();
    const auto offset = offsets[i];
    if (i == (length_ - 1)) {
      *out_length = static_cast<offset_type>(value_data_builder_.length()) - offset;
    } else {
      *out_length = offsets[i + 1] - offset;
    }
    return value_data_builder_.data() + offset;
  }

  offset_type offset(int64_t i) const { return offsets_data()[i]; }

  /// Temporary access to a value.
  ///
  /// This view becomes invalid on the next modifying operation.
  std::string_view GetView(int64_t i) const {
    offset_type value_length;
    const uint8_t* value_data = GetValue(i, &value_length);
    return std::string_view(reinterpret_cast<const char*>(value_data), value_length);
  }

  // Cannot make this a static attribute because of linking issues
  static constexpr int64_t memory_limit() {
    return std::numeric_limits<offset_type>::max() - 1;
  }

 protected:
  TypedBufferBuilder<offset_type> offsets_builder_;
  TypedBufferBuilder<uint8_t> value_data_builder_;

  Status AppendNextOffset() {
    const int64_t num_bytes = value_data_builder_.length();
    return offsets_builder_.Append(static_cast<offset_type>(num_bytes));
  }

  void UnsafeAppendNextOffset() {
    const int64_t num_bytes = value_data_builder_.length();
    offsets_builder_.UnsafeAppend(static_cast<offset_type>(num_bytes));
  }
};

/// \class BinaryBuilder
/// \brief Builder class for variable-length binary data
class ARROW_EXPORT BinaryBuilder : public BaseBinaryBuilder<BinaryType> {
 public:
  using BaseBinaryBuilder::BaseBinaryBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<BinaryArray>* out) { return FinishTyped(out); }

  std::shared_ptr<DataType> type() const override { return binary(); }
};

/// \class StringBuilder
/// \brief Builder class for UTF8 strings
class ARROW_EXPORT StringBuilder : public BinaryBuilder {
 public:
  using BinaryBuilder::BinaryBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<StringArray>* out) { return FinishTyped(out); }

  std::shared_ptr<DataType> type() const override { return utf8(); }
};

/// \class LargeBinaryBuilder
/// \brief Builder class for large variable-length binary data
class ARROW_EXPORT LargeBinaryBuilder : public BaseBinaryBuilder<LargeBinaryType> {
 public:
  using BaseBinaryBuilder::BaseBinaryBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<LargeBinaryArray>* out) { return FinishTyped(out); }

  std::shared_ptr<DataType> type() const override { return large_binary(); }
};

/// \class LargeStringBuilder
/// \brief Builder class for large UTF8 strings
class ARROW_EXPORT LargeStringBuilder : public LargeBinaryBuilder {
 public:
  using LargeBinaryBuilder::LargeBinaryBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<LargeStringArray>* out) { return FinishTyped(out); }

  std::shared_ptr<DataType> type() const override { return large_utf8(); }
};

// ----------------------------------------------------------------------
// BinaryViewBuilder, StringViewBuilder
//
// These builders do not support building raw pointer view arrays.

namespace internal {

// We allocate medium-sized memory chunks and accumulate data in those, which
// may result in some waste if there are many large-ish strings. If a string
// comes along that does not fit into a block, we allocate a new block and
// write into that.
//
// Later we can implement optimizations to continuing filling underfull blocks
// after encountering a large string that required allocating a new block.
class ARROW_EXPORT StringHeapBuilder {
 public:
  static constexpr int64_t kDefaultBlocksize = 32 << 10;  // 32KB

  StringHeapBuilder(MemoryPool* pool, int64_t alignment)
      : pool_(pool), alignment_(alignment) {}

  void SetBlockSize(int64_t blocksize) { blocksize_ = blocksize; }

  using c_type = BinaryViewType::c_type;

  template <bool Safe>
  std::conditional_t<Safe, Result<c_type>, c_type> Append(const uint8_t* value,
                                                          int64_t length) {
    if (length <= BinaryViewType::kInlineSize) {
      return util::ToInlineBinaryView(value, static_cast<int32_t>(length));
    }

    if constexpr (Safe) {
      ARROW_RETURN_NOT_OK(Reserve(length));
    }

    auto v =
        util::ToBinaryView(value, static_cast<int32_t>(length),
                           static_cast<int32_t>(blocks_.size() - 1), current_offset_);

    memcpy(current_out_buffer_, value, static_cast<size_t>(length));
    current_out_buffer_ += length;
    current_remaining_bytes_ -= length;
    current_offset_ += static_cast<int32_t>(length);
    return v;
  }

  static constexpr int64_t ValueSizeLimit() {
    return std::numeric_limits<int32_t>::max();
  }

  /// \brief Ensure that the indicated number of bytes can be appended via
  /// UnsafeAppend operations without the need to allocate more memory
  Status Reserve(int64_t num_bytes) {
    if (ARROW_PREDICT_FALSE(num_bytes > ValueSizeLimit())) {
      return Status::CapacityError(
          "BinaryView or StringView elements cannot reference "
          "strings larger than 2GB");
    }
    if (num_bytes > current_remaining_bytes_) {
      ARROW_RETURN_NOT_OK(FinishLastBlock());
      current_remaining_bytes_ = num_bytes > blocksize_ ? num_bytes : blocksize_;
      ARROW_ASSIGN_OR_RAISE(
          std::shared_ptr<ResizableBuffer> new_block,
          AllocateResizableBuffer(current_remaining_bytes_, alignment_, pool_));
      current_offset_ = 0;
      current_out_buffer_ = new_block->mutable_data();
      blocks_.emplace_back(std::move(new_block));
    }
    return Status::OK();
  }

  void Reset() {
    current_offset_ = 0;
    current_out_buffer_ = NULLPTR;
    current_remaining_bytes_ = 0;
    blocks_.clear();
  }

  int64_t current_remaining_bytes() const { return current_remaining_bytes_; }

  Result<std::vector<std::shared_ptr<ResizableBuffer>>> Finish() {
    if (!blocks_.empty()) {
      ARROW_RETURN_NOT_OK(FinishLastBlock());
    }
    current_offset_ = 0;
    current_out_buffer_ = NULLPTR;
    current_remaining_bytes_ = 0;
    return std::move(blocks_);
  }

 private:
  Status FinishLastBlock() {
    if (current_remaining_bytes_ > 0) {
      // Avoid leaking uninitialized bytes from the allocator
      ARROW_RETURN_NOT_OK(
          blocks_.back()->Resize(blocks_.back()->size() - current_remaining_bytes_,
                                 /*shrink_to_fit=*/true));
      blocks_.back()->ZeroPadding();
    }
    return Status::OK();
  }

  MemoryPool* pool_;
  int64_t alignment_;
  int64_t blocksize_ = kDefaultBlocksize;
  std::vector<std::shared_ptr<ResizableBuffer>> blocks_;

  int32_t current_offset_ = 0;
  uint8_t* current_out_buffer_ = NULLPTR;
  int64_t current_remaining_bytes_ = 0;
};

}  // namespace internal

class ARROW_EXPORT BinaryViewBuilder : public ArrayBuilder {
 public:
  using TypeClass = BinaryViewType;

  // this constructor provided for MakeBuilder compatibility
  BinaryViewBuilder(const std::shared_ptr<DataType>&, MemoryPool* pool);

  explicit BinaryViewBuilder(MemoryPool* pool = default_memory_pool(),
                             int64_t alignment = kDefaultBufferAlignment)
      : ArrayBuilder(pool, alignment),
        data_builder_(pool, alignment),
        data_heap_builder_(pool, alignment) {}

  /// Set the size for future preallocated data buffers.
  ///
  /// The default size is 32KB, so after each 32KB of string data appended to the builder
  /// a new data buffer will be allocated. Adjust this to a larger value to decrease the
  /// frequency of allocation, or to a smaller value to lower the overhead of each
  /// allocation.
  void SetBlockSize(int64_t blocksize) { data_heap_builder_.SetBlockSize(blocksize); }

  /// The number of bytes which can be appended to this builder without allocating another
  /// data buffer.
  int64_t current_block_bytes_remaining() const {
    return data_heap_builder_.current_remaining_bytes();
  }

  Status Append(const uint8_t* value, int64_t length) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppendToBitmap(true);
    ARROW_ASSIGN_OR_RAISE(auto v,
                          data_heap_builder_.Append</*Safe=*/true>(value, length));
    data_builder_.UnsafeAppend(v);
    return Status::OK();
  }

  Status Append(const char* value, int64_t length) {
    return Append(reinterpret_cast<const uint8_t*>(value), length);
  }

  Status Append(std::string_view value) {
    return Append(value.data(), static_cast<int64_t>(value.size()));
  }

  /// \brief Append without checking capacity
  ///
  /// Builder should have been presized using Reserve() and ReserveData(),
  /// respectively, and the value must not be larger than 2GB
  void UnsafeAppend(const uint8_t* value, int64_t length) {
    UnsafeAppendToBitmap(true);
    auto v = data_heap_builder_.Append</*Safe=*/false>(value, length);
    data_builder_.UnsafeAppend(v);
  }

  void UnsafeAppend(const char* value, int64_t length) {
    UnsafeAppend(reinterpret_cast<const uint8_t*>(value), length);
  }

  void UnsafeAppend(const std::string& value) {
    UnsafeAppend(value.c_str(), static_cast<int64_t>(value.size()));
  }

  void UnsafeAppend(std::string_view value) {
    UnsafeAppend(value.data(), static_cast<int64_t>(value.size()));
  }

  /// \brief Ensures there is enough allocated available capacity in the
  /// out-of-line data heap to append the indicated number of bytes without
  /// additional allocations
  Status ReserveData(int64_t length);

  Status AppendNulls(int64_t length) final {
    ARROW_RETURN_NOT_OK(Reserve(length));
    data_builder_.UnsafeAppend(length, BinaryViewType::c_type{});
    UnsafeSetNull(length);
    return Status::OK();
  }

  /// \brief Append a single null element
  Status AppendNull() final {
    ARROW_RETURN_NOT_OK(Reserve(1));
    data_builder_.UnsafeAppend(BinaryViewType::c_type{});
    UnsafeAppendToBitmap(false);
    return Status::OK();
  }

  /// \brief Append a empty element (length-0 inline string)
  Status AppendEmptyValue() final {
    ARROW_RETURN_NOT_OK(Reserve(1));
    data_builder_.UnsafeAppend(BinaryViewType::c_type{});
    UnsafeAppendToBitmap(true);
    return Status::OK();
  }

  /// \brief Append several empty elements
  Status AppendEmptyValues(int64_t length) final {
    ARROW_RETURN_NOT_OK(Reserve(length));
    data_builder_.UnsafeAppend(length, BinaryViewType::c_type{});
    UnsafeSetNotNull(length);
    return Status::OK();
  }

  void UnsafeAppendNull() {
    data_builder_.UnsafeAppend(BinaryViewType::c_type{});
    UnsafeAppendToBitmap(false);
  }

  void UnsafeAppendEmptyValue() {
    data_builder_.UnsafeAppend(BinaryViewType::c_type{});
    UnsafeAppendToBitmap(true);
  }

  /// \brief Append a slice of a BinaryViewArray passed as an ArraySpan. Copies
  /// the underlying out-of-line string memory to avoid memory lifetime issues
  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override;

  void Reset() override;

  Status Resize(int64_t capacity) override {
    ARROW_RETURN_NOT_OK(CheckCapacity(capacity));
    capacity = std::max(capacity, kMinBuilderCapacity);
    ARROW_RETURN_NOT_OK(data_builder_.Resize(capacity));
    return ArrayBuilder::Resize(capacity);
  }

  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  std::shared_ptr<DataType> type() const override { return binary_view(); }

 protected:
  TypedBufferBuilder<BinaryViewType::c_type> data_builder_;

  // Accumulates out-of-line data in fixed-size chunks which are then attached
  // to the resulting ArrayData
  internal::StringHeapBuilder data_heap_builder_;
};

class ARROW_EXPORT StringViewBuilder : public BinaryViewBuilder {
 public:
  using BinaryViewBuilder::BinaryViewBuilder;
  std::shared_ptr<DataType> type() const override { return utf8_view(); }
};

// ----------------------------------------------------------------------
// FixedSizeBinaryBuilder

class ARROW_EXPORT FixedSizeBinaryBuilder : public ArrayBuilder {
 public:
  using TypeClass = FixedSizeBinaryType;

  explicit FixedSizeBinaryBuilder(const std::shared_ptr<DataType>& type,
                                  MemoryPool* pool = default_memory_pool(),
                                  int64_t alignment = kDefaultBufferAlignment);

  Status Append(const uint8_t* value) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppend(value);
    return Status::OK();
  }

  Status Append(const char* value) {
    return Append(reinterpret_cast<const uint8_t*>(value));
  }

  Status Append(std::string_view view) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppend(view);
    return Status::OK();
  }

  Status Append(const std::string& s) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppend(s);
    return Status::OK();
  }

  Status Append(const Buffer& s) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppend(s);
    return Status::OK();
  }

  Status Append(const std::shared_ptr<Buffer>& s) { return Append(*s); }

  template <size_t NBYTES>
  Status Append(const std::array<uint8_t, NBYTES>& value) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppend(
        std::string_view(reinterpret_cast<const char*>(value.data()), value.size()));
    return Status::OK();
  }

  Status AppendValues(const uint8_t* data, int64_t length,
                      const uint8_t* valid_bytes = NULLPTR);

  Status AppendValues(const uint8_t* data, int64_t length, const uint8_t* validity,
                      int64_t bitmap_offset);

  Status AppendNull() final;
  Status AppendNulls(int64_t length) final;

  Status AppendEmptyValue() final;
  Status AppendEmptyValues(int64_t length) final;

  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override {
    return AppendValues(
        array.GetValues<uint8_t>(1, 0) + ((array.offset + offset) * byte_width_), length,
        array.GetValues<uint8_t>(0, 0), array.offset + offset);
  }

  void UnsafeAppend(const uint8_t* value) {
    UnsafeAppendToBitmap(true);
    if (ARROW_PREDICT_TRUE(byte_width_ > 0)) {
      byte_builder_.UnsafeAppend(value, byte_width_);
    }
  }

  void UnsafeAppend(const char* value) {
    UnsafeAppend(reinterpret_cast<const uint8_t*>(value));
  }

  void UnsafeAppend(std::string_view value) {
#ifndef NDEBUG
    CheckValueSize(static_cast<size_t>(value.size()));
#endif
    UnsafeAppend(reinterpret_cast<const uint8_t*>(value.data()));
  }

  void UnsafeAppend(const Buffer& s) { UnsafeAppend(std::string_view{s}); }

  void UnsafeAppend(const std::shared_ptr<Buffer>& s) { UnsafeAppend(*s); }

  void UnsafeAppendNull() {
    UnsafeAppendToBitmap(false);
    byte_builder_.UnsafeAppend(/*num_copies=*/byte_width_, 0);
  }

  Status ValidateOverflow(int64_t new_bytes) const {
    auto new_size = byte_builder_.length() + new_bytes;
    if (ARROW_PREDICT_FALSE(new_size > memory_limit())) {
      return Status::CapacityError("array cannot contain more than ", memory_limit(),
                                   " bytes, have ", new_size);
    } else {
      return Status::OK();
    }
  }

  /// \brief Ensures there is enough allocated capacity to append the indicated
  /// number of bytes to the value data buffer without additional allocations
  Status ReserveData(int64_t elements) {
    ARROW_RETURN_NOT_OK(ValidateOverflow(elements));
    return byte_builder_.Reserve(elements);
  }

  void Reset() override;
  Status Resize(int64_t capacity) override;
  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<FixedSizeBinaryArray>* out) { return FinishTyped(out); }

  /// \return size of values buffer so far
  int64_t value_data_length() const { return byte_builder_.length(); }

  int32_t byte_width() const { return byte_width_; }

  /// Temporary access to a value.
  ///
  /// This pointer becomes invalid on the next modifying operation.
  const uint8_t* GetValue(int64_t i) const;

  /// Temporary access to a value.
  ///
  /// This view becomes invalid on the next modifying operation.
  std::string_view GetView(int64_t i) const;

  static constexpr int64_t memory_limit() {
    return std::numeric_limits<int64_t>::max() - 1;
  }

  std::shared_ptr<DataType> type() const override {
    return fixed_size_binary(byte_width_);
  }

 protected:
  int32_t byte_width_;
  BufferBuilder byte_builder_;

  /// Temporary access to a value.
  ///
  /// This pointer becomes invalid on the next modifying operation.
  uint8_t* GetMutableValue(int64_t i) {
    uint8_t* data_ptr = byte_builder_.mutable_data();
    return data_ptr + i * byte_width_;
  }

  void CheckValueSize(int64_t size);
};

/// @}

// ----------------------------------------------------------------------
// Chunked builders: build a sequence of BinaryArray or StringArray that are
// limited to a particular size (to the upper limit of 2GB)

namespace internal {

class ARROW_EXPORT ChunkedBinaryBuilder {
 public:
  explicit ChunkedBinaryBuilder(int32_t max_chunk_value_length,
                                MemoryPool* pool = default_memory_pool());

  ChunkedBinaryBuilder(int32_t max_chunk_value_length, int32_t max_chunk_length,
                       MemoryPool* pool = default_memory_pool());

  virtual ~ChunkedBinaryBuilder() = default;

  Status Append(const uint8_t* value, int32_t length) {
    if (ARROW_PREDICT_FALSE(length + builder_->value_data_length() >
                            max_chunk_value_length_)) {
      if (builder_->value_data_length() == 0) {
        // The current item is larger than max_chunk_size_;
        // this chunk will be oversize and hold *only* this item
        ARROW_RETURN_NOT_OK(builder_->Append(value, length));
        return NextChunk();
      }
      // The current item would cause builder_->value_data_length() to exceed
      // max_chunk_size_, so finish this chunk and append the current item to the next
      // chunk
      ARROW_RETURN_NOT_OK(NextChunk());
      return Append(value, length);
    }

    if (ARROW_PREDICT_FALSE(builder_->length() == max_chunk_length_)) {
      // The current item would cause builder_->length() to exceed max_chunk_length_, so
      // finish this chunk and append the current item to the next chunk
      ARROW_RETURN_NOT_OK(NextChunk());
    }

    return builder_->Append(value, length);
  }

  Status Append(std::string_view value) {
    return Append(reinterpret_cast<const uint8_t*>(value.data()),
                  static_cast<int32_t>(value.size()));
  }

  Status AppendNull() {
    if (ARROW_PREDICT_FALSE(builder_->length() == max_chunk_length_)) {
      ARROW_RETURN_NOT_OK(NextChunk());
    }
    return builder_->AppendNull();
  }

  Status Reserve(int64_t values);

  virtual Status Finish(ArrayVector* out);

 protected:
  Status NextChunk();

  // maximum total character data size per chunk
  int64_t max_chunk_value_length_;

  // maximum elements allowed per chunk
  int64_t max_chunk_length_ = kListMaximumElements;

  // when Reserve() would cause builder_ to exceed its max_chunk_length_,
  // add to extra_capacity_ instead and wait to reserve until the next chunk
  int64_t extra_capacity_ = 0;

  std::unique_ptr<BinaryBuilder> builder_;
  std::vector<std::shared_ptr<Array>> chunks_;
};

class ARROW_EXPORT ChunkedStringBuilder : public ChunkedBinaryBuilder {
 public:
  using ChunkedBinaryBuilder::ChunkedBinaryBuilder;

  Status Finish(ArrayVector* out) override;
};

}  // namespace internal

}  // namespace arrow
