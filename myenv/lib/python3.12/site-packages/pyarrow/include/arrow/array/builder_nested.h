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

#include "arrow/array/array_nested.h"
#include "arrow/array/builder_base.h"
#include "arrow/array/data.h"
#include "arrow/buffer.h"
#include "arrow/buffer_builder.h"
#include "arrow/status.h"
#include "arrow/type.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

/// \addtogroup nested-builders
///
/// @{

// ----------------------------------------------------------------------
// VarLengthListLikeBuilder

template <typename TYPE>
class ARROW_EXPORT VarLengthListLikeBuilder : public ArrayBuilder {
 public:
  using TypeClass = TYPE;
  using offset_type = typename TypeClass::offset_type;

  /// Use this constructor to incrementally build the value array along with offsets and
  /// null bitmap.
  VarLengthListLikeBuilder(MemoryPool* pool,
                           std::shared_ptr<ArrayBuilder> const& value_builder,
                           const std::shared_ptr<DataType>& type,
                           int64_t alignment = kDefaultBufferAlignment)
      : ArrayBuilder(pool, alignment),
        offsets_builder_(pool, alignment),
        value_builder_(value_builder),
        value_field_(type->field(0)->WithType(NULLPTR)) {}

  VarLengthListLikeBuilder(MemoryPool* pool,
                           std::shared_ptr<ArrayBuilder> const& value_builder,
                           int64_t alignment = kDefaultBufferAlignment)
      : VarLengthListLikeBuilder(pool, value_builder,
                                 std::make_shared<TYPE>(value_builder->type()),
                                 alignment) {}

  ~VarLengthListLikeBuilder() override = default;

  Status Resize(int64_t capacity) override {
    if (ARROW_PREDICT_FALSE(capacity > maximum_elements())) {
      return Status::CapacityError(type_name(),
                                   " array cannot reserve space for more than ",
                                   maximum_elements(), " got ", capacity);
    }
    ARROW_RETURN_NOT_OK(CheckCapacity(capacity));

    // One more than requested for list offsets
    const int64_t offsets_capacity =
        is_list_view(TYPE::type_id) ? capacity : capacity + 1;
    ARROW_RETURN_NOT_OK(offsets_builder_.Resize(offsets_capacity));
    return ArrayBuilder::Resize(capacity);
  }

  void Reset() override {
    ArrayBuilder::Reset();
    offsets_builder_.Reset();
    value_builder_->Reset();
  }

  /// \brief Start a new variable-length list slot
  ///
  /// This function should be called before appending elements to the
  /// value builder. Elements appended to the value builder before this function
  /// is called for the first time, will not be members of any list value.
  ///
  /// After this function is called, list_length elements SHOULD be appended to
  /// the values builder. If this contract is violated, the behavior is defined by
  /// the concrete builder implementation and SHOULD NOT be relied upon unless
  /// the caller is specifically building a [Large]List or [Large]ListView array.
  ///
  /// For [Large]List arrays, the list slot length will be the number of elements
  /// appended to the values builder before the next call to Append* or Finish. For
  /// [Large]ListView arrays, the list slot length will be exactly list_length, but if
  /// Append* is called before at least list_length elements are appended to the values
  /// builder, the current list slot will share elements with the next list
  /// slots or an invalid [Large]ListView array will be generated because there
  /// aren't enough elements in the values builder to fill the list slots.
  ///
  /// If you're building a [Large]List and don't need to be compatible
  /// with [Large]ListView, then `BaseListBuilder::Append(bool is_valid)`
  /// is a simpler API.
  ///
  /// \pre if is_valid is false, list_length MUST be 0
  /// \param is_valid Whether the new list slot is valid
  /// \param list_length The number of elements in the list
  Status Append(bool is_valid, int64_t list_length) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    assert(is_valid || list_length == 0);
    UnsafeAppendToBitmap(is_valid);
    UnsafeAppendDimensions(/*offset=*/value_builder_->length(), /*size=*/list_length);
    return Status::OK();
  }

  Status AppendNull() final {
    // Append() a null list slot with list_length=0.
    //
    // When building [Large]List arrays, elements being appended to the values builder
    // before the next call to Append* or Finish will extend the list slot length, but
    // that is totally fine because list arrays admit non-empty null list slots.
    //
    // In the case of [Large]ListViews that's not a problem either because the
    // list slot length remains zero.
    return Append(false, 0);
  }

  Status AppendNulls(int64_t length) final {
    ARROW_RETURN_NOT_OK(Reserve(length));
    UnsafeAppendToBitmap(length, false);
    UnsafeAppendEmptyDimensions(/*num_values=*/length);
    return Status::OK();
  }

  /// \brief Append an empty list slot
  ///
  /// \post Another call to Append* or Finish should be made before appending to
  /// the values builder to ensure list slot remains empty
  Status AppendEmptyValue() final { return Append(true, 0); }

  /// \brief Append an empty list slot
  ///
  /// \post Another call to Append* or Finish should be made before appending to
  /// the values builder to ensure the last list slot remains empty
  Status AppendEmptyValues(int64_t length) final {
    ARROW_RETURN_NOT_OK(Reserve(length));
    UnsafeAppendToBitmap(length, true);
    UnsafeAppendEmptyDimensions(/*num_values=*/length);
    return Status::OK();
  }

  /// \brief Vector append
  ///
  /// For list-array builders, the sizes are inferred from the offsets.
  /// BaseListBuilder<T> provides an implementation that doesn't take sizes, but
  /// this virtual function allows dispatching calls to both list-array and
  /// list-view-array builders (which need the sizes)
  ///
  /// \param offsets The offsets of the variable-length lists
  /// \param sizes The sizes of the variable-length lists
  /// \param length The number of offsets, sizes, and validity bits to append
  /// \param valid_bytes If passed, valid_bytes is of equal length to values,
  /// and any zero byte will be considered as a null for that slot
  virtual Status AppendValues(const offset_type* offsets, const offset_type* sizes,
                              int64_t length, const uint8_t* valid_bytes) = 0;

  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override {
    const offset_type* offsets = array.GetValues<offset_type>(1);
    [[maybe_unused]] const offset_type* sizes = NULLPTR;
    if constexpr (is_list_view(TYPE::type_id)) {
      sizes = array.GetValues<offset_type>(2);
    }
    const bool all_valid = !array.MayHaveLogicalNulls();
    const uint8_t* validity = array.HasValidityBitmap() ? array.buffers[0].data : NULLPTR;
    ARROW_RETURN_NOT_OK(Reserve(length));
    for (int64_t row = offset; row < offset + length; row++) {
      const bool is_valid =
          all_valid || (validity && bit_util::GetBit(validity, array.offset + row)) ||
          array.IsValid(row);
      int64_t size = 0;
      if (is_valid) {
        if constexpr (is_list_view(TYPE::type_id)) {
          size = sizes[row];
        } else {
          size = offsets[row + 1] - offsets[row];
        }
      }
      UnsafeAppendToBitmap(is_valid);
      UnsafeAppendDimensions(/*offset=*/value_builder_->length(), size);
      if (is_valid) {
        ARROW_RETURN_NOT_OK(
            value_builder_->AppendArraySlice(array.child_data[0], offsets[row], size));
      }
    }
    return Status::OK();
  }

  Status ValidateOverflow(int64_t new_elements) const {
    auto new_length = value_builder_->length() + new_elements;
    if (ARROW_PREDICT_FALSE(new_length > maximum_elements())) {
      return Status::CapacityError(type_name(), " array cannot contain more than ",
                                   maximum_elements(), " elements, have ", new_elements);
    } else {
      return Status::OK();
    }
  }

  ArrayBuilder* value_builder() const { return value_builder_.get(); }

  // Cannot make this a static attribute because of linking issues
  static constexpr int64_t maximum_elements() {
    return std::numeric_limits<offset_type>::max() - 1;
  }

  std::shared_ptr<DataType> type() const override {
    return std::make_shared<TYPE>(value_field_->WithType(value_builder_->type()));
  }

 private:
  static constexpr const char* type_name() {
    if constexpr (is_list_view(TYPE::type_id)) {
      return "ListView";
    } else {
      return "List";
    }
  }

 protected:
  /// \brief Append dimensions for num_values empty list slots.
  ///
  /// ListViewBuilder overrides this to also append the sizes.
  virtual void UnsafeAppendEmptyDimensions(int64_t num_values) {
    const int64_t offset = value_builder_->length();
    for (int64_t i = 0; i < num_values; ++i) {
      offsets_builder_.UnsafeAppend(static_cast<offset_type>(offset));
    }
  }

  /// \brief Append dimensions for a single list slot.
  ///
  /// ListViewBuilder overrides this to also append the size.
  virtual void UnsafeAppendDimensions(int64_t offset, int64_t size) {
    offsets_builder_.UnsafeAppend(static_cast<offset_type>(offset));
  }

  TypedBufferBuilder<offset_type> offsets_builder_;
  std::shared_ptr<ArrayBuilder> value_builder_;
  std::shared_ptr<Field> value_field_;
};

// ----------------------------------------------------------------------
// ListBuilder / LargeListBuilder

template <typename TYPE>
class ARROW_EXPORT BaseListBuilder : public VarLengthListLikeBuilder<TYPE> {
 private:
  using BASE = VarLengthListLikeBuilder<TYPE>;

 public:
  using TypeClass = TYPE;
  using offset_type = typename BASE::offset_type;

  using BASE::BASE;

  using BASE::Append;

  ~BaseListBuilder() override = default;

  /// \brief Start a new variable-length list slot
  ///
  /// This function should be called before beginning to append elements to the
  /// value builder
  Status Append(bool is_valid = true) {
    // The value_length parameter to BASE::Append(bool, int64_t) is ignored when
    // building a list array, so we can pass 0 here.
    return BASE::Append(is_valid, 0);
  }

  /// \brief Vector append
  ///
  /// If passed, valid_bytes is of equal length to values, and any zero byte
  /// will be considered as a null for that slot
  Status AppendValues(const offset_type* offsets, int64_t length,
                      const uint8_t* valid_bytes = NULLPTR) {
    ARROW_RETURN_NOT_OK(this->Reserve(length));
    this->UnsafeAppendToBitmap(valid_bytes, length);
    this->offsets_builder_.UnsafeAppend(offsets, length);
    return Status::OK();
  }

  Status AppendValues(const offset_type* offsets, const offset_type* sizes,
                      int64_t length, const uint8_t* valid_bytes) final {
    // Offsets are assumed to be valid, but the first length-1 sizes have to be
    // consistent with the offsets to partially rule out the possibility that the
    // caller is passing sizes that could work if building a list-view, but don't
    // work on building a list that requires offsets to be non-decreasing.
    //
    // CAUTION: the last size element (`sizes[length - 1]`) is not
    // validated and could be inconsistent with the offsets given in a
    // subsequent call to AppendValues.
#ifndef NDEBUG
    if (sizes) {
      for (int64_t i = 0; i < length - 1; ++i) {
        if (ARROW_PREDICT_FALSE(offsets[i] != offsets[i + 1] - sizes[i])) {
          if (!valid_bytes || valid_bytes[i]) {
            return Status::Invalid(
                "BaseListBuilder: sizes are inconsistent with offsets provided");
          }
        }
      }
    }
#endif
    return AppendValues(offsets, length, valid_bytes);
  }

  Status AppendValues(const offset_type* offsets, const offset_type* sizes,
                      int64_t length) {
    return AppendValues(offsets, sizes, length, /*valid_bytes=*/NULLPTR);
  }

  Status AppendNextOffset() {
    ARROW_RETURN_NOT_OK(this->ValidateOverflow(0));
    const int64_t num_values = this->value_builder_->length();
    return this->offsets_builder_.Append(static_cast<offset_type>(num_values));
  }

  Status FinishInternal(std::shared_ptr<ArrayData>* out) override {
    ARROW_RETURN_NOT_OK(AppendNextOffset());

    // Offset padding zeroed by BufferBuilder
    std::shared_ptr<Buffer> offsets;
    std::shared_ptr<Buffer> null_bitmap;
    ARROW_RETURN_NOT_OK(this->offsets_builder_.Finish(&offsets));
    ARROW_RETURN_NOT_OK(this->null_bitmap_builder_.Finish(&null_bitmap));

    if (this->value_builder_->length() == 0) {
      // Try to make sure we get a non-null values buffer (ARROW-2744)
      ARROW_RETURN_NOT_OK(this->value_builder_->Resize(0));
    }

    std::shared_ptr<ArrayData> items;
    ARROW_RETURN_NOT_OK(this->value_builder_->FinishInternal(&items));

    *out = ArrayData::Make(this->type(), this->length_,
                           {std::move(null_bitmap), std::move(offsets)},
                           {std::move(items)}, this->null_count_);
    this->Reset();
    return Status::OK();
  }
};

/// \class ListBuilder
/// \brief Builder class for variable-length list array value types
///
/// To use this class, you must append values to the child array builder and use
/// the Append function to delimit each distinct list value (once the values
/// have been appended to the child array) or use the bulk API to append
/// a sequence of offsets and null values.
///
/// A note on types.  Per arrow/type.h all types in the c++ implementation are
/// logical so even though this class always builds list array, this can
/// represent multiple different logical types.  If no logical type is provided
/// at construction time, the class defaults to List<T> where t is taken from the
/// value_builder/values that the object is constructed with.
class ARROW_EXPORT ListBuilder : public BaseListBuilder<ListType> {
 public:
  using BaseListBuilder::BaseListBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<ListArray>* out) { return FinishTyped(out); }
};

/// \class LargeListBuilder
/// \brief Builder class for large variable-length list array value types
///
/// Like ListBuilder, but to create large list arrays (with 64-bit offsets).
class ARROW_EXPORT LargeListBuilder : public BaseListBuilder<LargeListType> {
 public:
  using BaseListBuilder::BaseListBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<LargeListArray>* out) { return FinishTyped(out); }
};

// ----------------------------------------------------------------------
// ListViewBuilder / LargeListViewBuilder

template <typename TYPE>
class ARROW_EXPORT BaseListViewBuilder : public VarLengthListLikeBuilder<TYPE> {
 private:
  using BASE = VarLengthListLikeBuilder<TYPE>;

 public:
  using TypeClass = TYPE;
  using offset_type = typename BASE::offset_type;

  using BASE::BASE;

  ~BaseListViewBuilder() override = default;

  Status Resize(int64_t capacity) override {
    ARROW_RETURN_NOT_OK(BASE::Resize(capacity));
    return sizes_builder_.Resize(capacity);
  }

  void Reset() override {
    BASE::Reset();
    sizes_builder_.Reset();
  }

  /// \brief Vector append
  ///
  /// If passed, valid_bytes is of equal length to values, and any zero byte
  /// will be considered as a null for that slot
  Status AppendValues(const offset_type* offsets, const offset_type* sizes,
                      int64_t length, const uint8_t* valid_bytes) final {
    ARROW_RETURN_NOT_OK(this->Reserve(length));
    this->UnsafeAppendToBitmap(valid_bytes, length);
    this->offsets_builder_.UnsafeAppend(offsets, length);
    this->sizes_builder_.UnsafeAppend(sizes, length);
    return Status::OK();
  }

  Status AppendValues(const offset_type* offsets, const offset_type* sizes,
                      int64_t length) {
    return AppendValues(offsets, sizes, length, /*valid_bytes=*/NULLPTR);
  }

  Status FinishInternal(std::shared_ptr<ArrayData>* out) override {
    // Offset and sizes padding zeroed by BufferBuilder
    std::shared_ptr<Buffer> null_bitmap;
    std::shared_ptr<Buffer> offsets;
    std::shared_ptr<Buffer> sizes;
    ARROW_RETURN_NOT_OK(this->null_bitmap_builder_.Finish(&null_bitmap));
    ARROW_RETURN_NOT_OK(this->offsets_builder_.Finish(&offsets));
    ARROW_RETURN_NOT_OK(this->sizes_builder_.Finish(&sizes));

    if (this->value_builder_->length() == 0) {
      // Try to make sure we get a non-null values buffer (ARROW-2744)
      ARROW_RETURN_NOT_OK(this->value_builder_->Resize(0));
    }

    std::shared_ptr<ArrayData> items;
    ARROW_RETURN_NOT_OK(this->value_builder_->FinishInternal(&items));

    *out = ArrayData::Make(this->type(), this->length_,
                           {std::move(null_bitmap), std::move(offsets), std::move(sizes)},
                           {std::move(items)}, this->null_count_);
    this->Reset();
    return Status::OK();
  }

 protected:
  void UnsafeAppendEmptyDimensions(int64_t num_values) override {
    for (int64_t i = 0; i < num_values; ++i) {
      this->offsets_builder_.UnsafeAppend(0);
    }
    for (int64_t i = 0; i < num_values; ++i) {
      this->sizes_builder_.UnsafeAppend(0);
    }
  }

  void UnsafeAppendDimensions(int64_t offset, int64_t size) override {
    this->offsets_builder_.UnsafeAppend(static_cast<offset_type>(offset));
    this->sizes_builder_.UnsafeAppend(static_cast<offset_type>(size));
  }

 private:
  TypedBufferBuilder<offset_type> sizes_builder_;
};

class ARROW_EXPORT ListViewBuilder final : public BaseListViewBuilder<ListViewType> {
 public:
  using BaseListViewBuilder::BaseListViewBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<ListViewArray>* out) { return FinishTyped(out); }
};

class ARROW_EXPORT LargeListViewBuilder final
    : public BaseListViewBuilder<LargeListViewType> {
 public:
  using BaseListViewBuilder::BaseListViewBuilder;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<LargeListViewArray>* out) { return FinishTyped(out); }
};

// ----------------------------------------------------------------------
// Map builder

/// \class MapBuilder
/// \brief Builder class for arrays of variable-size maps
///
/// To use this class, you must use the Append function to delimit each distinct
/// map before appending values to the key and item array builders, or use the
/// bulk API to append a sequence of offsets and null maps.
///
/// Key uniqueness and ordering are not validated.
class ARROW_EXPORT MapBuilder : public ArrayBuilder {
 public:
  /// Use this constructor to define the built array's type explicitly. If key_builder
  /// or item_builder has indeterminate type, this builder will also.
  MapBuilder(MemoryPool* pool, const std::shared_ptr<ArrayBuilder>& key_builder,
             const std::shared_ptr<ArrayBuilder>& item_builder,
             const std::shared_ptr<DataType>& type);

  /// Use this constructor to infer the built array's type. If key_builder or
  /// item_builder has indeterminate type, this builder will also.
  MapBuilder(MemoryPool* pool, const std::shared_ptr<ArrayBuilder>& key_builder,
             const std::shared_ptr<ArrayBuilder>& item_builder, bool keys_sorted = false);

  MapBuilder(MemoryPool* pool, const std::shared_ptr<ArrayBuilder>& item_builder,
             const std::shared_ptr<DataType>& type);

  Status Resize(int64_t capacity) override;
  void Reset() override;
  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<MapArray>* out) { return FinishTyped(out); }

  /// \brief Vector append
  ///
  /// If passed, valid_bytes is of equal length to values, and any zero byte
  /// will be considered as a null for that slot
  Status AppendValues(const int32_t* offsets, int64_t length,
                      const uint8_t* valid_bytes = NULLPTR);

  /// \brief Start a new variable-length map slot
  ///
  /// This function should be called before beginning to append elements to the
  /// key and item builders
  Status Append();

  Status AppendNull() final;

  Status AppendNulls(int64_t length) final;

  Status AppendEmptyValue() final;

  Status AppendEmptyValues(int64_t length) final;

  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override {
    const int32_t* offsets = array.GetValues<int32_t>(1);
    const bool all_valid = !array.MayHaveLogicalNulls();
    const uint8_t* validity = array.HasValidityBitmap() ? array.buffers[0].data : NULLPTR;
    for (int64_t row = offset; row < offset + length; row++) {
      const bool is_valid =
          all_valid || (validity && bit_util::GetBit(validity, array.offset + row)) ||
          array.IsValid(row);
      if (is_valid) {
        ARROW_RETURN_NOT_OK(Append());
        const int64_t slot_length = offsets[row + 1] - offsets[row];
        // Add together the inner StructArray offset to the Map/List offset
        int64_t key_value_offset = array.child_data[0].offset + offsets[row];
        ARROW_RETURN_NOT_OK(key_builder_->AppendArraySlice(
            array.child_data[0].child_data[0], key_value_offset, slot_length));
        ARROW_RETURN_NOT_OK(item_builder_->AppendArraySlice(
            array.child_data[0].child_data[1], key_value_offset, slot_length));
      } else {
        ARROW_RETURN_NOT_OK(AppendNull());
      }
    }
    return Status::OK();
  }

  /// \brief Get builder to append keys.
  ///
  /// Append a key with this builder should be followed by appending
  /// an item or null value with item_builder().
  ArrayBuilder* key_builder() const { return key_builder_.get(); }

  /// \brief Get builder to append items
  ///
  /// Appending an item with this builder should have been preceded
  /// by appending a key with key_builder().
  ArrayBuilder* item_builder() const { return item_builder_.get(); }

  /// \brief Get builder to add Map entries as struct values.
  ///
  /// This is used instead of key_builder()/item_builder() and allows
  /// the Map to be built as a list of struct values.
  ArrayBuilder* value_builder() const { return list_builder_->value_builder(); }

  std::shared_ptr<DataType> type() const override {
    // Key and Item builder may update types, but they don't contain the field names,
    // so we need to reconstruct the type. (See ARROW-13735.)
    return std::make_shared<MapType>(
        field(entries_name_,
              struct_({field(key_name_, key_builder_->type(), false),
                       field(item_name_, item_builder_->type(), item_nullable_)}),
              false),
        keys_sorted_);
  }

  Status ValidateOverflow(int64_t new_elements) {
    return list_builder_->ValidateOverflow(new_elements);
  }

 protected:
  inline Status AdjustStructBuilderLength();

 protected:
  bool keys_sorted_ = false;
  bool item_nullable_ = false;
  std::string entries_name_;
  std::string key_name_;
  std::string item_name_;
  std::shared_ptr<ListBuilder> list_builder_;
  std::shared_ptr<ArrayBuilder> key_builder_;
  std::shared_ptr<ArrayBuilder> item_builder_;
};

// ----------------------------------------------------------------------
// FixedSizeList builder

/// \class FixedSizeListBuilder
/// \brief Builder class for fixed-length list array value types
class ARROW_EXPORT FixedSizeListBuilder : public ArrayBuilder {
 public:
  /// Use this constructor to define the built array's type explicitly. If value_builder
  /// has indeterminate type, this builder will also.
  FixedSizeListBuilder(MemoryPool* pool,
                       std::shared_ptr<ArrayBuilder> const& value_builder,
                       int32_t list_size);

  /// Use this constructor to infer the built array's type. If value_builder has
  /// indeterminate type, this builder will also.
  FixedSizeListBuilder(MemoryPool* pool,
                       std::shared_ptr<ArrayBuilder> const& value_builder,
                       const std::shared_ptr<DataType>& type);

  Status Resize(int64_t capacity) override;
  void Reset() override;
  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<FixedSizeListArray>* out) { return FinishTyped(out); }

  /// \brief Append a valid fixed length list.
  ///
  /// This function affects only the validity bitmap; the child values must be appended
  /// using the child array builder.
  Status Append();

  /// \brief Vector append
  ///
  /// If passed, valid_bytes will be read and any zero byte
  /// will cause the corresponding slot to be null
  ///
  /// This function affects only the validity bitmap; the child values must be appended
  /// using the child array builder. This includes appending nulls for null lists.
  /// XXX this restriction is confusing, should this method be omitted?
  Status AppendValues(int64_t length, const uint8_t* valid_bytes = NULLPTR);

  /// \brief Append a null fixed length list.
  ///
  /// The child array builder will have the appropriate number of nulls appended
  /// automatically.
  Status AppendNull() final;

  /// \brief Append length null fixed length lists.
  ///
  /// The child array builder will have the appropriate number of nulls appended
  /// automatically.
  Status AppendNulls(int64_t length) final;

  Status ValidateOverflow(int64_t new_elements);

  Status AppendEmptyValue() final;

  Status AppendEmptyValues(int64_t length) final;

  Status AppendArraySlice(const ArraySpan& array, int64_t offset, int64_t length) final {
    const uint8_t* validity = array.MayHaveNulls() ? array.buffers[0].data : NULLPTR;
    for (int64_t row = offset; row < offset + length; row++) {
      if (!validity || bit_util::GetBit(validity, array.offset + row)) {
        ARROW_RETURN_NOT_OK(value_builder_->AppendArraySlice(
            array.child_data[0], list_size_ * (array.offset + row), list_size_));
        ARROW_RETURN_NOT_OK(Append());
      } else {
        ARROW_RETURN_NOT_OK(AppendNull());
      }
    }
    return Status::OK();
  }

  ArrayBuilder* value_builder() const { return value_builder_.get(); }

  std::shared_ptr<DataType> type() const override {
    return fixed_size_list(value_field_->WithType(value_builder_->type()), list_size_);
  }

  // Cannot make this a static attribute because of linking issues
  static constexpr int64_t maximum_elements() {
    return std::numeric_limits<FixedSizeListType::offset_type>::max() - 1;
  }

 protected:
  std::shared_ptr<Field> value_field_;
  const int32_t list_size_;
  std::shared_ptr<ArrayBuilder> value_builder_;
};

// ----------------------------------------------------------------------
// Struct

// ---------------------------------------------------------------------------------
// StructArray builder
/// Append, Resize and Reserve methods are acting on StructBuilder.
/// Please make sure all these methods of all child-builders' are consistently
/// called to maintain data-structure consistency.
class ARROW_EXPORT StructBuilder : public ArrayBuilder {
 public:
  /// If any of field_builders has indeterminate type, this builder will also
  StructBuilder(const std::shared_ptr<DataType>& type, MemoryPool* pool,
                std::vector<std::shared_ptr<ArrayBuilder>> field_builders);

  Status FinishInternal(std::shared_ptr<ArrayData>* out) override;

  /// \cond FALSE
  using ArrayBuilder::Finish;
  /// \endcond

  Status Finish(std::shared_ptr<StructArray>* out) { return FinishTyped(out); }

  /// Null bitmap is of equal length to every child field, and any zero byte
  /// will be considered as a null for that field, but users must using app-
  /// end methods or advance methods of the child builders' independently to
  /// insert data.
  Status AppendValues(int64_t length, const uint8_t* valid_bytes) {
    ARROW_RETURN_NOT_OK(Reserve(length));
    UnsafeAppendToBitmap(valid_bytes, length);
    return Status::OK();
  }

  /// Append an element to the Struct. All child-builders' Append method must
  /// be called independently to maintain data-structure consistency.
  Status Append(bool is_valid = true) {
    ARROW_RETURN_NOT_OK(Reserve(1));
    UnsafeAppendToBitmap(is_valid);
    return Status::OK();
  }

  /// \brief Append a null value. Automatically appends an empty value to each child
  /// builder.
  Status AppendNull() final {
    for (const auto& field : children_) {
      ARROW_RETURN_NOT_OK(field->AppendEmptyValue());
    }
    return Append(false);
  }

  /// \brief Append multiple null values. Automatically appends empty values to each
  /// child builder.
  Status AppendNulls(int64_t length) final {
    for (const auto& field : children_) {
      ARROW_RETURN_NOT_OK(field->AppendEmptyValues(length));
    }
    ARROW_RETURN_NOT_OK(Reserve(length));
    UnsafeAppendToBitmap(length, false);
    return Status::OK();
  }

  Status AppendEmptyValue() final {
    for (const auto& field : children_) {
      ARROW_RETURN_NOT_OK(field->AppendEmptyValue());
    }
    return Append(true);
  }

  Status AppendEmptyValues(int64_t length) final {
    for (const auto& field : children_) {
      ARROW_RETURN_NOT_OK(field->AppendEmptyValues(length));
    }
    ARROW_RETURN_NOT_OK(Reserve(length));
    UnsafeAppendToBitmap(length, true);
    return Status::OK();
  }

  Status AppendArraySlice(const ArraySpan& array, int64_t offset,
                          int64_t length) override {
    for (int i = 0; static_cast<size_t>(i) < children_.size(); i++) {
      ARROW_RETURN_NOT_OK(children_[i]->AppendArraySlice(array.child_data[i],
                                                         array.offset + offset, length));
    }
    const uint8_t* validity = array.MayHaveNulls() ? array.buffers[0].data : NULLPTR;
    ARROW_RETURN_NOT_OK(Reserve(length));
    UnsafeAppendToBitmap(validity, array.offset + offset, length);
    return Status::OK();
  }

  void Reset() override;

  ArrayBuilder* field_builder(int i) const { return children_[i].get(); }

  int num_fields() const { return static_cast<int>(children_.size()); }

  std::shared_ptr<DataType> type() const override;

 private:
  std::shared_ptr<DataType> type_;
};

/// @}

}  // namespace arrow
