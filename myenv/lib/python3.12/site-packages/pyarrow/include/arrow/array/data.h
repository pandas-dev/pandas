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

#include <atomic>  // IWYU pragma: export
#include <cassert>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "arrow/buffer.h"
#include "arrow/result.h"
#include "arrow/type.h"
#include "arrow/type_fwd.h"
#include "arrow/util/bit_util.h"
#include "arrow/util/macros.h"
#include "arrow/util/span.h"
#include "arrow/util/visibility.h"

namespace arrow {

namespace internal {
// ----------------------------------------------------------------------
// Null handling for types without a validity bitmap and the dictionary type

ARROW_EXPORT bool IsNullSparseUnion(const ArrayData& data, int64_t i);
ARROW_EXPORT bool IsNullDenseUnion(const ArrayData& data, int64_t i);
ARROW_EXPORT bool IsNullRunEndEncoded(const ArrayData& data, int64_t i);

ARROW_EXPORT bool UnionMayHaveLogicalNulls(const ArrayData& data);
ARROW_EXPORT bool RunEndEncodedMayHaveLogicalNulls(const ArrayData& data);
ARROW_EXPORT bool DictionaryMayHaveLogicalNulls(const ArrayData& data);
}  // namespace internal

// When slicing, we do not know the null count of the sliced range without
// doing some computation. To avoid doing this eagerly, we set the null count
// to -1 (any negative number will do). When Array::null_count is called the
// first time, the null count will be computed. See ARROW-33
constexpr int64_t kUnknownNullCount = -1;

// ----------------------------------------------------------------------
// Generic array data container

/// \class ArrayData
/// \brief Mutable container for generic Arrow array data
///
/// This data structure is a self-contained representation of the memory and
/// metadata inside an Arrow array data structure (called vectors in Java). The
/// classes arrow::Array and its subclasses provide strongly-typed accessors
/// with support for the visitor pattern and other affordances.
///
/// This class is designed for easy internal data manipulation, analytical data
/// processing, and data transport to and from IPC messages. For example, we
/// could cast from int64 to float64 like so:
///
/// Int64Array arr = GetMyData();
/// auto new_data = arr.data()->Copy();
/// new_data->type = arrow::float64();
/// DoubleArray double_arr(new_data);
///
/// This object is also useful in an analytics setting where memory may be
/// reused. For example, if we had a group of operations all returning doubles,
/// say:
///
/// Log(Sqrt(Expr(arr)))
///
/// Then the low-level implementations of each of these functions could have
/// the signatures
///
/// void Log(const ArrayData& values, ArrayData* out);
///
/// As another example a function may consume one or more memory buffers in an
/// input array and replace them with newly-allocated data, changing the output
/// data type as well.
struct ARROW_EXPORT ArrayData {
  ArrayData() = default;

  ArrayData(std::shared_ptr<DataType> type, int64_t length,
            int64_t null_count = kUnknownNullCount, int64_t offset = 0)
      : type(std::move(type)), length(length), null_count(null_count), offset(offset) {}

  ArrayData(std::shared_ptr<DataType> type, int64_t length,
            std::vector<std::shared_ptr<Buffer>> buffers,
            int64_t null_count = kUnknownNullCount, int64_t offset = 0)
      : ArrayData(std::move(type), length, null_count, offset) {
    this->buffers = std::move(buffers);
  }

  ArrayData(std::shared_ptr<DataType> type, int64_t length,
            std::vector<std::shared_ptr<Buffer>> buffers,
            std::vector<std::shared_ptr<ArrayData>> child_data,
            int64_t null_count = kUnknownNullCount, int64_t offset = 0)
      : ArrayData(std::move(type), length, null_count, offset) {
    this->buffers = std::move(buffers);
    this->child_data = std::move(child_data);
  }

  static std::shared_ptr<ArrayData> Make(std::shared_ptr<DataType> type, int64_t length,
                                         std::vector<std::shared_ptr<Buffer>> buffers,
                                         int64_t null_count = kUnknownNullCount,
                                         int64_t offset = 0);

  static std::shared_ptr<ArrayData> Make(
      std::shared_ptr<DataType> type, int64_t length,
      std::vector<std::shared_ptr<Buffer>> buffers,
      std::vector<std::shared_ptr<ArrayData>> child_data,
      int64_t null_count = kUnknownNullCount, int64_t offset = 0);

  static std::shared_ptr<ArrayData> Make(
      std::shared_ptr<DataType> type, int64_t length,
      std::vector<std::shared_ptr<Buffer>> buffers,
      std::vector<std::shared_ptr<ArrayData>> child_data,
      std::shared_ptr<ArrayData> dictionary, int64_t null_count = kUnknownNullCount,
      int64_t offset = 0);

  static std::shared_ptr<ArrayData> Make(std::shared_ptr<DataType> type, int64_t length,
                                         int64_t null_count = kUnknownNullCount,
                                         int64_t offset = 0);

  // Move constructor
  ArrayData(ArrayData&& other) noexcept
      : type(std::move(other.type)),
        length(other.length),
        offset(other.offset),
        buffers(std::move(other.buffers)),
        child_data(std::move(other.child_data)),
        dictionary(std::move(other.dictionary)) {
    SetNullCount(other.null_count);
  }

  // Copy constructor
  ArrayData(const ArrayData& other) noexcept
      : type(other.type),
        length(other.length),
        offset(other.offset),
        buffers(other.buffers),
        child_data(other.child_data),
        dictionary(other.dictionary) {
    SetNullCount(other.null_count);
  }

  // Move assignment
  ArrayData& operator=(ArrayData&& other) {
    type = std::move(other.type);
    length = other.length;
    SetNullCount(other.null_count);
    offset = other.offset;
    buffers = std::move(other.buffers);
    child_data = std::move(other.child_data);
    dictionary = std::move(other.dictionary);
    return *this;
  }

  // Copy assignment
  ArrayData& operator=(const ArrayData& other) {
    type = other.type;
    length = other.length;
    SetNullCount(other.null_count);
    offset = other.offset;
    buffers = other.buffers;
    child_data = other.child_data;
    dictionary = other.dictionary;
    return *this;
  }

  std::shared_ptr<ArrayData> Copy() const { return std::make_shared<ArrayData>(*this); }

  /// \brief Copy all buffers and children recursively to destination MemoryManager
  ///
  /// This utilizes MemoryManager::CopyBuffer to create a new ArrayData object
  /// recursively copying the buffers and all child buffers to the destination
  /// memory manager. This includes dictionaries if applicable.
  Result<std::shared_ptr<ArrayData>> CopyTo(
      const std::shared_ptr<MemoryManager>& to) const;
  /// \brief View or Copy this ArrayData to destination memory manager.
  ///
  /// Tries to view the buffer contents on the given memory manager's device
  /// if possible (to avoid a copy) but falls back to copying if a no-copy view
  /// isn't supported.
  Result<std::shared_ptr<ArrayData>> ViewOrCopyTo(
      const std::shared_ptr<MemoryManager>& to) const;

  bool IsNull(int64_t i) const { return !IsValid(i); }

  bool IsValid(int64_t i) const {
    if (buffers[0] != NULLPTR) {
      return bit_util::GetBit(buffers[0]->data(), i + offset);
    }
    const auto type = this->type->id();
    if (type == Type::SPARSE_UNION) {
      return !internal::IsNullSparseUnion(*this, i);
    }
    if (type == Type::DENSE_UNION) {
      return !internal::IsNullDenseUnion(*this, i);
    }
    if (type == Type::RUN_END_ENCODED) {
      return !internal::IsNullRunEndEncoded(*this, i);
    }
    return null_count.load() != length;
  }

  // Access a buffer's data as a typed C pointer
  template <typename T>
  inline const T* GetValues(int i, int64_t absolute_offset) const {
    if (buffers[i]) {
      return reinterpret_cast<const T*>(buffers[i]->data()) + absolute_offset;
    } else {
      return NULLPTR;
    }
  }

  template <typename T>
  inline const T* GetValues(int i) const {
    return GetValues<T>(i, offset);
  }

  // Like GetValues, but returns NULLPTR instead of aborting if the underlying
  // buffer is not a CPU buffer.
  template <typename T>
  inline const T* GetValuesSafe(int i, int64_t absolute_offset) const {
    if (buffers[i] && buffers[i]->is_cpu()) {
      return reinterpret_cast<const T*>(buffers[i]->data()) + absolute_offset;
    } else {
      return NULLPTR;
    }
  }

  template <typename T>
  inline const T* GetValuesSafe(int i) const {
    return GetValuesSafe<T>(i, offset);
  }

  // Access a buffer's data as a typed C pointer
  template <typename T>
  inline T* GetMutableValues(int i, int64_t absolute_offset) {
    if (buffers[i]) {
      return reinterpret_cast<T*>(buffers[i]->mutable_data()) + absolute_offset;
    } else {
      return NULLPTR;
    }
  }

  template <typename T>
  inline T* GetMutableValues(int i) {
    return GetMutableValues<T>(i, offset);
  }

  /// \brief Construct a zero-copy slice of the data with the given offset and length
  std::shared_ptr<ArrayData> Slice(int64_t offset, int64_t length) const;

  /// \brief Input-checking variant of Slice
  ///
  /// An Invalid Status is returned if the requested slice falls out of bounds.
  /// Note that unlike Slice, `length` isn't clamped to the available buffer size.
  Result<std::shared_ptr<ArrayData>> SliceSafe(int64_t offset, int64_t length) const;

  void SetNullCount(int64_t v) { null_count.store(v); }

  /// \brief Return physical null count, or compute and set it if it's not known
  int64_t GetNullCount() const;

  /// \brief Return true if the data has a validity bitmap and the physical null
  /// count is known to be non-zero or not yet known.
  ///
  /// Note that this is not the same as MayHaveLogicalNulls, which also checks
  /// for the presence of nulls in child data for types like unions and run-end
  /// encoded types.
  ///
  /// \see HasValidityBitmap
  /// \see MayHaveLogicalNulls
  bool MayHaveNulls() const {
    // If an ArrayData is slightly malformed it may have kUnknownNullCount set
    // but no buffer
    return null_count.load() != 0 && buffers[0] != NULLPTR;
  }

  /// \brief Return true if the data has a validity bitmap
  bool HasValidityBitmap() const { return buffers[0] != NULLPTR; }

  /// \brief Return true if the validity bitmap may have 0's in it, or if the
  /// child arrays (in the case of types without a validity bitmap) may have
  /// nulls, or if the dictionary of dictionay array may have nulls.
  ///
  /// This is not a drop-in replacement for MayHaveNulls, as historically
  /// MayHaveNulls() has been used to check for the presence of a validity
  /// bitmap that needs to be checked.
  ///
  /// Code that previously used MayHaveNulls() and then dealt with the validity
  /// bitmap directly can be fixed to handle all types correctly without
  /// performance degradation when handling most types by adopting
  /// HasValidityBitmap and MayHaveLogicalNulls.
  ///
  /// Before:
  ///
  ///     uint8_t* validity = array.MayHaveNulls() ? array.buffers[0].data : NULLPTR;
  ///     for (int64_t i = 0; i < array.length; ++i) {
  ///       if (validity && !bit_util::GetBit(validity, i)) {
  ///         continue;  // skip a NULL
  ///       }
  ///       ...
  ///     }
  ///
  /// After:
  ///
  ///     bool all_valid = !array.MayHaveLogicalNulls();
  ///     uint8_t* validity = array.HasValidityBitmap() ? array.buffers[0].data : NULLPTR;
  ///     for (int64_t i = 0; i < array.length; ++i) {
  ///       bool is_valid = all_valid ||
  ///                       (validity && bit_util::GetBit(validity, i)) ||
  ///                       array.IsValid(i);
  ///       if (!is_valid) {
  ///         continue;  // skip a NULL
  ///       }
  ///       ...
  ///     }
  bool MayHaveLogicalNulls() const {
    if (buffers[0] != NULLPTR) {
      return null_count.load() != 0;
    }
    const auto t = type->id();
    if (t == Type::SPARSE_UNION || t == Type::DENSE_UNION) {
      return internal::UnionMayHaveLogicalNulls(*this);
    }
    if (t == Type::RUN_END_ENCODED) {
      return internal::RunEndEncodedMayHaveLogicalNulls(*this);
    }
    if (t == Type::DICTIONARY) {
      return internal::DictionaryMayHaveLogicalNulls(*this);
    }
    return null_count.load() != 0;
  }

  /// \brief Computes the logical null count for arrays of all types including
  /// those that do not have a validity bitmap like union and run-end encoded
  /// arrays
  ///
  /// If the array has a validity bitmap, this function behaves the same as
  /// GetNullCount. For types that have no validity bitmap, this function will
  /// recompute the null count every time it is called.
  ///
  /// \see GetNullCount
  int64_t ComputeLogicalNullCount() const;

  std::shared_ptr<DataType> type;
  int64_t length = 0;
  mutable std::atomic<int64_t> null_count{0};
  // The logical start point into the physical buffers (in values, not bytes).
  // Note that, for child data, this must be *added* to the child data's own offset.
  int64_t offset = 0;
  std::vector<std::shared_ptr<Buffer>> buffers;
  std::vector<std::shared_ptr<ArrayData>> child_data;

  // The dictionary for this Array, if any. Only used for dictionary type
  std::shared_ptr<ArrayData> dictionary;
};

/// \brief A non-owning Buffer reference
struct ARROW_EXPORT BufferSpan {
  // It is the user of this class's responsibility to ensure that
  // buffers that were const originally are not written to
  // accidentally.
  uint8_t* data = NULLPTR;
  int64_t size = 0;
  // Pointer back to buffer that owns this memory
  const std::shared_ptr<Buffer>* owner = NULLPTR;

  template <typename T>
  const T* data_as() const {
    return reinterpret_cast<const T*>(data);
  }
  template <typename T>
  T* mutable_data_as() {
    return reinterpret_cast<T*>(data);
  }
};

/// \brief EXPERIMENTAL: A non-owning ArrayData reference that is cheaply
/// copyable and does not contain any shared_ptr objects. Do not use in public
/// APIs aside from compute kernels for now
struct ARROW_EXPORT ArraySpan {
  const DataType* type = NULLPTR;
  int64_t length = 0;
  mutable int64_t null_count = kUnknownNullCount;
  int64_t offset = 0;
  BufferSpan buffers[3];

  ArraySpan() = default;

  explicit ArraySpan(const DataType* type, int64_t length) : type(type), length(length) {}

  ArraySpan(const ArrayData& data) {  // NOLINT implicit conversion
    SetMembers(data);
  }
  explicit ArraySpan(const Scalar& data) { FillFromScalar(data); }

  /// If dictionary-encoded, put dictionary in the first entry
  std::vector<ArraySpan> child_data;

  /// \brief Populate ArraySpan to look like an array of length 1 pointing at
  /// the data members of a Scalar value
  void FillFromScalar(const Scalar& value);

  void SetMembers(const ArrayData& data);

  void SetBuffer(int index, const std::shared_ptr<Buffer>& buffer) {
    this->buffers[index].data = const_cast<uint8_t*>(buffer->data());
    this->buffers[index].size = buffer->size();
    this->buffers[index].owner = &buffer;
  }

  const ArraySpan& dictionary() const { return child_data[0]; }

  /// \brief Return the number of buffers (out of 3) that are used to
  /// constitute this array
  int num_buffers() const;

  // Access a buffer's data as a typed C pointer
  template <typename T>
  inline T* GetValues(int i, int64_t absolute_offset) {
    return reinterpret_cast<T*>(buffers[i].data) + absolute_offset;
  }

  template <typename T>
  inline T* GetValues(int i) {
    return GetValues<T>(i, this->offset);
  }

  // Access a buffer's data as a typed C pointer
  template <typename T>
  inline const T* GetValues(int i, int64_t absolute_offset) const {
    return reinterpret_cast<const T*>(buffers[i].data) + absolute_offset;
  }

  template <typename T>
  inline const T* GetValues(int i) const {
    return GetValues<T>(i, this->offset);
  }

  /// \brief Access a buffer's data as a span
  ///
  /// \param i The buffer index
  /// \param length The required length (in number of typed values) of the requested span
  /// \pre i > 0
  /// \pre length <= the length of the buffer (in number of values) that's expected for
  /// this array type
  /// \return A span<const T> of the requested length
  template <typename T>
  util::span<const T> GetSpan(int i, int64_t length) const {
    const int64_t buffer_length = buffers[i].size / static_cast<int64_t>(sizeof(T));
    assert(i > 0 && length + offset <= buffer_length);
    ARROW_UNUSED(buffer_length);
    return util::span<const T>(buffers[i].data_as<T>() + this->offset, length);
  }

  /// \brief Access a buffer's data as a span
  ///
  /// \param i The buffer index
  /// \param length The required length (in number of typed values) of the requested span
  /// \pre i > 0
  /// \pre length <= the length of the buffer (in number of values) that's expected for
  /// this array type
  /// \return A span<T> of the requested length
  template <typename T>
  util::span<T> GetSpan(int i, int64_t length) {
    const int64_t buffer_length = buffers[i].size / static_cast<int64_t>(sizeof(T));
    assert(i > 0 && length + offset <= buffer_length);
    ARROW_UNUSED(buffer_length);
    return util::span<T>(buffers[i].mutable_data_as<T>() + this->offset, length);
  }

  inline bool IsNull(int64_t i) const { return !IsValid(i); }

  inline bool IsValid(int64_t i) const {
    if (this->buffers[0].data != NULLPTR) {
      return bit_util::GetBit(this->buffers[0].data, i + this->offset);
    } else {
      const auto type = this->type->id();
      if (type == Type::SPARSE_UNION) {
        return !IsNullSparseUnion(i);
      }
      if (type == Type::DENSE_UNION) {
        return !IsNullDenseUnion(i);
      }
      if (type == Type::RUN_END_ENCODED) {
        return !IsNullRunEndEncoded(i);
      }
      return this->null_count != this->length;
    }
  }

  std::shared_ptr<ArrayData> ToArrayData() const;

  std::shared_ptr<Array> ToArray() const;

  std::shared_ptr<Buffer> GetBuffer(int index) const {
    const BufferSpan& buf = this->buffers[index];
    if (buf.owner) {
      return *buf.owner;
    } else if (buf.data != NULLPTR) {
      // Buffer points to some memory without an owning buffer
      return std::make_shared<Buffer>(buf.data, buf.size);
    } else {
      return NULLPTR;
    }
  }

  void SetSlice(int64_t offset, int64_t length) {
    this->offset = offset;
    this->length = length;
    if (this->type->id() == Type::NA) {
      this->null_count = this->length;
    } else if (this->MayHaveNulls()) {
      this->null_count = kUnknownNullCount;
    } else {
      this->null_count = 0;
    }
  }

  /// \brief Return physical null count, or compute and set it if it's not known
  int64_t GetNullCount() const;

  /// \brief Return true if the array has a validity bitmap and the physical null
  /// count is known to be non-zero or not yet known
  ///
  /// Note that this is not the same as MayHaveLogicalNulls, which also checks
  /// for the presence of nulls in child data for types like unions and run-end
  /// encoded types.
  ///
  /// \see HasValidityBitmap
  /// \see MayHaveLogicalNulls
  bool MayHaveNulls() const {
    // If an ArrayData is slightly malformed it may have kUnknownNullCount set
    // but no buffer
    return null_count != 0 && buffers[0].data != NULLPTR;
  }

  /// \brief Return true if the array has a validity bitmap
  bool HasValidityBitmap() const { return buffers[0].data != NULLPTR; }

  /// \brief Return true if the validity bitmap may have 0's in it, or if the
  /// child arrays (in the case of types without a validity bitmap) may have
  /// nulls, or if the dictionary of dictionay array may have nulls.
  ///
  /// \see ArrayData::MayHaveLogicalNulls
  bool MayHaveLogicalNulls() const {
    if (buffers[0].data != NULLPTR) {
      return null_count != 0;
    }
    const auto t = type->id();
    if (t == Type::SPARSE_UNION || t == Type::DENSE_UNION) {
      return UnionMayHaveLogicalNulls();
    }
    if (t == Type::RUN_END_ENCODED) {
      return RunEndEncodedMayHaveLogicalNulls();
    }
    if (t == Type::DICTIONARY) {
      return DictionaryMayHaveLogicalNulls();
    }
    return null_count != 0;
  }

  /// \brief Compute the logical null count for arrays of all types including
  /// those that do not have a validity bitmap like union and run-end encoded
  /// arrays
  ///
  /// If the array has a validity bitmap, this function behaves the same as
  /// GetNullCount. For types that have no validity bitmap, this function will
  /// recompute the logical null count every time it is called.
  ///
  /// \see GetNullCount
  int64_t ComputeLogicalNullCount() const;

  /// Some DataTypes (StringView, BinaryView) may have an arbitrary number of variadic
  /// buffers. Since ArraySpan only has 3 buffers, we pack the variadic buffers into
  /// buffers[2]; IE buffers[2].data points to the first shared_ptr<Buffer> of the
  /// variadic set and buffers[2].size is the number of variadic buffers times
  /// sizeof(shared_ptr<Buffer>).
  ///
  /// \see HasVariadicBuffers
  util::span<const std::shared_ptr<Buffer>> GetVariadicBuffers() const;
  bool HasVariadicBuffers() const;

 private:
  ARROW_FRIEND_EXPORT friend bool internal::IsNullRunEndEncoded(const ArrayData& span,
                                                                int64_t i);

  bool IsNullSparseUnion(int64_t i) const;
  bool IsNullDenseUnion(int64_t i) const;

  /// \brief Return true if the value at logical index i is null
  ///
  /// This function uses binary-search, so it has a O(log N) cost.
  /// Iterating over the whole array and calling IsNull is O(N log N), so
  /// for better performance it is recommended to use a
  /// ree_util::RunEndEncodedArraySpan to iterate run by run instead.
  bool IsNullRunEndEncoded(int64_t i) const;

  bool UnionMayHaveLogicalNulls() const;
  bool RunEndEncodedMayHaveLogicalNulls() const;
  bool DictionaryMayHaveLogicalNulls() const;
};

namespace internal {

void FillZeroLengthArray(const DataType* type, ArraySpan* span);

/// Construct a zero-copy view of this ArrayData with the given type.
///
/// This method checks if the types are layout-compatible.
/// Nested types are traversed in depth-first order. Data buffers must have
/// the same item sizes, even though the logical types may be different.
/// An error is returned if the types are not layout-compatible.
ARROW_EXPORT
Result<std::shared_ptr<ArrayData>> GetArrayView(const std::shared_ptr<ArrayData>& data,
                                                const std::shared_ptr<DataType>& type);

}  // namespace internal
}  // namespace arrow
