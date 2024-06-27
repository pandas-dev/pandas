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

#include <optional>
#include <vector>
#include "arrow/array/builder_base.h"
#include "arrow/array/builder_binary.h"
#include "arrow/array/builder_primitive.h"
#include "arrow/memory_pool.h"
#include "arrow/record_batch.h"
#include "arrow/type_traits.h"
#include "arrow/util/logging.h"

namespace arrow::acero {

/// Lightweight representation of a cell of an unmaterialized table.
///
struct CompositeEntry {
  RecordBatch* batch;
  uint64_t start;
  uint64_t end;
};

// Forward declare the builder
template <size_t MAX_COMPOSITE_TABLES>
class UnmaterializedSliceBuilder;

/// A table of composite reference rows.  Rows maintain pointers to the
/// constituent record batches, but the overall table retains shared_ptr
/// references to ensure memory remains resident while the table is live.
///
/// The main reason for this is that, especially for wide tables, some operations
/// such as sorted_merge or asof_join are effectively row-oriented, rather than
/// column-oriented.  Separating the join part from the columnar materialization
/// part simplifies the logic around data types and increases efficiency.
///
/// We don't put the shared_ptr's into the rows for efficiency reasons. Use
/// UnmaterializedSliceBuilder to add ranges of record batches to this table
template <size_t MAX_COMPOSITE_TABLES>
class UnmaterializedCompositeTable {
 public:
  UnmaterializedCompositeTable(
      const std::shared_ptr<arrow::Schema>& output_schema, size_t num_composite_tables,
      std::unordered_map<int, std::pair<int, int>> output_col_to_src_,
      arrow::MemoryPool* pool_ = arrow::default_memory_pool())
      : schema(output_schema),
        num_composite_tables(num_composite_tables),
        output_col_to_src(std::move(output_col_to_src_)),
        pool{pool_} {}

  // Shallow wrappers around std::vector for performance
  inline size_t capacity() { return slices.capacity(); }
  inline void reserve(size_t num_slices) { slices.reserve(num_slices); }

  inline size_t Size() const { return num_rows; }
  inline size_t Empty() const { return num_rows == 0; }

  Result<std::optional<std::shared_ptr<RecordBatch>>> Materialize() {
    // Don't build empty batches
    if (Empty()) {
      return std::nullopt;
    }
    DCHECK_LE(Size(), (uint64_t)std::numeric_limits<int64_t>::max());
    std::vector<std::shared_ptr<arrow::Array>> arrays(schema->num_fields());

#define MATERIALIZE_CASE(id)                                                          \
  case arrow::Type::id: {                                                             \
    using T = typename arrow::TypeIdTraits<arrow::Type::id>::Type;                    \
    ARROW_ASSIGN_OR_RAISE(arrays.at(i_col), materializeColumn<T>(field_type, i_col)); \
    break;                                                                            \
  }

    // Build the arrays column-by-column from the rows
    for (int i_col = 0; i_col < schema->num_fields(); ++i_col) {
      const std::shared_ptr<arrow::Field>& field = schema->field(i_col);
      const auto& field_type = field->type();

      switch (field_type->id()) {
        MATERIALIZE_CASE(BOOL)
        MATERIALIZE_CASE(INT8)
        MATERIALIZE_CASE(INT16)
        MATERIALIZE_CASE(INT32)
        MATERIALIZE_CASE(INT64)
        MATERIALIZE_CASE(UINT8)
        MATERIALIZE_CASE(UINT16)
        MATERIALIZE_CASE(UINT32)
        MATERIALIZE_CASE(UINT64)
        MATERIALIZE_CASE(FLOAT)
        MATERIALIZE_CASE(DOUBLE)
        MATERIALIZE_CASE(DATE32)
        MATERIALIZE_CASE(DATE64)
        MATERIALIZE_CASE(TIME32)
        MATERIALIZE_CASE(TIME64)
        MATERIALIZE_CASE(TIMESTAMP)
        MATERIALIZE_CASE(STRING)
        MATERIALIZE_CASE(LARGE_STRING)
        MATERIALIZE_CASE(BINARY)
        MATERIALIZE_CASE(LARGE_BINARY)
        default:
          return arrow::Status::Invalid("Unsupported data type ",
                                        field->type()->ToString(), " for field ",
                                        field->name());
      }
    }

#undef MATERIALIZE_CASE

    std::shared_ptr<arrow::RecordBatch> r =
        arrow::RecordBatch::Make(schema, (int64_t)num_rows, arrays);
    return r;
  }

 private:
  struct UnmaterializedSlice {
    CompositeEntry components[MAX_COMPOSITE_TABLES];
    size_t num_components;

    inline int64_t Size() const {
      if (num_components == 0) {
        return 0;
      }
      return components[0].end - components[0].start;
    }
  };

  // Mapping from an output column ID to a source table ID and column ID
  std::shared_ptr<arrow::Schema> schema;
  size_t num_composite_tables;
  std::unordered_map<int, std::pair<int, int>> output_col_to_src;

  arrow::MemoryPool* pool;

  /// A map from address of a record batch to the record batch. Used to
  /// maintain the lifetime of the record batch in case it goes out of scope
  /// by the main exec node thread
  std::unordered_map<uintptr_t, std::shared_ptr<arrow::RecordBatch>> ptr2Ref = {};
  std::vector<UnmaterializedSlice> slices;

  size_t num_rows = 0;

  // for AddRecordBatchRef/AddSlice and access to UnmaterializedSlice
  friend class UnmaterializedSliceBuilder<MAX_COMPOSITE_TABLES>;

  void AddRecordBatchRef(const std::shared_ptr<arrow::RecordBatch>& ref) {
    ptr2Ref[(uintptr_t)ref.get()] = ref;
  }
  void AddSlice(const UnmaterializedSlice& slice) {
    slices.push_back(slice);
    num_rows += slice.Size();
  }

  template <class Type, class Builder = typename TypeTraits<Type>::BuilderType>
  enable_if_boolean<Type, Status> static BuilderAppend(
      Builder& builder, const std::shared_ptr<ArrayData>& source, uint64_t row) {
    if (source->IsNull(row)) {
      builder.UnsafeAppendNull();
      return Status::OK();
    }
    builder.UnsafeAppend(bit_util::GetBit(source->template GetValues<uint8_t>(1), row));
    return Status::OK();
  }

  template <class Type, class Builder = typename TypeTraits<Type>::BuilderType>
  enable_if_t<is_fixed_width_type<Type>::value && !is_boolean_type<Type>::value,
              Status> static BuilderAppend(Builder& builder,
                                           const std::shared_ptr<ArrayData>& source,
                                           uint64_t row) {
    if (source->IsNull(row)) {
      builder.UnsafeAppendNull();
      return Status::OK();
    }
    using CType = typename TypeTraits<Type>::CType;
    builder.UnsafeAppend(source->template GetValues<CType>(1)[row]);
    return Status::OK();
  }

  template <class Type, class Builder = typename TypeTraits<Type>::BuilderType>
  enable_if_base_binary<Type, Status> static BuilderAppend(
      Builder& builder, const std::shared_ptr<ArrayData>& source, uint64_t row) {
    if (source->IsNull(row)) {
      return builder.AppendNull();
    }
    using offset_type = typename Type::offset_type;
    const uint8_t* data = source->buffers[2]->data();
    const offset_type* offsets = source->GetValues<offset_type>(1);
    const offset_type offset0 = offsets[row];
    const offset_type offset1 = offsets[row + 1];
    return builder.Append(data + offset0, offset1 - offset0);
  }

  template <class Type, class Builder = typename arrow::TypeTraits<Type>::BuilderType>
  arrow::Result<std::shared_ptr<arrow::Array>> materializeColumn(
      const std::shared_ptr<arrow::DataType>& type, int i_col) {
    ARROW_ASSIGN_OR_RAISE(auto builderPtr, arrow::MakeBuilder(type, pool));
    Builder& builder = *arrow::internal::checked_cast<Builder*>(builderPtr.get());
    ARROW_RETURN_NOT_OK(builder.Reserve(num_rows));

    const auto& [table_index, column_index] = output_col_to_src[i_col];

    for (const auto& unmaterialized_slice : slices) {
      const auto& [batch, start, end] = unmaterialized_slice.components[table_index];
      if (batch) {
        for (uint64_t rowNum = start; rowNum < end; ++rowNum) {
          arrow::Status st = BuilderAppend<Type, Builder>(
              builder, batch->column_data(column_index), rowNum);
          ARROW_RETURN_NOT_OK(st);
        }
      } else {
        for (uint64_t rowNum = start; rowNum < end; ++rowNum) {
          ARROW_RETURN_NOT_OK(builder.AppendNull());
        }
      }
    }
    std::shared_ptr<arrow::Array> result;
    ARROW_RETURN_NOT_OK(builder.Finish(&result));
    return Result{std::move(result)};
  }
};

/// A builder class that can append blocks of data to a row. A "slice"
/// is built by horizontally concatenating record batches.
template <size_t MAX_COMPOSITE_TABLES>
class UnmaterializedSliceBuilder {
 public:
  explicit UnmaterializedSliceBuilder(
      UnmaterializedCompositeTable<MAX_COMPOSITE_TABLES>* table_)
      : table(table_) {}

  void AddEntry(std::shared_ptr<RecordBatch> rb, uint64_t start, uint64_t end) {
    if (rb) {
      table->AddRecordBatchRef(rb);
    }
    if (slice.num_components) {
      size_t last_index = slice.num_components - 1;
      DCHECK_EQ(slice.components[last_index].end - slice.components[last_index].start,
                end - start)
          << "Slices should be the same length. ";
    }
    slice.components[slice.num_components++] = CompositeEntry{rb.get(), start, end};
  }

  void Finalize() { table->AddSlice(slice); }
  int64_t Size() { return slice.Size(); }

 private:
  using TUnmaterializedCompositeTable =
      UnmaterializedCompositeTable<MAX_COMPOSITE_TABLES>;
  using TUnmaterializedSlice =
      typename TUnmaterializedCompositeTable::UnmaterializedSlice;

  TUnmaterializedCompositeTable* table;
  TUnmaterializedSlice slice{};
};

}  // namespace arrow::acero
