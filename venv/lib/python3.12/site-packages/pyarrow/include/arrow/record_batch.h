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
#include <memory>
#include <string>
#include <vector>

#include "arrow/compare.h"
#include "arrow/device.h"
#include "arrow/result.h"
#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/iterator.h"
#include "arrow/util/macros.h"
#include "arrow/util/visibility.h"

namespace arrow {

/// \class RecordBatch
/// \brief Collection of equal-length arrays matching a particular Schema
///
/// A record batch is table-like data structure that is semantically a sequence
/// of fields, each a contiguous Arrow array
class ARROW_EXPORT RecordBatch {
 public:
  virtual ~RecordBatch() = default;

  /// \param[in] schema The record batch schema
  /// \param[in] num_rows length of fields in the record batch. Each array
  /// should have the same length as num_rows
  /// \param[in] columns the record batch fields as vector of arrays
  /// \param[in] sync_event optional synchronization event for non-CPU device
  /// memory used by buffers
  static std::shared_ptr<RecordBatch> Make(
      std::shared_ptr<Schema> schema, int64_t num_rows,
      std::vector<std::shared_ptr<Array>> columns,
      std::shared_ptr<Device::SyncEvent> sync_event = NULLPTR);

  /// \brief Construct record batch from vector of internal data structures
  /// \since 0.5.0
  ///
  /// This class is intended for internal use, or advanced users.
  ///
  /// \param schema the record batch schema
  /// \param num_rows the number of semantic rows in the record batch. This
  /// should be equal to the length of each field
  /// \param columns the data for the batch's columns
  /// \param device_type the type of the device that the Arrow columns are
  /// allocated on
  /// \param sync_event optional synchronization event for non-CPU device
  /// memory used by buffers
  static std::shared_ptr<RecordBatch> Make(
      std::shared_ptr<Schema> schema, int64_t num_rows,
      std::vector<std::shared_ptr<ArrayData>> columns,
      DeviceAllocationType device_type = DeviceAllocationType::kCPU,
      std::shared_ptr<Device::SyncEvent> sync_event = NULLPTR);

  /// \brief Create an empty RecordBatch of a given schema
  ///
  /// The output RecordBatch will be created with DataTypes from
  /// the given schema.
  ///
  /// \param[in] schema the schema of the empty RecordBatch
  /// \param[in] pool the memory pool to allocate memory from
  /// \return the resulting RecordBatch
  static Result<std::shared_ptr<RecordBatch>> MakeEmpty(
      std::shared_ptr<Schema> schema, MemoryPool* pool = default_memory_pool());

  /// \brief Convert record batch to struct array
  ///
  /// Create a struct array whose child arrays are the record batch's columns.
  /// Note that the record batch's top-level field metadata cannot be reflected
  /// in the resulting struct array.
  Result<std::shared_ptr<StructArray>> ToStructArray() const;

  /// \brief Convert record batch with one data type to Tensor
  ///
  /// Create a Tensor object with shape (number of rows, number of columns) and
  /// strides (type size in bytes, type size in bytes * number of rows).
  /// Generated Tensor will have column-major layout.
  ///
  /// \param[in] null_to_nan if true, convert nulls to NaN
  /// \param[in] row_major if true, create row-major Tensor else column-major Tensor
  /// \param[in] pool the memory pool to allocate the tensor buffer
  /// \return the resulting Tensor
  Result<std::shared_ptr<Tensor>> ToTensor(
      bool null_to_nan = false, bool row_major = true,
      MemoryPool* pool = default_memory_pool()) const;

  /// \brief Construct record batch from struct array
  ///
  /// This constructs a record batch using the child arrays of the given
  /// array, which must be a struct array.
  ///
  /// \param[in] array the source array, must be a StructArray
  /// \param[in] pool the memory pool to allocate new validity bitmaps
  ///
  /// This operation will usually be zero-copy.  However, if the struct array has an
  /// offset or a validity bitmap then these will need to be pushed into the child arrays.
  /// Pushing the offset is zero-copy but pushing the validity bitmap is not.
  static Result<std::shared_ptr<RecordBatch>> FromStructArray(
      const std::shared_ptr<Array>& array, MemoryPool* pool = default_memory_pool());

  /// \brief Determine if two record batches are equal
  ///
  /// \param[in] other the RecordBatch to compare with
  /// \param[in] check_metadata if true, the schema metadata will be compared,
  ///            regardless of the value set in \ref EqualOptions::use_metadata
  /// \param[in] opts the options for equality comparisons
  /// \return true if batches are equal
  bool Equals(const RecordBatch& other, bool check_metadata = false,
              const EqualOptions& opts = EqualOptions::Defaults()) const;

  /// \brief Determine if two record batches are equal
  ///
  /// \param[in] other the RecordBatch to compare with
  /// \param[in] opts the options for equality comparisons
  /// \return true if batches are equal
  bool Equals(const RecordBatch& other, const EqualOptions& opts) const;

  /// \brief Determine if two record batches are approximately equal
  ///
  /// \param[in] other the RecordBatch to compare with
  /// \param[in] opts the options for equality comparisons
  /// \return true if batches are approximately equal
  bool ApproxEquals(const RecordBatch& other,
                    const EqualOptions& opts = EqualOptions::Defaults()) const {
    return Equals(other, opts.use_schema(false).use_atol(true));
  }

  /// \return the record batch's schema
  const std::shared_ptr<Schema>& schema() const { return schema_; }

  /// \brief Replace the schema with another schema with the same types, but potentially
  /// different field names and/or metadata.
  Result<std::shared_ptr<RecordBatch>> ReplaceSchema(
      std::shared_ptr<Schema> schema) const;

  /// \brief Retrieve all columns at once
  virtual const std::vector<std::shared_ptr<Array>>& columns() const = 0;

  /// \brief Retrieve an array from the record batch
  /// \param[in] i field index, does not boundscheck
  /// \return an Array object
  virtual std::shared_ptr<Array> column(int i) const = 0;

  /// \brief Retrieve an array from the record batch
  /// \param[in] name field name
  /// \return an Array or null if no field was found
  std::shared_ptr<Array> GetColumnByName(const std::string& name) const;

  /// \brief Retrieve an array's internal data from the record batch
  /// \param[in] i field index, does not boundscheck
  /// \return an internal ArrayData object
  virtual std::shared_ptr<ArrayData> column_data(int i) const = 0;

  /// \brief Retrieve all arrays' internal data from the record batch.
  virtual const ArrayDataVector& column_data() const = 0;

  /// \brief Add column to the record batch, producing a new RecordBatch
  ///
  /// \param[in] i field index, which will be boundschecked
  /// \param[in] field field to be added
  /// \param[in] column column to be added
  virtual Result<std::shared_ptr<RecordBatch>> AddColumn(
      int i, const std::shared_ptr<Field>& field,
      const std::shared_ptr<Array>& column) const = 0;

  /// \brief Add new nullable column to the record batch, producing a new
  /// RecordBatch.
  ///
  /// For non-nullable columns, use the Field-based version of this method.
  ///
  /// \param[in] i field index, which will be boundschecked
  /// \param[in] field_name name of field to be added
  /// \param[in] column column to be added
  virtual Result<std::shared_ptr<RecordBatch>> AddColumn(
      int i, std::string field_name, const std::shared_ptr<Array>& column) const;

  /// \brief Replace a column in the record batch, producing a new RecordBatch
  ///
  /// \param[in] i field index, does boundscheck
  /// \param[in] field field to be replaced
  /// \param[in] column column to be replaced
  virtual Result<std::shared_ptr<RecordBatch>> SetColumn(
      int i, const std::shared_ptr<Field>& field,
      const std::shared_ptr<Array>& column) const = 0;

  /// \brief Remove column from the record batch, producing a new RecordBatch
  ///
  /// \param[in] i field index, does boundscheck
  virtual Result<std::shared_ptr<RecordBatch>> RemoveColumn(int i) const = 0;

  virtual std::shared_ptr<RecordBatch> ReplaceSchemaMetadata(
      const std::shared_ptr<const KeyValueMetadata>& metadata) const = 0;

  /// \brief Name in i-th column
  const std::string& column_name(int i) const;

  /// \return the number of columns in the table
  int num_columns() const;

  /// \return the number of rows (the corresponding length of each column)
  int64_t num_rows() const { return num_rows_; }

  /// \brief Copy the entire RecordBatch to destination MemoryManager
  ///
  /// This uses Array::CopyTo on each column of the record batch to create
  /// a new record batch where all underlying buffers for the columns have
  /// been copied to the destination MemoryManager. This uses
  /// MemoryManager::CopyBuffer under the hood.
  Result<std::shared_ptr<RecordBatch>> CopyTo(
      const std::shared_ptr<MemoryManager>& to) const;

  /// \brief View or Copy the entire RecordBatch to destination MemoryManager
  ///
  /// This uses Array::ViewOrCopyTo on each column of the record batch to create
  /// a new record batch where all underlying buffers for the columns have
  /// been zero-copy viewed on the destination MemoryManager, falling back
  /// to performing a copy if it can't be viewed as a zero-copy buffer. This uses
  /// Buffer::ViewOrCopy under the hood.
  Result<std::shared_ptr<RecordBatch>> ViewOrCopyTo(
      const std::shared_ptr<MemoryManager>& to) const;

  /// \brief Slice each of the arrays in the record batch
  /// \param[in] offset the starting offset to slice, through end of batch
  /// \return new record batch
  virtual std::shared_ptr<RecordBatch> Slice(int64_t offset) const;

  /// \brief Slice each of the arrays in the record batch
  /// \param[in] offset the starting offset to slice
  /// \param[in] length the number of elements to slice from offset
  /// \return new record batch
  virtual std::shared_ptr<RecordBatch> Slice(int64_t offset, int64_t length) const = 0;

  /// \return PrettyPrint representation suitable for debugging
  std::string ToString() const;

  /// \brief Return names of all columns
  std::vector<std::string> ColumnNames() const;

  /// \brief Rename columns with provided names
  Result<std::shared_ptr<RecordBatch>> RenameColumns(
      const std::vector<std::string>& names) const;

  /// \brief Return new record batch with specified columns
  Result<std::shared_ptr<RecordBatch>> SelectColumns(
      const std::vector<int>& indices) const;

  /// \brief Perform cheap validation checks to determine obvious inconsistencies
  /// within the record batch's schema and internal data.
  ///
  /// This is O(k) where k is the total number of fields and array descendents.
  ///
  /// \return Status
  virtual Status Validate() const;

  /// \brief Perform extensive validation checks to determine inconsistencies
  /// within the record batch's schema and internal data.
  ///
  /// This is potentially O(k*n) where n is the number of rows.
  ///
  /// \return Status
  virtual Status ValidateFull() const;

  /// \brief EXPERIMENTAL: Return a top-level sync event object for this record batch
  ///
  /// If all of the data for this record batch is in CPU memory, then this
  /// will return null. If the data for this batch is
  /// on a device, then if synchronization is needed before accessing the
  /// data the returned sync event will allow for it.
  ///
  /// \return null or a Device::SyncEvent
  virtual const std::shared_ptr<Device::SyncEvent>& GetSyncEvent() const = 0;

  virtual DeviceAllocationType device_type() const = 0;

  /// \brief Create a statistics array of this record batch
  ///
  /// The created array follows the C data interface statistics
  /// specification. See
  /// https://arrow.apache.org/docs/format/StatisticsSchema.html
  /// for details.
  ///
  /// \param[in] pool the memory pool to allocate memory from
  /// \return the statistics array of this record batch
  Result<std::shared_ptr<Array>> MakeStatisticsArray(
      MemoryPool* pool = default_memory_pool()) const;

 protected:
  RecordBatch(std::shared_ptr<Schema> schema, int64_t num_rows);

  std::shared_ptr<Schema> schema_;
  int64_t num_rows_;

 private:
  ARROW_DISALLOW_COPY_AND_ASSIGN(RecordBatch);
};

struct ARROW_EXPORT RecordBatchWithMetadata {
  std::shared_ptr<RecordBatch> batch;
  std::shared_ptr<KeyValueMetadata> custom_metadata;
};

template <>
struct IterationTraits<RecordBatchWithMetadata> {
  static RecordBatchWithMetadata End() { return {NULLPTR, NULLPTR}; }
  static bool IsEnd(const RecordBatchWithMetadata& val) { return val.batch == NULLPTR; }
};

/// \brief Abstract interface for reading stream of record batches
class ARROW_EXPORT RecordBatchReader {
 public:
  using ValueType = std::shared_ptr<RecordBatch>;

  virtual ~RecordBatchReader();

  /// \return the shared schema of the record batches in the stream
  virtual std::shared_ptr<Schema> schema() const = 0;

  /// \brief Read the next record batch in the stream. Return null for batch
  /// when reaching end of stream
  ///
  /// Example:
  ///
  /// ```
  /// while (true) {
  ///   std::shared_ptr<RecordBatch> batch;
  ///   ARROW_RETURN_NOT_OK(reader->ReadNext(&batch));
  ///   if (!batch) {
  ///     break;
  ///   }
  ///   // handling the `batch`, the `batch->num_rows()`
  ///   // might be 0.
  /// }
  /// ```
  ///
  /// \param[out] batch the next loaded batch, null at end of stream. Returning
  /// an empty batch doesn't mean the end of stream because it is valid data.
  /// \return Status
  virtual Status ReadNext(std::shared_ptr<RecordBatch>* batch) = 0;

  virtual Result<RecordBatchWithMetadata> ReadNext() {
    return Status::NotImplemented("ReadNext with custom metadata");
  }

  /// \brief Iterator interface
  Result<std::shared_ptr<RecordBatch>> Next() {
    std::shared_ptr<RecordBatch> batch;
    ARROW_RETURN_NOT_OK(ReadNext(&batch));
    return batch;
  }

  /// \brief finalize reader
  virtual Status Close() { return Status::OK(); }

  /// \brief EXPERIMENTAL: Get the device type for record batches this reader produces
  ///
  /// default implementation is to return DeviceAllocationType::kCPU
  virtual DeviceAllocationType device_type() const { return DeviceAllocationType::kCPU; }

  class RecordBatchReaderIterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = std::shared_ptr<RecordBatch>;
    using pointer = const value_type*;
    using reference = const value_type&;

    RecordBatchReaderIterator() : batch_(RecordBatchEnd()), reader_(NULLPTR) {}

    explicit RecordBatchReaderIterator(RecordBatchReader* reader)
        : batch_(RecordBatchEnd()), reader_(reader) {
      Next();
    }

    bool operator==(const RecordBatchReaderIterator& other) const {
      return batch_ == other.batch_;
    }

    bool operator!=(const RecordBatchReaderIterator& other) const {
      return !(*this == other);
    }

    Result<std::shared_ptr<RecordBatch>> operator*() {
      ARROW_RETURN_NOT_OK(batch_);

      return batch_;
    }

    RecordBatchReaderIterator& operator++() {
      Next();
      return *this;
    }

    RecordBatchReaderIterator operator++(int) {
      RecordBatchReaderIterator tmp(*this);
      Next();
      return tmp;
    }

   private:
    std::shared_ptr<RecordBatch> RecordBatchEnd() {
      return std::shared_ptr<RecordBatch>(NULLPTR);
    }

    void Next() {
      if (reader_ == NULLPTR) {
        batch_ = RecordBatchEnd();
        return;
      }
      batch_ = reader_->Next();
    }

    Result<std::shared_ptr<RecordBatch>> batch_;
    RecordBatchReader* reader_;
  };
  /// \brief Return an iterator to the first record batch in the stream
  RecordBatchReaderIterator begin() { return RecordBatchReaderIterator(this); }

  /// \brief Return an iterator to the end of the stream
  RecordBatchReaderIterator end() { return RecordBatchReaderIterator(); }

  /// \brief Consume entire stream as a vector of record batches
  Result<RecordBatchVector> ToRecordBatches();

  /// \brief Read all batches and concatenate as arrow::Table
  Result<std::shared_ptr<Table>> ToTable();

  /// \brief Create a RecordBatchReader from a vector of RecordBatch.
  ///
  /// \param[in] batches the vector of RecordBatch to read from
  /// \param[in] schema schema to conform to. Will be inferred from the first
  ///            element if not provided.
  /// \param[in] device_type the type of device that the batches are allocated on
  static Result<std::shared_ptr<RecordBatchReader>> Make(
      RecordBatchVector batches, std::shared_ptr<Schema> schema = NULLPTR,
      DeviceAllocationType device_type = DeviceAllocationType::kCPU);

  /// \brief Create a RecordBatchReader from an Iterator of RecordBatch.
  ///
  /// \param[in] batches an iterator of RecordBatch to read from.
  /// \param[in] schema schema that each record batch in iterator will conform to.
  /// \param[in] device_type the type of device that the batches are allocated on
  static Result<std::shared_ptr<RecordBatchReader>> MakeFromIterator(
      Iterator<std::shared_ptr<RecordBatch>> batches, std::shared_ptr<Schema> schema,
      DeviceAllocationType device_type = DeviceAllocationType::kCPU);
};

/// \brief Concatenate record batches
///
/// The columns of the new batch are formed by concatenate the same columns of each input
/// batch. Concatenate multiple batches into a new batch requires that the schema must be
/// consistent. It supports merging batches without columns (only length, scenarios such
/// as count(*)).
///
/// \param[in] batches a vector of record batches to be concatenated
/// \param[in] pool memory to store the result will be allocated from this memory pool
/// \return the concatenated record batch
ARROW_EXPORT
Result<std::shared_ptr<RecordBatch>> ConcatenateRecordBatches(
    const RecordBatchVector& batches, MemoryPool* pool = default_memory_pool());

}  // namespace arrow
