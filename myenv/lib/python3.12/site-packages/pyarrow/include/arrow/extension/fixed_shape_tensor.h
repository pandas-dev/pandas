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

#include "arrow/extension_type.h"

namespace arrow {
namespace extension {

class ARROW_EXPORT FixedShapeTensorArray : public ExtensionArray {
 public:
  using ExtensionArray::ExtensionArray;

  /// \brief Create a FixedShapeTensorArray from a Tensor
  ///
  /// This method will create a FixedShapeTensorArray from a Tensor, taking its first
  /// dimension as the number of elements in the resulting array and the remaining
  /// dimensions as the shape of the individual tensors. If Tensor provides strides,
  /// they will be used to determine dimension permutation. Otherwise, row-major layout
  /// (i.e. no permutation) will be assumed.
  ///
  /// \param[in] tensor The Tensor to convert to a FixedShapeTensorArray
  static Result<std::shared_ptr<FixedShapeTensorArray>> FromTensor(
      const std::shared_ptr<Tensor>& tensor);

  /// \brief Create a Tensor from FixedShapeTensorArray
  ///
  /// This method will create a Tensor from a FixedShapeTensorArray, setting its first
  /// dimension as length equal to the FixedShapeTensorArray's length and the remaining
  /// dimensions as the FixedShapeTensorType's shape. Shape and dim_names will be
  /// permuted according to permutation stored in the FixedShapeTensorType metadata.
  const Result<std::shared_ptr<Tensor>> ToTensor() const;
};

/// \brief Concrete type class for constant-size Tensor data.
/// This is a canonical arrow extension type.
/// See: https://arrow.apache.org/docs/format/CanonicalExtensions.html
class ARROW_EXPORT FixedShapeTensorType : public ExtensionType {
 public:
  FixedShapeTensorType(const std::shared_ptr<DataType>& value_type, const int32_t& size,
                       const std::vector<int64_t>& shape,
                       const std::vector<int64_t>& permutation = {},
                       const std::vector<std::string>& dim_names = {})
      : ExtensionType(fixed_size_list(value_type, size)),
        value_type_(value_type),
        shape_(shape),
        permutation_(permutation),
        dim_names_(dim_names) {}

  std::string extension_name() const override { return "arrow.fixed_shape_tensor"; }
  std::string ToString(bool show_metadata = false) const override;

  /// Number of dimensions of tensor elements
  size_t ndim() const { return shape_.size(); }

  /// Shape of tensor elements
  const std::vector<int64_t> shape() const { return shape_; }

  /// Value type of tensor elements
  const std::shared_ptr<DataType> value_type() const { return value_type_; }

  /// Strides of tensor elements. Strides state offset in bytes between adjacent
  /// elements along each dimension. In case permutation is non-empty strides are
  /// computed from permuted tensor element's shape.
  const std::vector<int64_t>& strides();

  /// Permutation mapping from logical to physical memory layout of tensor elements
  const std::vector<int64_t>& permutation() const { return permutation_; }

  /// Dimension names of tensor elements. Dimensions are ordered physically.
  const std::vector<std::string>& dim_names() const { return dim_names_; }

  bool ExtensionEquals(const ExtensionType& other) const override;

  std::string Serialize() const override;

  Result<std::shared_ptr<DataType>> Deserialize(
      std::shared_ptr<DataType> storage_type,
      const std::string& serialized_data) const override;

  /// Create a FixedShapeTensorArray from ArrayData
  std::shared_ptr<Array> MakeArray(std::shared_ptr<ArrayData> data) const override;

  /// \brief Create a Tensor from an ExtensionScalar from a FixedShapeTensorArray
  ///
  /// This method will return a Tensor from ExtensionScalar with strides
  /// derived from shape and permutation of FixedShapeTensorType. Shape and
  /// dim_names will be permuted according to permutation stored in the
  /// FixedShapeTensorType metadata.
  static Result<std::shared_ptr<Tensor>> MakeTensor(
      const std::shared_ptr<ExtensionScalar>& scalar);

  /// \brief Create a FixedShapeTensorType instance
  static Result<std::shared_ptr<DataType>> Make(
      const std::shared_ptr<DataType>& value_type, const std::vector<int64_t>& shape,
      const std::vector<int64_t>& permutation = {},
      const std::vector<std::string>& dim_names = {});

 private:
  std::shared_ptr<DataType> storage_type_;
  std::shared_ptr<DataType> value_type_;
  std::vector<int64_t> shape_;
  std::vector<int64_t> strides_;
  std::vector<int64_t> permutation_;
  std::vector<std::string> dim_names_;
};

/// \brief Return a FixedShapeTensorType instance.
ARROW_EXPORT std::shared_ptr<DataType> fixed_shape_tensor(
    const std::shared_ptr<DataType>& storage_type, const std::vector<int64_t>& shape,
    const std::vector<int64_t>& permutation = {},
    const std::vector<std::string>& dim_names = {});

}  // namespace extension
}  // namespace arrow
