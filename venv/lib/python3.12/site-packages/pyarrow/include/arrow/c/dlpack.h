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

#include "arrow/array/array_base.h"
#include "arrow/c/dlpack_abi.h"

namespace arrow::dlpack {

/// \brief Export Arrow array as DLPack tensor.
///
/// DLMangedTensor is produced as defined by the DLPack protocol,
/// see https://dmlc.github.io/dlpack/latest/.
///
/// Data types for which the protocol is supported are
/// integer and floating-point data types.
///
/// DLPack protocol only supports arrays with one contiguous
/// memory region which means Arrow Arrays with validity buffers
/// are not supported.
///
/// \param[in] arr Arrow array
/// \return DLManagedTensor struct
ARROW_EXPORT
Result<DLManagedTensor*> ExportArray(const std::shared_ptr<Array>& arr);

ARROW_EXPORT
Result<DLManagedTensor*> ExportTensor(const std::shared_ptr<Tensor>& t);

/// \brief Get DLDevice with enumerator specifying the
/// type of the device data is stored on and index of the
/// device which is 0 by default for CPU.
///
/// \param[in] arr Arrow array
/// \return DLDevice struct
ARROW_EXPORT
Result<DLDevice> ExportDevice(const std::shared_ptr<Array>& arr);

ARROW_EXPORT
Result<DLDevice> ExportDevice(const std::shared_ptr<Tensor>& t);

}  // namespace arrow::dlpack
