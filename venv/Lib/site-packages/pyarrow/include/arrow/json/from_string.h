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

// Implement a simple JSON representation format for arrays

#pragma once

#include <memory>
#include <string>
#include <string_view>

#include "arrow/status.h"
#include "arrow/type_fwd.h"
#include "arrow/util/visibility.h"

namespace arrow {

class Array;
class DataType;

namespace json {

/// \defgroup array-from-json-string FromJSONString Helpers
///
/// These helpers are intended to be used in examples, tests, or for quick
/// prototyping and are not intended to be used where performance matters.
///
/// See the <a href="../arrays.html#fromjsonstring-helpers">User Guide</a> for
/// more information.
///
/// @{

/// \brief Create an Array from a JSON string
///
/// \code {.cpp}
/// Result<std::shared_ptr<Array>> maybe_array =
///     ArrayFromJSONString(int64(), "[2, 3, null, 7, 11]");
/// \endcode
ARROW_EXPORT
Result<std::shared_ptr<Array>> ArrayFromJSONString(const std::shared_ptr<DataType>&,
                                                   const std::string& json);

/// \copydoc ArrayFromJSONString(const std::shared_ptr<DataType>&, const std::string&)
ARROW_EXPORT
Result<std::shared_ptr<Array>> ArrayFromJSONString(const std::shared_ptr<DataType>&,
                                                   std::string_view json);

/// \copydoc ArrayFromJSONString(const std::shared_ptr<DataType>&, const std::string&)
ARROW_EXPORT
Result<std::shared_ptr<Array>> ArrayFromJSONString(const std::shared_ptr<DataType>&,
                                                   const char* json);

/// \brief Create a ChunkedArray from a JSON string
///
/// \code {.cpp}
/// Result<std::shared_ptr<ChunkedArray>> maybe_chunked_array =
///     ChunkedArrayFromJSONString(int64(), {R"([5, 10])", R"([null])", R"([16])"});
/// \endcode
ARROW_EXPORT
Result<std::shared_ptr<ChunkedArray>> ChunkedArrayFromJSONString(
    const std::shared_ptr<DataType>& type, const std::vector<std::string>& json_strings);

/// \brief Create a DictionaryArray from a JSON string
///
/// \code {.cpp}
/// Result<std::shared_ptr<Array>> maybe_dict_array =
///     DictArrayFromJSONString(dictionary(int32(), utf8()), "[0, 1, 0, 2, 0, 3]",
///     R"(["k1", "k2", "k3", "k4"])");
/// \endcode
ARROW_EXPORT
Result<std::shared_ptr<Array>> DictArrayFromJSONString(const std::shared_ptr<DataType>&,
                                                       std::string_view indices_json,
                                                       std::string_view dictionary_json);

/// \brief Create a Scalar from a JSON string
/// \code {.cpp}
/// Result<std::shared_ptr<Scalar>> maybe_scalar =
///     ScalarFromJSONString(float64(), "42", &scalar);
/// \endcode
ARROW_EXPORT
Result<std::shared_ptr<Scalar>> ScalarFromJSONString(const std::shared_ptr<DataType>&,
                                                     std::string_view json);

/// \brief Create a DictionaryScalar from a JSON string
/// \code {.cpp}
/// Result<std::shared_ptr<Scalar>> maybe_dict_scalar =
///     DictScalarFromJSONString(dictionary(int32(), utf8()), "3", R"(["k1", "k2", "k3",
///     "k4"])", &scalar);
/// \endcode
ARROW_EXPORT
Result<std::shared_ptr<Scalar>> DictScalarFromJSONString(
    const std::shared_ptr<DataType>&, std::string_view index_json,
    std::string_view dictionary_json);

/// @}

}  // namespace json
}  // namespace arrow
