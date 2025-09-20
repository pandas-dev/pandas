# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

# distutils: language=c++
# cython: language_level = 3

from pyarrow.lib cimport *

cdef extern from * namespace "arrow::py" nogil:
    """
    #include "arrow/status.h"
    #include "arrow/extension_type.h"
    #include "arrow/json/from_string.h"

    namespace arrow {
    namespace py {

    class UuidArray : public ExtensionArray {
    public:
        using ExtensionArray::ExtensionArray;
    };

    class UuidType : public ExtensionType {
    public:
        UuidType() : ExtensionType(fixed_size_binary(16)) {}
        std::string extension_name() const override { return "example-uuid"; }

        bool ExtensionEquals(const ExtensionType& other) const override {
            return other.extension_name() == this->extension_name();
        }

        std::shared_ptr<Array> MakeArray(std::shared_ptr<ArrayData> data) const override {
            return std::make_shared<ExtensionArray>(data);
        }

        Result<std::shared_ptr<DataType>> Deserialize(
            std::shared_ptr<DataType> storage_type,
            const std::string& serialized) const override {
            return std::make_shared<UuidType>();
        }

        std::string Serialize() const override { return ""; }
    };


    std::shared_ptr<DataType> MakeUuidType() {
        return std::make_shared<UuidType>();
    }

    std::shared_ptr<Array> MakeUuidArray() {
        auto uuid_type = MakeUuidType();
        auto json = "[\\"abcdefghijklmno0\\", \\"0onmlkjihgfedcba\\"]";
        auto result = json::ArrayFromJSONString(fixed_size_binary(16), json);
        return ExtensionType::WrapArray(uuid_type, result.ValueOrDie());
    }

    std::once_flag uuid_registered;

    static bool RegisterUuidType() {
        std::call_once(uuid_registered, RegisterExtensionType,
                       std::make_shared<UuidType>());
        return true;
    }

    static auto uuid_type_registered = RegisterUuidType();

    }  // namespace py
    }  // namespace arrow
    """

    cdef shared_ptr[CDataType] CMakeUuidType" arrow::py::MakeUuidType"()
    cdef shared_ptr[CArray] CMakeUuidArray" arrow::py::MakeUuidArray"()


def _make_uuid_type():
    return pyarrow_wrap_data_type(CMakeUuidType())


def _make_uuid_array():
    return pyarrow_wrap_array(CMakeUuidArray())
