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

# distutils: language = c++

from pyarrow.includes.common cimport *
from pyarrow.includes.libarrow cimport CSecureString
from pyarrow.includes.libarrow_fs cimport CFileSystem
from pyarrow._parquet cimport (ParquetCipher,
                               CFileEncryptionProperties,
                               CFileDecryptionProperties,
                               ParquetCipher_AES_GCM_V1,
                               ParquetCipher_AES_GCM_CTR_V1)


cdef extern from "parquet/encryption/kms_client.h" \
        namespace "parquet::encryption" nogil:
    cdef cppclass CKmsClient" parquet::encryption::KmsClient":
        c_string WrapKey(const CSecureString& key,
                         const c_string& master_key_identifier) except +
        CSecureString UnwrapKey(const c_string& wrapped_key,
                                const c_string& master_key_identifier) except +

    cdef cppclass CKeyAccessToken" parquet::encryption::KeyAccessToken":
        CKeyAccessToken(const c_string value)
        void Refresh(const c_string& new_value)
        const c_string& value() const

    cdef cppclass CKmsConnectionConfig \
            " parquet::encryption::KmsConnectionConfig":
        CKmsConnectionConfig()
        c_string kms_instance_id
        c_string kms_instance_url
        shared_ptr[CKeyAccessToken] refreshable_key_access_token
        unordered_map[c_string, c_string] custom_kms_conf

# Callbacks for implementing Python kms clients
# Use typedef to emulate syntax for std::function<void(..)>
ctypedef void CallbackWrapKey(
    object, const CSecureString&, const c_string&, c_string*)
ctypedef void CallbackUnwrapKey(
    object, const c_string&, const c_string&, CSecureString*)

cdef extern from "parquet/encryption/kms_client_factory.h" \
        namespace "parquet::encryption" nogil:
    cdef cppclass CKmsClientFactory" parquet::encryption::KmsClientFactory":
        shared_ptr[CKmsClient] CreateKmsClient(
            const CKmsConnectionConfig& kms_connection_config) except +

# Callbacks for implementing Python kms client factories
# Use typedef to emulate syntax for std::function<void(..)>
ctypedef void CallbackCreateKmsClient(
    object,
    const CKmsConnectionConfig&, shared_ptr[CKmsClient]*)

cdef extern from "parquet/encryption/crypto_factory.h" \
        namespace "parquet::encryption" nogil:
    cdef cppclass CEncryptionConfiguration\
            " parquet::encryption::EncryptionConfiguration":
        CEncryptionConfiguration(const c_string& footer_key) except +
        c_string footer_key
        c_string column_keys
        c_bool uniform_encryption
        ParquetCipher encryption_algorithm
        c_bool plaintext_footer
        c_bool double_wrapping
        double cache_lifetime_seconds
        c_bool internal_key_material
        int32_t data_key_length_bits

    cdef cppclass CDecryptionConfiguration\
            " parquet::encryption::DecryptionConfiguration":
        CDecryptionConfiguration() except +
        double cache_lifetime_seconds

    cdef cppclass CCryptoFactory" parquet::encryption::CryptoFactory":
        void RegisterKmsClientFactory(
            shared_ptr[CKmsClientFactory] kms_client_factory) except +
        shared_ptr[CFileEncryptionProperties] GetFileEncryptionProperties(
            const CKmsConnectionConfig& kms_connection_config,
            const CEncryptionConfiguration& encryption_config,
            const c_string parquet_file_path,
            const shared_ptr[CFileSystem] file_system) except +*
        shared_ptr[CFileDecryptionProperties] GetFileDecryptionProperties(
            const CKmsConnectionConfig& kms_connection_config,
            const CDecryptionConfiguration& decryption_config,
            const c_string parquet_file_path,
            const shared_ptr[CFileSystem] file_system) except +*
        void RemoveCacheEntriesForToken(const c_string& access_token) except +
        void RemoveCacheEntriesForAllTokens() except +
        void RotateMasterKeys(const CKmsConnectionConfig& kms_connection_config,
                              const c_string parquet_file_path,
                              const shared_ptr[CFileSystem] file_system,
                              c_bool double_wrapping,
                              double cache_lifetime_seconds)

cdef extern from "parquet/encryption/file_key_material_store.h" \
        namespace "parquet::encryption" nogil:
    cdef cppclass CFileKeyMaterialStore\
            "parquet::encryption::FileKeyMaterialStore":
        @staticmethod
        c_string GetKeyMaterial(c_string key_id_in_file) except +
        vector[c_string] GetKeyIDSet() except +

cdef extern from "parquet/encryption/file_system_key_material_store.h" \
        namespace "parquet::encryption" nogil:
    cdef cppclass CFileSystemKeyMaterialStore\
            "parquet::encryption::FileSystemKeyMaterialStore":

        @staticmethod
        shared_ptr[CFileSystemKeyMaterialStore] Make(c_string parquet_file_path,
                                                     shared_ptr[CFileSystem] file_system,
                                                     c_bool use_tmp_prefix) except +

        c_string GetKeyMaterial(c_string key_id_in_file) except +

        vector[c_string] GetKeyIDSet() except +

cdef extern from "parquet/encryption/key_material.h" \
        namespace "parquet::encryption" nogil:
    cdef cppclass CKeyMaterial "parquet::encryption::KeyMaterial":
        @staticmethod
        CKeyMaterial Parse(const c_string& key_material_string)
        c_bool is_footer_key()
        c_bool is_double_wrapped()
        const c_string& master_key_id()
        const c_string& wrapped_dek()
        const c_string& kek_id()
        const c_string& wrapped_kek()
        const c_string& kms_instance_id()
        const c_string& kms_instance_url()

cdef extern from "arrow/python/parquet_encryption.h" \
        namespace "arrow::py::parquet::encryption" nogil:
    cdef cppclass CPyKmsClientVtable \
            " arrow::py::parquet::encryption::PyKmsClientVtable":
        CPyKmsClientVtable()
        function[CallbackWrapKey] wrap_key
        function[CallbackUnwrapKey] unwrap_key

    cdef cppclass CPyKmsClient\
            " arrow::py::parquet::encryption::PyKmsClient"(CKmsClient):
        CPyKmsClient(object handler, CPyKmsClientVtable vtable)

    cdef cppclass CPyKmsClientFactoryVtable\
            " arrow::py::parquet::encryption::PyKmsClientFactoryVtable":
        CPyKmsClientFactoryVtable()
        function[CallbackCreateKmsClient] create_kms_client

    cdef cppclass CPyKmsClientFactory\
            " arrow::py::parquet::encryption::PyKmsClientFactory"(
                CKmsClientFactory):
        CPyKmsClientFactory(object handler, CPyKmsClientFactoryVtable vtable)

    cdef cppclass CPyCryptoFactory\
            " arrow::py::parquet::encryption::PyCryptoFactory"(CCryptoFactory):
        CResult[shared_ptr[CFileEncryptionProperties]] \
            SafeGetFileEncryptionProperties(
            const CKmsConnectionConfig& kms_connection_config,
            const CEncryptionConfiguration& encryption_config,
            const c_string parquet_file_path,
            const shared_ptr[CFileSystem] filesystem)
        CResult[shared_ptr[CFileDecryptionProperties]] \
            SafeGetFileDecryptionProperties(
            const CKmsConnectionConfig& kms_connection_config,
            const CDecryptionConfiguration& decryption_config,
            const c_string parquet_file_path,
            const shared_ptr[CFileSystem] filesystem)
        CStatus SafeRotateMasterKeys(const CKmsConnectionConfig& kms_connection_config,
                                     const c_string parquet_file_path,
                                     const shared_ptr[CFileSystem] filesystem,
                                     c_bool double_wrapping,
                                     double cache_lifetime_seconds)
