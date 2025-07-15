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
# cython: language_level = 3

from pyarrow.includes.common cimport *
from pyarrow.includes.libparquet_encryption cimport *
from pyarrow._parquet cimport (ParquetCipher,
                               CFileEncryptionProperties,
                               CFileDecryptionProperties,
                               FileEncryptionProperties,
                               FileDecryptionProperties,
                               ParquetCipher_AES_GCM_V1,
                               ParquetCipher_AES_GCM_CTR_V1)
from pyarrow.lib cimport _Weakrefable

cdef class CryptoFactory(_Weakrefable):
    cdef shared_ptr[CPyCryptoFactory] factory
    cdef init(self, callable_client_factory)
    cdef inline shared_ptr[CPyCryptoFactory] unwrap(self)

cdef class EncryptionConfiguration(_Weakrefable):
    cdef shared_ptr[CEncryptionConfiguration] configuration
    cdef inline shared_ptr[CEncryptionConfiguration] unwrap(self) nogil

cdef class DecryptionConfiguration(_Weakrefable):
    cdef shared_ptr[CDecryptionConfiguration] configuration
    cdef inline shared_ptr[CDecryptionConfiguration] unwrap(self) nogil

cdef class KmsConnectionConfig(_Weakrefable):
    cdef shared_ptr[CKmsConnectionConfig] configuration
    cdef inline shared_ptr[CKmsConnectionConfig] unwrap(self) nogil

    @staticmethod
    cdef wrap(const CKmsConnectionConfig& config)


cdef shared_ptr[CCryptoFactory] pyarrow_unwrap_cryptofactory(object crypto_factory) except *
cdef shared_ptr[CKmsConnectionConfig] pyarrow_unwrap_kmsconnectionconfig(object kmsconnectionconfig) except *
cdef shared_ptr[CEncryptionConfiguration] pyarrow_unwrap_encryptionconfig(object encryptionconfig) except *
cdef shared_ptr[CDecryptionConfiguration] pyarrow_unwrap_decryptionconfig(object decryptionconfig) except *
