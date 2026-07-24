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

# cython: language_level = 3

"""Dataset support for Parquet encryption."""

from pyarrow.includes.libarrow_dataset_parquet cimport *
from pyarrow._parquet_encryption cimport *
from pyarrow._dataset_parquet cimport ParquetFragmentScanOptions, ParquetFileWriteOptions


cdef class ParquetEncryptionConfig(_Weakrefable):
    """
    Core configuration class encapsulating parameters for high-level encryption
    within the Parquet framework.

    The ParquetEncryptionConfig class serves as a bridge for passing encryption-related
    parameters to the appropriate components within the Parquet library. It maintains references
    to objects that define the encryption strategy, Key Management Service (KMS) configuration,
    and specific encryption configurations for Parquet data.

    Parameters
    ----------
    crypto_factory : pyarrow.parquet.encryption.CryptoFactory
        Shared pointer to a `CryptoFactory` object. The `CryptoFactory` is responsible for
        creating cryptographic components, such as encryptors and decryptors.
    kms_connection_config : pyarrow.parquet.encryption.KmsConnectionConfig
        Shared pointer to a `KmsConnectionConfig` object. This object holds the configuration
        parameters necessary for connecting to a Key Management Service (KMS).
    encryption_config : pyarrow.parquet.encryption.EncryptionConfiguration
        Shared pointer to an `EncryptionConfiguration` object. This object defines specific
        encryption settings for Parquet data, including the keys assigned to different columns.

    Raises
    ------
    ValueError
        Raised if `encryption_config` is None.
    """
    cdef:
        shared_ptr[CParquetEncryptionConfig] c_config

    # Avoid mistakenly creating attributes
    __slots__ = ()

    def __cinit__(self, CryptoFactory crypto_factory, KmsConnectionConfig kms_connection_config,
                  EncryptionConfiguration encryption_config):

        cdef shared_ptr[CEncryptionConfiguration] c_encryption_config

        if crypto_factory is None:
            raise ValueError("crypto_factory cannot be None")

        if kms_connection_config is None:
            raise ValueError("kms_connection_config cannot be None")

        if encryption_config is None:
            raise ValueError("encryption_config cannot be None")

        self.c_config.reset(new CParquetEncryptionConfig())

        c_encryption_config = pyarrow_unwrap_encryptionconfig(
            encryption_config)

        self.c_config.get().crypto_factory = pyarrow_unwrap_cryptofactory(crypto_factory)
        self.c_config.get().kms_connection_config = pyarrow_unwrap_kmsconnectionconfig(
            kms_connection_config)
        self.c_config.get().encryption_config = c_encryption_config

    @staticmethod
    cdef wrap(shared_ptr[CParquetEncryptionConfig] c_config):
        cdef ParquetEncryptionConfig python_config = ParquetEncryptionConfig.__new__(ParquetEncryptionConfig)
        python_config.c_config = c_config
        return python_config

    cdef shared_ptr[CParquetEncryptionConfig] unwrap(self):
        return self.c_config


cdef class ParquetDecryptionConfig(_Weakrefable):
    """
    Core configuration class encapsulating parameters for high-level decryption
    within the Parquet framework.

    ParquetDecryptionConfig is designed to pass decryption-related parameters to
    the appropriate decryption components within the Parquet library. It holds references to
    objects that define the decryption strategy, Key Management Service (KMS) configuration,
    and specific decryption configurations for reading encrypted Parquet data.

    Parameters
    ----------
    crypto_factory : pyarrow.parquet.encryption.CryptoFactory
        Shared pointer to a `CryptoFactory` object, pivotal in creating cryptographic
        components for the decryption process.
    kms_connection_config : pyarrow.parquet.encryption.KmsConnectionConfig
        Shared pointer to a `KmsConnectionConfig` object, containing parameters necessary
        for connecting to a Key Management Service (KMS) during decryption.
    decryption_config : pyarrow.parquet.encryption.DecryptionConfiguration
        Shared pointer to a `DecryptionConfiguration` object, specifying decryption settings
        for reading encrypted Parquet data.

    Raises
    ------
    ValueError
        Raised if `decryption_config` is None.
    """

    cdef:
        shared_ptr[CParquetDecryptionConfig] c_config

    # Avoid mistakingly creating attributes
    __slots__ = ()

    def __cinit__(self, CryptoFactory crypto_factory, KmsConnectionConfig kms_connection_config,
                  DecryptionConfiguration decryption_config):

        cdef shared_ptr[CDecryptionConfiguration] c_decryption_config

        if decryption_config is None:
            raise ValueError(
                "decryption_config cannot be None")

        self.c_config.reset(new CParquetDecryptionConfig())

        c_decryption_config = pyarrow_unwrap_decryptionconfig(
            decryption_config)

        self.c_config.get().crypto_factory = pyarrow_unwrap_cryptofactory(crypto_factory)
        self.c_config.get().kms_connection_config = pyarrow_unwrap_kmsconnectionconfig(
            kms_connection_config)
        self.c_config.get().decryption_config = c_decryption_config

    @staticmethod
    cdef wrap(shared_ptr[CParquetDecryptionConfig] c_config):
        cdef ParquetDecryptionConfig python_config = ParquetDecryptionConfig.__new__(ParquetDecryptionConfig)
        python_config.c_config = c_config
        return python_config

    cdef shared_ptr[CParquetDecryptionConfig] unwrap(self):
        return self.c_config


def set_encryption_config(
    ParquetFileWriteOptions opts not None,
    ParquetEncryptionConfig config not None
):
    cdef shared_ptr[CParquetEncryptionConfig] c_config = config.unwrap()
    opts.parquet_options.parquet_encryption_config = c_config


def set_decryption_properties(
    ParquetFragmentScanOptions opts not None,
    FileDecryptionProperties config not None
):
    cdef CReaderProperties* reader_props = opts.reader_properties()
    reader_props.file_decryption_properties(config.unwrap())


def set_decryption_config(
    ParquetFragmentScanOptions opts not None,
    ParquetDecryptionConfig config not None
):
    cdef shared_ptr[CParquetDecryptionConfig] c_config = config.unwrap()
    opts.parquet_options.parquet_decryption_config = c_config
