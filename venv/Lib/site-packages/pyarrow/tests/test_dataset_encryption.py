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

import base64
from datetime import timedelta
import random
import pyarrow.fs as fs
import pyarrow as pa

import pytest

encryption_unavailable = False

try:
    import pyarrow.parquet as pq
    import pyarrow.dataset as ds
except ImportError:
    pq = None
    ds = None

try:
    from pyarrow.tests.parquet.encryption import InMemoryKmsClient
    import pyarrow.parquet.encryption as pe
except ImportError:
    encryption_unavailable = True


# Marks all of the tests in this module
pytestmark = pytest.mark.dataset


FOOTER_KEY = b"0123456789112345"
FOOTER_KEY_NAME = "footer_key"
COL_KEY = b"1234567890123450"
COL_KEY_NAME = "col_key"


def create_sample_table():
    return pa.table(
        {
            "year": [2020, 2022, 2021, 2022, 2019, 2021],
            "n_legs": [2, 2, 4, 4, 5, 100],
            "animal": [
                "Flamingo",
                "Parrot",
                "Dog",
                "Horse",
                "Brittle stars",
                "Centipede",
            ],
        }
    )


def create_encryption_config():
    return pe.EncryptionConfiguration(
        footer_key=FOOTER_KEY_NAME,
        plaintext_footer=False,
        column_keys={COL_KEY_NAME: ["n_legs", "animal"]},
        encryption_algorithm="AES_GCM_V1",
        # requires timedelta or an assertion is raised
        cache_lifetime=timedelta(minutes=5.0),
        data_key_length_bits=256,
    )


def create_decryption_config():
    return pe.DecryptionConfiguration(cache_lifetime=300)


def create_kms_connection_config():
    return pe.KmsConnectionConfig(
        custom_kms_conf={
            FOOTER_KEY_NAME: FOOTER_KEY.decode("UTF-8"),
            COL_KEY_NAME: COL_KEY.decode("UTF-8"),
        }
    )


def kms_factory(kms_connection_configuration):
    return InMemoryKmsClient(kms_connection_configuration)


@pytest.mark.skipif(
    encryption_unavailable, reason="Parquet Encryption is not currently enabled"
)
def test_dataset_encryption_decryption():
    table = create_sample_table()

    encryption_config = create_encryption_config()
    decryption_config = create_decryption_config()
    kms_connection_config = create_kms_connection_config()

    crypto_factory = pe.CryptoFactory(kms_factory)
    parquet_encryption_cfg = ds.ParquetEncryptionConfig(
        crypto_factory, kms_connection_config, encryption_config
    )
    parquet_decryption_cfg = ds.ParquetDecryptionConfig(
        crypto_factory, kms_connection_config, decryption_config
    )

    # create write_options with dataset encryption config
    pformat = pa.dataset.ParquetFileFormat()
    write_options = pformat.make_write_options(encryption_config=parquet_encryption_cfg)

    mockfs = fs._MockFileSystem()
    mockfs.create_dir("/")

    ds.write_dataset(
        data=table,
        base_dir="sample_dataset",
        format=pformat,
        file_options=write_options,
        filesystem=mockfs,
    )

    # read without decryption config -> should error is dataset was properly encrypted
    pformat = pa.dataset.ParquetFileFormat()
    with pytest.raises(IOError, match=r"no decryption"):
        ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)

    # set decryption config for parquet fragment scan options
    pq_scan_opts = ds.ParquetFragmentScanOptions(
        decryption_config=parquet_decryption_cfg
    )
    pformat = pa.dataset.ParquetFileFormat(default_fragment_scan_options=pq_scan_opts)
    dataset = ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)

    assert table.equals(dataset.to_table())

    # set decryption properties for parquet fragment scan options
    decryption_properties = crypto_factory.file_decryption_properties(
        kms_connection_config, decryption_config)
    pq_scan_opts = ds.ParquetFragmentScanOptions(
        decryption_properties=decryption_properties
    )

    pformat = pa.dataset.ParquetFileFormat(default_fragment_scan_options=pq_scan_opts)
    dataset = ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)

    assert table.equals(dataset.to_table())


@pytest.mark.skipif(
    not encryption_unavailable, reason="Parquet Encryption is currently enabled"
)
def test_write_dataset_parquet_without_encryption():
    """Test write_dataset with ParquetFileFormat and test if an exception is thrown
    if you try to set encryption_config using make_write_options"""

    # Set the encryption configuration using ParquetFileFormat
    # and make_write_options
    pformat = pa.dataset.ParquetFileFormat()

    with pytest.raises(NotImplementedError):
        _ = pformat.make_write_options(encryption_config="some value")


@pytest.mark.skipif(
    encryption_unavailable, reason="Parquet Encryption is not currently enabled"
)
def test_large_row_encryption_decryption():
    """Test encryption and decryption of a large number of rows."""

    class NoOpKmsClient(pe.KmsClient):
        def wrap_key(self, key_bytes: bytes, _: str) -> bytes:
            b = base64.b64encode(key_bytes)
            return b

        def unwrap_key(self, wrapped_key: bytes, _: str) -> bytes:
            b = base64.b64decode(wrapped_key)
            return b

    row_count = 2**15 + 1
    table = pa.Table.from_arrays(
        [pa.array(
            [random.random() for _ in range(row_count)],
            type=pa.float32()
        )], names=["foo"]
    )

    kms_config = pe.KmsConnectionConfig()
    crypto_factory = pe.CryptoFactory(lambda _: NoOpKmsClient())
    encryption_config = pe.EncryptionConfiguration(
        footer_key="UNIMPORTANT_KEY",
        column_keys={"UNIMPORTANT_KEY": ["foo"]},
        double_wrapping=True,
        plaintext_footer=False,
        data_key_length_bits=128,
    )
    pqe_config = ds.ParquetEncryptionConfig(
        crypto_factory, kms_config, encryption_config
    )
    pqd_config = ds.ParquetDecryptionConfig(
        crypto_factory, kms_config, pe.DecryptionConfiguration()
    )
    scan_options = ds.ParquetFragmentScanOptions(decryption_config=pqd_config)
    file_format = ds.ParquetFileFormat(default_fragment_scan_options=scan_options)
    write_options = file_format.make_write_options(encryption_config=pqe_config)
    file_decryption_properties = crypto_factory.file_decryption_properties(kms_config)

    mockfs = fs._MockFileSystem()
    mockfs.create_dir("/")

    path = "large-row-test-dataset"
    ds.write_dataset(table, path, format=file_format,
                     file_options=write_options, filesystem=mockfs)

    file_path = path + "/part-0.parquet"
    new_table = pq.ParquetFile(
        file_path, decryption_properties=file_decryption_properties,
        filesystem=mockfs
    ).read()
    assert table == new_table

    dataset = ds.dataset(path, format=file_format, filesystem=mockfs)
    new_table = dataset.to_table()
    assert table == new_table
