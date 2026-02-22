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
from contextlib import contextmanager
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
KEYS = {FOOTER_KEY_NAME: FOOTER_KEY, COL_KEY_NAME: COL_KEY}
EXTRA_COL_KEY = b"2345678901234501"
EXTRA_COL_KEY_NAME = "col2_key"
COLUMNS = ["year", "n_legs", "animal"]
COLUMN_KEYS = {COL_KEY_NAME: ["n_legs", "animal"]}


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


def create_encryption_config(footer_key=FOOTER_KEY_NAME, column_keys=COLUMN_KEYS):
    return pe.EncryptionConfiguration(
        footer_key=footer_key,
        plaintext_footer=False,
        column_keys=column_keys,
        encryption_algorithm="AES_GCM_V1",
        # requires timedelta or an assertion is raised
        cache_lifetime=timedelta(minutes=5.0),
        data_key_length_bits=256,
    )


def create_decryption_config():
    return pe.DecryptionConfiguration(cache_lifetime=300)


def create_kms_connection_config(keys=KEYS):
    return pe.KmsConnectionConfig(
        custom_kms_conf={
            key_name: key.decode("UTF-8")
            for key_name, key in keys.items()
        }
    )


def kms_factory(kms_connection_configuration):
    return InMemoryKmsClient(kms_connection_configuration)


@contextmanager
def cond_raises(success, error_type, match):
    if success:
        yield
    else:
        with pytest.raises(error_type, match=match):
            yield


def do_test_dataset_encryption_decryption(table, extra_column_path=None):
    # use extra column key for column extra_column_path, if given
    if extra_column_path:
        keys = dict(**KEYS, **{EXTRA_COL_KEY_NAME: EXTRA_COL_KEY})
        column_keys = dict(**COLUMN_KEYS, **{EXTRA_COL_KEY_NAME: [extra_column_path]})
        extra_column_name = extra_column_path.split(".")[0]
    else:
        keys = KEYS
        column_keys = COLUMN_KEYS
        extra_column_name = None

    # define the actual test
    def assert_decrypts(
        read_keys,
        read_columns,
        to_table_success,
        dataset_success=True,
    ):
        # use all keys for writing
        write_keys = keys
        encryption_config = create_encryption_config(FOOTER_KEY_NAME, column_keys)
        decryption_config = create_decryption_config()
        encrypt_kms_connection_config = create_kms_connection_config(write_keys)
        decrypt_kms_connection_config = create_kms_connection_config(read_keys)

        crypto_factory = pe.CryptoFactory(kms_factory)
        parquet_encryption_cfg = ds.ParquetEncryptionConfig(
            crypto_factory, encrypt_kms_connection_config, encryption_config
        )
        parquet_decryption_cfg = ds.ParquetDecryptionConfig(
            crypto_factory, decrypt_kms_connection_config, decryption_config
        )

        # create write_options with dataset encryption config
        pformat = pa.dataset.ParquetFileFormat()
        write_options = pformat.make_write_options(
            encryption_config=parquet_encryption_cfg
        )

        mockfs = fs._MockFileSystem()
        mockfs.create_dir("/")

        ds.write_dataset(
            data=table,
            base_dir="sample_dataset",
            format=pformat,
            file_options=write_options,
            filesystem=mockfs,
        )

        # read without decryption config -> errors if dataset was properly encrypted
        pformat = pa.dataset.ParquetFileFormat()
        with pytest.raises(IOError, match=r"no decryption"):
            ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)

        # set decryption config for parquet fragment scan options
        pq_scan_opts = ds.ParquetFragmentScanOptions(
            decryption_config=parquet_decryption_cfg
        )
        pformat = pa.dataset.ParquetFileFormat(
            default_fragment_scan_options=pq_scan_opts
        )
        with cond_raises(dataset_success, ValueError, match="Unknown master key"):
            dataset = ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)
            with cond_raises(to_table_success, ValueError, match="Unknown master key"):
                assert table.select(read_columns).equals(dataset.to_table(read_columns))

        # set decryption properties for parquet fragment scan options
        decryption_properties = crypto_factory.file_decryption_properties(
            decrypt_kms_connection_config, decryption_config)
        pq_scan_opts = ds.ParquetFragmentScanOptions(
            decryption_properties=decryption_properties
        )

        pformat = pa.dataset.ParquetFileFormat(
            default_fragment_scan_options=pq_scan_opts
        )
        with cond_raises(dataset_success, ValueError, match="Unknown master key"):
            dataset = ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)
            with cond_raises(to_table_success, ValueError, match="Unknown master key"):
                assert table.select(read_columns).equals(dataset.to_table(read_columns))

    # some notable column names and keys
    all_column_names = table.column_names
    encrypted_column_names = [column_name.split(".")[0]
                              for key_name, column_names in column_keys.items()
                              for column_name in column_names]
    plaintext_column_names = [column_name
                              for column_name in all_column_names
                              if column_name not in encrypted_column_names and
                              (extra_column_path is None or
                               not extra_column_path.startswith(f"{column_name}."))]
    assert len(encrypted_column_names) > 0
    assert len(plaintext_column_names) > 0
    footer_key_only = {FOOTER_KEY_NAME: FOOTER_KEY}
    column_keys_only = {key_name: key
                        for key_name, key in keys.items()
                        if key_name != FOOTER_KEY_NAME}

    # the test scenarios

    # read with footer key only, can only read plaintext columns
    assert_decrypts(footer_key_only, plaintext_column_names, True)
    assert_decrypts(footer_key_only, encrypted_column_names, False)
    assert_decrypts(footer_key_only, all_column_names, False)

    # read with all but footer key, cannot read any columns
    assert_decrypts(column_keys_only, plaintext_column_names, False, False)
    assert_decrypts(column_keys_only, encrypted_column_names, False, False)
    assert_decrypts(column_keys_only, all_column_names, False, False)

    # with footer key and one column key, all plaintext and
    # those encrypted columns that use that key, can be read
    if len(column_keys) > 1:
        for column_key_name, column_key_column_names in column_keys.items():
            for encrypted_column_name in column_key_column_names:
                # if one nested field of a column is encrypted,
                # the entire column is considered encrypted
                encrypted_column_name = encrypted_column_name.split(".")[0]

                # decrypt with footer key and one column key
                read_keys = {key_name: key
                             for key_name, key in keys.items()
                             if key_name in [FOOTER_KEY_NAME, column_key_name]}

                # that one encrypted column can only be read
                # if it is not a column path / nested field
                plaintext_and_one_success = encrypted_column_name != extra_column_name
                plaintext_and_one = plaintext_column_names + [encrypted_column_name]

                assert_decrypts(read_keys, plaintext_column_names, True)
                assert_decrypts(read_keys, plaintext_and_one, plaintext_and_one_success)
                assert_decrypts(read_keys, encrypted_column_names, False)
                assert_decrypts(read_keys, all_column_names, False)

    # with all column keys, all columns can be read
    assert_decrypts(keys, plaintext_column_names, True)
    assert_decrypts(keys, encrypted_column_names, True)
    assert_decrypts(keys, all_column_names, True)


@pytest.mark.skipif(
    encryption_unavailable, reason="Parquet Encryption is not currently enabled"
)
def test_dataset_encryption_decryption():
    do_test_dataset_encryption_decryption(create_sample_table())


@pytest.mark.skipif(
    encryption_unavailable, reason="Parquet Encryption is not currently enabled"
)
@pytest.mark.parametrize("column_name", ["list", "list.list.element"])
def test_list_encryption_decryption(column_name):
    list_data = pa.array(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9], [-1], [-2], [-3]],
        type=pa.list_(pa.int32()),
    )
    table = create_sample_table().append_column("list", list_data)

    do_test_dataset_encryption_decryption(table, column_name)


@pytest.mark.skipif(
    encryption_unavailable,
    reason="Parquet Encryption is not currently enabled"
)
@pytest.mark.parametrize(
    "column_name", ["map", "map.key_value.key", "map.key_value.value"]
)
def test_map_encryption_decryption(column_name):
    map_type = pa.map_(pa.string(), pa.int32())
    map_data = pa.array(
        [
            [("k1", 1), ("k2", 2)], [("k1", 3), ("k3", 4)], [("k2", 5), ("k3", 6)],
            [("k4", 7)], [], []
        ],
        type=map_type
    )
    table = create_sample_table().append_column("map", map_data)

    do_test_dataset_encryption_decryption(table, column_name)


@pytest.mark.skipif(
    encryption_unavailable, reason="Parquet Encryption is not currently enabled"
)
@pytest.mark.parametrize("column_name", ["struct", "struct.f1", "struct.f2"])
def test_struct_encryption_decryption(column_name):
    struct_fields = [("f1", pa.int32()), ("f2", pa.string())]
    struct_type = pa.struct(struct_fields)
    struct_data = pa.array(
        [(1, "one"), (2, "two"), (3, "three"), (4, "four"), (5, "five"), (6, "six")],
        type=struct_type
    )
    table = create_sample_table().append_column("struct", struct_data)

    do_test_dataset_encryption_decryption(table, column_name)


@pytest.mark.skipif(
    encryption_unavailable,
    reason="Parquet Encryption is not currently enabled"
)
@pytest.mark.parametrize(
    "column_name",
    [
        "col",
        "col.list.element.key_value.key",
        "col.list.element.key_value.value.f1",
        "col.list.element.key_value.value.f2"
    ]
)
def test_deep_nested_encryption_decryption(column_name):
    struct_fields = [("f1", pa.int32()), ("f2", pa.string())]
    struct_type = pa.struct(struct_fields)
    struct1 = (1, "one")
    struct2 = (2, "two")
    struct3 = (3, "three")
    struct4 = (4, "four")
    struct5 = (5, "five")
    struct6 = (6, "six")

    map_type = pa.map_(pa.int32(), struct_type)
    map1 = {1: struct1, 2: struct2}
    map2 = {3: struct3}
    map3 = {4: struct4}
    map4 = {5: struct5, 6: struct6}

    list_type = pa.list_(map_type)
    list1 = [map1, map2]
    list2 = [map3]
    list3 = [map4]
    list_data = [pa.array([list1, list2, None, list3, None, None], type=list_type)]
    table = create_sample_table().append_column("col", list_data)

    do_test_dataset_encryption_decryption(table, column_name)


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


@pytest.mark.skipif(
    encryption_unavailable, reason="Parquet Encryption is not currently enabled"
)
def test_dataset_encryption_with_selected_column_statistics():
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
    # and specify that statistics should be enabled for a subset of columns.
    pformat = pa.dataset.ParquetFileFormat()
    write_options = pformat.make_write_options(
        encryption_config=parquet_encryption_cfg,
        write_statistics=["year", "n_legs"]
    )

    mockfs = fs._MockFileSystem()
    mockfs.create_dir("/")

    ds.write_dataset(
        data=table,
        base_dir="sample_dataset",
        format=pformat,
        file_options=write_options,
        filesystem=mockfs,
    )

    # Open Parquet files directly and check that statistics are present
    # for the expected columns.
    pq_scan_opts = ds.ParquetFragmentScanOptions(
        decryption_config=parquet_decryption_cfg
    )
    pformat = pa.dataset.ParquetFileFormat(default_fragment_scan_options=pq_scan_opts)
    dataset = ds.dataset("sample_dataset", format=pformat, filesystem=mockfs)

    for fragment in dataset.get_fragments():
        decryption_properties = crypto_factory.file_decryption_properties(
            kms_connection_config, decryption_config, fragment.path, mockfs)
        with pq.ParquetFile(
            fragment.path,
            decryption_properties=decryption_properties,
            filesystem=mockfs,
        ) as parquet_file:
            for rg_idx in range(parquet_file.metadata.num_row_groups):
                row_group = parquet_file.metadata.row_group(rg_idx)

                assert row_group.column(0).statistics is not None
                assert row_group.column(0).statistics.min == 2019
                assert row_group.column(0).statistics.max == 2022

                assert row_group.column(1).statistics is not None
                assert row_group.column(1).statistics.min == 2
                assert row_group.column(1).statistics.max == 100

                assert row_group.column(2).statistics is None
