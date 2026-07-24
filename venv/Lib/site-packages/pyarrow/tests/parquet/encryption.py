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
import pyarrow.parquet.encryption as pe
from pyarrow._parquet_encryption import FileSystemKeyMaterialStore
import re


class InMemoryKmsClient(pe.KmsClient):
    """This is a mock class implementation of KmsClient, built for testing
    only.
    """

    def __init__(self, config):
        """Create an InMemoryKmsClient instance."""
        pe.KmsClient.__init__(self)
        self.master_keys_map = config.custom_kms_conf

    def wrap_key(self, key_bytes, master_key_identifier):
        """Not a secure cipher - the wrapped key
        is just the master key concatenated with key bytes"""
        master_key_bytes = self.master_keys_map[master_key_identifier].encode(
            'utf-8')
        wrapped_key = b"".join([master_key_bytes, key_bytes])
        result = base64.b64encode(wrapped_key)
        return result

    def unwrap_key(self, wrapped_key, master_key_identifier):
        """Not a secure cipher - just extract the key from
        the wrapped key"""
        if master_key_identifier not in self.master_keys_map:
            raise ValueError("Unknown master key", master_key_identifier)
        expected_master_key = self.master_keys_map[master_key_identifier]
        decoded_wrapped_key = base64.b64decode(wrapped_key)
        master_key_bytes = decoded_wrapped_key[:16]
        decrypted_key = decoded_wrapped_key[16:]
        if (expected_master_key == master_key_bytes.decode('utf-8')):
            return decrypted_key
        raise ValueError("Incorrect master key used",
                         master_key_bytes, decrypted_key)


def parse_wrapped_key(wrapped_key: str) -> tuple[str, int, bytes]:
    """Parses a wrapped key string into a tuple: (key id, version, key) given
    input in the form: <key id>:v<version>:<bas64 encoded key>"""
    ptn = re.compile("(.+?):v([0-9]+?):(.+)")
    if m := ptn.fullmatch(wrapped_key):
        id, version, b64key = m.groups()
        version = int(version)
        key = base64.b64decode(b64key)
        return (id, version, key)
    else:
        raise ValueError("Cannot parse wrapped key", wrapped_key)


MASTER_KEY_VERSION = "master_key_version"


class MockVersioningKmsClient(pe.KmsClient):
    """This is a mock class implementation of KmsClient, built for testing
    only.

    During tests that involve CryptoFactory.rotate_master_keys, separate
    instances of this client will be created when writing, rotating keys, and
    reading back parquet data. To help unit tests verify that external key
    material was stored correctly at each step, this client wraps keys with a
    master_key_identifier and a version number. To ensure each client wraps
    with the correct version, the current version is persisted in the
    key_access_token attribute of the KmsConnectionConfig shared by all clients
    """

    def __init__(self, connection_config) -> None:
        pe.KmsClient.__init__(self)
        self.connection_config = connection_config

    @property
    def master_key_version(self) -> int:
        return int(self.connection_config.key_access_token)

    def wrap_key(self, key_bytes: bytes, master_key_identifier: str) -> str:
        b64key = base64.b64encode(key_bytes).decode('utf-8')
        return f"{master_key_identifier}:v{self.master_key_version}:{b64key}"

    def unwrap_key(
            self,
            wrapped_key: str,
            master_key_identifier: str) -> bytes:
        key_id, _, key = parse_wrapped_key(wrapped_key)
        if key_id != master_key_identifier:
            raise ValueError("Mismatched master key identifiers:",
                             key_id, master_key_identifier)
        return key


def verify_file_encrypted(path):
    """Verify that the file is encrypted by looking at its first 4 bytes.
    If it's the magic string PARE
    then this is a parquet with encrypted footer."""
    with open(path, "rb") as file:
        magic_str = file.read(4)
        # Verify magic string for parquet with encrypted footer is PARE
        assert magic_str == b'PARE'


def read_external_keys_to_dict(path):
    """Reads an external key material store given a parquet file path and
    returns a dict mapping master_key_id to KeyMaterial objects"""
    store = FileSystemKeyMaterialStore.for_file(path)
    keys = dict()
    for id in store.get_key_id_set():
        key_material = store.get_key_material(id)
        keys[key_material.master_key_id] = key_material
    return keys
