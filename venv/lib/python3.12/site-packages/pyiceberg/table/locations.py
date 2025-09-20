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
import importlib
import logging
import uuid
from abc import ABC, abstractmethod
from typing import Optional

import mmh3

from pyiceberg.partitioning import PartitionKey
from pyiceberg.typedef import Properties
from pyiceberg.utils.properties import property_as_bool

logger = logging.getLogger(__name__)


class LocationProvider(ABC):
    """A base class for location providers, that provide file locations for a table's write tasks.

    Args:
        table_location (str): The table's base storage location.
        table_properties (Properties): The table's properties.
    """

    table_location: str
    table_properties: Properties

    data_path: str
    metadata_path: str

    def __init__(self, table_location: str, table_properties: Properties):
        self.table_location = table_location
        self.table_properties = table_properties

        from pyiceberg.table import TableProperties

        if path := table_properties.get(TableProperties.WRITE_DATA_PATH):
            self.data_path = path.rstrip("/")
        else:
            self.data_path = f"{self.table_location.rstrip('/')}/data"

        if path := table_properties.get(TableProperties.WRITE_METADATA_PATH):
            self.metadata_path = path.rstrip("/")
        else:
            self.metadata_path = f"{self.table_location.rstrip('/')}/metadata"

    @abstractmethod
    def new_data_location(self, data_file_name: str, partition_key: Optional[PartitionKey] = None) -> str:
        """Return a fully-qualified data file location for the given filename.

        Args:
            data_file_name (str): The name of the data file.
            partition_key (Optional[PartitionKey]): The data file's partition key. If None, the data is not partitioned.

        Returns:
            str: A fully-qualified location URI for the data file.
        """

    def new_table_metadata_file_location(self, new_version: int = 0) -> str:
        """Return a fully-qualified metadata file location for a new table version.

        Args:
            new_version (int): Version number of the metadata file.

        Returns:
            str: fully-qualified URI for the new table metadata file.

        Raises:
            ValueError: If the version is negative.
        """
        if new_version < 0:
            raise ValueError(f"Table metadata version: `{new_version}` must be a non-negative integer")

        file_name = f"{new_version:05d}-{uuid.uuid4()}.metadata.json"
        return self.new_metadata_location(file_name)

    def new_metadata_location(self, metadata_file_name: str) -> str:
        """Return a fully-qualified metadata file location for the given filename.

        Args:
            metadata_file_name (str): Name of the metadata file.

        Returns:
            str: A fully-qualified location URI for the metadata file.
        """
        return f"{self.metadata_path}/{metadata_file_name}"


class SimpleLocationProvider(LocationProvider):
    def __init__(self, table_location: str, table_properties: Properties):
        super().__init__(table_location, table_properties)

    def new_data_location(self, data_file_name: str, partition_key: Optional[PartitionKey] = None) -> str:
        return (
            f"{self.data_path}/{partition_key.to_path()}/{data_file_name}"
            if partition_key
            else f"{self.data_path}/{data_file_name}"
        )


class ObjectStoreLocationProvider(LocationProvider):
    HASH_BINARY_STRING_BITS = 20
    ENTROPY_DIR_LENGTH = 4
    ENTROPY_DIR_DEPTH = 3

    _include_partition_paths: bool

    def __init__(self, table_location: str, table_properties: Properties):
        super().__init__(table_location, table_properties)
        from pyiceberg.table import TableProperties

        self._include_partition_paths = property_as_bool(
            self.table_properties,
            TableProperties.WRITE_OBJECT_STORE_PARTITIONED_PATHS,
            TableProperties.WRITE_OBJECT_STORE_PARTITIONED_PATHS_DEFAULT,
        )

    def new_data_location(self, data_file_name: str, partition_key: Optional[PartitionKey] = None) -> str:
        if self._include_partition_paths and partition_key:
            return self.new_data_location(f"{partition_key.to_path()}/{data_file_name}")

        hashed_path = self._compute_hash(data_file_name)

        return (
            f"{self.data_path}/{hashed_path}/{data_file_name}"
            if self._include_partition_paths
            else f"{self.data_path}/{hashed_path}-{data_file_name}"
        )

    @staticmethod
    def _compute_hash(data_file_name: str) -> str:
        # Bitwise AND to combat sign-extension; bitwise OR to preserve leading zeroes that `bin` would otherwise strip.
        top_mask = 1 << ObjectStoreLocationProvider.HASH_BINARY_STRING_BITS
        hash_code = mmh3.hash(data_file_name) & (top_mask - 1) | top_mask
        return ObjectStoreLocationProvider._dirs_from_hash(bin(hash_code)[-ObjectStoreLocationProvider.HASH_BINARY_STRING_BITS :])

    @staticmethod
    def _dirs_from_hash(file_hash: str) -> str:
        """Divides hash into directories for optimized orphan removal operation using ENTROPY_DIR_DEPTH and ENTROPY_DIR_LENGTH."""
        total_entropy_length = ObjectStoreLocationProvider.ENTROPY_DIR_DEPTH * ObjectStoreLocationProvider.ENTROPY_DIR_LENGTH

        hash_with_dirs = []
        for i in range(0, total_entropy_length, ObjectStoreLocationProvider.ENTROPY_DIR_LENGTH):
            hash_with_dirs.append(file_hash[i : i + ObjectStoreLocationProvider.ENTROPY_DIR_LENGTH])

        if len(file_hash) > total_entropy_length:
            hash_with_dirs.append(file_hash[total_entropy_length:])

        return "/".join(hash_with_dirs)


def _import_location_provider(
    location_provider_impl: str, table_location: str, table_properties: Properties
) -> Optional[LocationProvider]:
    try:
        path_parts = location_provider_impl.split(".")
        if len(path_parts) < 2:
            from pyiceberg.table import TableProperties

            raise ValueError(
                f"{TableProperties.WRITE_PY_LOCATION_PROVIDER_IMPL} should be full path (module.CustomLocationProvider), got: {location_provider_impl}"
            )
        module_name, class_name = ".".join(path_parts[:-1]), path_parts[-1]
        module = importlib.import_module(module_name)
        class_ = getattr(module, class_name)
        return class_(table_location, table_properties)
    except ModuleNotFoundError as exc:
        logger.warning(f"Could not initialize LocationProvider: {location_provider_impl}", exc_info=exc)
        return None


def load_location_provider(table_location: str, table_properties: Properties) -> LocationProvider:
    from pyiceberg.table import TableProperties

    table_location = table_location.rstrip("/")

    if location_provider_impl := table_properties.get(TableProperties.WRITE_PY_LOCATION_PROVIDER_IMPL):
        if location_provider := _import_location_provider(location_provider_impl, table_location, table_properties):
            logger.info("Loaded LocationProvider: %s", location_provider_impl)
            return location_provider
        else:
            raise ValueError(f"Could not initialize LocationProvider: {location_provider_impl}")

    if property_as_bool(table_properties, TableProperties.OBJECT_STORE_ENABLED, TableProperties.OBJECT_STORE_ENABLED_DEFAULT):
        return ObjectStoreLocationProvider(table_location, table_properties)
    else:
        return SimpleLocationProvider(table_location, table_properties)
