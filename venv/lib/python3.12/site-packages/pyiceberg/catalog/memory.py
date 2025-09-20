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

from pyiceberg.catalog import URI
from pyiceberg.catalog.sql import SqlCatalog


class InMemoryCatalog(SqlCatalog):
    """
    An in-memory catalog implementation that uses SqlCatalog with SQLite in-memory database.

    This is useful for test, demo, and playground but not in production as it does not support concurrent access.
    """

    def __init__(self, name: str, warehouse: str = "file:///tmp/iceberg/warehouse", **kwargs: str) -> None:
        self._warehouse_location = warehouse
        if URI not in kwargs:
            kwargs[URI] = "sqlite:///:memory:"
        super().__init__(name=name, warehouse=warehouse, **kwargs)
