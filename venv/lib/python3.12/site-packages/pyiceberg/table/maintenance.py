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
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

logger = logging.getLogger(__name__)


if TYPE_CHECKING:
    from pyiceberg.table import Table
    from pyiceberg.table.update.snapshot import ExpireSnapshots


class MaintenanceTable:
    tbl: Table

    def __init__(self, tbl: Table) -> None:
        self.tbl = tbl

    def expire_snapshots(self) -> ExpireSnapshots:
        """Return an ExpireSnapshots builder for snapshot expiration operations.

        Returns:
            ExpireSnapshots builder for configuring and executing snapshot expiration.
        """
        from pyiceberg.table import Transaction
        from pyiceberg.table.update.snapshot import ExpireSnapshots

        return ExpireSnapshots(transaction=Transaction(self.tbl, autocommit=True))
