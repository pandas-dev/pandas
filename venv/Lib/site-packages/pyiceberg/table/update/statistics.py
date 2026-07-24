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
from typing import TYPE_CHECKING

from pyiceberg.table.statistics import StatisticsFile
from pyiceberg.table.update import (
    RemoveStatisticsUpdate,
    SetStatisticsUpdate,
    TableUpdate,
    UpdatesAndRequirements,
    UpdateTableMetadata,
)

if TYPE_CHECKING:
    from pyiceberg.table import Transaction


class UpdateStatistics(UpdateTableMetadata["UpdateStatistics"]):
    """
    Run statistics management operations using APIs.

    APIs include set_statistics and remove statistics operations.

    Use table.update_statistics().<operation>().commit() to run a specific operation.
    Use table.update_statistics().<operation-one>().<operation-two>().commit() to run multiple operations.

    Pending changes are applied on commit.

    We can also use context managers to make more changes. For example:

    with table.update_statistics() as update:
        update.set_statistics(statistics_file=statistics_file)
        update.remove_statistics(snapshot_id=2)
    """

    _updates: tuple[TableUpdate, ...] = ()

    def __init__(self, transaction: "Transaction") -> None:
        super().__init__(transaction)

    def set_statistics(self, statistics_file: StatisticsFile) -> "UpdateStatistics":
        self._updates += (
            SetStatisticsUpdate(
                statistics=statistics_file,
            ),
        )

        return self

    def remove_statistics(self, snapshot_id: int) -> "UpdateStatistics":
        self._updates = (
            RemoveStatisticsUpdate(
                snapshot_id=snapshot_id,
            ),
        )

        return self

    def _commit(self) -> UpdatesAndRequirements:
        return self._updates, ()
