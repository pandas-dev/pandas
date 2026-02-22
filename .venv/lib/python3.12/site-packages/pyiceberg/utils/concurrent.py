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
"""Concurrency concepts that support efficient multi-threading."""

import os
from concurrent.futures import Executor, ThreadPoolExecutor

from pyiceberg.utils.config import Config


class ExecutorFactory:
    _instance: Executor | None = None
    _instance_pid: int | None = None

    @staticmethod
    def max_workers() -> int | None:
        """Return the max number of workers configured."""
        return Config().get_int("max-workers")

    @staticmethod
    def get_or_create() -> Executor:
        """Return the same executor in each call."""
        # ThreadPoolExecutor cannot be shared across processes.  If a new pid is found it means
        # there is a new process so a new executor is needed.  Otherwise, the executor may be in
        # an invalid state and tasks submitted will not be started.
        if ExecutorFactory._instance_pid != os.getpid():
            ExecutorFactory._instance_pid = os.getpid()
            ExecutorFactory._instance = None

        if ExecutorFactory._instance is None:
            max_workers = ExecutorFactory.max_workers()
            ExecutorFactory._instance = ThreadPoolExecutor(max_workers=max_workers)

        return ExecutorFactory._instance
