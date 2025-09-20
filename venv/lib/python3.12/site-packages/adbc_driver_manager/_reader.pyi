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

import typing

import pandas
import pyarrow

class AdbcRecordBatchReader(pyarrow.RecordBatchReader):
    def close(self) -> None: ...
    def read_all(self) -> pyarrow.Table: ...
    def read_next_batch(self) -> pyarrow.RecordBatch: ...
    def read_pandas(self, **kwargs) -> pandas.DataFrame: ...
    @property
    def schema(self) -> pyarrow.Schema: ...
    @classmethod
    def _import_from_c(cls, address: int) -> AdbcRecordBatchReader: ...
    def __enter__(self) -> AdbcRecordBatchReader: ...
    def __exit__(self, type, value, traceback) -> None: ...
    def __iter__(self) -> typing.Iterator[pyarrow.RecordBatch]: ...
