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

"""Dataset support for Parquet file format."""

from pyarrow.includes.libarrow_dataset cimport *
from pyarrow.includes.libarrow_dataset_parquet cimport *

from pyarrow._dataset cimport FragmentScanOptions, FileWriteOptions


cdef class ParquetFragmentScanOptions(FragmentScanOptions):
    cdef:
        CParquetFragmentScanOptions* parquet_options
        object _parquet_decryption_config

    cdef void init(self, const shared_ptr[CFragmentScanOptions]& sp)
    cdef CReaderProperties* reader_properties(self)
    cdef ArrowReaderProperties* arrow_reader_properties(self)


cdef class ParquetFileWriteOptions(FileWriteOptions):

    cdef:
        CParquetFileWriteOptions* parquet_options
        object _properties
