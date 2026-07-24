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

# distutils: language = c++
# cython: language_level = 3

from pyarrow.includes.libparquet cimport *
from pyarrow.lib cimport _Weakrefable


cdef class FileEncryptionProperties:
    """File-level encryption properties for the low-level API"""
    cdef:
        shared_ptr[CFileEncryptionProperties] properties

    @staticmethod
    cdef inline FileEncryptionProperties wrap(
            shared_ptr[CFileEncryptionProperties] properties):

        result = FileEncryptionProperties()
        result.properties = properties
        return result

    cdef inline shared_ptr[CFileEncryptionProperties] unwrap(self):
        return self.properties

cdef shared_ptr[WriterProperties] _create_writer_properties(
    use_dictionary=*,
    compression=*,
    version=*,
    write_statistics=*,
    data_page_size=*,
    max_rows_per_page=*,
    compression_level=*,
    use_byte_stream_split=*,
    column_encoding=*,
    data_page_version=*,
    FileEncryptionProperties encryption_properties=*,
    write_batch_size=*,
    dictionary_pagesize_limit=*,
    write_page_index=*,
    write_page_checksum=*,
    sorting_columns=*,
    store_decimal_as_integer=*,
    use_content_defined_chunking=*
) except *


cdef shared_ptr[ArrowWriterProperties] _create_arrow_writer_properties(
    use_deprecated_int96_timestamps=*,
    coerce_timestamps=*,
    allow_truncated_timestamps=*,
    writer_engine_version=*,
    use_compliant_nested_type=*,
    store_schema=*,
    write_time_adjusted_to_utc=*,
) except *


# Unwrap the "list_type" argument for ArrowReaderProperties
cdef Type _unwrap_list_type(obj) except *


cdef class ParquetSchema(_Weakrefable):
    cdef:
        FileMetaData parent  # the FileMetaData owning the SchemaDescriptor
        const SchemaDescriptor* schema

cdef class FileMetaData(_Weakrefable):
    cdef:
        shared_ptr[CFileMetaData] sp_metadata
        CFileMetaData* _metadata
        ParquetSchema _schema

    cdef inline init(self, const shared_ptr[CFileMetaData]& metadata):
        self.sp_metadata = metadata
        self._metadata = metadata.get()

cdef class RowGroupMetaData(_Weakrefable):
    cdef:
        int index  # for pickling support
        unique_ptr[CRowGroupMetaData] up_metadata
        CRowGroupMetaData* metadata
        FileMetaData parent

    cdef inline init(self, FileMetaData parent, int index):
        if index < 0 or index >= parent.num_row_groups:
            raise IndexError('{0} out of bounds'.format(index))
        self.up_metadata = parent._metadata.RowGroup(index)
        self.metadata = self.up_metadata.get()
        self.parent = parent
        self.index = index


cdef class ColumnChunkMetaData(_Weakrefable):
    cdef:
        unique_ptr[CColumnChunkMetaData] up_metadata
        CColumnChunkMetaData* metadata
        RowGroupMetaData parent

    cdef inline init(self, RowGroupMetaData parent, int i):
        self.up_metadata = parent.metadata.ColumnChunk(i)
        self.metadata = self.up_metadata.get()
        self.parent = parent

cdef class Statistics(_Weakrefable):
    cdef:
        shared_ptr[CStatistics] statistics
        ColumnChunkMetaData parent

    cdef inline init(self, const shared_ptr[CStatistics]& statistics,
                     ColumnChunkMetaData parent):
        self.statistics = statistics
        self.parent = parent

cdef class GeoStatistics(_Weakrefable):
    cdef:
        shared_ptr[CParquetGeoStatistics] statistics
        ColumnChunkMetaData parent

    cdef inline init(self, const shared_ptr[CParquetGeoStatistics]& statistics,
                     ColumnChunkMetaData parent):
        self.statistics = statistics
        self.parent = parent

cdef class FileDecryptionProperties:
    """File-level decryption properties for the low-level API"""
    cdef:
        shared_ptr[CFileDecryptionProperties] properties

    @staticmethod
    cdef inline FileDecryptionProperties wrap(
            shared_ptr[CFileDecryptionProperties] properties):

        result = FileDecryptionProperties()
        result.properties = properties
        return result

    cdef inline shared_ptr[CFileDecryptionProperties] unwrap(self):
        return self.properties
