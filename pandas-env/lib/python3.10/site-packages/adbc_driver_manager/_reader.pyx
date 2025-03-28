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

import pyarrow
from cython.operator cimport dereference as deref
from libc.stdint cimport uintptr_t

from ._lib cimport *


cdef class _AdbcErrorHelper:
    cdef:
        CArrowArrayStream c_stream

    def check_error(self, exception):
        cdef:
            CAdbcStatusCode c_status = ADBC_STATUS_OK
            const CAdbcError* error = \
                PyAdbcErrorFromArrayStream(&self.c_stream, &c_status)

        exc = convert_error(c_status, <CAdbcError*> error)
        if exc is not None:
            # Suppress "During handling of the above exception, another
            # exception occurred"
            raise exc from None

        raise exception


# Can't directly inherit pyarrow.RecordBatchReader, but we want to in order to
# pass isinstance checks
class AdbcRecordBatchReader(pyarrow.RecordBatchReader):
    def __init__(self, reader, helper):
        self._reader = reader
        self._helper = helper

    def close(self):
        self._reader.close()

    @classmethod
    def _import_from_c(cls, address) -> AdbcRecordBatchReader:
        cdef:
            CArrowArrayStream* c_stream = <CArrowArrayStream*> <uintptr_t> address
            _AdbcErrorHelper helper = _AdbcErrorHelper.__new__(_AdbcErrorHelper)
        # Save a copy of the stream to use
        helper.c_stream = deref(c_stream)
        try:
            reader = pyarrow.RecordBatchReader._import_from_c(int(address))
        except Exception as e:
            helper.check_error(e)
        return cls(reader, helper)

    def __enter__(self):
        return self

    def __exit__(self, exc_info, exc_val, exc_tb):
        pass

    def __iter__(self):
        try:
            yield from self._reader
        except Exception as e:
            self._helper.check_error(e)

    def iter_batches_with_custom_metadata(self):
        try:
            return self._reader.iter_batches_with_custom_metadata()
        except Exception as e:
            self._helper.check_error(e)

    def read_all(self):
        try:
            return self._reader.read_all()
        except Exception as e:
            self._helper.check_error(e)

    def read_next_batch(self, *args, **kwargs):
        try:
            return self._reader.read_next_batch(*args, **kwargs)
        except Exception as e:
            self._helper.check_error(e)

    def read_next_batch_with_custom_metadata(self):
        try:
            return self._reader.read_next_batch_with_custom_metadata()
        except Exception as e:
            self._helper.check_error(e)

    def read_pandas(self, *args, **kwargs):
        try:
            return self._reader.read_pandas(*args, **kwargs)
        except Exception as e:
            self._helper.check_error(e)

    @property
    def schema(self):
        try:
            return self._reader.schema
        except Exception as e:
            self._helper.check_error(e)
