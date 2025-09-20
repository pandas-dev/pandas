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

from libc.stdint cimport int32_t, int64_t, uint8_t, uint32_t


cdef extern from "arrow-adbc/adbc.h" nogil:
    # C ABI

    ctypedef void (*CArrowSchemaRelease)(void*)
    ctypedef void (*CArrowArrayRelease)(void*)

    cdef struct CArrowSchema"ArrowSchema":
        CArrowSchemaRelease release

    cdef struct CArrowArray"ArrowArray":
        CArrowArrayRelease release

    ctypedef int (*CArrowArrayStreamGetLastError)(void*)
    ctypedef int (*CArrowArrayStreamGetNext)(void*, CArrowArray*)
    ctypedef char* (*CArrowArrayStreamGetSchema)(void*, CArrowSchema*)
    ctypedef void (*CArrowArrayStreamRelease)(void*)

    cdef struct CArrowArrayStream"ArrowArrayStream":
        CArrowArrayStreamGetLastError get_last_error
        CArrowArrayStreamGetNext get_next
        CArrowArrayStreamGetSchema get_schema
        CArrowArrayStreamRelease release

    # ADBC
    ctypedef uint8_t CAdbcStatusCode"AdbcStatusCode"
    cdef CAdbcStatusCode ADBC_STATUS_OK
    cdef CAdbcStatusCode ADBC_STATUS_UNKNOWN
    cdef CAdbcStatusCode ADBC_STATUS_NOT_IMPLEMENTED
    cdef CAdbcStatusCode ADBC_STATUS_NOT_FOUND
    cdef CAdbcStatusCode ADBC_STATUS_ALREADY_EXISTS
    cdef CAdbcStatusCode ADBC_STATUS_INVALID_ARGUMENT
    cdef CAdbcStatusCode ADBC_STATUS_INVALID_STATE
    cdef CAdbcStatusCode ADBC_STATUS_INVALID_DATA
    cdef CAdbcStatusCode ADBC_STATUS_INTEGRITY
    cdef CAdbcStatusCode ADBC_STATUS_INTERNAL
    cdef CAdbcStatusCode ADBC_STATUS_IO
    cdef CAdbcStatusCode ADBC_STATUS_CANCELLED
    cdef CAdbcStatusCode ADBC_STATUS_TIMEOUT
    cdef CAdbcStatusCode ADBC_STATUS_UNAUTHENTICATED
    cdef CAdbcStatusCode ADBC_STATUS_UNAUTHORIZED

    cdef const char* ADBC_OPTION_VALUE_DISABLED
    cdef const char* ADBC_OPTION_VALUE_ENABLED

    cdef const char* ADBC_CONNECTION_OPTION_AUTOCOMMIT
    cdef const char* ADBC_INGEST_OPTION_TARGET_TABLE
    cdef const char* ADBC_INGEST_OPTION_MODE
    cdef const char* ADBC_INGEST_OPTION_MODE_APPEND
    cdef const char* ADBC_INGEST_OPTION_MODE_CREATE
    cdef const char* ADBC_INGEST_OPTION_MODE_REPLACE
    cdef const char* ADBC_INGEST_OPTION_MODE_CREATE_APPEND

    cdef int ADBC_OBJECT_DEPTH_ALL
    cdef int ADBC_OBJECT_DEPTH_CATALOGS
    cdef int ADBC_OBJECT_DEPTH_DB_SCHEMAS
    cdef int ADBC_OBJECT_DEPTH_TABLES
    cdef int ADBC_OBJECT_DEPTH_COLUMNS

    cdef uint32_t ADBC_INFO_VENDOR_NAME
    cdef uint32_t ADBC_INFO_VENDOR_VERSION
    cdef uint32_t ADBC_INFO_VENDOR_ARROW_VERSION
    cdef uint32_t ADBC_INFO_DRIVER_NAME
    cdef uint32_t ADBC_INFO_DRIVER_VERSION
    cdef uint32_t ADBC_INFO_DRIVER_ARROW_VERSION

    ctypedef void (*CAdbcErrorRelease)(CAdbcError*)
    ctypedef void (*CAdbcPartitionsRelease)(CAdbcPartitions*)

    cdef struct CAdbcError"AdbcError":
        char* message
        int32_t vendor_code
        char[5] sqlstate
        CAdbcErrorRelease release

    cdef struct CAdbcErrorDetail"AdbcErrorDetail":
        char* key
        uint8_t* value
        size_t value_length

    int AdbcErrorGetDetailCount(CAdbcError* error)
    CAdbcErrorDetail AdbcErrorGetDetail(CAdbcError* error, int index)
    CAdbcError* AdbcErrorFromArrayStream(CArrowArrayStream*, CAdbcStatusCode*)

    cdef int ADBC_ERROR_VENDOR_CODE_PRIVATE_DATA

    cdef struct CAdbcDriver"AdbcDriver":
        pass

    cdef struct CAdbcDatabase"AdbcDatabase":
        void* private_data

    cdef struct CAdbcConnection"AdbcConnection":
        void* private_data

    cdef struct CAdbcStatement"AdbcStatement":
        void* private_data

    cdef struct CAdbcPartitions"AdbcPartitions":
        size_t num_partitions
        const uint8_t** partitions
        const size_t* partition_lengths
        void* private_data
        CAdbcPartitionsRelease release

    CAdbcStatusCode AdbcDatabaseNew(CAdbcDatabase* database, CAdbcError* error)
    CAdbcStatusCode AdbcDatabaseGetOption(
        CAdbcDatabase*, const char*, char*, size_t*, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseGetOptionBytes(
        CAdbcDatabase*, const char*, uint8_t*, size_t*, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseGetOptionDouble(
        CAdbcDatabase*, const char*, double*, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseGetOptionInt(
        CAdbcDatabase*, const char*, int64_t*, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseSetOption(
        CAdbcDatabase*, const char*, const char*, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseSetOptionBytes(
        CAdbcDatabase*, const char*, const uint8_t*, size_t, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseSetOptionDouble(
        CAdbcDatabase*, const char*, double, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseSetOptionInt(
        CAdbcDatabase*, const char*, int64_t, CAdbcError*)
    CAdbcStatusCode AdbcDatabaseInit(CAdbcDatabase* database, CAdbcError* error)
    CAdbcStatusCode AdbcDatabaseRelease(CAdbcDatabase* database, CAdbcError* error)

    ctypedef void (*CAdbcDriverInitFunc "AdbcDriverInitFunc")(int, void*, CAdbcError*)
    CAdbcStatusCode AdbcDriverManagerDatabaseSetInitFunc(
        CAdbcDatabase* database,
        CAdbcDriverInitFunc init_func,
        CAdbcError* error)

    CAdbcStatusCode AdbcConnectionCancel(CAdbcConnection*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionCommit(
        CAdbcConnection* connection,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionRollback(
        CAdbcConnection* connection,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionReadPartition(
        CAdbcConnection* connection,
        const uint8_t* serialized_partition,
        size_t serialized_length,
        CArrowArrayStream* out,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionGetInfo(
        CAdbcConnection* connection,
        const uint32_t* info_codes,
        size_t info_codes_length,
        CArrowArrayStream* stream,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionGetObjects(
        CAdbcConnection* connection,
        int depth,
        const char* catalog,
        const char* db_schema,
        const char* table_name,
        const char** table_type,
        const char* column_name,
        CArrowArrayStream* stream,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionGetOption(
        CAdbcConnection*, const char*, char*, size_t*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionGetOptionBytes(
        CAdbcConnection*, const char*, uint8_t*, size_t*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionGetOptionDouble(
        CAdbcConnection*, const char*, double*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionGetOptionInt(
        CAdbcConnection*, const char*, int64_t*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionGetStatistics(
        CAdbcConnection*, const char*, const char*, const char*,
        char, CArrowArrayStream*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionGetStatisticNames(
        CAdbcConnection*, CArrowArrayStream*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionGetTableSchema(
        CAdbcConnection* connection,
        const char* catalog,
        const char* db_schema,
        const char* table_name,
        CArrowSchema* schema,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionGetTableTypes(
        CAdbcConnection* connection,
        CArrowArrayStream* stream,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionInit(
        CAdbcConnection* connection,
        CAdbcDatabase* database,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionNew(
        CAdbcConnection* connection,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionRelease(
        CAdbcConnection* connection,
        CAdbcError* error)
    CAdbcStatusCode AdbcConnectionSetOption(
        CAdbcConnection*, const char*, const char*, CAdbcError*)
    CAdbcStatusCode AdbcConnectionSetOptionBytes(
        CAdbcConnection*, const char*, const uint8_t*, size_t, CAdbcError*)
    CAdbcStatusCode AdbcConnectionSetOptionDouble(
        CAdbcConnection*, const char*, double, CAdbcError*)
    CAdbcStatusCode AdbcConnectionSetOptionInt(
        CAdbcConnection*, const char*, int64_t, CAdbcError*)
    CAdbcStatusCode AdbcStatementBind(
        CAdbcStatement* statement,
        CArrowArray*,
        CArrowSchema*,
        CAdbcError* error)

    CAdbcStatusCode AdbcStatementCancel(CAdbcStatement*, CAdbcError*)
    CAdbcStatusCode AdbcStatementBindStream(
        CAdbcStatement* statement,
        CArrowArrayStream*,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementExecutePartitions(
        CAdbcStatement* statement,
        CArrowSchema* schema, CAdbcPartitions* partitions,
        int64_t* rows_affected,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementExecuteQuery(
        CAdbcStatement* statement,
        CArrowArrayStream* out, int64_t* rows_affected,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementExecuteSchema(
        CAdbcStatement*, CArrowSchema*, CAdbcError*)
    CAdbcStatusCode AdbcStatementGetOption(
        CAdbcStatement*, const char*, char*, size_t*, CAdbcError*)
    CAdbcStatusCode AdbcStatementGetOptionBytes(
        CAdbcStatement*, const char*, uint8_t*, size_t*, CAdbcError*)
    CAdbcStatusCode AdbcStatementGetOptionDouble(
        CAdbcStatement*, const char*, double*, CAdbcError*)
    CAdbcStatusCode AdbcStatementGetOptionInt(
        CAdbcStatement*, const char*, int64_t*, CAdbcError*)
    CAdbcStatusCode AdbcStatementGetParameterSchema(
        CAdbcStatement* statement,
        CArrowSchema* schema,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementNew(
        CAdbcConnection* connection,
        CAdbcStatement* statement,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementPrepare(
        CAdbcStatement* statement,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementSetOption(
        CAdbcStatement*, const char*, const char*, CAdbcError*)
    CAdbcStatusCode AdbcStatementSetOptionBytes(
        CAdbcStatement*, const char*, const uint8_t*, size_t, CAdbcError*)
    CAdbcStatusCode AdbcStatementSetOptionDouble(
        CAdbcStatement*, const char*, double, CAdbcError*)
    CAdbcStatusCode AdbcStatementSetOptionInt(
        CAdbcStatement*, const char*, int64_t, CAdbcError*)
    CAdbcStatusCode AdbcStatementSetSqlQuery(
        CAdbcStatement* statement,
        const char* query,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementSetSubstraitPlan(
        CAdbcStatement* statement,
        const uint8_t* plan,
        size_t length,
        CAdbcError* error)
    CAdbcStatusCode AdbcStatementRelease(
        CAdbcStatement* statement,
        CAdbcError* error)

cdef const CAdbcError* PyAdbcErrorFromArrayStream(
    CArrowArrayStream* stream, CAdbcStatusCode* status)

cdef void check_error(CAdbcStatusCode status, CAdbcError* error) except *
cdef object convert_error(CAdbcStatusCode status, CAdbcError* error)

cdef extern from "arrow-adbc/adbc_driver_manager.h":
    const char* CAdbcStatusCodeMessage"AdbcStatusCodeMessage"(CAdbcStatusCode code)

    ctypedef uint32_t CAdbcLoadFlags"AdbcLoadFlags"
    cdef CAdbcLoadFlags CAdbcLoadFlagDefault"ADBC_LOAD_FLAG_DEFAULT"
    CAdbcStatusCode AdbcDriverManagerDatabaseSetLoadFlags(
        CAdbcDatabase* database,
        CAdbcLoadFlags flags,
        CAdbcError* error)
    CAdbcStatusCode AdbcDriverManagerDatabaseSetAdditionalSearchPathList(
        CAdbcDatabase* database,
        const char* path_list,
        CAdbcError* error)
