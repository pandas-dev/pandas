from typing import Any, Final, final, overload
from typing_extensions import Self

__version__: Final[str]
ATTR_CASE: Final = 3271982
CASE_LOWER: Final = 1
CASE_NATURAL: Final = 0
CASE_UPPER: Final = 2
PARAM_FILE: Final = 11
QUOTED_LITERAL_REPLACEMENT_OFF: Final = 0
QUOTED_LITERAL_REPLACEMENT_ON: Final = 1
SQL_API_SQLROWCOUNT: Final[int]
SQL_ATTR_AUTOCOMMIT: Final[int]
SQL_ATTR_CALL_RETURN: Final[int]
SQL_ATTR_CURRENT_SCHEMA: Final[int]
SQL_ATTR_CURSOR_TYPE: Final[int]
SQL_ATTR_INFO_ACCTSTR: Final[int]
SQL_ATTR_INFO_APPLNAME: Final[int]
SQL_ATTR_INFO_PROGRAMNAME: Final[int]
SQL_ATTR_INFO_USERID: Final[int]
SQL_ATTR_INFO_WRKSTNNAME: Final[int]
SQL_ATTR_PARAMSET_SIZE: Final[int]
SQL_ATTR_PARAM_BIND_TYPE: Final[int]
SQL_ATTR_QUERY_TIMEOUT: Final[int]
SQL_ATTR_ROWCOUNT_PREFETCH: Final[int]
SQL_ATTR_TRUSTED_CONTEXT_PASSWORD: Final[int]
SQL_ATTR_TRUSTED_CONTEXT_USERID: Final[int]
SQL_ATTR_TXN_ISOLATION: Final[int]
SQL_ATTR_USE_TRUSTED_CONTEXT: Final[int]
SQL_ATTR_XML_DECLARATION: Final[int]
SQL_AUTOCOMMIT_OFF: Final[int]
SQL_AUTOCOMMIT_ON: Final[int]
SQL_BIGINT: Final[int]
SQL_BINARY: Final[int]
SQL_BIT: Final[int]
SQL_BLOB: Final[int]
SQL_BLOB_LOCATOR: Final[int]
SQL_BOOLEAN: Final[int]
SQL_CHAR: Final[int]
SQL_CLOB: Final[int]
SQL_CLOB_LOCATOR: Final[int]
SQL_CURSOR_DYNAMIC: Final[int]
SQL_CURSOR_FORWARD_ONLY: Final[int]
SQL_CURSOR_KEYSET_DRIVEN: Final[int]
SQL_CURSOR_STATIC: Final[int]
SQL_DBCLOB: Final[int]
SQL_DBCLOB_LOCATOR: Final[int]
SQL_DBMS_NAME: Final[int]
SQL_DBMS_VER: Final[int]
SQL_DECFLOAT: Final[int]
SQL_DECIMAL: Final[int]
SQL_DOUBLE: Final[int]
SQL_FALSE: Final[int]
SQL_FLOAT: Final[int]
SQL_GRAPHIC: Final[int]
SQL_INDEX_CLUSTERED: Final[int]
SQL_INDEX_OTHER: Final[int]
SQL_INTEGER: Final[int]
SQL_LONGVARBINARY: Final[int]
SQL_LONGVARCHAR: Final[int]
SQL_LONGVARGRAPHIC: Final[int]
SQL_NUMERIC: Final[int]
SQL_PARAM_BIND_BY_COLUMN: Final[int]
SQL_PARAM_INPUT: Final[int]
SQL_PARAM_INPUT_OUTPUT: Final[int]
SQL_PARAM_OUTPUT: Final[int]
SQL_REAL: Final[int]
SQL_ROWCOUNT_PREFETCH_OFF: Final[int]
SQL_ROWCOUNT_PREFETCH_ON: Final[int]
SQL_SMALLINT: Final[int]
SQL_TABLE_STAT: Final[int]
SQL_TINYINT: Final[int]
SQL_TRUE: Final[int]
SQL_TXN_NO_COMMIT: Final[int]
SQL_TXN_READ_COMMITTED: Final[int]
SQL_TXN_READ_UNCOMMITTED: Final[int]
SQL_TXN_REPEATABLE_READ: Final[int]
SQL_TXN_SERIALIZABLE: Final[int]
SQL_TYPE_DATE: Final[int]
SQL_TYPE_TIME: Final[int]
SQL_TYPE_TIMESTAMP: Final[int]
SQL_VARBINARY: Final[int]
SQL_VARCHAR: Final[int]
SQL_VARGRAPHIC: Final[int]
SQL_WCHAR: Final[int]
SQL_WLONGVARCHAR: Final[int]
SQL_WVARCHAR: Final[int]
SQL_XML: Final[int]
USE_WCHAR: Final = 100
WCHAR_NO: Final = 0
WCHAR_YES: Final = 1

SQL_ATTR_ACCESS_MODE: Final[int]
SQL_ATTR_ALLOW_INTERLEAVED_GETDATA: Final[int]
SQL_ATTR_ANSI_APP: Final[int]
SQL_ATTR_APPEND_FOR_FETCH_ONLY: Final[int]
SQL_ATTR_APP_USES_LOB_LOCATOR: Final[int]
SQL_ATTR_ASYNC_ENABLE: Final[int]
SQL_ATTR_AUTO_IPD: Final[int]
SQL_ATTR_CACHE_USRLIBL: Final[int]
SQL_ATTR_CLIENT_APPLCOMPAT: Final[int]
SQL_ATTR_CLIENT_CODEPAGE: Final[int]
SQL_ATTR_COLUMNWISE_MRI: Final[int]
SQL_ATTR_COMMITONEOF: Final[int]
SQL_ATTR_CONCURRENT_ACCESS_RESOLUTION: Final[int]
SQL_ATTR_CONFIG_KEYWORDS_ARRAY_SIZE: Final[int]
SQL_ATTR_CONFIG_KEYWORDS_MAXLEN: Final[int]
SQL_ATTR_CONNECTION_DEAD: Final[int]
SQL_ATTR_CONNECTTYPE: Final[int]
SQL_ATTR_CONNECT_NODE: Final[int]
SQL_ATTR_CONNECT_PASSIVE: Final[int]
SQL_ATTR_CONN_CONTEXT: Final[int]
SQL_ATTR_CURRENT_CATALOG: Final[int]
SQL_ATTR_CURRENT_IMPLICIT_XMLPARSE_OPTION: Final[int]
SQL_ATTR_CURRENT_PACKAGE_PATH: Final[int]
SQL_ATTR_CURRENT_PACKAGE_SET: Final[int]
SQL_ATTR_DATE_FMT: Final[int]
SQL_ATTR_DATE_SEP: Final[int]
SQL_ATTR_DB2EXPLAIN: Final[int]
SQL_ATTR_DB2_APPLICATION_HANDLE: Final[int]
SQL_ATTR_DB2_APPLICATION_ID: Final[int]
SQL_ATTR_DB2_SQLERRP: Final[int]
SQL_ATTR_DECFLOAT_ROUNDING_MODE: Final[int]
SQL_ATTR_DECIMAL_SEP: Final[int]
SQL_ATTR_DESCRIBE_CALL: Final[int]
SQL_ATTR_DESCRIBE_OUTPUT_LEVEL: Final[int]
SQL_ATTR_DETECT_READ_ONLY_TXN: Final[int]
SQL_ATTR_ENLIST_IN_DTC: Final[int]
SQL_ATTR_EXTENDED_INDICATORS: Final[int]
SQL_ATTR_FET_BUF_SIZE: Final[int]
SQL_ATTR_FORCE_ROLLBACK: Final[int]
SQL_ATTR_FREE_LOCATORS_ON_FETCH: Final[int]
SQL_ATTR_GET_LATEST_MEMBER: Final[int]
SQL_ATTR_GET_LATEST_MEMBER_NAME: Final[int]
SQL_ATTR_IGNORE_SERVER_LIST: Final[int]
SQL_ATTR_INFO_CRRTKN: Final[int]
SQL_ATTR_INFO_PROGRAMID: Final[int]
SQL_ATTR_KEEP_DYNAMIC: Final[int]
SQL_ATTR_LOB_CACHE_SIZE: Final[int]
SQL_ATTR_LOB_FILE_THRESHOLD: Final[int]
SQL_ATTR_LOGIN_TIMEOUT: Final[int]
SQL_ATTR_LONGDATA_COMPAT: Final[int]
SQL_ATTR_MAPCHAR: Final[int]
SQL_ATTR_MAXBLKEXT: Final[int]
SQL_ATTR_MAX_LOB_BLOCK_SIZE: Final[int]
SQL_ATTR_NETWORK_STATISTICS: Final[int]
SQL_ATTR_OVERRIDE_CHARACTER_CODEPAGE: Final[int]
SQL_ATTR_OVERRIDE_CODEPAGE: Final[int]
SQL_ATTR_OVERRIDE_PRIMARY_AFFINITY: Final[int]
SQL_ATTR_PARC_BATCH: Final[int]
SQL_ATTR_PING_DB: Final[int]
SQL_ATTR_PING_NTIMES: Final[int]
SQL_ATTR_PING_REQUEST_PACKET_SIZE: Final[int]
SQL_ATTR_QUERY_PREFETCH: Final[int]
SQL_ATTR_QUIET_MODE: Final[int]
SQL_ATTR_READ_ONLY_CONNECTION: Final[int]
SQL_ATTR_RECEIVE_TIMEOUT: Final[int]
SQL_ATTR_REOPT: Final[int]
SQL_ATTR_REPORT_ISLONG_FOR_LONGTYPES_OLEDB: Final[int]
SQL_ATTR_REPORT_SEAMLESSFAILOVER_WARNING: Final[int]
SQL_ATTR_REPORT_TIMESTAMP_TRUNC_AS_WARN: Final[int]
SQL_ATTR_RETRYONERROR: Final[int]
SQL_ATTR_RETRY_ON_MERGE: Final[int]
SQL_ATTR_SERVER_MSGTXT_MASK: Final[int]
SQL_ATTR_SERVER_MSGTXT_SP: Final[int]
SQL_ATTR_SESSION_GLOBAL_VAR: Final[int]
SQL_ATTR_SESSION_TIME_ZONE: Final[int]
SQL_ATTR_SPECIAL_REGISTER: Final[int]
SQL_ATTR_SQLCOLUMNS_SORT_BY_ORDINAL_OLEDB: Final[int]
SQL_ATTR_STMT_CONCENTRATOR: Final[int]
SQL_ATTR_STREAM_GETDATA: Final[int]
SQL_ATTR_STREAM_OUTPUTLOB_ON_CALL: Final[int]
SQL_ATTR_TIME_FMT: Final[int]
SQL_ATTR_TIME_SEP: Final[int]
SQL_ATTR_TRUSTED_CONTEXT_ACCESSTOKEN: Final[int]
SQL_ATTR_USER_REGISTRY_NAME: Final[int]
SQL_ATTR_WCHARTYPE: Final[int]

@final
class IBM_DBClientInfo:
    def __new__(cls, *args: object, **kwargs: object) -> Self: ...
    APPL_CODEPAGE: int
    CONN_CODEPAGE: int
    DATA_SOURCE_NAME: str
    DRIVER_NAME: str
    DRIVER_ODBC_VER: str
    DRIVER_VER: str
    ODBC_SQL_CONFORMANCE: str
    ODBC_VER: str

@final
class IBM_DBConnection:
    def __new__(cls, *args: object, **kwargs: object) -> Self: ...

@final
class IBM_DBServerInfo:
    def __new__(cls, *args: object, **kwargs: object) -> Self: ...
    DBMS_NAME: str
    DBMS_VER: str
    DB_CODEPAGE: int
    DB_NAME: str
    DFT_ISOLATION: str
    IDENTIFIER_QUOTE_CHAR: str
    INST_NAME: str
    ISOLATION_OPTION: tuple[str, str, str, str, str]
    KEYWORDS: str
    LIKE_ESCAPE_CLAUSE: bool
    MAX_COL_NAME_LEN: int
    MAX_IDENTIFIER_LEN: int
    MAX_INDEX_SIZE: int
    MAX_PROC_NAME_LEN: int
    MAX_ROW_SIZE: int
    MAX_SCHEMA_NAME_LEN: int
    MAX_STATEMENT_LEN: int
    MAX_TABLE_NAME_LEN: int
    NON_NULLABLE_COLUMNS: bool
    PROCEDURES: bool
    SPECIAL_CHARS: str
    SQL_CONFORMANCE: str

@final
class IBM_DBStatement:
    def __new__(cls, *args: object, **kwargs: object) -> Self: ...

def active(connection: IBM_DBConnection | None, /) -> bool: ...
def autocommit(connection: IBM_DBConnection, value: int = ..., /) -> int | bool: ...
def bind_param(
    stmt: IBM_DBStatement,
    parameter_number: int,
    variable: str,
    parameter_type: int | None = ...,
    data_type: int | None = ...,
    precision: int | None = ...,
    scale: int | None = ...,
    size: int | None = ...,
    /,
) -> bool: ...
@overload
def callproc(connection: IBM_DBConnection, procname: str, /) -> IBM_DBStatement | None: ...
@overload
def callproc(connection: IBM_DBConnection, procname: str, parameters: tuple[object, ...], /) -> tuple[object, ...] | None: ...
def check_function_support(connection: IBM_DBConnection, function_id: int, /) -> bool: ...
def client_info(connection: IBM_DBConnection, /) -> IBM_DBClientInfo | bool: ...
def close(connection: IBM_DBConnection, /) -> bool: ...
def column_privileges(
    connection: IBM_DBConnection,
    qualifier: str | None = ...,
    schema: str | None = ...,
    table_name: str | None = ...,
    column_name: str | None = ...,
    /,
) -> IBM_DBStatement: ...
def columns(
    connection: IBM_DBConnection,
    qualifier: str | None = ...,
    schema: str | None = ...,
    table_name: str | None = ...,
    column_name: str | None = ...,
    /,
) -> IBM_DBStatement: ...
def commit(connection: IBM_DBConnection, /) -> bool: ...
def conn_error(connection: IBM_DBConnection | None = ..., /) -> str: ...
def conn_errormsg(connection: IBM_DBConnection | None = ..., /) -> str: ...
def conn_warn(connection: IBM_DBConnection | None = ..., /) -> str: ...
def connect(
    database: str, user: str, password: str, options: dict[int, int | str] | None = ..., replace_quoted_literal: int = ..., /
) -> IBM_DBConnection | None: ...
def createdb(connection: IBM_DBConnection, dbName: str, codeSet: str = ..., mode: str = ..., /) -> bool: ...
def createdbNX(connection: IBM_DBConnection, dbName: str, codeSet: str = ..., mode: str = ..., /) -> bool: ...
def cursor_type(stmt: IBM_DBStatement, /) -> int: ...
def debug(option: str | bool) -> None: ...
def dropdb(connection: IBM_DBConnection, dbName: str, /) -> bool: ...
def exec_immediate(
    connection: IBM_DBConnection, statement: str | None, options: dict[int, int] = ..., /
) -> IBM_DBStatement | bool: ...
def execute(stmt: IBM_DBStatement, parameters: tuple[object, ...] | None = ..., /) -> bool: ...
def execute_many(
    stmt: IBM_DBStatement, seq_of_parameters: tuple[object, ...], options: dict[int, int] = ..., /
) -> int | None: ...
def fetchall(stmt: IBM_DBStatement, /) -> list[tuple[object, ...]]: ...
def fetchmany(stmt: IBM_DBStatement, numberOfRows: int, /) -> list[tuple[object, ...]]: ...
def fetchone(stmt: IBM_DBStatement, /) -> tuple[object, ...]: ...
def fetch_assoc(stmt: IBM_DBStatement, row_number: int = ..., /) -> dict[str, object] | bool: ...
def fetch_both(stmt: IBM_DBStatement, row_number: int = ..., /) -> dict[int | str, object] | bool: ...
def fetch_row(stmt: IBM_DBStatement, row_number: int = ..., /) -> bool: ...
def fetch_tuple(stmt: IBM_DBStatement, row_number: int = ..., /) -> tuple[object, ...]: ...
def field_display_size(stmt: IBM_DBStatement, column: int | str, /) -> int | bool: ...
def field_name(stmt: IBM_DBStatement, column: int | str, /) -> str | bool: ...
def field_nullable(stmt: IBM_DBStatement, column: int | str, /) -> bool: ...
def field_num(stmt: IBM_DBStatement, column: int | str, /) -> int | bool: ...
def field_precision(stmt: IBM_DBStatement, column: int | str, /) -> int | bool: ...
def field_scale(stmt: IBM_DBStatement, column: int | str, /) -> int | bool: ...
def field_type(stmt: IBM_DBStatement, column: int | str, /) -> str | bool: ...
def field_width(stmt: IBM_DBStatement, column: int | str, /) -> int | bool: ...
def foreign_keys(
    connection: IBM_DBConnection,
    pk_qualifier: str | None,
    pk_schema: str | None,
    pk_table_name: str | None,
    fk_qualifier: str | None = ...,
    fk_schema: str | None = ...,
    fk_table_name: str | None = ...,
    /,
) -> IBM_DBStatement: ...
def free_result(stmt: IBM_DBStatement, /) -> bool: ...
def free_stmt(stmt: IBM_DBStatement, /) -> bool: ...
def get_db_info(connection: IBM_DBConnection, option: int, /) -> str | bool: ...
def get_last_serial_value(stmt: IBM_DBStatement, /) -> str | bool: ...
def get_num_result(stmt: IBM_DBStatement, /) -> int | bool: ...
def get_option(resc: IBM_DBConnection | IBM_DBStatement, options: int, type: int, /) -> Any: ...
def get_sqlcode(connection_or_stmt: IBM_DBConnection | IBM_DBStatement | None = None, /) -> str: ...
def next_result(stmt: IBM_DBStatement, /) -> IBM_DBStatement | bool: ...
def num_fields(stmt: IBM_DBStatement, /) -> int | bool: ...
def num_rows(stmt: IBM_DBStatement, /) -> int: ...
def pconnect(
    database: str, username: str, password: str, options: dict[int, int | str] | None = ..., /
) -> IBM_DBConnection | None: ...
def prepare(
    connection: IBM_DBConnection, statement: str, options: dict[int, int | str] | None = ..., /
) -> IBM_DBStatement | bool: ...
def primary_keys(
    connection: IBM_DBConnection, qualifier: str | None, schema: str | None, table_name: str | None, /
) -> IBM_DBStatement: ...
def procedure_columns(
    connection: IBM_DBConnection, qualifier: str | None, schema: str | None, procedure: str | None, parameter: str | None, /
) -> IBM_DBStatement | bool: ...
def procedures(
    connection: IBM_DBConnection, qualifier: str | None, schema: str | None, procedure: str | None, /
) -> IBM_DBStatement | bool: ...
def recreatedb(connection: IBM_DBConnection, dbName: str, codeSet: str | None = ..., mode: str | None = ..., /) -> bool: ...
def result(stmt: IBM_DBStatement, column: int | str, /) -> Any: ...
def rollback(connection: IBM_DBConnection, /) -> bool: ...
def server_info(connection: IBM_DBConnection, /) -> IBM_DBServerInfo | bool: ...
def set_option(resc: IBM_DBConnection | IBM_DBStatement, options: dict[int, int | str], type: int, /) -> bool: ...
def special_columns(
    connection: IBM_DBConnection, qualifier: str | None, schema: str | None, table_name: str | None, scope: int, /
) -> IBM_DBStatement: ...
def statistics(
    connection: IBM_DBConnection, qualifier: str | None, schema: str | None, table_name: str | None, unique: bool | None, /
) -> IBM_DBStatement: ...
def stmt_error(stmt: IBM_DBStatement = ..., /) -> str: ...
def stmt_errormsg(stmt: IBM_DBStatement = ..., /) -> str: ...
def stmt_warn(connection: IBM_DBConnection = ..., /) -> IBM_DBStatement: ...
def table_privileges(
    connection: IBM_DBConnection, qualifier: str | None = ..., schema: str | None = ..., table_name: str | None = ..., /
) -> IBM_DBStatement | bool: ...
def tables(
    connection: IBM_DBConnection,
    qualifier: str | None = ...,
    schema: str | None = ...,
    table_name: str | None = ...,
    table_type: str | None = ...,
    /,
) -> IBM_DBStatement | bool: ...
