import struct

from .constants import ER


class MySQLError(Exception):
    """Exception related to operation with MySQL."""


class Warning(Warning, MySQLError):
    """Exception raised for important warnings like data truncations
    while inserting, etc."""


class Error(MySQLError):
    """Exception that is the base class of all other error exceptions
    (not Warning)."""


class InterfaceError(Error):
    """Exception raised for errors that are related to the database
    interface rather than the database itself."""


class DatabaseError(Error):
    """Exception raised for errors that are related to the
    database."""


class DataError(DatabaseError):
    """Exception raised for errors that are due to problems with the
    processed data like division by zero, numeric value out of range,
    etc."""


class OperationalError(DatabaseError):
    """Exception raised for errors that are related to the database's
    operation and not necessarily under the control of the programmer,
    e.g. an unexpected disconnect occurs, the data source name is not
    found, a transaction could not be processed, a memory allocation
    error occurred during processing, etc."""


class IntegrityError(DatabaseError):
    """Exception raised when the relational integrity of the database
    is affected, e.g. a foreign key check fails, duplicate key,
    etc."""


class InternalError(DatabaseError):
    """Exception raised when the database encounters an internal
    error, e.g. the cursor is not valid anymore, the transaction is
    out of sync, etc."""


class ProgrammingError(DatabaseError):
    """Exception raised for programming errors, e.g. table not found
    or already exists, syntax error in the SQL statement, wrong number
    of parameters specified, etc."""


class NotSupportedError(DatabaseError):
    """Exception raised in case a method or database API was used
    which is not supported by the database, e.g. requesting a
    .rollback() on a connection that does not support transaction or
    has transactions turned off."""


error_map = {}


def _map_error(exc, *errors):
    for error in errors:
        error_map[error] = exc


_map_error(
    ProgrammingError,
    ER.DB_CREATE_EXISTS,
    ER.SYNTAX_ERROR,
    ER.PARSE_ERROR,
    ER.NO_SUCH_TABLE,
    ER.WRONG_DB_NAME,
    ER.WRONG_TABLE_NAME,
    ER.FIELD_SPECIFIED_TWICE,
    ER.INVALID_GROUP_FUNC_USE,
    ER.UNSUPPORTED_EXTENSION,
    ER.TABLE_MUST_HAVE_COLUMNS,
    ER.CANT_DO_THIS_DURING_AN_TRANSACTION,
    ER.WRONG_DB_NAME,
    ER.WRONG_COLUMN_NAME,
)
_map_error(
    DataError,
    ER.WARN_DATA_TRUNCATED,
    ER.WARN_NULL_TO_NOTNULL,
    ER.WARN_DATA_OUT_OF_RANGE,
    ER.NO_DEFAULT,
    ER.PRIMARY_CANT_HAVE_NULL,
    ER.DATA_TOO_LONG,
    ER.DATETIME_FUNCTION_OVERFLOW,
    ER.TRUNCATED_WRONG_VALUE_FOR_FIELD,
    ER.ILLEGAL_VALUE_FOR_TYPE,
)
_map_error(
    IntegrityError,
    ER.DUP_ENTRY,
    ER.NO_REFERENCED_ROW,
    ER.NO_REFERENCED_ROW_2,
    ER.ROW_IS_REFERENCED,
    ER.ROW_IS_REFERENCED_2,
    ER.CANNOT_ADD_FOREIGN,
    ER.BAD_NULL_ERROR,
)
_map_error(
    NotSupportedError,
    ER.WARNING_NOT_COMPLETE_ROLLBACK,
    ER.NOT_SUPPORTED_YET,
    ER.FEATURE_DISABLED,
    ER.UNKNOWN_STORAGE_ENGINE,
)
_map_error(
    OperationalError,
    ER.DBACCESS_DENIED_ERROR,
    ER.ACCESS_DENIED_ERROR,
    ER.CON_COUNT_ERROR,
    ER.TABLEACCESS_DENIED_ERROR,
    ER.COLUMNACCESS_DENIED_ERROR,
    ER.CONSTRAINT_FAILED,
    ER.LOCK_DEADLOCK,
)


del _map_error, ER


def raise_mysql_exception(data):
    errno = struct.unpack("<h", data[1:3])[0]
    # https://dev.mysql.com/doc/dev/mysql-server/latest/page_protocol_basic_err_packet.html
    # Error packet has optional sqlstate that is 5 bytes and starts with '#'.
    if data[3] == 0x23:  # '#'
        # sqlstate = data[4:9].decode()
        # TODO: Append (sqlstate) in the error message. This will be come in next minor release.
        errval = data[9:].decode("utf-8", "replace")
    else:
        errval = data[3:].decode("utf-8", "replace")
    errorclass = error_map.get(errno)
    if errorclass is None:
        errorclass = InternalError if errno < 1000 else OperationalError
    raise errorclass(errno, errval)
