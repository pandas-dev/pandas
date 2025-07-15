# dialects/mysql/mysqlconnector.py
# Copyright (C) 2005-2025 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: ignore-errors


r"""
.. dialect:: mysql+mysqlconnector
    :name: MySQL Connector/Python
    :dbapi: myconnpy
    :connectstring: mysql+mysqlconnector://<user>:<password>@<host>[:<port>]/<dbname>
    :url: https://pypi.org/project/mysql-connector-python/

Driver Status
-------------

MySQL Connector/Python is supported as of SQLAlchemy 2.0.39 to the
degree which the driver is functional.   There are still ongoing issues
with features such as server side cursors which remain disabled until
upstream issues are repaired.

.. warning:: The MySQL Connector/Python driver published by Oracle is subject
   to frequent, major regressions of essential functionality such as being able
   to correctly persist simple binary strings which indicate it is not well
   tested.  The SQLAlchemy project is not able to maintain this dialect fully as
   regressions in the driver prevent it from being included in continuous
   integration.

.. versionchanged:: 2.0.39

    The MySQL Connector/Python dialect has been updated to support the
    latest version of this DBAPI.   Previously, MySQL Connector/Python
    was not fully supported.  However, support remains limited due to ongoing
    regressions introduced in this driver.

Connecting to MariaDB with MySQL Connector/Python
--------------------------------------------------

MySQL Connector/Python may attempt to pass an incompatible collation to the
database when connecting to MariaDB.  Experimentation has shown that using
``?charset=utf8mb4&collation=utfmb4_general_ci`` or similar MariaDB-compatible
charset/collation will allow connectivity.


"""  # noqa

import re

from .base import BIT
from .base import MariaDBIdentifierPreparer
from .base import MySQLCompiler
from .base import MySQLDialect
from .base import MySQLExecutionContext
from .base import MySQLIdentifierPreparer
from .mariadb import MariaDBDialect
from ... import util


class MySQLExecutionContext_mysqlconnector(MySQLExecutionContext):
    def create_server_side_cursor(self):
        return self._dbapi_connection.cursor(buffered=False)

    def create_default_cursor(self):
        return self._dbapi_connection.cursor(buffered=True)


class MySQLCompiler_mysqlconnector(MySQLCompiler):
    def visit_mod_binary(self, binary, operator, **kw):
        return (
            self.process(binary.left, **kw)
            + " % "
            + self.process(binary.right, **kw)
        )


class IdentifierPreparerCommon_mysqlconnector:
    @property
    def _double_percents(self):
        return False

    @_double_percents.setter
    def _double_percents(self, value):
        pass

    def _escape_identifier(self, value):
        value = value.replace(self.escape_quote, self.escape_to_quote)
        return value


class MySQLIdentifierPreparer_mysqlconnector(
    IdentifierPreparerCommon_mysqlconnector, MySQLIdentifierPreparer
):
    pass


class MariaDBIdentifierPreparer_mysqlconnector(
    IdentifierPreparerCommon_mysqlconnector, MariaDBIdentifierPreparer
):
    pass


class _myconnpyBIT(BIT):
    def result_processor(self, dialect, coltype):
        """MySQL-connector already converts mysql bits, so."""

        return None


class MySQLDialect_mysqlconnector(MySQLDialect):
    driver = "mysqlconnector"
    supports_statement_cache = True

    supports_sane_rowcount = True
    supports_sane_multi_rowcount = True

    supports_native_decimal = True

    supports_native_bit = True

    # not until https://bugs.mysql.com/bug.php?id=117548
    supports_server_side_cursors = False

    default_paramstyle = "format"
    statement_compiler = MySQLCompiler_mysqlconnector

    execution_ctx_cls = MySQLExecutionContext_mysqlconnector

    preparer = MySQLIdentifierPreparer_mysqlconnector

    colspecs = util.update_copy(MySQLDialect.colspecs, {BIT: _myconnpyBIT})

    @classmethod
    def import_dbapi(cls):
        from mysql import connector

        return connector

    def do_ping(self, dbapi_connection):
        dbapi_connection.ping(False)
        return True

    def create_connect_args(self, url):
        opts = url.translate_connect_args(username="user")

        opts.update(url.query)

        util.coerce_kw_type(opts, "allow_local_infile", bool)
        util.coerce_kw_type(opts, "autocommit", bool)
        util.coerce_kw_type(opts, "buffered", bool)
        util.coerce_kw_type(opts, "client_flag", int)
        util.coerce_kw_type(opts, "compress", bool)
        util.coerce_kw_type(opts, "connection_timeout", int)
        util.coerce_kw_type(opts, "connect_timeout", int)
        util.coerce_kw_type(opts, "consume_results", bool)
        util.coerce_kw_type(opts, "force_ipv6", bool)
        util.coerce_kw_type(opts, "get_warnings", bool)
        util.coerce_kw_type(opts, "pool_reset_session", bool)
        util.coerce_kw_type(opts, "pool_size", int)
        util.coerce_kw_type(opts, "raise_on_warnings", bool)
        util.coerce_kw_type(opts, "raw", bool)
        util.coerce_kw_type(opts, "ssl_verify_cert", bool)
        util.coerce_kw_type(opts, "use_pure", bool)
        util.coerce_kw_type(opts, "use_unicode", bool)

        # note that "buffered" is set to False by default in MySQL/connector
        # python.  If you set it to True, then there is no way to get a server
        # side cursor because the logic is written to disallow that.

        # leaving this at True until
        # https://bugs.mysql.com/bug.php?id=117548 can be fixed
        opts["buffered"] = True

        # FOUND_ROWS must be set in ClientFlag to enable
        # supports_sane_rowcount.
        if self.dbapi is not None:
            try:
                from mysql.connector.constants import ClientFlag

                client_flags = opts.get(
                    "client_flags", ClientFlag.get_default()
                )
                client_flags |= ClientFlag.FOUND_ROWS
                opts["client_flags"] = client_flags
            except Exception:
                pass

        return [[], opts]

    @util.memoized_property
    def _mysqlconnector_version_info(self):
        if self.dbapi and hasattr(self.dbapi, "__version__"):
            m = re.match(r"(\d+)\.(\d+)(?:\.(\d+))?", self.dbapi.__version__)
            if m:
                return tuple(int(x) for x in m.group(1, 2, 3) if x is not None)

    def _detect_charset(self, connection):
        return connection.connection.charset

    def _extract_error_code(self, exception):
        return exception.errno

    def is_disconnect(self, e, connection, cursor):
        errnos = (2006, 2013, 2014, 2045, 2055, 2048)
        exceptions = (
            self.dbapi.OperationalError,
            self.dbapi.InterfaceError,
            self.dbapi.ProgrammingError,
        )
        if isinstance(e, exceptions):
            return (
                e.errno in errnos
                or "MySQL Connection not available." in str(e)
                or "Connection to MySQL is not available" in str(e)
            )
        else:
            return False

    def _compat_fetchall(self, rp, charset=None):
        return rp.fetchall()

    def _compat_fetchone(self, rp, charset=None):
        return rp.fetchone()

    def get_isolation_level_values(self, dbapi_connection):
        return (
            "SERIALIZABLE",
            "READ UNCOMMITTED",
            "READ COMMITTED",
            "REPEATABLE READ",
            "AUTOCOMMIT",
        )

    def set_isolation_level(self, connection, level):
        if level == "AUTOCOMMIT":
            connection.autocommit = True
        else:
            connection.autocommit = False
            super().set_isolation_level(connection, level)


class MariaDBDialect_mysqlconnector(
    MariaDBDialect, MySQLDialect_mysqlconnector
):
    supports_statement_cache = True
    _allows_uuid_binds = False
    preparer = MariaDBIdentifierPreparer_mysqlconnector


dialect = MySQLDialect_mysqlconnector
mariadb_dialect = MariaDBDialect_mysqlconnector
