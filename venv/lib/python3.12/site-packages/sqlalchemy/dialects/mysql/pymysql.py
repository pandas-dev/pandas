# dialects/mysql/pymysql.py
# Copyright (C) 2005-2025 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php

r"""

.. dialect:: mysql+pymysql
    :name: PyMySQL
    :dbapi: pymysql
    :connectstring: mysql+pymysql://<username>:<password>@<host>/<dbname>[?<options>]
    :url: https://pymysql.readthedocs.io/

Unicode
-------

Please see :ref:`mysql_unicode` for current recommendations on unicode
handling.

.. _pymysql_ssl:

SSL Connections
------------------

The PyMySQL DBAPI accepts the same SSL arguments as that of MySQLdb,
described at :ref:`mysqldb_ssl`.   See that section for additional examples.

If the server uses an automatically-generated certificate that is self-signed
or does not match the host name (as seen from the client), it may also be
necessary to indicate ``ssl_check_hostname=false`` in PyMySQL::

    connection_uri = (
        "mysql+pymysql://scott:tiger@192.168.0.134/test"
        "?ssl_ca=/home/gord/client-ssl/ca.pem"
        "&ssl_cert=/home/gord/client-ssl/client-cert.pem"
        "&ssl_key=/home/gord/client-ssl/client-key.pem"
        "&ssl_check_hostname=false"
    )

MySQL-Python Compatibility
--------------------------

The pymysql DBAPI is a pure Python port of the MySQL-python (MySQLdb) driver,
and targets 100% compatibility.   Most behavioral notes for MySQL-python apply
to the pymysql driver as well.

"""  # noqa
from __future__ import annotations

from typing import Any
from typing import Dict
from typing import Optional
from typing import TYPE_CHECKING
from typing import Union

from .mysqldb import MySQLDialect_mysqldb
from ...util import langhelpers
from ...util.typing import Literal

if TYPE_CHECKING:

    from ...engine.interfaces import ConnectArgsType
    from ...engine.interfaces import DBAPIConnection
    from ...engine.interfaces import DBAPICursor
    from ...engine.interfaces import DBAPIModule
    from ...engine.interfaces import PoolProxiedConnection
    from ...engine.url import URL


class MySQLDialect_pymysql(MySQLDialect_mysqldb):
    driver = "pymysql"
    supports_statement_cache = True

    description_encoding = None

    @langhelpers.memoized_property
    def supports_server_side_cursors(self) -> bool:
        try:
            cursors = __import__("pymysql.cursors").cursors
            self._sscursor = cursors.SSCursor
            return True
        except (ImportError, AttributeError):
            return False

    @classmethod
    def import_dbapi(cls) -> DBAPIModule:
        return __import__("pymysql")

    @langhelpers.memoized_property
    def _send_false_to_ping(self) -> bool:
        """determine if pymysql has deprecated, changed the default of,
        or removed the 'reconnect' argument of connection.ping().

        See #10492 and
        https://github.com/PyMySQL/mysqlclient/discussions/651#discussioncomment-7308971
        for background.

        """  # noqa: E501

        try:
            Connection = __import__(
                "pymysql.connections"
            ).connections.Connection
        except (ImportError, AttributeError):
            return True
        else:
            insp = langhelpers.get_callable_argspec(Connection.ping)
            try:
                reconnect_arg = insp.args[1]
            except IndexError:
                return False
            else:
                return reconnect_arg == "reconnect" and (
                    not insp.defaults or insp.defaults[0] is not False
                )

    def do_ping(self, dbapi_connection: DBAPIConnection) -> Literal[True]:
        if self._send_false_to_ping:
            dbapi_connection.ping(False)
        else:
            dbapi_connection.ping()

        return True

    def create_connect_args(
        self, url: URL, _translate_args: Optional[Dict[str, Any]] = None
    ) -> ConnectArgsType:
        if _translate_args is None:
            _translate_args = dict(username="user")
        return super().create_connect_args(
            url, _translate_args=_translate_args
        )

    def is_disconnect(
        self,
        e: DBAPIModule.Error,
        connection: Optional[Union[PoolProxiedConnection, DBAPIConnection]],
        cursor: Optional[DBAPICursor],
    ) -> bool:
        if super().is_disconnect(e, connection, cursor):
            return True
        elif isinstance(e, self.loaded_dbapi.Error):
            str_e = str(e).lower()
            return (
                "already closed" in str_e or "connection was killed" in str_e
            )
        else:
            return False

    def _extract_error_code(self, exception: BaseException) -> Any:
        if isinstance(exception.args[0], Exception):
            exception = exception.args[0]
        return exception.args[0]


dialect = MySQLDialect_pymysql
