# dialects/oracle/oracledb.py
# Copyright (C) 2005-2026 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: ignore-errors

r""".. dialect:: oracle+oracledb
    :name: python-oracledb
    :dbapi: oracledb
    :connectstring: oracle+oracledb://user:pass@hostname:port[/dbname][?service_name=<service>[&key=value&key=value...]]
    :url: https://oracle.github.io/python-oracledb/

Description
-----------

Python-oracledb is the Oracle Database driver for Python. It features a default
"thin" client mode that requires no dependencies, and an optional "thick" mode
that uses Oracle Client libraries.  It supports SQLAlchemy features including
two phase transactions and Asyncio.

Python-oracle is the renamed, updated cx_Oracle driver. Oracle is no longer
doing any releases in the cx_Oracle namespace.

The SQLAlchemy ``oracledb`` dialect provides both a sync and an async
implementation under the same dialect name. The proper version is
selected depending on how the engine is created:

* calling :func:`_sa.create_engine` with ``oracle+oracledb://...`` will
  automatically select the sync version::

    from sqlalchemy import create_engine

    sync_engine = create_engine(
        "oracle+oracledb://scott:tiger@localhost?service_name=FREEPDB1"
    )

* calling :func:`_asyncio.create_async_engine` with ``oracle+oracledb://...``
  will automatically select the async version::

    from sqlalchemy.ext.asyncio import create_async_engine

    asyncio_engine = create_async_engine(
        "oracle+oracledb://scott:tiger@localhost?service_name=FREEPDB1"
    )

  The asyncio version of the dialect may also be specified explicitly using the
  ``oracledb_async`` suffix::

      from sqlalchemy.ext.asyncio import create_async_engine

      asyncio_engine = create_async_engine(
          "oracle+oracledb_async://scott:tiger@localhost?service_name=FREEPDB1"
      )

.. versionadded:: 2.0.25 added support for the async version of oracledb.

Thick mode support
------------------

By default, the python-oracledb driver runs in a "thin" mode that does not
require Oracle Client libraries to be installed. The driver also supports a
"thick" mode that uses Oracle Client libraries to get functionality such as
Oracle Application Continuity.

To enable thick mode, call `oracledb.init_oracle_client()
<https://python-oracledb.readthedocs.io/en/latest/api_manual/module.html#oracledb.init_oracle_client>`_
explicitly, or pass the parameter ``thick_mode=True`` to
:func:`_sa.create_engine`. To pass custom arguments to
``init_oracle_client()``, like the ``lib_dir`` path, a dict may be passed, for
example::

    engine = sa.create_engine(
        "oracle+oracledb://...",
        thick_mode={
            "lib_dir": "/path/to/oracle/client/lib",
            "config_dir": "/path/to/network_config_file_directory",
            "driver_name": "my-app : 1.0.0",
        },
    )

Note that passing a ``lib_dir`` path should only be done on macOS or
Windows. On Linux it does not behave as you might expect.

.. seealso::

    python-oracledb documentation `Enabling python-oracledb Thick mode
    <https://python-oracledb.readthedocs.io/en/latest/user_guide/initialization.html#enabling-python-oracledb-thick-mode>`_

Connecting to Oracle Database
-----------------------------

python-oracledb provides several methods of indicating the target database.
The dialect translates from a series of different URL forms.

Given the hostname, port and service name of the target database, you can
connect in SQLAlchemy using the ``service_name`` query string parameter::

    engine = create_engine(
        "oracle+oracledb://scott:tiger@hostname:port?service_name=myservice"
    )

Connecting with Easy Connect strings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can pass any valid python-oracledb connection string as the ``dsn`` key
value in a :paramref:`_sa.create_engine.connect_args` dictionary.  See
python-oracledb documentation `Oracle Net Services Connection Strings
<https://python-oracledb.readthedocs.io/en/latest/user_guide/connection_handling.html#oracle-net-services-connection-strings>`_.

For example to use an `Easy Connect string
<https://download.oracle.com/ocomdocs/global/Oracle-Net-Easy-Connect-Plus.pdf>`_
with a timeout to prevent connection establishment from hanging if the network
transport to the database cannot be established in 30 seconds, and also setting
a keep-alive time of 60 seconds to stop idle network connections from being
terminated by a firewall::

    e = create_engine(
        "oracle+oracledb://@",
        connect_args={
            "user": "scott",
            "password": "tiger",
            "dsn": "hostname:port/myservice?transport_connect_timeout=30&expire_time=60",
        },
    )

The Easy Connect syntax has been enhanced during the life of Oracle Database.
Review the documentation for your database version.  The current documentation
is at `Understanding the Easy Connect Naming Method
<https://www.oracle.com/pls/topic/lookup?ctx=dblatest&id=GUID-B0437826-43C1-49EC-A94D-B650B6A4A6EE>`_.

The general syntax is similar to:

.. sourcecode:: text

    [[protocol:]//]host[:port][/[service_name]][?parameter_name=value{&parameter_name=value}]

Note that although the SQLAlchemy URL syntax ``hostname:port/dbname`` looks
like Oracle's Easy Connect syntax, it is different. SQLAlchemy's URL requires a
system identifier (SID) for the ``dbname`` component::

    engine = create_engine("oracle+oracledb://scott:tiger@hostname:port/sid")

Easy Connect syntax does not support SIDs. It uses services names, which are
the preferred choice for connecting to Oracle Database.

Passing python-oracledb connect arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Other python-oracledb driver `connection options
<https://python-oracledb.readthedocs.io/en/latest/api_manual/module.html#oracledb.connect>`_
can be passed in ``connect_args``.  For example::

    e = create_engine(
        "oracle+oracledb://@",
        connect_args={
            "user": "scott",
            "password": "tiger",
            "dsn": "hostname:port/myservice",
            "events": True,
            "mode": oracledb.AUTH_MODE_SYSDBA,
        },
    )

Connecting with tnsnames.ora TNS aliases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If no port, database name, or service name is provided, the dialect will use an
Oracle Database DSN "connection string".  This takes the "hostname" portion of
the URL as the data source name.  For example, if the ``tnsnames.ora`` file
contains a `TNS Alias
<https://python-oracledb.readthedocs.io/en/latest/user_guide/connection_handling.html#tns-aliases-for-connection-strings>`_
of ``myalias`` as below:

.. sourcecode:: text

    myalias =
      (DESCRIPTION =
        (ADDRESS = (PROTOCOL = TCP)(HOST = mymachine.example.com)(PORT = 1521))
        (CONNECT_DATA =
          (SERVER = DEDICATED)
          (SERVICE_NAME = orclpdb1)
        )
      )

The python-oracledb dialect connects to this database service when ``myalias`` is the
hostname portion of the URL, without specifying a port, database name or
``service_name``::

    engine = create_engine("oracle+oracledb://scott:tiger@myalias")

Connecting to Oracle Autonomous Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users of Oracle Autonomous Database should use either use the TNS Alias URL
shown above, or pass the TNS Alias as the ``dsn`` key value in a
:paramref:`_sa.create_engine.connect_args` dictionary.

If Oracle Autonomous Database is configured for mutual TLS ("mTLS")
connections, then additional configuration is required as shown in `Connecting
to Oracle Cloud Autonomous Databases
<https://python-oracledb.readthedocs.io/en/latest/user_guide/connection_handling.html#connecting-to-oracle-cloud-autonomous-databases>`_. In
summary, Thick mode users should configure file locations and set the wallet
path in ``sqlnet.ora`` appropriately::

    e = create_engine(
        "oracle+oracledb://@",
        thick_mode={
            # directory containing tnsnames.ora and cwallet.so
            "config_dir": "/opt/oracle/wallet_dir",
        },
        connect_args={
            "user": "scott",
            "password": "tiger",
            "dsn": "mydb_high",
        },
    )

Thin mode users of mTLS should pass the appropriate directories and PEM wallet
password when creating the engine, similar to::

    e = create_engine(
        "oracle+oracledb://@",
        connect_args={
            "user": "scott",
            "password": "tiger",
            "dsn": "mydb_high",
            "config_dir": "/opt/oracle/wallet_dir",  # directory containing tnsnames.ora
            "wallet_location": "/opt/oracle/wallet_dir",  # directory containing ewallet.pem
            "wallet_password": "top secret",  # password for the PEM file
        },
    )

Typically ``config_dir`` and ``wallet_location`` are the same directory, which
is where the Oracle Autonomous Database wallet zip file was extracted.  Note
this directory should be protected.

Using python-oracledb Connection Pooling
----------------------------------------

The python-oracledb driver provides its own connection pool implementation that
may be used in place of SQLAlchemy's pooling functionality.  The driver pool
gives support for high availability features such as dead connection detection,
connection draining for planned database downtime, support for Oracle
Application Continuity and Transparent Application Continuity, and gives
support for `Database Resident Connection Pooling (DRCP)
<https://python-oracledb.readthedocs.io/en/latest/user_guide/connection_handling.html#database-resident-connection-pooling-drcp>`_.

To take advantage of python-oracledb's pool, use the
:paramref:`_sa.create_engine.creator` parameter to provide a function that
returns a new connection, along with setting
:paramref:`_sa.create_engine.pool_class` to ``NullPool`` to disable
SQLAlchemy's pooling::

    import oracledb
    from sqlalchemy import create_engine
    from sqlalchemy import text
    from sqlalchemy.pool import NullPool

    # Uncomment to use the optional python-oracledb Thick mode.
    # Review the python-oracledb doc for the appropriate parameters
    # oracledb.init_oracle_client(<your parameters>)

    pool = oracledb.create_pool(
        user="scott",
        password="tiger",
        dsn="localhost:1521/freepdb1",
        min=1,
        max=4,
        increment=1,
    )
    engine = create_engine(
        "oracle+oracledb://", creator=pool.acquire, poolclass=NullPool
    )

The above engine may then be used normally. Internally, python-oracledb handles
connection pooling::

    with engine.connect() as conn:
        print(conn.scalar(text("select 1 from dual")))

Refer to the python-oracledb documentation for `oracledb.create_pool()
<https://python-oracledb.readthedocs.io/en/latest/api_manual/module.html#oracledb.create_pool>`_
for the arguments that can be used when creating a connection pool.

.. _drcp:

Using Oracle Database Resident Connection Pooling (DRCP)
--------------------------------------------------------

When using Oracle Database's Database Resident Connection Pooling (DRCP), the
best practice is to specify a connection class and "purity". Refer to the
`python-oracledb documentation on DRCP
<https://python-oracledb.readthedocs.io/en/latest/user_guide/connection_handling.html#database-resident-connection-pooling-drcp>`_.
For example::

    import oracledb
    from sqlalchemy import create_engine
    from sqlalchemy import text
    from sqlalchemy.pool import NullPool

    # Uncomment to use the optional python-oracledb Thick mode.
    # Review the python-oracledb doc for the appropriate parameters
    # oracledb.init_oracle_client(<your parameters>)

    pool = oracledb.create_pool(
        user="scott",
        password="tiger",
        dsn="localhost:1521/freepdb1",
        min=1,
        max=4,
        increment=1,
        cclass="MYCLASS",
        purity=oracledb.PURITY_SELF,
    )
    engine = create_engine(
        "oracle+oracledb://", creator=pool.acquire, poolclass=NullPool
    )

The above engine may then be used normally where python-oracledb handles
application connection pooling and Oracle Database additionally uses DRCP::

    with engine.connect() as conn:
        print(conn.scalar(text("select 1 from dual")))

If you wish to use different connection classes or purities for different
connections, then wrap ``pool.acquire()``::

    import oracledb
    from sqlalchemy import create_engine
    from sqlalchemy import text
    from sqlalchemy.pool import NullPool

    # Uncomment to use python-oracledb Thick mode.
    # Review the python-oracledb doc for the appropriate parameters
    # oracledb.init_oracle_client(<your parameters>)

    pool = oracledb.create_pool(
        user="scott",
        password="tiger",
        dsn="localhost:1521/freepdb1",
        min=1,
        max=4,
        increment=1,
        cclass="MYCLASS",
        purity=oracledb.PURITY_SELF,
    )


    def creator():
        return pool.acquire(cclass="MYOTHERCLASS", purity=oracledb.PURITY_NEW)


    engine = create_engine(
        "oracle+oracledb://", creator=creator, poolclass=NullPool
    )

Engine Options consumed by the SQLAlchemy oracledb dialect outside of the driver
--------------------------------------------------------------------------------

There are also options that are consumed by the SQLAlchemy oracledb dialect
itself.  These options are always passed directly to :func:`_sa.create_engine`,
such as::

    e = create_engine("oracle+oracledb://user:pass@tnsalias", arraysize=500)

The parameters accepted by the oracledb dialect are as follows:

* ``arraysize`` - set the driver cursor.arraysize value. It defaults to
  ``None``, indicating that the driver default value of 100 should be used.
  This setting controls how many rows are buffered when fetching rows, and can
  have a significant effect on performance if increased for queries that return
  large numbers of rows.

  .. versionchanged:: 2.0.26 - changed the default value from 50 to None,
    to use the default value of the driver itself.

* ``auto_convert_lobs`` - defaults to True; See :ref:`oracledb_lob`.

* ``coerce_to_decimal`` - see :ref:`oracledb_numeric` for detail.

* ``encoding_errors`` - see :ref:`oracledb_unicode_encoding_errors` for detail.

.. _oracledb_unicode:

Unicode
-------

As is the case for all DBAPIs under Python 3, all strings are inherently
Unicode strings.

Ensuring the Correct Client Encoding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In python-oracledb, the encoding used for all character data is "UTF-8".

Unicode-specific Column datatypes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Core expression language handles unicode data by use of the
:class:`.Unicode` and :class:`.UnicodeText` datatypes.  These types correspond
to the VARCHAR2 and CLOB Oracle Database datatypes by default.  When using
these datatypes with Unicode data, it is expected that the database is
configured with a Unicode-aware character set so that the VARCHAR2 and CLOB
datatypes can accommodate the data.

In the case that Oracle Database is not configured with a Unicode character
set, the two options are to use the :class:`_types.NCHAR` and
:class:`_oracle.NCLOB` datatypes explicitly, or to pass the flag
``use_nchar_for_unicode=True`` to :func:`_sa.create_engine`, which will cause
the SQLAlchemy dialect to use NCHAR/NCLOB for the :class:`.Unicode` /
:class:`.UnicodeText` datatypes instead of VARCHAR/CLOB.

.. versionchanged:: 1.3 The :class:`.Unicode` and :class:`.UnicodeText`
   datatypes now correspond to the ``VARCHAR2`` and ``CLOB`` Oracle Database
   datatypes unless the ``use_nchar_for_unicode=True`` is passed to the dialect
   when :func:`_sa.create_engine` is called.


.. _oracledb_unicode_encoding_errors:

Encoding Errors
^^^^^^^^^^^^^^^

For the unusual case that data in Oracle Database is present with a broken
encoding, the dialect accepts a parameter ``encoding_errors`` which will be
passed to Unicode decoding functions in order to affect how decoding errors are
handled.  The value is ultimately consumed by the Python `decode
<https://docs.python.org/3/library/stdtypes.html#bytes.decode>`_ function, and
is passed both via python-oracledb's ``encodingErrors`` parameter consumed by
``Cursor.var()``, as well as SQLAlchemy's own decoding function, as the
python-oracledb dialect makes use of both under different circumstances.

.. versionadded:: 1.3.11


.. _oracledb_setinputsizes:

Fine grained control over python-oracledb data binding with setinputsizes
-------------------------------------------------------------------------

The python-oracle DBAPI has a deep and fundamental reliance upon the usage of
the DBAPI ``setinputsizes()`` call.  The purpose of this call is to establish
the datatypes that are bound to a SQL statement for Python values being passed
as parameters.  While virtually no other DBAPI assigns any use to the
``setinputsizes()`` call, the python-oracledb DBAPI relies upon it heavily in
its interactions with the Oracle Database, and in some scenarios it is not
possible for SQLAlchemy to know exactly how data should be bound, as some
settings can cause profoundly different performance characteristics, while
altering the type coercion behavior at the same time.

Users of the oracledb dialect are **strongly encouraged** to read through
python-oracledb's list of built-in datatype symbols at `Database Types
<https://python-oracledb.readthedocs.io/en/latest/api_manual/module.html#database-types>`_
Note that in some cases, significant performance degradation can occur when
using these types vs. not.

On the SQLAlchemy side, the :meth:`.DialectEvents.do_setinputsizes` event can
be used both for runtime visibility (e.g. logging) of the setinputsizes step as
well as to fully control how ``setinputsizes()`` is used on a per-statement
basis.

.. versionadded:: 1.2.9 Added :meth:`.DialectEvents.setinputsizes`


Example 1 - logging all setinputsizes calls
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following example illustrates how to log the intermediary values from a
SQLAlchemy perspective before they are converted to the raw ``setinputsizes()``
parameter dictionary.  The keys of the dictionary are :class:`.BindParameter`
objects which have a ``.key`` and a ``.type`` attribute::

    from sqlalchemy import create_engine, event

    engine = create_engine(
        "oracle+oracledb://scott:tiger@localhost:1521?service_name=freepdb1"
    )


    @event.listens_for(engine, "do_setinputsizes")
    def _log_setinputsizes(inputsizes, cursor, statement, parameters, context):
        for bindparam, dbapitype in inputsizes.items():
            log.info(
                "Bound parameter name: %s  SQLAlchemy type: %r DBAPI object: %s",
                bindparam.key,
                bindparam.type,
                dbapitype,
            )

Example 2 - remove all bindings to CLOB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For performance, fetching LOB datatypes from Oracle Database is set by default
for the ``Text`` type within SQLAlchemy.  This setting can be modified as
follows::


    from sqlalchemy import create_engine, event
    from oracledb import CLOB

    engine = create_engine(
        "oracle+oracledb://scott:tiger@localhost:1521?service_name=freepdb1"
    )


    @event.listens_for(engine, "do_setinputsizes")
    def _remove_clob(inputsizes, cursor, statement, parameters, context):
        for bindparam, dbapitype in list(inputsizes.items()):
            if dbapitype is CLOB:
                del inputsizes[bindparam]

.. _oracledb_lob:

LOB Datatypes
--------------

LOB datatypes refer to the "large object" datatypes such as CLOB, NCLOB and
BLOB. Oracle Database can efficiently return these datatypes as a single
buffer. SQLAlchemy makes use of type handlers to do this by default.

To disable the use of the type handlers and deliver LOB objects as classic
buffered objects with a ``read()`` method, the parameter
``auto_convert_lobs=False`` may be passed to :func:`_sa.create_engine`.

.. _oracledb_returning:

RETURNING Support
-----------------

The oracledb dialect implements RETURNING using OUT parameters.  The dialect
supports RETURNING fully.

Two Phase Transaction Support
-----------------------------

Two phase transactions are fully supported with python-oracledb. (Thin mode
requires python-oracledb 2.3).  APIs for two phase transactions are provided at
the Core level via :meth:`_engine.Connection.begin_twophase` and
:paramref:`_orm.Session.twophase` for transparent ORM use.

.. versionchanged:: 2.0.32 added support for two phase transactions

.. _oracledb_numeric:

Precision Numerics
------------------

SQLAlchemy's numeric types can handle receiving and returning values as Python
``Decimal`` objects or float objects.  When a :class:`.Numeric` object, or a
subclass such as :class:`.Float`, :class:`_oracle.DOUBLE_PRECISION` etc. is in
use, the :paramref:`.Numeric.asdecimal` flag determines if values should be
coerced to ``Decimal`` upon return, or returned as float objects.  To make
matters more complicated under Oracle Database, the ``NUMBER`` type can also
represent integer values if the "scale" is zero, so the Oracle
Database-specific :class:`_oracle.NUMBER` type takes this into account as well.

The oracledb dialect makes extensive use of connection- and cursor-level
"outputtypehandler" callables in order to coerce numeric values as requested.
These callables are specific to the specific flavor of :class:`.Numeric` in
use, as well as if no SQLAlchemy typing objects are present.  There are
observed scenarios where Oracle Database may send incomplete or ambiguous
information about the numeric types being returned, such as a query where the
numeric types are buried under multiple levels of subquery.  The type handlers
do their best to make the right decision in all cases, deferring to the
underlying python-oracledb DBAPI for all those cases where the driver can make
the best decision.

When no typing objects are present, as when executing plain SQL strings, a
default "outputtypehandler" is present which will generally return numeric
values which specify precision and scale as Python ``Decimal`` objects.  To
disable this coercion to decimal for performance reasons, pass the flag
``coerce_to_decimal=False`` to :func:`_sa.create_engine`::

    engine = create_engine(
        "oracle+oracledb://scott:tiger@tnsalias", coerce_to_decimal=False
    )

The ``coerce_to_decimal`` flag only impacts the results of plain string
SQL statements that are not otherwise associated with a :class:`.Numeric`
SQLAlchemy type (or a subclass of such).

.. versionchanged:: 1.2 The numeric handling system for the oracle dialects has
   been reworked to take advantage of newer driver features as well as better
   integration of outputtypehandlers.

.. versionadded:: 2.0.0 added support for the python-oracledb driver.

"""  # noqa
from __future__ import annotations

import collections
import re
from typing import Any
from typing import TYPE_CHECKING

from . import cx_oracle as _cx_oracle
from ... import exc
from ... import pool
from ...connectors.asyncio import AsyncAdapt_dbapi_connection
from ...connectors.asyncio import AsyncAdapt_dbapi_cursor
from ...connectors.asyncio import AsyncAdapt_dbapi_ss_cursor
from ...connectors.asyncio import AsyncAdaptFallback_dbapi_connection
from ...engine import default
from ...util import asbool
from ...util import await_fallback
from ...util import await_only

if TYPE_CHECKING:
    from oracledb import AsyncConnection
    from oracledb import AsyncCursor


class OracleExecutionContext_oracledb(
    _cx_oracle.OracleExecutionContext_cx_oracle
):
    pass


class OracleDialect_oracledb(_cx_oracle.OracleDialect_cx_oracle):
    supports_statement_cache = True
    execution_ctx_cls = OracleExecutionContext_oracledb

    driver = "oracledb"
    _min_version = (1,)

    def __init__(
        self,
        auto_convert_lobs=True,
        coerce_to_decimal=True,
        arraysize=None,
        encoding_errors=None,
        thick_mode=None,
        **kwargs,
    ):
        super().__init__(
            auto_convert_lobs,
            coerce_to_decimal,
            arraysize,
            encoding_errors,
            **kwargs,
        )

        if self.dbapi is not None and (
            thick_mode or isinstance(thick_mode, dict)
        ):
            kw = thick_mode if isinstance(thick_mode, dict) else {}
            self.dbapi.init_oracle_client(**kw)

    @classmethod
    def import_dbapi(cls):
        import oracledb

        return oracledb

    @classmethod
    def is_thin_mode(cls, connection):
        return connection.connection.dbapi_connection.thin

    @classmethod
    def get_async_dialect_cls(cls, url):
        return OracleDialectAsync_oracledb

    def _load_version(self, dbapi_module):
        version = (0, 0, 0)
        if dbapi_module is not None:
            m = re.match(r"(\d+)\.(\d+)(?:\.(\d+))?", dbapi_module.version)
            if m:
                version = tuple(
                    int(x) for x in m.group(1, 2, 3) if x is not None
                )
        self.oracledb_ver = version
        if (
            self.oracledb_ver > (0, 0, 0)
            and self.oracledb_ver < self._min_version
        ):
            raise exc.InvalidRequestError(
                f"oracledb version {self._min_version} and above are supported"
            )

    def do_begin_twophase(self, connection, xid):
        conn_xis = connection.connection.xid(*xid)
        connection.connection.tpc_begin(conn_xis)
        connection.connection.info["oracledb_xid"] = conn_xis

    def do_prepare_twophase(self, connection, xid):
        should_commit = connection.connection.tpc_prepare()
        connection.info["oracledb_should_commit"] = should_commit

    def do_rollback_twophase(
        self, connection, xid, is_prepared=True, recover=False
    ):
        if recover:
            conn_xid = connection.connection.xid(*xid)
        else:
            conn_xid = None
        connection.connection.tpc_rollback(conn_xid)

    def do_commit_twophase(
        self, connection, xid, is_prepared=True, recover=False
    ):
        conn_xid = None
        if not is_prepared:
            should_commit = connection.connection.tpc_prepare()
        elif recover:
            conn_xid = connection.connection.xid(*xid)
            should_commit = True
        else:
            should_commit = connection.info["oracledb_should_commit"]
        if should_commit:
            connection.connection.tpc_commit(conn_xid)

    def do_recover_twophase(self, connection):
        return [
            # oracledb seems to return bytes
            (
                fi,
                gti.decode() if isinstance(gti, bytes) else gti,
                bq.decode() if isinstance(bq, bytes) else bq,
            )
            for fi, gti, bq in connection.connection.tpc_recover()
        ]

    def _check_max_identifier_length(self, connection):
        if self.oracledb_ver >= (2, 5):
            max_len = connection.connection.max_identifier_length
            if max_len is not None:
                return max_len
        return super()._check_max_identifier_length(connection)


class AsyncAdapt_oracledb_cursor(AsyncAdapt_dbapi_cursor):
    _cursor: AsyncCursor
    _awaitable_cursor_close: bool = False

    __slots__ = ()

    @property
    def outputtypehandler(self):
        return self._cursor.outputtypehandler

    @outputtypehandler.setter
    def outputtypehandler(self, value):
        self._cursor.outputtypehandler = value

    def var(self, *args, **kwargs):
        return self._cursor.var(*args, **kwargs)

    def setinputsizes(self, *args: Any, **kwargs: Any) -> Any:
        return self._cursor.setinputsizes(*args, **kwargs)

    def _aenter_cursor(self, cursor: AsyncCursor) -> AsyncCursor:
        try:
            return cursor.__enter__()
        except Exception as error:
            self._adapt_connection._handle_exception(error)

    async def _execute_async(self, operation, parameters):
        # override to not use mutex, oracledb already has a mutex

        if parameters is None:
            result = await self._cursor.execute(operation)
        else:
            result = await self._cursor.execute(operation, parameters)

        if self._cursor.description and not self.server_side:
            self._rows = collections.deque(await self._cursor.fetchall())
        return result

    async def _executemany_async(
        self,
        operation,
        seq_of_parameters,
    ):
        # override to not use mutex, oracledb already has a mutex
        return await self._cursor.executemany(operation, seq_of_parameters)

    def __enter__(self):
        return self

    def __exit__(self, type_: Any, value: Any, traceback: Any) -> None:
        self.close()


class AsyncAdapt_oracledb_ss_cursor(
    AsyncAdapt_dbapi_ss_cursor, AsyncAdapt_oracledb_cursor
):
    __slots__ = ()

    def close(self) -> None:
        if self._cursor is not None:
            self._cursor.close()
            self._cursor = None  # type: ignore


class AsyncAdapt_oracledb_connection(AsyncAdapt_dbapi_connection):
    _connection: AsyncConnection
    __slots__ = ()

    thin = True

    _cursor_cls = AsyncAdapt_oracledb_cursor
    _ss_cursor_cls = None

    @property
    def autocommit(self):
        return self._connection.autocommit

    @autocommit.setter
    def autocommit(self, value):
        self._connection.autocommit = value

    @property
    def outputtypehandler(self):
        return self._connection.outputtypehandler

    @outputtypehandler.setter
    def outputtypehandler(self, value):
        self._connection.outputtypehandler = value

    @property
    def version(self):
        return self._connection.version

    @property
    def stmtcachesize(self):
        return self._connection.stmtcachesize

    @stmtcachesize.setter
    def stmtcachesize(self, value):
        self._connection.stmtcachesize = value

    @property
    def max_identifier_length(self):
        return self._connection.max_identifier_length

    def cursor(self):
        return AsyncAdapt_oracledb_cursor(self)

    def ss_cursor(self):
        return AsyncAdapt_oracledb_ss_cursor(self)

    def xid(self, *args: Any, **kwargs: Any) -> Any:
        return self._connection.xid(*args, **kwargs)

    def tpc_begin(self, *args: Any, **kwargs: Any) -> Any:
        return self.await_(self._connection.tpc_begin(*args, **kwargs))

    def tpc_commit(self, *args: Any, **kwargs: Any) -> Any:
        return self.await_(self._connection.tpc_commit(*args, **kwargs))

    def tpc_prepare(self, *args: Any, **kwargs: Any) -> Any:
        return self.await_(self._connection.tpc_prepare(*args, **kwargs))

    def tpc_recover(self, *args: Any, **kwargs: Any) -> Any:
        return self.await_(self._connection.tpc_recover(*args, **kwargs))

    def tpc_rollback(self, *args: Any, **kwargs: Any) -> Any:
        return self.await_(self._connection.tpc_rollback(*args, **kwargs))


class AsyncAdaptFallback_oracledb_connection(
    AsyncAdaptFallback_dbapi_connection, AsyncAdapt_oracledb_connection
):
    __slots__ = ()


class OracledbAdaptDBAPI:
    def __init__(self, oracledb) -> None:
        self.oracledb = oracledb

        for k, v in self.oracledb.__dict__.items():
            if k != "connect":
                self.__dict__[k] = v

    def connect(self, *arg, **kw):
        async_fallback = kw.pop("async_fallback", False)
        creator_fn = kw.pop("async_creator_fn", self.oracledb.connect_async)

        if asbool(async_fallback):
            return AsyncAdaptFallback_oracledb_connection(
                self, await_fallback(creator_fn(*arg, **kw))
            )

        else:
            return AsyncAdapt_oracledb_connection(
                self, await_only(creator_fn(*arg, **kw))
            )


class OracleExecutionContextAsync_oracledb(OracleExecutionContext_oracledb):
    # restore default create cursor
    create_cursor = default.DefaultExecutionContext.create_cursor

    def create_default_cursor(self):
        # copy of OracleExecutionContext_cx_oracle.create_cursor
        c = self._dbapi_connection.cursor()
        if self.dialect.arraysize:
            c.arraysize = self.dialect.arraysize

        return c

    def create_server_side_cursor(self):
        c = self._dbapi_connection.ss_cursor()
        if self.dialect.arraysize:
            c.arraysize = self.dialect.arraysize

        return c


class OracleDialectAsync_oracledb(OracleDialect_oracledb):
    is_async = True
    supports_server_side_cursors = True
    supports_statement_cache = True
    execution_ctx_cls = OracleExecutionContextAsync_oracledb

    _min_version = (2,)

    # thick_mode mode is not supported by asyncio, oracledb will raise
    @classmethod
    def import_dbapi(cls):
        import oracledb

        return OracledbAdaptDBAPI(oracledb)

    @classmethod
    def get_pool_class(cls, url):
        async_fallback = url.query.get("async_fallback", False)

        if asbool(async_fallback):
            return pool.FallbackAsyncAdaptedQueuePool
        else:
            return pool.AsyncAdaptedQueuePool

    def get_driver_connection(self, connection):
        return connection._connection


dialect = OracleDialect_oracledb
dialect_async = OracleDialectAsync_oracledb
