.. _io.sql:

===========
SQL queries
===========

The :mod:`pandas.io.sql` module provides a collection of query wrappers to both
facilitate data retrieval and to reduce dependency on DB-specific API.

Where available, users may first want to opt for `Apache Arrow ADBC
<https://arrow.apache.org/adbc/current/index.html>`_ drivers. These drivers
should provide the best performance, null handling, and type detection.

  .. versionadded:: 2.2.0

     Added native support for ADBC drivers

For a full list of ADBC drivers and their development status, see the `ADBC Driver
Implementation Status <https://arrow.apache.org/adbc/current/driver/status.html>`_
documentation.

Where an ADBC driver is not available or may be missing functionality,
users should opt for installing SQLAlchemy alongside their database driver library.
Examples of such drivers are `psycopg2 <https://www.psycopg.org/>`__
for PostgreSQL or `pymysql <https://github.com/PyMySQL/PyMySQL>`__ for MySQL.
For `SQLite <https://docs.python.org/3/library/sqlite3.html>`__ this is
included in Python's standard library by default.
You can find an overview of supported drivers for each SQL dialect in the
`SQLAlchemy docs <https://docs.sqlalchemy.org/en/latest/dialects/index.html>`__.

If SQLAlchemy is not installed, you can use a :class:`sqlite3.Connection` in place of
a SQLAlchemy engine, connection, or URI string.

See also some :ref:`cookbook examples <cookbook.sql>` for some advanced strategies.

The key functions are:

.. currentmodule:: pandas
.. autosummary::

    read_sql_table
    read_sql_query
    read_sql
    DataFrame.to_sql

.. note::

    The function :func:`~pandas.read_sql` is a convenience wrapper around
    :func:`~pandas.read_sql_table` and :func:`~pandas.read_sql_query` (and for
    backward compatibility) and will delegate to specific function depending on
    the provided input (database table name or sql query).
    Table names do not need to be quoted if they have special characters.

In the following example, we use the `SQlite <https://www.sqlite.org/index.html>`__ SQL database
engine. You can use a temporary SQLite database where data are stored in
"memory".

To connect using an ADBC driver you will want to install the ``adbc_driver_sqlite`` using your
package manager. Once installed, you can use the DBAPI interface provided by the ADBC driver
to connect to your database.

.. code-block:: python

   import adbc_driver_sqlite.dbapi as sqlite_dbapi

   # Create the connection
   with sqlite_dbapi.connect("sqlite:///:memory:") as conn:
        df = pd.read_sql_table("data", conn)

To connect with SQLAlchemy you use the :func:`create_engine` function to create an engine
object from database URI. You only need to create the engine once per database you are
connecting to.
For more information on :func:`create_engine` and the URI formatting, see the examples
below and the SQLAlchemy `documentation <https://docs.sqlalchemy.org/en/latest/core/engines.html>`__

.. ipython:: python

   from sqlalchemy import create_engine

   # Create your engine.
   engine = create_engine("sqlite:///:memory:")

If you want to manage your own connections you can pass one of those instead. The example below opens a
connection to the database using a Python context manager that automatically closes the connection after
the block has completed.
See the `SQLAlchemy docs <https://docs.sqlalchemy.org/en/latest/core/connections.html#basic-usage>`__
for an explanation of how the database connection is handled.

.. code-block:: python

   with engine.connect() as conn, conn.begin():
       data = pd.read_sql_table("data", conn)

.. warning::

        When you open a connection to a database you are also responsible for closing it.
        Side effects of leaving a connection open may include locking the database or
        other breaking behaviour.

Writing DataFrames
''''''''''''''''''

Assuming the following data is in a ``DataFrame`` ``data``, we can insert it into
the database using :func:`~pandas.DataFrame.to_sql`.

+-----+------------+-------+-------+-------+
| id  |    Date    | Col_1 | Col_2 | Col_3 |
+=====+============+=======+=======+=======+
| 26  | 2012-10-18 |   X   |  25.7 | True  |
+-----+------------+-------+-------+-------+
| 42  | 2012-10-19 |   Y   | -12.4 | False |
+-----+------------+-------+-------+-------+
| 63  | 2012-10-20 |   Z   |  5.73 | True  |
+-----+------------+-------+-------+-------+


.. ipython:: python

   import datetime

   c = ["id", "Date", "Col_1", "Col_2", "Col_3"]
   d = [
       (26, datetime.datetime(2010, 10, 18), "X", 27.5, True),
       (42, datetime.datetime(2010, 10, 19), "Y", -12.5, False),
       (63, datetime.datetime(2010, 10, 20), "Z", 5.73, True),
   ]

   data = pd.DataFrame(d, columns=c)

   data
   data.to_sql("data", con=engine)

With some databases, writing large DataFrames can result in errors due to
packet size limitations being exceeded. This can be avoided by setting the
``chunksize`` parameter when calling ``to_sql``.  For example, the following
writes ``data`` to the database in batches of 1000 rows at a time:

.. ipython:: python

    data.to_sql("data_chunked", con=engine, chunksize=1000)

SQL data types
++++++++++++++

Ensuring consistent data type management across SQL databases is challenging.
Not every SQL database offers the same types, and even when they do the implementation
of a given type can vary in ways that have subtle effects on how types can be
preserved.

For the best odds at preserving database types users are advised to use
ADBC drivers when available. The Arrow type system offers a wider array of
types that more closely match database types than the historical pandas/NumPy
type system. To illustrate, note this (non-exhaustive) listing of types
available in different databases and pandas backends:

+-----------------+-----------------------+----------------+---------+
|numpy/pandas     |arrow                  |postgres        |sqlite   |
+=================+=======================+================+=========+
|int16/Int16      |int16                  |SMALLINT        |INTEGER  |
+-----------------+-----------------------+----------------+---------+
|int32/Int32      |int32                  |INTEGER         |INTEGER  |
+-----------------+-----------------------+----------------+---------+
|int64/Int64      |int64                  |BIGINT          |INTEGER  |
+-----------------+-----------------------+----------------+---------+
|float32          |float32                |REAL            |REAL     |
+-----------------+-----------------------+----------------+---------+
|float64          |float64                |DOUBLE PRECISION|REAL     |
+-----------------+-----------------------+----------------+---------+
|object           |string                 |TEXT            |TEXT     |
+-----------------+-----------------------+----------------+---------+
|bool             |``bool_``              |BOOLEAN         |         |
+-----------------+-----------------------+----------------+---------+
|datetime64[ns]   |timestamp(us)          |TIMESTAMP       |         |
+-----------------+-----------------------+----------------+---------+
|datetime64[ns,tz]|timestamp(us,tz)       |TIMESTAMPTZ     |         |
+-----------------+-----------------------+----------------+---------+
|                 |date32                 |DATE            |         |
+-----------------+-----------------------+----------------+---------+
|                 |month_day_nano_interval|INTERVAL        |         |
+-----------------+-----------------------+----------------+---------+
|                 |binary                 |BINARY          |BLOB     |
+-----------------+-----------------------+----------------+---------+
|                 |decimal128             |DECIMAL [#f1]_  |         |
+-----------------+-----------------------+----------------+---------+
|                 |list                   |ARRAY [#f1]_    |         |
+-----------------+-----------------------+----------------+---------+
|                 |struct                 |COMPOSITE TYPE  |         |
|                 |                       | [#f1]_         |         |
+-----------------+-----------------------+----------------+---------+

.. rubric:: Footnotes

.. [#f1] Not implemented as of writing, but theoretically possible

If you are interested in preserving database types as best as possible
throughout the lifecycle of your DataFrame, users are encouraged to
leverage the ``dtype_backend="pyarrow"`` argument of :func:`~pandas.read_sql`

.. code-block:: ipython

   # for roundtripping
   with pg_dbapi.connect(uri) as conn:
       df2 = pd.read_sql("pandas_table", conn, dtype_backend="pyarrow")

This will prevent your data from being converted to the traditional pandas/NumPy
type system, which often converts SQL types in ways that make them impossible to
round-trip.

In case an ADBC driver is not available, :func:`~pandas.DataFrame.to_sql`
will try to map your data to an appropriate SQL data type based on the dtype of
the data. When you have columns of dtype ``object``, pandas will try to infer
the data type.

You can always override the default type by specifying the desired SQL type of
any of the columns by using the ``dtype`` argument. This argument needs a
dictionary mapping column names to SQLAlchemy types (or strings for the sqlite3
fallback mode).
For example, specifying to use the sqlalchemy ``String`` type instead of the
default ``Text`` type for string columns:

.. ipython:: python

    from sqlalchemy.types import String

    data.to_sql("data_dtype", con=engine, dtype={"Col_1": String})

.. note::

    Due to the limited support for timedelta's in the different database
    flavors, columns with type ``timedelta64`` will be written as integer
    values as nanoseconds to the database and a warning will be raised. The only
    exception to this is when using the ADBC PostgreSQL driver in which case a
    timedelta will be written to the database as an ``INTERVAL``

.. note::

    Columns of ``category`` dtype will be converted to the dense representation
    as you would get with ``np.asarray(categorical)`` (e.g. for string categories
    this gives an array of strings).
    Because of this, reading the database table back in does **not** generate
    a categorical.

.. _io.sql_datetime_data:

Datetime data types
'''''''''''''''''''

Using ADBC or SQLAlchemy, :func:`~pandas.DataFrame.to_sql` is capable of writing
datetime data that is timezone naive or timezone aware. However, the resulting
data stored in the database ultimately depends on the supported data type
for datetime data of the database system being used.

The following table lists supported data types for datetime data for some
common databases. Other database dialects may have different data types for
datetime data.

===========   =============================================  ===================
Database      SQL Datetime Types                             Timezone Support
===========   =============================================  ===================
SQLite        ``TEXT``                                       No
MySQL         ``TIMESTAMP`` or ``DATETIME``                  No
PostgreSQL    ``TIMESTAMP`` or ``TIMESTAMP WITH TIME ZONE``  Yes
===========   =============================================  ===================

When writing timezone aware data to databases that do not support timezones,
the data will be written as timezone naive timestamps that are in local time
with respect to the timezone.

:func:`~pandas.read_sql_table` is also capable of reading datetime data that is
timezone aware or naive. When reading ``TIMESTAMP WITH TIME ZONE`` types, pandas
will convert the data to UTC.

.. _io.sql.method:

Insertion method
++++++++++++++++

The parameter ``method`` controls the SQL insertion clause used.
Possible values are:

- ``None``: Uses standard SQL ``INSERT`` clause (one per row).
- ``'multi'``: Pass multiple values in a single ``INSERT`` clause.
  It uses a *special* SQL syntax not supported by all backends.
  This usually provides better performance for analytic databases
  like *Presto* and *Redshift*, but has worse performance for
  traditional SQL backend if the table contains many columns.
  For more information check the SQLAlchemy `documentation
  <https://docs.sqlalchemy.org/en/latest/core/dml.html#sqlalchemy.sql.expression.Insert.values.params.*args>`__.
- callable with signature ``(pd_table, conn, keys, data_iter)``:
  This can be used to implement a more performant insertion method based on
  specific backend dialect features.

Example of a callable using PostgreSQL `COPY clause
<https://www.postgresql.org/docs/current/sql-copy.html>`__::

  # Alternative to_sql() *method* for DBs that support COPY FROM
  import csv
  from io import StringIO

  def psql_insert_copy(table, conn, keys, data_iter):
      """
      Execute SQL statement inserting data

      Parameters
      ----------
      table : pandas.io.sql.SQLTable
      conn : sqlalchemy.engine.Engine or sqlalchemy.engine.Connection
      keys : list of str
          Column names
      data_iter : Iterable that iterates the values to be inserted
      """
      # gets a DBAPI connection that can provide a cursor
      dbapi_conn = conn.connection
      with dbapi_conn.cursor() as cur:
          s_buf = StringIO()
          writer = csv.writer(s_buf)
          writer.writerows(data_iter)
          s_buf.seek(0)

          columns = ', '.join(['"{}"'.format(k) for k in keys])
          if table.schema:
              table_name = '{}.{}'.format(table.schema, table.name)
          else:
              table_name = table.name

          sql = 'COPY {} ({}) FROM STDIN WITH CSV'.format(
              table_name, columns)
          cur.copy_expert(sql=sql, file=s_buf)

Reading tables
''''''''''''''

:func:`~pandas.read_sql_table` will read a database table given the
table name and optionally a subset of columns to read.

.. note::

    In order to use :func:`~pandas.read_sql_table`, you **must** have the
    ADBC driver or SQLAlchemy optional dependency installed.

.. ipython:: python

   pd.read_sql_table("data", engine)

.. note::

  ADBC drivers will map database types directly back to arrow types. For other drivers
  note that pandas infers column dtypes from query outputs, and not by looking
  up data types in the physical database schema. For example, assume ``userid``
  is an integer column in a table. Then, intuitively, ``select userid ...`` will
  return integer-valued series, while ``select cast(userid as text) ...`` will
  return object-valued (str) series. Accordingly, if the query output is empty,
  then all resulting columns will be returned as object-valued (since they are
  most general). If you foresee that your query will sometimes generate an empty
  result, you may want to explicitly typecast afterwards to ensure dtype
  integrity.

You can also specify the name of the column as the ``DataFrame`` index,
and specify a subset of columns to be read.

.. ipython:: python

   pd.read_sql_table("data", engine, index_col="id")
   pd.read_sql_table("data", engine, columns=["Col_1", "Col_2"])

And you can explicitly force columns to be parsed as dates:

.. ipython:: python

   pd.read_sql_table("data", engine, parse_dates=["Date"])

If needed you can explicitly specify a format string, or a dict of arguments
to pass to :func:`pandas.to_datetime`:

.. code-block:: python

   pd.read_sql_table("data", engine, parse_dates={"Date": "%Y-%m-%d"})
   pd.read_sql_table(
       "data",
       engine,
       parse_dates={"Date": {"format": "%Y-%m-%d %H:%M:%S"}},
   )


You can check if a table exists using :func:`~pandas.io.sql.has_table`

Schema support
''''''''''''''

Reading from and writing to different schemas is supported through the ``schema``
keyword in the :func:`~pandas.read_sql_table` and :func:`~pandas.DataFrame.to_sql`
functions. Note however that this depends on the database flavor (sqlite does not
have schemas). For example:

.. code-block:: python

   df.to_sql(name="table", con=engine, schema="other_schema")
   pd.read_sql_table("table", engine, schema="other_schema")

Querying
''''''''

You can query using raw SQL in the :func:`~pandas.read_sql_query` function.
In this case you must use the SQL variant appropriate for your database.
When using SQLAlchemy, you can also pass SQLAlchemy Expression language constructs,
which are database-agnostic.

.. ipython:: python

   pd.read_sql_query("SELECT * FROM data", engine)

Of course, you can specify a more "complex" query.

.. ipython:: python

   pd.read_sql_query("SELECT id, Col_1, Col_2 FROM data WHERE id = 42;", engine)

The :func:`~pandas.read_sql_query` function supports a ``chunksize`` argument.
Specifying this will return an iterator through chunks of the query result:

.. ipython:: python

    df = pd.DataFrame(np.random.randn(20, 3), columns=list("abc"))
    df.to_sql(name="data_chunks", con=engine, index=False)

.. ipython:: python

    for chunk in pd.read_sql_query("SELECT * FROM data_chunks", engine, chunksize=5):
        print(chunk)


Engine connection examples
''''''''''''''''''''''''''

To connect with SQLAlchemy you use the :func:`create_engine` function to create an engine
object from database URI. You only need to create the engine once per database you are
connecting to.

.. code-block:: python

   from sqlalchemy import create_engine

   engine = create_engine("postgresql://scott:tiger@localhost:5432/mydatabase")

   engine = create_engine("mysql+mysqldb://scott:tiger@localhost/foo")

   engine = create_engine("oracle://scott:tiger@127.0.0.1:1521/sidname")

   engine = create_engine("mssql+pyodbc://mydsn")

   # sqlite://<nohostname>/<path>
   # where <path> is relative:
   engine = create_engine("sqlite:///foo.db")

   # or absolute, starting with a slash:
   engine = create_engine("sqlite:////absolute/path/to/foo.db")

For more information see the examples the SQLAlchemy `documentation <https://docs.sqlalchemy.org/en/latest/core/engines.html>`__


Advanced SQLAlchemy queries
'''''''''''''''''''''''''''

You can use SQLAlchemy constructs to describe your query.

Use :func:`sqlalchemy.text` to specify query parameters in a backend-neutral way

.. ipython:: python

   import sqlalchemy as sa

   pd.read_sql(
       sa.text("SELECT * FROM data where Col_1=:col1"), engine, params={"col1": "X"}
   )

If you have an SQLAlchemy description of your database you can express where conditions using SQLAlchemy expressions

.. ipython:: python

   metadata = sa.MetaData()
   data_table = sa.Table(
       "data",
       metadata,
       sa.Column("index", sa.Integer),
       sa.Column("Date", sa.DateTime),
       sa.Column("Col_1", sa.String),
       sa.Column("Col_2", sa.Float),
       sa.Column("Col_3", sa.Boolean),
   )

   pd.read_sql(sa.select(data_table).where(data_table.c.Col_3 is True), engine)

You can combine SQLAlchemy expressions with parameters passed to :func:`read_sql` using :func:`sqlalchemy.bindparam`

.. ipython:: python

    import datetime as dt

    expr = sa.select(data_table).where(data_table.c.Date > sa.bindparam("date"))
    pd.read_sql(expr, engine, params={"date": dt.datetime(2010, 10, 18)})


Sqlite fallback
'''''''''''''''

The use of sqlite is supported without using SQLAlchemy.
This mode requires a Python database adapter which respect the `Python
DB-API <https://www.python.org/dev/peps/pep-0249/>`__.

You can create connections like so:

.. code-block:: python

   import sqlite3

   con = sqlite3.connect(":memory:")

And then issue the following queries:

.. code-block:: python

   data.to_sql("data", con)
   pd.read_sql_query("SELECT * FROM data", con)
