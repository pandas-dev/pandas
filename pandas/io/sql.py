"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""
from __future__ import print_function
from datetime import datetime, date
import warnings
from pandas.compat import range, lzip, map, zip, raise_with_traceback
import numpy as np


from pandas.core.api import DataFrame
from pandas.core.base import PandasObject
from pandas.tseries.tools import to_datetime


class SQLAlchemyRequired(ImportError):
    pass


class DatabaseError(IOError):
    pass


#------------------------------------------------------------------------------
# Helper execution functions

def _convert_params(sql, params):
    """convert sql and params args to DBAPI2.0 compliant format"""
    args = [sql]
    if params is not None:
        args += list(params)
    return args


def execute(sql, con, cur=None, params=None, flavor='sqlite'):
    """
    Execute the given SQL query using the provided connection object.

    Parameters
    ----------
    sql: string
        Query to be executed
    con: SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object is given, a supported SQL flavor must also be provided
    cur: depreciated, cursor is obtained from connection
    params: list or tuple, optional
        List of parameters to pass to execute method.
    flavor : string "sqlite", "mysql"
        Specifies the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2 connection.
    Returns
    -------
    Results Iterable
    """
    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    args = _convert_params(sql, params)
    return pandas_sql.execute(*args)


def tquery(sql, con, cur=None, params=None, flavor='sqlite'):
    """
    Returns list of tuples corresponding to each row in given sql
    query.

    If only one column selected, then plain list is returned.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object is given, a supported SQL flavor must also be provided
    cur: depreciated, cursor is obtained from connection
    params: list or tuple, optional
        List of parameters to pass to execute method.
    flavor : string "sqlite", "mysql"
        Specifies the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2
        connection.
    Returns
    -------
    Results Iterable
    """
    warnings.warn(
        "tquery is depreciated, and will be removed in future versions",
        DeprecationWarning)

    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    args = _convert_params(sql, params)
    return pandas_sql.tquery(*args)


def uquery(sql, con, cur=None, params=None, engine=None, flavor='sqlite'):
    """
    Does the same thing as tquery, but instead of returning results, it
    returns the number of rows affected.  Good for update queries.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object is given, a supported SQL flavor must also be provided
    cur: depreciated, cursor is obtained from connection
    params: list or tuple, optional
        List of parameters to pass to execute method.
    flavor : string "sqlite", "mysql"
        Specifies the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2
        connection.
    Returns
    -------
    Number of affected rows
    """
    warnings.warn(
        "uquery is depreciated, and will be removed in future versions",
        DeprecationWarning)
    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    args = _convert_params(sql, params)
    return pandas_sql.uquery(*args)


#------------------------------------------------------------------------------
# Read and write to DataFrames


def read_sql(sql, con, index_col=None, flavor='sqlite', coerce_float=True,
             params=None, parse_dates=None):
    """
    Returns a DataFrame corresponding to the result set of the query
    string.

    Optionally provide an index_col parameter to use one of the
    columns as the index, otherwise default integer index will be used

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object is given, a supported SQL flavor must also be provided
    index_col: string, optional
        column name to use for the returned DataFrame object.
    flavor : string specifying the flavor of SQL to use. Ignored when using
        SQLAlchemy engine. Required when using DBAPI2 connection.
    coerce_float : boolean, default True
        Attempt to convert values to non-string, non-numeric objects (like
        decimal.Decimal) to floating point, useful for SQL result sets
    cur: depreciated, cursor is obtained from connection
    params: list or tuple, optional
        List of parameters to pass to execute method.
    parse_dates: list or dict
        List of column names to parse as dates
        Or
        Dict of {column_name: format string} where format string is
        strftime compatible in case of parsing string times or is one of
        (D, s, ns, ms, us) in case of parsing integer timestamps
        Or
        Dict of {column_name: arg dict}, where the arg dict corresponds
        to the keyword arguments of :func:`pandas.tseries.tools.to_datetime`
        Especially useful with databases without native Datetime support,
        such as SQLite

    Returns
    -------
    DataFrame
    """
    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    return pandas_sql.read_sql(sql,
                               index_col=index_col,
                               params=params,
                               coerce_float=coerce_float,
                               parse_dates=parse_dates)


def to_sql(frame, name, con, flavor='sqlite', if_exists='fail', index=True):
    """
    Write records stored in a DataFrame to a SQL database.

    Parameters
    ----------
    frame: DataFrame
    name: name of SQL table
    con: SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object is given, a supported SQL flavor must also be provided
    flavor: {'sqlite', 'mysql', 'postgres'}, default 'sqlite'
        ignored when SQLAlchemy engine. Required when using DBAPI2 connection.
    if_exists: {'fail', 'replace', 'append'}, default 'fail'
        fail: If table exists, do nothing.
        replace: If table exists, drop it, recreate it, and insert data.
        append: If table exists, insert data. Create if does not exist.
    index : boolean, default True
        Write DataFrame index as an column
    """
    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    pandas_sql.to_sql(frame, name, if_exists=if_exists, index=index)


def has_table(table_name, con, meta=None, flavor='sqlite'):
    """
    Check if DB has named table

    Parameters
    ----------
    frame: DataFrame
    name: name of SQL table
    con: SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object is given, a supported SQL flavor name must also be provided
    flavor: {'sqlite', 'mysql'}, default 'sqlite', ignored when using engine
    Returns
    -------
    boolean
    """
    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    return pandas_sql.has_table(table_name)


def read_table(table_name, con, meta=None, index_col=None, coerce_float=True,
               parse_dates=None, columns=None):
    """Given a table name and SQLAlchemy engine, return a DataFrame.
    Type convertions will be done automatically

    Parameters
    ----------
    table_name: name of SQL table in database
    con: SQLAlchemy engine. Legacy mode not supported
    meta: SQLAlchemy meta, optional. If omitted MetaData is reflected from engine
    index_col: column to set as index, optional
    coerce_float : boolean, default True
        Attempt to convert values to non-string, non-numeric objects (like
        decimal.Decimal) to floating point. Can result in loss of Precision.
    parse_dates: list or dict
        List of column names to parse as dates
        Or
        Dict of {column_name: format string} where format string is
        strftime compatible in case of parsing string times or is one of
        (D, s, ns, ms, us) in case of parsing integer timestamps
        Or
        Dict of {column_name: arg dict}, where the arg dict corresponds
        to the keyword arguments of :func:`pandas.tseries.tools.to_datetime`
        Especially useful with databases without native Datetime support,
        such as SQLite
    columns: list
        List of column names to select from sql table
    Returns
    -------
    DataFrame
    """
    pandas_sql = PandasSQLAlchemy(con, meta=meta)
    table = pandas_sql.read_table(table_name,
                                  index_col=index_col,
                                  coerce_float=coerce_float,
                                  parse_dates=parse_dates)

    if table is not None:
        return table
    else:
        raise ValueError("Table %s not found" % table_name, con)


def pandasSQL_builder(con, flavor=None, meta=None):
    """
    Convenience function to return the correct PandasSQL subclass based on the
    provided parameters
    """
    try:
        import sqlalchemy

        if isinstance(con, sqlalchemy.engine.Engine):
            return PandasSQLAlchemy(con, meta=meta)
        else:
            warnings.warn(
                "Not an SQLAlchemy engine, attempting to use as legacy DBAPI connection")
            if flavor is None:
                raise ValueError("""PandasSQL must be created with an SQLAlchemy engine
                    or a DBAPI2 connection and SQL flavour""")
            else:
                return PandasSQLLegacy(con, flavor)

    except ImportError:
        warnings.warn("SQLAlchemy not installed, using legacy mode")
        if flavor is None:
            raise SQLAlchemyRequired
        else:
            return PandasSQLLegacy(con, flavor)


def _safe_col_name(col_name):
    return col_name.strip().replace(' ', '_')


def _parse_date_column(col, format=None):
    if isinstance(format, dict):
        return to_datetime(col, **format)
    else:
        if format in ['D', 's', 'ms', 'us', 'ns']:
            return to_datetime(col, coerce=True, unit=format)
        elif issubclass(col.dtype.type, np.floating) or issubclass(col.dtype.type, np.integer):
            # parse dates as timestamp
            format = 's' if format is None else format
            return to_datetime(col, coerce=True, unit=format)
        else:
            return to_datetime(col, coerce=True, format=format)


def _frame_from_data_and_columns(data, columns, index_col=None,
                                 coerce_float=True):
    df = DataFrame.from_records(
        data, columns=columns, coerce_float=coerce_float)
    if index_col is not None:
        df.set_index(index_col, inplace=True)
    return df


class PandasSQLTable(PandasObject):

    def __init__(self, name, pandas_sql_engine, frame=None, index=True, if_exists='fail', prefix='pandas'):
        self.name = name
        self.pd_sql = pandas_sql_engine
        self.prefix = prefix
        self.frame = frame
        self.index = self._index_name(index)

        if frame is not None:
            # We want to write a frame
            if self.pd_sql.has_table(self.name):
                if if_exists == 'fail':
                    raise ValueError("Table '%s' already exists." % name)
                elif if_exists == 'replace':
                    self.pd_sql.drop_table(self.name)
                    self.table = self._create_table_statement()
                    self.create()
                elif if_exists == 'append':
                    self.table = self.pd_sql.get_table(self.name)
                    if self.table is None:
                        self.table = self._create_table_statement()
            else:
                self.table = self._create_table_statement()
                self.create()
        else:
            # no data provided, read-only mode
            self.table = self.pd_sql.get_table(self.name)

        if self.table is None:
            raise ValueError("Could not init table '%s'" % name)

    def exists(self):
        return self.pd_sql.has_table(self.name)

    def sql_schema(self):
        return str(self.table.compile())

    def create(self):
        self.table.create()

    def insert_statement(self):
        return self.table.insert()

    def maybe_asscalar(self, i):
        try:
            return np.asscalar(i)
        except AttributeError:
            return i

    def insert(self):
        ins = self.insert_statement()

        for t in self.frame.iterrows():
            data = dict((k, self.maybe_asscalar(v))
                        for k, v in t[1].iteritems())
            if self.index is not None:
                data[self.index] = self.maybe_asscalar(t[0])
            self.pd_sql.execute(ins, **data)

    def read(self, coerce_float=True, parse_dates=None, columns=None):

        if columns is not None and len(columns) > 0:
            from sqlalchemy import select
            cols = [self.table.c[n] for n in columns]
            if self.index is not None:
                cols.insert(0, self.table.c[self.index])
            sql_select = select(cols)
        else:
            sql_select = self.table.select()

        result = self.pd_sql.execute(sql_select)
        data = result.fetchall()
        column_names = result.keys()

        self.frame = _frame_from_data_and_columns(data, column_names,
                                                  index_col=self.index,
                                                  coerce_float=coerce_float)

        self._harmonize_columns(parse_dates=parse_dates)

        # Assume that if the index was in prefix_index format, we gave it a name
        # and should return it nameless
        if self.index == self.prefix + '_index':
            self.frame.index.name = None

        return self.frame

    def _index_name(self, index):
        if index is True:
            if self.frame.index.name is not None:
                return _safe_col_name(self.frame.index.name)
            else:
                return self.prefix + '_index'
        elif isinstance(index, basestring):
            return index
        else:
            return None

    def _create_table_statement(self):
        from sqlalchemy import Table, Column

        safe_columns = map(_safe_col_name, self.frame.dtypes.index)
        column_types = map(self._sqlalchemy_type, self.frame.dtypes)

        columns = [Column(name, typ)
                   for name, typ in zip(safe_columns, column_types)]

        if self.index is not None:
            columns.insert(0, Column(self.index,
                                     self._sqlalchemy_type(
                                         self.frame.index.dtype),
                                     index=True))

        return Table(self.name, self.pd_sql.meta, *columns)

    def _harmonize_columns(self, parse_dates=None):
        """ Make a data_frame's column type align with an sql_table column types
            Need to work around limited NA value support.
            Floats are always fine, ints must always
            be floats if there are Null values.
            Booleans are hard because converting bool column with None replaces
            all Nones with false. Therefore only convert bool if there are no NA
            values.
            Datetimes should already be converted
            to np.datetime if supported, but here we also force conversion
            if required
        """
        # handle non-list entries for parse_dates gracefully
        if parse_dates is True or parse_dates is None or parse_dates is False:
            parse_dates = []

        if not hasattr(parse_dates, '__iter__'):
            parse_dates = [parse_dates]

        for sql_col in self.table.columns:
            col_name = sql_col.name
            try:
                df_col = self.frame[col_name]
                # the type the dataframe column should have
                col_type = self._numpy_type(sql_col.type)

                if col_type is datetime or col_type is date:
                    if not issubclass(df_col.dtype.type, np.datetime64):
                        self.frame[col_name] = _parse_date_column(df_col)

                elif col_type is float:
                    # floats support NA, can always convert!
                    self.frame[col_name].astype(col_type, copy=False)

                elif len(df_col) == df_col.count():
                    # No NA values, can convert ints and bools
                    if col_type is int or col_type is bool:
                        self.frame[col_name].astype(col_type, copy=False)

                # Handle date parsing
                if col_name in parse_dates:
                    try:
                        fmt = parse_dates[col_name]
                    except TypeError:
                        fmt = None
                    self.frame[col_name] = _parse_date_column(
                        df_col, format=fmt)

            except KeyError:
                pass  # this column not in results

    def _sqlalchemy_type(self, dtype):
        from sqlalchemy.types import Integer, Float, Text, Boolean, DateTime, Date

        pytype = dtype.type

        if pytype is date:
            return Date
        if issubclass(pytype, np.datetime64) or pytype is datetime:
            # Caution: np.datetime64 is also a subclass of np.number.
            return DateTime
        if issubclass(pytype, np.floating):
            return Float
        if issubclass(pytype, np.integer):
            # TODO: Refine integer size.
            return Integer
        if issubclass(pytype, np.bool_):
            return Boolean
        return Text

    def _numpy_type(self, sqltype):
        from sqlalchemy.types import Integer, Float, Boolean, DateTime, Date

        if isinstance(sqltype, Float):
            return float
        if isinstance(sqltype, Integer):
            # TODO: Refine integer size.
            return int
        if isinstance(sqltype, DateTime):
            # Caution: np.datetime64 is also a subclass of np.number.
            return datetime
        if isinstance(sqltype, Date):
            return date
        if isinstance(sqltype, Boolean):
            return bool
        return object


class PandasSQL(PandasObject):

    """
    Subclasses Should define read_sql and to_sql
    """

    def read_sql(self, *args, **kwargs):
        raise ValueError(
            "PandasSQL must be created with an SQLAlchemy engine or connection+sql flavor")

    def to_sql(self, *args, **kwargs):
        raise ValueError(
            "PandasSQL must be created with an SQLAlchemy engine or connection+sql flavor")

    def _parse_date_columns(self, data_frame, parse_dates):
        """ Force non-datetime columns to be read as such.
            Supports both string formatted and integer timestamp columns
        """
        # handle non-list entries for parse_dates gracefully
        if parse_dates is True or parse_dates is None or parse_dates is False:
            parse_dates = []

        if not hasattr(parse_dates, '__iter__'):
            parse_dates = [parse_dates]

        for col_name in parse_dates:
            df_col = data_frame[col_name]
            try:
                fmt = parse_dates[col_name]
            except TypeError:
                fmt = None
            data_frame[col_name] = _parse_date_column(df_col, format=fmt)

        return data_frame


class PandasSQLAlchemy(PandasSQL):

    """
    This class enables convertion between DataFrame and SQL databases
    using SQLAlchemy to handle DataBase abstraction
    """

    def __init__(self, engine, meta=None):
        self.engine = engine
        if not meta:
            from sqlalchemy.schema import MetaData
            meta = MetaData(self.engine)
            meta.reflect(self.engine)

        self.meta = meta

    def execute(self, *args, **kwargs):
        """Simple passthrough to SQLAlchemy engine"""
        return self.engine.execute(*args, **kwargs)

    def tquery(self, *args, **kwargs):
        result = self.execute(*args, **kwargs)
        return result.fetchall()

    def uquery(self, *args, **kwargs):
        result = self.execute(*args, **kwargs)
        return result.rowcount

    def read_sql(self, sql, index_col=None, coerce_float=True, parse_dates=None, params=None):
        args = _convert_params(sql, params)
        result = self.execute(*args)
        data = result.fetchall()
        columns = result.keys()

        data_frame = _frame_from_data_and_columns(data, columns,
                                                  index_col=index_col,
                                                  coerce_float=coerce_float)

        return self._parse_date_columns(data_frame, parse_dates)

    def to_sql(self, frame, name, if_exists='fail', index=True):
        table = PandasSQLTable(
            name, self, frame=frame, index=index, if_exists=if_exists)
        table.insert()

    def has_table(self, name):
        return self.engine.has_table(name)

    def get_table(self, table_name):
        if self.engine.has_table(table_name):
            return self.meta.tables[table_name]
        else:
            return None

    def read_table(self, table_name, index_col=None, coerce_float=True,
                   parse_dates=None, columns=None):

        table = PandasSQLTable(table_name, self, index=index_col)
        return table.read(coerce_float=parse_dates,
                          parse_dates=parse_dates, columns=columns)

    def drop_table(self, table_name):
        if self.engine.has_table(table_name):
            self.get_table(table_name).drop()
            self.meta.clear()
            self.meta.reflect()

    def _create_sql_schema(self, frame, table_name):
        table = PandasSQLTable(table_name, self, frame=frame)
        return str(table.compile())


# ---- SQL without SQLAlchemy ---
# Flavour specific sql strings and handler class for access to DBs without
# SQLAlchemy installed
# SQL type convertions for each DB
_SQL_TYPES = {
    'text': {
        'mysql': 'VARCHAR (63)',
        'sqlite': 'TEXT',
    },
    'float': {
        'mysql': 'FLOAT',
        'sqlite': 'REAL',
    },
    'int': {
        'mysql': 'BIGINT',
        'sqlite': 'INTEGER',
    },
    'datetime': {
        'mysql': 'DATETIME',
        'sqlite': 'TIMESTAMP',
    },
    'date': {
        'mysql': 'DATE',
        'sqlite': 'TIMESTAMP',
    },
    'bool': {
        'mysql': 'BOOLEAN',
        'sqlite': 'INTEGER',
    }
}

# SQL enquote and wildcard symbols
_SQL_SYMB = {
    'mysql': {
        'br_l': '`',
        'br_r': '`',
        'wld': '%s'
    },
    'sqlite': {
        'br_l': '[',
        'br_r': ']',
        'wld': '?'
    }
}


class PandasSQLTableLegacy(PandasSQLTable):
    """Patch the PandasSQLTable for legacy support.
        Instead of a table variable just use the Create Table
        statement"""
    def sql_schema(self):
        return str(self.table)

    def create(self):
        self.pd_sql.execute(self.table)

    def insert_statement(self):
        # Replace spaces in DataFrame column names with _.
        safe_names = [_safe_col_name(n) for n in self.frame.dtypes.index]
        flv = self.pd_sql.flavor
        br_l = _SQL_SYMB[flv]['br_l']  # left val quote char
        br_r = _SQL_SYMB[flv]['br_r']  # right val quote char
        wld = _SQL_SYMB[flv]['wld']  # wildcard char

        if self.index is not None:
            safe_names.insert(0, self.index)

        bracketed_names = [br_l + column + br_r for column in safe_names]
        col_names = ','.join(bracketed_names)
        wildcards = ','.join([wld] * len(safe_names))
        insert_statement = 'INSERT INTO %s (%s) VALUES (%s)' % (
            self.name, col_names, wildcards)
        return insert_statement

    def insert(self):
        ins = self.insert_statement()
        cur = self.pd_sql.con.cursor()
        for r in self.frame.iterrows():
            data = [self.maybe_asscalar(v) for v in r[1].values]
            if self.index is not None:
                data.insert(0, self.maybe_asscalar(r[0]))
            print(type(data[2]))
            print(type(r[0]))
            cur.execute(ins, tuple(data))
        cur.close()

    def _create_table_statement(self):
        "Return a CREATE TABLE statement to suit the contents of a DataFrame."

        # Replace spaces in DataFrame column names with _.
        safe_columns = [_safe_col_name(n) for n in self.frame.dtypes.index]
        column_types = [self._sql_type_name(typ) for typ in self.frame.dtypes]

        if self.index is not None:
            safe_columns.insert(0, self.index)
            column_types.insert(0, self._sql_type_name(self.frame.index.dtype))
        flv = self.pd_sql.flavor

        br_l = _SQL_SYMB[flv]['br_l']  # left val quote char
        br_r = _SQL_SYMB[flv]['br_r']  # right val quote char

        col_template = br_l + '%s' + br_r + ' %s'

        columns = ',\n  '.join(col_template %
                               x for x in zip(safe_columns, column_types))
        template = """CREATE TABLE %(name)s (
                      %(columns)s
                      );"""
        create_statement = template % {'name': self.name, 'columns': columns}
        return create_statement

    def _sql_type_name(self, dtype):
        pytype = dtype.type
        pytype_name = "text"
        if issubclass(pytype, np.floating):
            pytype_name = "float"
        elif issubclass(pytype, np.integer):
            pytype_name = "int"
        elif issubclass(pytype, np.datetime64) or pytype is datetime:
            # Caution: np.datetime64 is also a subclass of np.number.
            pytype_name = "datetime"
        elif pytype is datetime.date:
            pytype_name = "date"
        elif issubclass(pytype, np.bool_):
            pytype_name = "bool"

        return _SQL_TYPES[pytype_name][self.pd_sql.flavor]


class PandasSQLLegacy(PandasSQL):

    def __init__(self, con, flavor):
        self.con = con
        if flavor not in ['sqlite', 'mysql']:
            raise NotImplementedError
        else:
            self.flavor = flavor

    def execute(self, *args, **kwargs):
        try:
            cur = self.con.cursor()
            if kwargs:
                cur.execute(*args, **kwargs)
            else:
                cur.execute(*args)
            return cur
        except Exception as e:
            try:
                self.con.rollback()
            except Exception:  # pragma: no cover
                ex = DatabaseError(
                    "Execution failed on sql: %s\n%s\nunable to rollback" % (args[0], e))
                raise_with_traceback(ex)

            ex = DatabaseError("Execution failed on sql: %s" % args[0])
            raise_with_traceback(ex)

    def tquery(self, *args):
        cur = self.execute(*args)
        result = self._fetchall_as_list(cur)

        # This makes into tuples
        if result and len(result[0]) == 1:
            # python 3 compat
            result = list(lzip(*result)[0])
        elif result is None:  # pragma: no cover
            result = []
        return result

    def uquery(self, *args):
        cur = self.execute(*args)
        return cur.rowcount

    def read_sql(self, sql, index_col=None, coerce_float=True, params=None,
                 parse_dates=None):
        args = _convert_params(sql, params)
        cursor = self.execute(*args)
        columns = [col_desc[0] for col_desc in cursor.description]
        data = self._fetchall_as_list(cursor)
        cursor.close()

        data_frame = _frame_from_data_and_columns(data, columns,
                                                  index_col=index_col,
                                                  coerce_float=coerce_float)
        return self._parse_date_columns(data_frame, parse_dates=parse_dates)

    def _fetchall_as_list(self, cur):
        result = cur.fetchall()
        if not isinstance(result, list):
            result = list(result)
        return result

    def to_sql(self, frame, name, if_exists='fail', index=True):
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        frame: DataFrame
        name: name of SQL table
        flavor: {'sqlite', 'mysql', 'postgres'}, default 'sqlite'
        if_exists: {'fail', 'replace', 'append'}, default 'fail'
            fail: If table exists, do nothing.
            replace: If table exists, drop it, recreate it, and insert data.
            append: If table exists, insert data. Create if does not exist.
        """
        table = PandasSQLTableLegacy(
            name, self, frame=frame, index=index, if_exists=if_exists)
        table.insert()

    def has_table(self, name):
        flavor_map = {
            'sqlite': ("SELECT name FROM sqlite_master "
                       "WHERE type='table' AND name='%s';") % name,
            'mysql': "SHOW TABLES LIKE '%s'" % name}
        query = flavor_map.get(self.flavor)

        return len(self.tquery(query)) > 0

    def get_table(self, table_name):
        return None  # not supported in Legacy mode

    def drop_table(self, name):
        drop_sql = "DROP TABLE %s" % name
        self.execute(drop_sql)


# legacy names, with depreciation warnings and copied docs
def get_schema(frame, name, con, flavor='sqlite'):
    """
    Get the SQL db table schema for the given frame

    Parameters
    ----------
    frame: DataFrame
    name: name of SQL table
    con: an open SQL database connection object
    engine: an SQLAlchemy engine - replaces connection and flavor
    flavor: {'sqlite', 'mysql', 'postgres'}, default 'sqlite'

    """
    warnings.warn(
        "get_schema is depreciated", DeprecationWarning)
    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    return pandas_sql._create_sql_schema(frame, name)


def read_frame(*args, **kwargs):
    """DEPRECIATED - use read_sql
    """
    warnings.warn(
        "read_frame is depreciated, use read_sql", DeprecationWarning)
    return read_sql(*args, **kwargs)


def write_frame(*args, **kwargs):
    """DEPRECIATED - use to_sql
    """
    warnings.warn("write_frame is depreciated, use to_sql", DeprecationWarning)
    return to_sql(*args, **kwargs)


# Append wrapped function docstrings
read_frame.__doc__ += read_sql.__doc__
write_frame.__doc__ += to_sql.__doc__
