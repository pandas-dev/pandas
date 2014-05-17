"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""
from __future__ import print_function, division
from datetime import datetime, date, timedelta

import warnings
import traceback
import itertools
import re
import numpy as np

import pandas.core.common as com
from pandas.compat import lzip, map, zip, raise_with_traceback, string_types
from pandas.core.api import DataFrame, Series
from pandas.core.base import PandasObject
from pandas.tseries.tools import to_datetime


class SQLAlchemyRequired(ImportError):
    pass


class DatabaseError(IOError):
    pass


#------------------------------------------------------------------------------
# Helper functions

def _convert_params(sql, params):
    """convert sql and params args to DBAPI2.0 compliant format"""
    args = [sql]
    if params is not None:
        if hasattr(params, 'keys'):  # test if params is a mapping
            args += [params]
        else:
            args += [list(params)]
    return args


def _handle_date_column(col, format=None):
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


def _parse_date_columns(data_frame, parse_dates):
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
        data_frame[col_name] = _handle_date_column(df_col, format=fmt)

    return data_frame


def execute(sql, con, cur=None, params=None):
    """
    Execute the given SQL query using the provided connection object.

    Parameters
    ----------
    sql : string
        Query to be executed
    con : SQLAlchemy engine or sqlite3 DBAPI2 connection
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object, only sqlite3 is supported.
    cur : depreciated, cursor is obtained from connection
    params : list or tuple, optional
        List of parameters to pass to execute method.

    Returns
    -------
    Results Iterable
    """
    if cur is None:
        pandas_sql = pandasSQL_builder(con)
    else:
        pandas_sql = pandasSQL_builder(cur, is_cursor=True)
    args = _convert_params(sql, params)
    return pandas_sql.execute(*args)


#------------------------------------------------------------------------------
#--- Deprecated tquery and uquery

def _safe_fetch(cur):
    try:
        result = cur.fetchall()
        if not isinstance(result, list):
            result = list(result)
        return result
    except Exception as e: # pragma: no cover
        excName = e.__class__.__name__
        if excName == 'OperationalError':
            return []

def tquery(sql, con=None, cur=None, retry=True):
    """
    DEPRECATED. Returns list of tuples corresponding to each row in given sql
    query.

    If only one column selected, then plain list is returned.

    To obtain the same result in the future, you can use the following:

    >>> execute(sql, con, params).fetchall()

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: DBAPI2 connection
    cur: depreciated, cursor is obtained from connection

    Returns
    -------
    Results Iterable

    """
    warnings.warn(
        "tquery is depreciated, and will be removed in future versions. "
        "You can use ``execute(...).fetchall()`` instead.",
        FutureWarning)

    cur = execute(sql, con, cur=cur)
    result = _safe_fetch(cur)

    if con is not None:
        try:
            cur.close()
            con.commit()
        except Exception as e:
            excName = e.__class__.__name__
            if excName == 'OperationalError': # pragma: no cover
                print('Failed to commit, may need to restart interpreter')
            else:
                raise

            traceback.print_exc()
            if retry:
                return tquery(sql, con=con, retry=False)

    if result and len(result[0]) == 1:
        # python 3 compat
        result = list(lzip(*result)[0])
    elif result is None: # pragma: no cover
        result = []

    return result


def uquery(sql, con=None, cur=None, retry=True, params=None):
    """
    DEPRECATED. Does the same thing as tquery, but instead of returning results, it
    returns the number of rows affected.  Good for update queries.

    To obtain the same result in the future, you can use the following:

    >>> execute(sql, con).rowcount

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: DBAPI2 connection
    cur: depreciated, cursor is obtained from connection
    params: list or tuple, optional
        List of parameters to pass to execute method.

    Returns
    -------
    Number of affected rows

    """
    warnings.warn(
        "uquery is depreciated, and will be removed in future versions. "
        "You can use ``execute(...).rowcount`` instead.",
        FutureWarning)

    cur = execute(sql, con, cur=cur, params=params)

    result = cur.rowcount
    try:
        con.commit()
    except Exception as e:
        excName = e.__class__.__name__
        if excName != 'OperationalError':
            raise

        traceback.print_exc()
        if retry:
            print('Looks like your connection failed, reconnecting...')
            return uquery(sql, con, retry=False)
    return result


#------------------------------------------------------------------------------
#--- Read and write to DataFrames

def read_sql_table(table_name, con, index_col=None, coerce_float=True,
                   parse_dates=None, columns=None):
    """Read SQL database table into a DataFrame.

    Given a table name and an SQLAlchemy engine, returns a DataFrame.
    This function does not support DBAPI connections.

    Parameters
    ----------
    table_name : string
        Name of SQL table in database
    con : SQLAlchemy engine
        Sqlite DBAPI conncection mode not supported
    index_col : string, optional
        Column to set as index
    coerce_float : boolean, default True
        Attempt to convert values to non-string, non-numeric objects (like
        decimal.Decimal) to floating point. Can result in loss of Precision.
    parse_dates : list or dict
        - List of column names to parse as dates
        - Dict of ``{column_name: format string}`` where format string is
          strftime compatible in case of parsing string times or is one of
          (D, s, ns, ms, us) in case of parsing integer timestamps
        - Dict of ``{column_name: arg dict}``, where the arg dict corresponds
          to the keyword arguments of :func:`pandas.to_datetime`
          Especially useful with databases without native Datetime support,
          such as SQLite
    columns : list
        List of column names to select from sql table

    Returns
    -------
    DataFrame

    See also
    --------
    read_sql_query : Read SQL query into a DataFrame.
    read_sql


    """
    pandas_sql = PandasSQLAlchemy(con)
    table = pandas_sql.read_table(
        table_name, index_col=index_col, coerce_float=coerce_float,
        parse_dates=parse_dates, columns=columns)

    if table is not None:
        return table
    else:
        raise ValueError("Table %s not found" % table_name, con)


def read_sql_query(sql, con, index_col=None, coerce_float=True, params=None,
                   parse_dates=None):
    """Read SQL query into a DataFrame.

    Returns a DataFrame corresponding to the result set of the query
    string. Optionally provide an `index_col` parameter to use one of the
    columns as the index, otherwise default integer index will be used.

    Parameters
    ----------
    sql : string
        SQL query to be executed
    con : SQLAlchemy engine or sqlite3 DBAPI2 connection
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object, only sqlite3 is supported.
    index_col : string, optional
        Column name to use as index for the returned DataFrame object.
    coerce_float : boolean, default True
        Attempt to convert values to non-string, non-numeric objects (like
        decimal.Decimal) to floating point, useful for SQL result sets
    params : list, tuple or dict, optional
        List of parameters to pass to execute method.
    parse_dates : list or dict
        - List of column names to parse as dates
        - Dict of ``{column_name: format string}`` where format string is
          strftime compatible in case of parsing string times or is one of
          (D, s, ns, ms, us) in case of parsing integer timestamps
        - Dict of ``{column_name: arg dict}``, where the arg dict corresponds
          to the keyword arguments of :func:`pandas.to_datetime`
          Especially useful with databases without native Datetime support,
          such as SQLite

    Returns
    -------
    DataFrame

    See also
    --------
    read_sql_table : Read SQL database table into a DataFrame
    read_sql

    """
    pandas_sql = pandasSQL_builder(con)
    return pandas_sql.read_sql(
        sql, index_col=index_col, params=params, coerce_float=coerce_float,
        parse_dates=parse_dates)


def read_sql(sql, con, index_col=None, coerce_float=True, params=None,
             parse_dates=None, columns=None):
    """
    Read SQL query or database table into a DataFrame.

    Parameters
    ----------
    sql : string
        SQL query to be executed or database table name.
    con : SQLAlchemy engine or DBAPI2 connection (legacy mode)
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object, only sqlite3 is supported.
    index_col : string, optional
        column name to use as index for the returned DataFrame object.
    coerce_float : boolean, default True
        Attempt to convert values to non-string, non-numeric objects (like
        decimal.Decimal) to floating point, useful for SQL result sets
    params : list, tuple or dict, optional
        List of parameters to pass to execute method.
    parse_dates : list or dict
        - List of column names to parse as dates
        - Dict of ``{column_name: format string}`` where format string is
          strftime compatible in case of parsing string times or is one of
          (D, s, ns, ms, us) in case of parsing integer timestamps
        - Dict of ``{column_name: arg dict}``, where the arg dict corresponds
          to the keyword arguments of :func:`pandas.to_datetime`
          Especially useful with databases without native Datetime support,
          such as SQLite
    columns : list
        List of column names to select from sql table (only used when reading
        a table).

    Returns
    -------
    DataFrame

    Notes
    -----
    This function is a convenience wrapper around ``read_sql_table`` and
    ``read_sql_query`` (and for backward compatibility) and will delegate
    to the specific function depending on the provided input (database
    table name or sql query).

    See also
    --------
    read_sql_table : Read SQL database table into a DataFrame
    read_sql_query : Read SQL query into a DataFrame

    """
    pandas_sql = pandasSQL_builder(con)

    if 'select' in sql.lower():
        try:
            if pandas_sql.has_table(sql):
                return pandas_sql.read_table(
                    sql, index_col=index_col, coerce_float=coerce_float,
                    parse_dates=parse_dates, columns=columns)
        except:
            pass

        return pandas_sql.read_sql(
            sql, index_col=index_col, params=params,
            coerce_float=coerce_float, parse_dates=parse_dates)
    else:
        if isinstance(pandas_sql, PandasSQLLegacy):
            raise ValueError("Reading a table with read_sql is not supported "
                             "for a DBAPI2 connection. Use an SQLAlchemy "
                             "engine or specify an sql query")
        return pandas_sql.read_table(
            sql, index_col=index_col, coerce_float=coerce_float,
            parse_dates=parse_dates, columns=columns)


def to_sql(frame, name, con, flavor='sqlite', if_exists='fail', index=True,
           index_label=None):
    """
    Write records stored in a DataFrame to a SQL database.

    Parameters
    ----------
    frame : DataFrame
    name : string
        Name of SQL table
    con : SQLAlchemy engine or sqlite3 DBAPI2 connection
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object, only sqlite3 is supported.
    flavor : {'sqlite', 'mysql'}, default 'sqlite'
        The flavor of SQL to use. Ignored when using SQLAlchemy engine.
        'mysql' is deprecated and will be removed in future versions, but it
        will be further supported through SQLAlchemy engines.
    if_exists : {'fail', 'replace', 'append'}, default 'fail'
        - fail: If table exists, do nothing.
        - replace: If table exists, drop it, recreate it, and insert data.
        - append: If table exists, insert data. Create if does not exist.
    index : boolean, default True
        Write DataFrame index as a column
    index_label : string or sequence, default None
        Column label for index column(s). If None is given (default) and
        `index` is True, then the index names are used.
        A sequence should be given if the DataFrame uses MultiIndex.

    """
    if if_exists not in ('fail', 'replace', 'append'):
        raise ValueError("'{0}' is not valid for if_exists".format(if_exists))

    pandas_sql = pandasSQL_builder(con, flavor=flavor)

    if isinstance(frame, Series):
        frame = frame.to_frame()
    elif not isinstance(frame, DataFrame):
        raise NotImplementedError

    pandas_sql.to_sql(frame, name, if_exists=if_exists, index=index,
                      index_label=index_label)


def has_table(table_name, con, flavor='sqlite'):
    """
    Check if DataBase has named table.

    Parameters
    ----------
    table_name: string
        Name of SQL table
    con: SQLAlchemy engine or sqlite3 DBAPI2 connection
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object, only sqlite3 is supported.
    flavor: {'sqlite', 'mysql'}, default 'sqlite'
        The flavor of SQL to use. Ignored when using SQLAlchemy engine.
        'mysql' is deprecated and will be removed in future versions, but it
        will be further supported through SQLAlchemy engines.

    Returns
    -------
    boolean
    """
    pandas_sql = pandasSQL_builder(con, flavor=flavor)
    return pandas_sql.has_table(table_name)

table_exists = has_table


_MYSQL_WARNING = ("The 'mysql' flavor with DBAPI connection is deprecated "
                  "and will be removed in future versions. "
                  "MySQL will be further supported with SQLAlchemy engines.")

def pandasSQL_builder(con, flavor=None, meta=None, is_cursor=False):
    """
    Convenience function to return the correct PandasSQL subclass based on the
    provided parameters
    """
    # When support for DBAPI connections is removed,
    # is_cursor should not be necessary.
    try:
        import sqlalchemy

        if isinstance(con, sqlalchemy.engine.Engine):
            return PandasSQLAlchemy(con, meta=meta)
        else:
            if flavor == 'mysql':
                warnings.warn(_MYSQL_WARNING, FutureWarning)
            return PandasSQLLegacy(con, flavor, is_cursor=is_cursor)

    except ImportError:
        if flavor == 'mysql':
            warnings.warn(_MYSQL_WARNING, FutureWarning)
        return PandasSQLLegacy(con, flavor, is_cursor=is_cursor)


class PandasSQLTable(PandasObject):
    """
    For mapping Pandas tables to SQL tables.
    Uses fact that table is reflected by SQLAlchemy to
    do better type convertions.
    Also holds various flags needed to avoid having to
    pass them between functions all the time.
    """
    # TODO: support for multiIndex
    def __init__(self, name, pandas_sql_engine, frame=None, index=True,
                 if_exists='fail', prefix='pandas', index_label=None):
        self.name = name
        self.pd_sql = pandas_sql_engine
        self.prefix = prefix
        self.frame = frame
        self.index = self._index_name(index, index_label)

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
                    raise ValueError(
                        "'{0}' is not valid for if_exists".format(if_exists))
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
        from sqlalchemy.schema import CreateTable
        return str(CreateTable(self.table))

    def create(self):
        self.table.create()

    def insert_statement(self):
        return self.table.insert()

    def maybe_asscalar(self, i):
        try:
            return np.asscalar(i)
        except AttributeError:
            return i

    def insert_data(self):
        if self.index is not None:
            temp = self.frame.copy()
            temp.index.names = self.index
            try:
                temp.reset_index(inplace=True)
            except ValueError as err:
                raise ValueError(
                    "duplicate name in index/columns: {0}".format(err))
        else:
            temp = self.frame

        return temp

    def insert(self):
        ins = self.insert_statement()
        data_list = []
        temp = self.insert_data()
        keys = temp.columns

        for t in temp.itertuples():
            data = dict((k, self.maybe_asscalar(v))
                        for k, v in zip(keys, t[1:]))
            data_list.append(data)

        self.pd_sql.execute(ins, data_list)

    def read(self, coerce_float=True, parse_dates=None, columns=None):

        if columns is not None and len(columns) > 0:
            from sqlalchemy import select
            cols = [self.table.c[n] for n in columns]
            if self.index is not None:
                [cols.insert(0, self.table.c[idx]) for idx in self.index[::-1]]
            sql_select = select(cols)
        else:
            sql_select = self.table.select()

        result = self.pd_sql.execute(sql_select)
        data = result.fetchall()
        column_names = result.keys()

        self.frame = DataFrame.from_records(
            data, columns=column_names, coerce_float=coerce_float)

        self._harmonize_columns(parse_dates=parse_dates)

        if self.index is not None:
            self.frame.set_index(self.index, inplace=True)

        return self.frame

    def _index_name(self, index, index_label):
        # for writing: index=True to include index in sql table
        if index is True:
            nlevels = self.frame.index.nlevels
            # if index_label is specified, set this as index name(s)
            if index_label is not None:
                if not isinstance(index_label, list):
                    index_label = [index_label]
                if len(index_label) != nlevels:
                    raise ValueError(
                        "Length of 'index_label' should match number of "
                        "levels, which is {0}".format(nlevels))
                else:
                    return index_label
            # return the used column labels for the index columns
            if nlevels == 1 and 'index' not in self.frame.columns and self.frame.index.name is None:
                return ['index']
            else:
                return [l if l is not None else "level_{0}".format(i)
                        for i, l in enumerate(self.frame.index.names)]

        # for reading: index=(list of) string to specify column to set as index
        elif isinstance(index, string_types):
            return [index]
        elif isinstance(index, list):
            return index
        else:
            return None

    def _create_table_statement(self):
        from sqlalchemy import Table, Column

        columns = list(map(str, self.frame.columns))
        column_types = map(self._sqlalchemy_type, self.frame.dtypes)

        columns = [Column(name, typ)
                   for name, typ in zip(columns, column_types)]

        if self.index is not None:
            for i, idx_label in enumerate(self.index[::-1]):
                idx_type = self._sqlalchemy_type(
                    self.frame.index.get_level_values(i))
                columns.insert(0, Column(idx_label, idx_type, index=True))

        return Table(self.name, self.pd_sql.meta, *columns)

    def _harmonize_columns(self, parse_dates=None):
        """ Make a data_frame's column type align with an sql_table
            column types
            Need to work around limited NA value support.
            Floats are always fine, ints must always
            be floats if there are Null values.
            Booleans are hard because converting bool column with None replaces
            all Nones with false. Therefore only convert bool if there are no
            NA values.
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
                        self.frame[col_name] = _handle_date_column(df_col)

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
                    self.frame[col_name] = _handle_date_column(
                        df_col, format=fmt)

            except KeyError:
                pass  # this column not in results

    def _sqlalchemy_type(self, arr_or_dtype):
        from sqlalchemy.types import Integer, Float, Text, Boolean, DateTime, Date, Interval

        if arr_or_dtype is date:
            return Date
        if com.is_datetime64_dtype(arr_or_dtype):
            try:
                tz = arr_or_dtype.tzinfo
                return DateTime(timezone=True)
            except:
                return DateTime
        if com.is_timedelta64_dtype(arr_or_dtype):
            warnings.warn("the 'timedelta' type is not supported, and will be "
                          "written as integer values (ns frequency) to the "
                          "database.", UserWarning)
            return Integer
        elif com.is_float_dtype(arr_or_dtype):
            return Float
        elif com.is_integer_dtype(arr_or_dtype):
            # TODO: Refine integer size.
            return Integer
        elif com.is_bool(arr_or_dtype):
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

    def read_table(self, table_name, index_col=None, coerce_float=True,
                   parse_dates=None, columns=None):

        table = PandasSQLTable(table_name, self, index=index_col)
        return table.read(coerce_float=coerce_float,
                          parse_dates=parse_dates, columns=columns)

    def read_sql(self, sql, index_col=None, coerce_float=True,
                 parse_dates=None, params=None):
        args = _convert_params(sql, params)

        result = self.execute(*args)
        data = result.fetchall()
        columns = result.keys()

        data_frame = DataFrame.from_records(
            data, columns=columns, coerce_float=coerce_float)

        _parse_date_columns(data_frame, parse_dates)

        if index_col is not None:
            data_frame.set_index(index_col, inplace=True)

        return data_frame

    def to_sql(self, frame, name, if_exists='fail', index=True,
               index_label=None):
        table = PandasSQLTable(
            name, self, frame=frame, index=index, if_exists=if_exists,
            index_label=index_label)
        table.insert()

    @property
    def tables(self):
        return self.meta.tables

    def has_table(self, name):
        if self.meta.tables.get(name) is not None:
            return True
        else:
            return False

    def get_table(self, table_name):
        return self.meta.tables.get(table_name)

    def drop_table(self, table_name):
        if self.engine.has_table(table_name):
            self.get_table(table_name).drop()
            self.meta.clear()
            self.meta.reflect()

    def _create_sql_schema(self, frame, table_name):
        table = PandasSQLTable(table_name, self, frame=frame)
        return str(table.sql_schema())


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


_SAFE_NAMES_WARNING = ("The spaces in these column names will not be changed. "
                       "In pandas versions < 0.14, spaces were converted to "
                       "underscores.")


class PandasSQLTableLegacy(PandasSQLTable):
    """Patch the PandasSQLTable for legacy support.
        Instead of a table variable just use the Create Table
        statement"""
    def sql_schema(self):
        return str(self.table)

    def create(self):
        self.pd_sql.execute(self.table)

    def insert_statement(self):
        names = list(map(str, self.frame.columns))
        flv = self.pd_sql.flavor
        br_l = _SQL_SYMB[flv]['br_l']  # left val quote char
        br_r = _SQL_SYMB[flv]['br_r']  # right val quote char
        wld = _SQL_SYMB[flv]['wld']  # wildcard char

        if self.index is not None:
            [names.insert(0, idx) for idx in self.index[::-1]]

        bracketed_names = [br_l + column + br_r for column in names]
        col_names = ','.join(bracketed_names)
        wildcards = ','.join([wld] * len(names))
        insert_statement = 'INSERT INTO %s (%s) VALUES (%s)' % (
            self.name, col_names, wildcards)
        return insert_statement

    def insert(self):
        ins = self.insert_statement()
        temp = self.insert_data()
        data_list = []

        for t in temp.itertuples():
            data = tuple((self.maybe_asscalar(v) for v in t[1:]))
            data_list.append(data)

        cur = self.pd_sql.con.cursor()
        cur.executemany(ins, data_list)
        cur.close()
        self.pd_sql.con.commit()

    def _create_table_statement(self):
        "Return a CREATE TABLE statement to suit the contents of a DataFrame."

        columns = list(map(str, self.frame.columns))
        pat = re.compile('\s+')
        if any(map(pat.search, columns)):
            warnings.warn(_SAFE_NAMES_WARNING)
        column_types = [self._sql_type_name(typ) for typ in self.frame.dtypes]

        if self.index is not None:
            for i, idx_label in enumerate(self.index[::-1]):
                columns.insert(0, idx_label)
                column_types.insert(0, self._sql_type_name(self.frame.index.get_level_values(i).dtype))

        flv = self.pd_sql.flavor

        br_l = _SQL_SYMB[flv]['br_l']  # left val quote char
        br_r = _SQL_SYMB[flv]['br_r']  # right val quote char

        col_template = br_l + '%s' + br_r + ' %s'

        columns = ',\n  '.join(col_template %
                               x for x in zip(columns, column_types))
        template = """CREATE TABLE %(name)s (
                      %(columns)s
                      )"""
        create_statement = template % {'name': self.name, 'columns': columns}
        return create_statement

    def _sql_type_name(self, dtype):
        pytype = dtype.type
        pytype_name = "text"
        if issubclass(pytype, np.floating):
            pytype_name = "float"
        elif com.is_timedelta64_dtype(pytype):
            warnings.warn("the 'timedelta' type is not supported, and will be "
                          "written as integer values (ns frequency) to the "
                          "database.", UserWarning)
            pytype_name = "int"
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

    def __init__(self, con, flavor, is_cursor=False):
        self.is_cursor = is_cursor
        self.con = con
        if flavor is None:
            flavor = 'sqlite'
        if flavor not in ['sqlite', 'mysql']:
            raise NotImplementedError
        else:
            self.flavor = flavor

    def execute(self, *args, **kwargs):
        if self.is_cursor:
            cur = self.con
        else:
            cur = self.con.cursor()
        try:
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

    def read_sql(self, sql, index_col=None, coerce_float=True, params=None,
                 parse_dates=None):
        args = _convert_params(sql, params)
        cursor = self.execute(*args)
        columns = [col_desc[0] for col_desc in cursor.description]
        data = self._fetchall_as_list(cursor)
        cursor.close()

        data_frame = DataFrame.from_records(
            data, columns=columns, coerce_float=coerce_float)

        _parse_date_columns(data_frame, parse_dates)

        if index_col is not None:
            data_frame.set_index(index_col, inplace=True)
        return data_frame

    def _fetchall_as_list(self, cur):
        result = cur.fetchall()
        if not isinstance(result, list):
            result = list(result)
        return result

    def to_sql(self, frame, name, if_exists='fail', index=True,
               index_label=None):
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        frame: DataFrame
        name: name of SQL table
        flavor: {'sqlite', 'mysql'}, default 'sqlite'
        if_exists: {'fail', 'replace', 'append'}, default 'fail'
            fail: If table exists, do nothing.
            replace: If table exists, drop it, recreate it, and insert data.
            append: If table exists, insert data. Create if does not exist.

        """
        table = PandasSQLTableLegacy(
            name, self, frame=frame, index=index, if_exists=if_exists,
            index_label=index_label)
        table.insert()

    def has_table(self, name):
        flavor_map = {
            'sqlite': ("SELECT name FROM sqlite_master "
                       "WHERE type='table' AND name='%s';") % name,
            'mysql': "SHOW TABLES LIKE '%s'" % name}
        query = flavor_map.get(self.flavor)

        return len(self.execute(query).fetchall()) > 0

    def get_table(self, table_name):
        return None  # not supported in Legacy mode

    def drop_table(self, name):
        drop_sql = "DROP TABLE %s" % name
        self.execute(drop_sql)

    def _create_sql_schema(self, frame, table_name):
        table = PandasSQLTableLegacy(table_name, self, frame=frame)
        return str(table.sql_schema())


def get_schema(frame, name, flavor='sqlite', keys=None, con=None):
    """
    Get the SQL db table schema for the given frame.

    Parameters
    ----------
    frame : DataFrame
    name : string
        name of SQL table
    flavor : {'sqlite', 'mysql'}, default 'sqlite'
        The flavor of SQL to use. Ignored when using SQLAlchemy engine.
        'mysql' is deprecated and will be removed in future versions, but it
        will be further supported through SQLAlchemy engines.
    keys : string or sequence
        columns to use a primary key
    con: an open SQL database connection object or an SQLAlchemy engine
        Using SQLAlchemy makes it possible to use any DB supported by that
        library.
        If a DBAPI2 object, only sqlite3 is supported.

    """

    if con is None:
        if flavor == 'mysql':
            warnings.warn(_MYSQL_WARNING, FutureWarning)
        return _get_schema_legacy(frame, name, flavor, keys)

    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    return pandas_sql._create_sql_schema(frame, name)


def _get_schema_legacy(frame, name, flavor, keys=None):
    """Old function from 0.13.1. To keep backwards compatibility.
    When mysql legacy support is dropped, it should be possible to
    remove this code
    """

    def get_sqltype(dtype, flavor):
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

        return _SQL_TYPES[pytype_name][flavor]

    lookup_type = lambda dtype: get_sqltype(dtype, flavor)

    column_types = lzip(frame.dtypes.index, map(lookup_type, frame.dtypes))
    if flavor == 'sqlite':
        columns = ',\n  '.join('[%s] %s' % x for x in column_types)
    else:
        columns = ',\n  '.join('`%s` %s' % x for x in column_types)

    keystr = ''
    if keys is not None:
        if isinstance(keys, string_types):
            keys = (keys,)
        keystr = ', PRIMARY KEY (%s)' % ','.join(keys)
    template = """CREATE TABLE %(name)s (
                  %(columns)s
                  %(keystr)s
                  );"""
    create_statement = template % {'name': name, 'columns': columns,
                                   'keystr': keystr}
    return create_statement


# legacy names, with depreciation warnings and copied docs

def read_frame(*args, **kwargs):
    """DEPRECIATED - use read_sql
    """
    warnings.warn("read_frame is depreciated, use read_sql", FutureWarning)
    return read_sql(*args, **kwargs)


def frame_query(*args, **kwargs):
    """DEPRECIATED - use read_sql
    """
    warnings.warn("frame_query is depreciated, use read_sql", FutureWarning)
    return read_sql(*args, **kwargs)


def write_frame(frame, name, con, flavor='sqlite', if_exists='fail', **kwargs):
    """DEPRECIATED - use to_sql

    Write records stored in a DataFrame to a SQL database.

    Parameters
    ----------
    frame : DataFrame
    name : string
    con : DBAPI2 connection
    flavor : {'sqlite', 'mysql'}, default 'sqlite'
        The flavor of SQL to use.
    if_exists : {'fail', 'replace', 'append'}, default 'fail'
        - fail: If table exists, do nothing.
        - replace: If table exists, drop it, recreate it, and insert data.
        - append: If table exists, insert data. Create if does not exist.
    index : boolean, default False
        Write DataFrame index as a column

    Notes
    -----
    This function is deprecated in favor of ``to_sql``. There are however
    two differences:

    - With ``to_sql`` the index is written to the sql database by default. To
      keep the behaviour this function you need to specify ``index=False``.
    - The new ``to_sql`` function supports sqlalchemy engines to work with
      different sql flavors.

    See also
    --------
    pandas.DataFrame.to_sql

    """
    warnings.warn("write_frame is depreciated, use to_sql", FutureWarning)

    # for backwards compatibility, set index=False when not specified
    index = kwargs.pop('index', False)
    return to_sql(frame, name, con, flavor=flavor, if_exists=if_exists,
                  index=index, **kwargs)


# Append wrapped function docstrings
read_frame.__doc__ += read_sql.__doc__
frame_query.__doc__ += read_sql.__doc__
