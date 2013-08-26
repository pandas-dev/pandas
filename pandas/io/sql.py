"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""
from __future__ import print_function
from datetime import datetime, date

from pandas.compat import range, lzip, map, zip, raise_with_traceback
import pandas.compat as compat
import numpy as np
import traceback

import sqlite3
import warnings

from pandas.core.datetools import format as date_format

from pandas.core.api import DataFrame, isnull
from pandas.io import sql_legacy
from pandas.core.base import PandasObject


class SQLAlchemyRequired(ImportError):
    pass


class LegacyMySQLConnection(Exception):
    pass


class DatabaseError(IOError):
    pass


#------------------------------------------------------------------------------
# Helper execution functions

def execute(sql, con=None, retry=True, cur=None, params=None, engine=None):
    """
    Execute the given SQL query using the provided connection object.

    Parameters
    ----------
    sql: string
        Query to be executed
    con: database connection instance
        Database connection.  Must implement PEP249 (Database API v2.0).
    retry: bool
        Not currently implemented
    cur: database cursor, optional
        Must implement PEP249 (Datbase API v2.0).  If cursor is not provided,
        one will be obtained from the database connection.
    params: list or tuple, optional
        List of parameters to pass to execute method.

    Returns
    -------
    Cursor object
    """
    pandas_sql = PandasSQL(engine=engine, cur=cur, con=con)
    return pandas_sql.execute(sql=sql, retry=retry, params=params)


def tquery(sql, con=None, retry=True, cur=None, engine=None, params=None):
    """
    Returns list of tuples corresponding to each row in given sql
    query.

    If only one column selected, then plain list is returned.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: SQLConnection or DB API 2.0-compliant connection
    retry : bool
    cur: DB API 2.0 cursor

    Provide a specific connection or a specific cursor if you are executing a
    lot of sequential statements and want to commit outside.
    """
    pandas_sql = PandasSQL(engine=engine, cur=cur, con=con)
    return pandas_sql.tquery(sql=sql, retry=retry, params=params)


def uquery(sql, con=None, cur=None, retry=True, params=None, engine=None):
    """
    Does the same thing as tquery, but instead of returning results, it
    returns the number of rows affected.  Good for update queries.
    """
    pandas_sql = PandasSQL(engine=engine, cur=cur, con=con)
    return pandas_sql.uquery(sql, retry=retry, params=params)


#------------------------------------------------------------------------------
# Read and write to DataFrames


def read_sql(sql, con=None, index_col=None, flavor=None, driver=None,
             username=None, password=None, host=None, port=None,
             database=None, coerce_float=True, params=None, engine=None):
    """
    Returns a DataFrame corresponding to the result set of the query
    string.

    Optionally provide an index_col parameter to use one of the
    columns as the index. Otherwise will be 0 to len(results) - 1.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con : Connection object, SQLAlchemy Engine object, a filepath string
        (sqlite only) or the string ':memory:' (sqlite only). Alternatively,
        specify a user, passwd, host, and db below.
    index_col: string, optional
        column name to use for the returned DataFrame object.
    flavor : string specifying the flavor of SQL to use
    driver : string specifying SQL driver (e.g., MySQLdb), optional
    username: username for database authentication
        only needed if a Connection, Engine, or filepath are not given
    password: password for database authentication
        only needed if a Connection, Engine, or filepath are not given
    host: host for database connection
        only needed if a Connection, Engine, or filepath are not given
    database: database name
        only needed if a Connection, Engine, or filepath are not given
    coerce_float : boolean, default True
        Attempt to convert values to non-string, non-numeric objects (like
        decimal.Decimal) to floating point, useful for SQL result sets
    params: list or tuple, optional
        List of parameters to pass to execute method.
    engine : SQLAlchemy engine, optional

    """
    if params is None:
        params = []

    # TODO make sure con is a con and not a cursor etc.

    pandas_sql = PandasSQL(con=con, engine=engine)
    return pandas_sql.read_sql(sql, index_col=index_col, params=params,
                               flavor=flavor, coerce_float=coerce_float)


def write_frame(frame, name, con=None, flavor='sqlite', if_exists='fail', engine=None, **kwargs):
    """
    Write records stored in a DataFrame to a SQL database.

    Parameters
    ----------
    frame: DataFrame
    name: name of SQL table
    con: an open SQL database connection object
    flavor: {'sqlite', 'mysql', 'oracle'}, default 'sqlite'
    if_exists: {'fail', 'replace', 'append'}, default 'fail'
        fail: If table exists, do nothing.
        replace: If table exists, drop it, recreate it, and insert data.
        append: If table exists, insert data. Create if does not exist.
    """
    pandas_sql = PandasSQL(con=con, engine=engine)
    pandas_sql.write_frame(frame, name, flavor=flavor, if_exists=if_exists, **kwargs)

# This is an awesome function
def _engine_read_table_name(table_name, engine, meta=None, index_col=None):
    '''Given a table name and SQLAlchemy engine, return a DataFrame.

    Parameters
    ----------
    table_name: name of SQL table in database
    engine: SQLAlchemy engine
    meta: SQLAlchemy meta, optional
    index_col: columns to set as index, optional

    '''
    table = PandasSQLWithEngine(engine=engine)._get_table(table_name)

    if table is not None:
        sql_select = table.select()
        return PandasSQLWithEngine(engine=engine).read_sql(sql_select, index_col=index_col)
    else:
        raise ValueError("Table %s not found with %s." % table_name, engine)

#------------------------------------------------------------------------------
# Helper connection functions

# TODO improve this API
def get_connection(con, dialect, driver, username, password,
                   host, port, database):
    if isinstance(con, basestring):
        try:
            import sqlalchemy
            return _alchemy_connect_sqlite(con)
        except:
            return sqlite3.connect(con)
    if isinstance(con, sqlite3.Connection):
        return con
    try:
        import MySQLdb
    except ImportError:
        # If we don't have MySQLdb, this can't be a MySQLdb connection.
        pass
    else:
        if isinstance(con, MySQLdb.connection):
            raise LegacyMySQLConnection
    # If we reach here, SQLAlchemy will be needed.
    try:
        import sqlalchemy
    except ImportError:
        raise SQLAlchemyRequired
    if isinstance(con, sqlalchemy.engine.Engine):
        return con.connect()
    if isinstance(con, sqlalchemy.engine.Connection):
        return con
    if con is None:
        url_params = (dialect, driver, username,
                      password, host, port, database)
        url = _build_engine_url(*url_params)
        engine = sqlalchemy.create_engine(url)
        return engine.connect()
    if hasattr(con, 'cursor') and callable(con.cursor):
        # This looks like some Connection object from a driver module.
        raise NotImplementedError(
            """To ensure robust support of varied SQL dialects, pandas
           only support database connections from SQLAlchemy. See
           documentation.""")
    else:
        raise ValueError(
            """con must be a string, a Connection to a sqlite Database,
           or a SQLAlchemy Connection or Engine object.""")


def _alchemy_connect_sqlite(path):
    if path == ':memory:':
        return create_engine('sqlite://').connect()
    else:
        return create_engine('sqlite:///%s' % path).connect()


# I don't think this is used...
def sequence2dict(seq):
    """Helper function for cx_Oracle.

    For each element in the sequence, creates a dictionary item equal
    to the element and keyed by the position of the item in the list.
    >>> sequence2dict(("Matt", 1))
    {'1': 'Matt', '2': 1}

    Source:
    http://www.gingerandjohn.com/archives/2004/02/26/cx_oracle-executemany-example/
    """
    d = {}
    for k, v in zip(range(1, 1 + len(seq)), seq):
        d[str(k)] = v
    return d


class PandasSQLFactory(type):
    def __call__(cls, engine=None, con=None, cur=None):
        if cls is PandasSQL:
            if engine:
                return PandasSQLWithEngine(engine=engine)
            elif con:
                if cur is None:
                    cur = con.cursor()
                return PandasSQLWithCur(cur=cur)
            elif cur:
                return PandasSQLWithCur(cur=cur)
            else:
                raise ValueError("PandasSQL must be created with an engine,"
                                 " connection or cursor.")

        if cls is PandasSQLWithEngine:
            return type.__call__(cls, engine)
        elif cls is PandasSQLWithCon:
            return type.__call__(cls, con)
        elif cls is PandasSQLWithCur:
            return type.__call__(cls, cur)
        else:
            raise NotImplementedError


class PandasSQL(PandasObject):
    __metaclass__ = PandasSQLFactory

    def read_sql(self, sql, index_col=None, coerce_float=True, params=None, flavor='sqlite'):
        # Note coerce_float ignored in engine
        if params is None:
            params = []
        data, columns = self._get_data_and_columns(
            sql=sql, params=params, flavor=flavor)
        return self._frame_from_data_and_columns(data, columns,
                                                 index_col=index_col,
                                                 coerce_float=coerce_float)

    def write_frame(self, frame, name, con=None, flavor='sqlite', if_exists='fail', engine=None, **kwargs):
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        frame: DataFrame
        name: name of SQL table
        con: an open SQL database connection object
        flavor: {'sqlite', 'mysql', 'oracle'}, default 'sqlite'
        if_exists: {'fail', 'replace', 'append'}, default 'fail'
            fail: If table exists, do nothing.
            replace: If table exists, drop it, recreate it, and insert data.
            append: If table exists, insert data. Create if does not exist.
        """
        kwargs = self._pop_dep_append_kwarg(kwargs)

        exists = self._has_table(name, flavor=flavor)
        if exists:
            if if_exists == 'fail':
                raise ValueError("Table '%s' already exists." % name)
            elif if_exists == 'replace':
                self._drop_table(name)
                # do I have to create too?
        else:
            self._create_table(
                frame, name, flavor)  # flavor is ignored with engine

        self._write(frame, name, flavor)

    def _get_data_and_columns(self, *args, **kwargs):
        raise ValueError("PandasSQL must be created with an engine,"
                         " connection or cursor.")

    def _has_table(self, name, flavor='sqlite'):
        # Note: engine overrides this, to use engine.has_table
        flavor_map = {
            'sqlite': ("SELECT name FROM sqlite_master "
                       "WHERE type='table' AND name='%s';") % name,
            'mysql': "SHOW TABLES LIKE '%s'" % name}
        query = flavor_map.get(flavor, None)
        if query is None:
            raise NotImplementedError
        return len(self.tquery(query)) > 0

    def _drop_table(self, name):
        # Previously this worried about connection tp cursor then closing...
        drop_sql = "DROP TABLE %s" % name
        self.execute(drop_sql)

    @staticmethod
    def _create_schema(frame, name, flavor, keys=None):
        "Return a CREATE TABLE statement to suit the contents of a DataFrame."
        lookup_type = lambda dtype: PandasSQL._get_sqltype(dtype.type, flavor)
        # Replace spaces in DataFrame column names with _.
        safe_columns = [s.replace(
            ' ', '_').strip() for s in frame.dtypes.index]
        column_types = lzip(safe_columns, map(lookup_type, frame.dtypes))
        if flavor == 'sqlite':
            columns = ',\n  '.join('[%s] %s' % x for x in column_types)
        elif flavor == 'mysql':
             columns = ',\n  '.join('`%s` %s' % x for x in column_types)
        elif flavor == 'postgres':
            columns = ',\n  '.join('%s %s' % x for x in column_types)
        else:
            raise ValueError("Don't have a template for that database flavor.")

        keystr = ''
        if keys is not None:
            if isinstance(keys, compat.string_types):
                keys = (keys,)
            keystr = ', PRIMARY KEY (%s)' % ','.join(keys)
        template = """CREATE TABLE %(name)s (
                      %(columns)s
                      %(keystr)s
                      );"""
        create_statement = template % {'name': name, 'columns': columns,
                                       'keystr': keystr}
        return create_statement

    @staticmethod
    def _frame_from_data_and_columns(data, columns, index_col=None, coerce_float=True):
        df = DataFrame.from_records(data, columns=columns)
        if index_col is not None:
            df.set_index(index_col, inplace=True)
        return df

    @staticmethod
    def _pop_dep_append_kwarg(kwargs):
        if 'append' in kwargs:
            from warnings import warn
            warn("append is deprecated, use if_exists='append'", FutureWarning)
            if kwargs.pop('append'):
                kwargs['if_exists'] = 'append'
            else:
                kwargs['if_exists'] = 'fail'
            return kwargs

    @staticmethod
    def _safe_fetch(cur):
        '''ensures result of fetchall is a list'''
        try:
            result = cur.fetchall()
            if not isinstance(result, list):
                result = list(result)
            return result
        except Exception as e:  # pragma: no cover
            excName = e.__class__.__name__
            if excName == 'OperationalError':
                return []

    @staticmethod
    def _get_sqltype(pytype, flavor):
        sqltype = {'mysql': 'VARCHAR (63)',
                   'sqlite': 'TEXT',
                   'postgres': 'text'}

        if issubclass(pytype, np.floating):
            sqltype['mysql'] = 'FLOAT'
            sqltype['sqlite'] = 'REAL'
            sqltype['postgres'] = 'real'

        if issubclass(pytype, np.integer):
            # TODO: Refine integer size.
            sqltype['mysql'] = 'BIGINT'
            sqltype['sqlite'] = 'INTEGER'
            sqltype['postgres'] = 'integer'

        if issubclass(pytype, np.datetime64) or pytype is datetime:
            # Caution: np.datetime64 is also a subclass of np.number.
            sqltype['mysql'] = 'DATETIME'
            sqltype['sqlite'] = 'TIMESTAMP'
            sqltype['postgres'] = 'timestamp'

        if pytype is datetime.date:
            sqltype['mysql'] = 'DATE'
            sqltype['sqlite'] = 'TIMESTAMP'
            sqltype['postgres'] = 'date'

        if issubclass(pytype, np.bool_):
            sqltype['sqlite'] = 'INTEGER'
            sqltype['postgres'] = 'boolean'

        return sqltype[flavor]


class PandasSQLWithCur(PandasSQL):
    def __init__(self, cur):
        self.cur = cur

    def execute(self, sql, retry=False, params=None):
        if params is None:
            params = []
        try:
            if params:
                self.cur.execute(sql, params)
            else:
                self.cur.execute(sql)
            return self.cur
        except Exception as e:
            try:
                con.rollback()
            except Exception:  # pragma: no cover
                ex = DatabaseError(
                    "Execution failed on sql: %s\n%s\nunable to rollback" % (sql, e))
                raise_with_traceback(ex)

            ex = DatabaseError("Execution failed on sql: %s" % sql)
            raise_with_traceback(ex)

    def tquery(self, sql, retry=True, params=None, con=None):
        if params is None:
            params = []

        # result = _cur_tquery(sql, retry=retry, params=params)

        ### previously _cur_tquery
        cur = self.execute(sql, retry=retry, params=params)
        result = self._safe_fetch(self.cur)

        if con is not None:
            try:
                cur.close()
                con.commit()
            except Exception as e:
                excName = e.__class__.__name__
                if excName == 'OperationalError':  # pragma: no cover
                    print ('Failed to commit, may need to restart interpreter')
                else:
                    raise

                traceback.print_exc()
                if retry:
                    return self.tquery(sql, retry=False, params=params, con=con)
        ###

        # This makes into tuples?
        if result and len(result[0]) == 1:
            # python 3 compat
            result = list(lzip(*result)[0])
        elif result is None:  # pragma: no cover
            result = []
        return result

    def uquery(self, sql, retry=True, params=None, con=None):
        """
        Does the same thing as tquery, but instead of returning results, it
        returns the number of rows affected.  Good for update queries.
        """
        if params is None:
            params = []

        cur = self.execute(sql, retry=retry, params=params)
        row_count = cur.rowcount

        if con:  # Not sure if this is needed....
            try:
                con.commit()
            except Exception as e:
                excName = e.__class__.__name__
                if excName != 'OperationalError':
                    raise

                traceback.print_exc()
                if retry:
                    print (
                        'Looks like your connection failed, reconnecting...')
                    return self.uquery(sql, con, retry=False)
        return row_count

    def _get_data_and_columns(self, sql, con=None, flavor=None, driver=None,
                              username=None, password=None, host=None, port=None,
                              database=None, coerce_float=True, params=None):

        # dialect = flavor
        # try:
        #     connection = get_connection(con, dialect, driver, username, password,
        #                                 host, port, database)
        # except LegacyMySQLConnection:
        #     warnings.warn("For more robust support, connect using "
        #                   "SQLAlchemy. See documentation.")
        # return sql_legacy.read_frame(sql, con, index_col, coerce_float,
        # params)

        cursor = self.execute(sql, params=params)
        data = PandasSQL._safe_fetch(cursor)
        columns = [col_desc[0] for col_desc in cursor.description]
        cursor.close()
        return data, columns

    def _create_table(self, frame, name, flavor, keys=None):
        create_sql = get_schema(frame, name, flavor, keys)
        self.execute(create_sql)

    def _write(self, frame, name, flavor='sqlite'):
        # cur = con.cursor()
        # Replace spaces in DataFrame column names with _.
        safe_names = [s.replace(' ', '_').strip() for s in frame.columns]
        flavor_picker = {'sqlite': self._cur_write_sqlite,
                         'mysql': self._cur_write_mysql}

        func = flavor_picker.get(flavor, None)
        if func is None:
            raise NotImplementedError
        func(frame, name, safe_names, self.cur)
        self.cur.close()

    @staticmethod
    def _cur_write_sqlite(frame, table, names, cur):
        bracketed_names = ['[' + column + ']' for column in names]
        col_names = ','.join(bracketed_names)
        wildcards = ','.join(['?'] * len(names))
        insert_query = 'INSERT INTO %s (%s) VALUES (%s)' % (
            table, col_names, wildcards)
        # pandas types are badly handled if there is only 1 column ( Issue
        # #3628 )
        if len(frame.columns) != 1:
            data = [tuple(x) for x in frame.values]
        else:
            data = [tuple(x) for x in frame.values.tolist()]
        cur.executemany(insert_query, data)

    @staticmethod
    def _cur_write_mysql(frame, table, names, cur):
        bracketed_names = ['`' + column + '`' for column in names]
        col_names = ','.join(bracketed_names)
        wildcards = ','.join([r'%s'] * len(names))
        insert_query = "INSERT INTO %s (%s) VALUES (%s)" % (
            table, col_names, wildcards)
        # pandas types are badly handled if there is only 1 column ( Issue
        # #3628 )
        if len(frame.columns) != 1:
            data = [tuple(x) for x in frame.values]
        else:
            data = [tuple(x) for x in frame.values.tolist()]
        cur.executemany(insert_query, data)


class PandasSQLWithCon(PandasSQL):
    def __init__(self, con):
        raise NotImplementedError


class PandasSQLWithEngine(PandasSQL):
    def __init__(self, engine):
        self.engine = engine

    def execute(self, sql, retry=True, params=None):
        try:
            return self.engine.execute(sql, params=params)
        except Exception as e:
            ex = DatabaseError("Execution failed with: %s" % e)
            raise_with_traceback(ex)

    def tquery(self, sql, retry=False, params=None):
        if params is None:
            params = []
        result = self.execute(sql, params=params)
        return result.fetchall()

    def uquery(self, sql, retry=False, params=None):
        if params is None:
            params = []
        result = self.execute(sql, params=params)
        return result.rowcount

    def _get_data_and_columns(self, sql, params, **kwargs):
        result = self.execute(sql, params=params)
        data = result.fetchall()
        columns = result.keys()
        return data, columns

    def write_frame(self, frame, name, if_exists='fail', **kwargs):
        self._pop_dep_append_kwarg(kwargs)

        exists = self.engine.has_table(name)
        if exists:
            if if_exists == 'fail':
                raise ValueError("Table '%s' already exists." % name)
            elif if_exists == 'replace':
                self._drop_table(name)

        self._create_table(frame, name)
        self._write(frame, name)

    def _write(self, frame, table_name, **kwargs):
        table = self._get_table(table_name)
        ins = table.insert()
        # TODO: do this in one pass
        # engine.execute(ins, *(t[1:] for t in frame.itertuples())) # t[1:] doesn't include index
        # engine.execute(ins, *[tuple(x) for x in frame.values])

        # TODO this should be done globally first (or work out how to pass np
        # dtypes to sql)
        def maybe_asscalar(i):
            try:
                return np.asscalar(i)
            except AttributeError:
                return i

        for t in frame.iterrows():
            self.engine.execute(ins, **dict((k, maybe_asscalar(
                v)) for k, v in t[1].iteritems()))
            # TODO more efficient, I'm *sure* this was just working with tuples

    def _has_table(self, name):
        return self.engine.has_table(name)

    def _get_table(self, table_name, meta=None):
        if self.engine.has_table(table_name):
            if not meta:
                from sqlalchemy.schema import MetaData
                meta = MetaData(self.engine)
                meta.reflect(self.engine)
                return meta.tables[table_name]
        else:
            return None

    def _drop_table(self, table_name):
        if self.engine.has_table(table_name):
            table = self._get_table(table_name)
            table.drop

    def _create_table(self, frame, table_name, keys=None, meta=None):
        from sqlalchemy import Table, Column
        if keys is None:
            keys = []
        if not meta:  # Do we need meta param?
            from sqlalchemy.schema import MetaData
            meta = MetaData(self.engine)
            meta.reflect(self.engine)

        safe_columns = [s.replace(
            ' ', '_').strip() for s in frame.dtypes.index]  # may not be safe enough...
        column_types = map(self._lookup_type, frame.dtypes)

        columns = [(col_name, col_sqltype, col_name in keys)
                   for col_name, col_sqltype in zip(safe_columns, column_types)]
        columns = map(lambda (name, typ, pk): Column(
            name, typ, primary_key=pk), columns)

        table = Table(table_name, meta, *columns)
        table.create()

    def _lookup_type(self, dtype):
        from sqlalchemy import INT, FLOAT, TEXT, BOOLEAN

        pytype = dtype.type

        if issubclass(pytype, np.floating):
            return FLOAT
        if issubclass(pytype, np.integer):
            # TODO: Refine integer size.
            return INT
        if issubclass(pytype, np.datetime64) or pytype is datetime:
            # Caution: np.datetime64 is also a subclass of np.number.
            return DATETIME
        if pytype is datetime.date:
            return DATE
        if issubclass(pytype, np.bool_):
            return BOOLEAN
        return TEXT


# legacy names
get_schema = PandasSQL._create_schema
frame_query = read_sql
read_frame = read_sql
