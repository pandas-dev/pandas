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



class SQLAlchemyRequired(Exception):
    pass

class LegacyMySQLConnection(Exception):
    pass

class DatabaseError(Exception):
    pass


#------------------------------------------------------------------------------
# Helper execution function

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
    if engine is not None:
        try:
            return engine.execute(sql, params=params)
        except Exception as e:
            ex = DatabaseError("Execution failed with: %s" % e)
            raise_with_traceback(ex)

    try:
        if cur is None:
            cur = con.cursor()

        if params is None:
            cur.execute(sql)
        else:
            cur.execute(sql, params)
        return cur
    except Exception as e:
        try:
            con.rollback()
        except Exception:  # pragma: no cover
            ex = DatabaseError("Execution failed on sql: %s\n%s\nunable to rollback" % (sql, e))
            raise_with_traceback(ex)

        ex = DatabaseError("Execution failed on sql: %s\nunable to rollback" % sql)
        raise_with_traceback(ex)

def _safe_fetch(cur=None):
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
    if params is None:
        params = []
    if engine:
        result = execute(sql, *params, engine=engine)
        return result.fetchall()  # is this tuples?
    else:
        result = _cur_tquery(sql, con=con, retry=retry, cur=cur, params=params)

    # This makes into tuples?
    if result and len(result[0]) == 1:
        # python 3 compat
        result = list(lzip(*result)[0])
    elif result is None:  # pragma: no cover
        result = []
    return result


def _cur_tquery(sql, con=None, retry=True, cur=None, engine=None, params=None):

    cur = execute(sql, con, cur=cur, params=params)
    result = _safe_fetch(cur)

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
                return tquery(sql, con=con, retry=False)

    return result


def uquery(sql, con=None, cur=None, retry=True, params=None, engine=None):
    """
    Does the same thing as tquery, but instead of returning results, it
    returns the number of rows affected.  Good for update queries.
    """
    if params is None:
        params = []

    if engine:
        result = execute(sql, *params, engine=engine)
        return result.rowcount

    else:
        return _cur_uquery(sql, con=con, cur=cur, retry=retry, params=params)


def _cur_uquery(sql, con=None, cur=None, retry=True, params=None, engine=None):
    cur = execute(sql, con, cur=cur, retry=retry, params=params)
    row_count = cur.rowcount
    try:
        con.commit()
    except Exception as e:
        excName = e.__class__.__name__
        if excName != 'OperationalError':
            raise

        traceback.print_exc()
        if retry:
            print ('Looks like your connection failed, reconnecting...')
            return uquery(sql, con, retry=False)
    return row_count


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
        url_params = (dialect, driver, username, \
                      password, host, port, database)
        url = _build_url(*url_params)
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

def _build_url(dialect, driver, username, password, host, port, database):
    # Create an Engine and from that a Connection.
    # We use a string instead of sqlalchemy.engine.url.URL because
    # we do not necessarily know the driver; we know the dialect.
    required_params = [dialect, username, password, host, database]
    for p in required_params:
        if not isinstance(p, basestring):
            raise ValueError("Insufficient information to connect to a database; see docstring.")
    url = dialect
    if driver is not None:
        url += "+%s" % driver
    url += "://%s:%s@%s" % (username, password, host)
    if port is not None:
        url += ":%d" % port
    url += "/%s" % database
    return url

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

    if engine:
        result = engine.execute(sql, *params)
        data = result.fetchall()
        columns = result.keys()

    else:
        dialect = flavor
        try:
            connection = get_connection(con, dialect, driver, username, password, 
                                        host, port, database)
        except LegacyMySQLConnection:
            warnings.warn("For more robust support, connect using " \
                          "SQLAlchemy. See documentation.")
            return sql_legacy.read_frame(sql, con, index_col, coerce_float, params)

        cursor = connection.execute(sql, *params)
        data = _safe_fetch(cursor)
        columns = [col_desc[0] for col_desc in cursor.description]
        cursor.close()

    result = DataFrame.from_records(data, columns=columns)

    if index_col is not None:
        result = result.set_index(index_col)

    return result

frame_query = read_sql
read_frame = read_sql


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

    if 'append' in kwargs:
        import warnings
        warnings.warn("append is deprecated, use if_exists instead",
                      FutureWarning)
        if kwargs['append']:
            if_exists = 'append'
        else:
            if_exists='fail'

    if engine:
        exists = engine.has_table(name)
    else:
        exists = table_exists(name, con, flavor)

    create = None  #create or drop-recreate if necessary
    if exists:
        if if_exists == 'fail':
            raise ValueError("Table '%s' already exists." % name)
        elif if_exists == 'replace':
            if engine:
                _engine_drop_table(name)
            else:
                create = "DROP TABLE %s" % name
    else:
        if engine:
            _engine_create_table(frame, name, engine=engine)
        else:
            create = get_schema(frame, name, flavor)

    if create is not None:
        cur = con.cursor()
        cur.execute(create)
        cur.close()

    if engine:
        _engine_write(frame, name, engine)
    else:
        cur = con.cursor()
        # Replace spaces in DataFrame column names with _.
        safe_names = [s.replace(' ', '_').strip() for s in frame.columns]
        flavor_picker = {'sqlite' : _cur_write_sqlite,
                         'mysql' : _cur_write_mysql}

        func = flavor_picker.get(flavor, None)
        if func is None:
            raise NotImplementedError
        func(frame, name, safe_names, cur)
        cur.close()
        con.commit()


def _cur_write_sqlite(frame, table, names, cur):
    bracketed_names = ['[' + column + ']' for column in names]
    col_names = ','.join(bracketed_names)
    wildcards = ','.join(['?'] * len(names))
    insert_query = 'INSERT INTO %s (%s) VALUES (%s)' % (
        table, col_names, wildcards)
    # pandas types are badly handled if there is only 1 column ( Issue #3628 )
    if not len(frame.columns) == 1:
        data = [tuple(x) for x in frame.values]
    else:
        data = [tuple(x) for x in frame.values.tolist()]
    cur.executemany(insert_query, data)

def _cur_write_mysql(frame, table, names, cur):
    bracketed_names = ['`' + column + '`' for column in names]
    col_names = ','.join(bracketed_names)
    wildcards = ','.join([r'%s'] * len(names))
    insert_query = "INSERT INTO %s (%s) VALUES (%s)" % (
        table, col_names, wildcards)
    data = [tuple(x) for x in frame.values]
    cur.executemany(insert_query, data)

def _engine_write(frame, table_name, engine):
    table = _engine_get_table(table_name, engine)
    ins = table.insert()
    # TODO: do this in one pass
    # engine.execute(ins, *(t[1:] for t in frame.itertuples())) # t[1:] doesn't include index
    # engine.execute(ins, *[tuple(x) for x in frame.values])

    # TODO this should be done globally first (or work out how to pass np dtypes to sql)
    def maybe_asscalar(i):
        try:
            return np.asscalar(i)
        except AttributeError:
            return i

    for t in frame.iterrows():
        engine.execute(ins, **dict((k, maybe_asscalar(v)) for k, v in t[1].iteritems()))
        # TODO more efficient, I'm *sure* this was just working with tuples


def table_exists(name, con=None, flavor=None, engine=None):
    if engine:
        return engine.has_table(name)

    else:
        flavor_map = {
            'sqlite': ("SELECT name FROM sqlite_master "
                       "WHERE type='table' AND name='%s';") % name,
            'mysql' : "SHOW TABLES LIKE '%s'" % name}
        query = flavor_map.get(flavor, None)
        if query is None:
            raise NotImplementedError
        return len(tquery(query, con)) > 0


def get_sqltype(pytype, flavor):
    sqltype = {'mysql': 'VARCHAR (63)',
               'sqlite': 'TEXT'}

    if issubclass(pytype, np.floating):
        sqltype['mysql'] = 'FLOAT'
        sqltype['sqlite'] = 'REAL'

    if issubclass(pytype, np.integer):
        #TODO: Refine integer size.
        sqltype['mysql'] = 'BIGINT'
        sqltype['sqlite'] = 'INTEGER'

    if issubclass(pytype, np.datetime64) or pytype is datetime:
        # Caution: np.datetime64 is also a subclass of np.number.
        sqltype['mysql'] = 'DATETIME'
        sqltype['sqlite'] = 'TIMESTAMP'

    if pytype is datetime.date:
        sqltype['mysql'] = 'DATE'
        sqltype['sqlite'] = 'TIMESTAMP'

    if issubclass(pytype, np.bool_):
        sqltype['sqlite'] = 'INTEGER'

    return sqltype[flavor]


def get_schema(frame, name, flavor, keys=None):
    "Return a CREATE TABLE statement to suit the contents of a DataFrame."
    lookup_type = lambda dtype: get_sqltype(dtype.type, flavor)
    # Replace spaces in DataFrame column names with _.
    safe_columns = [s.replace(' ', '_').strip() for s in frame.dtypes.index]
    column_types = lzip(safe_columns, map(lookup_type, frame.dtypes))
    if flavor == 'sqlite':
        columns = ',\n  '.join('[%s] %s' % x for x in column_types)
    else:
        columns = ',\n  '.join('`%s` %s' % x for x in column_types)

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


def _engine_drop_table(table_name, engine):
    if engine.has_table(table_name):
        table = _engine_get_table(table_name, engine=engine)
        table.drop()

def _engine_lookup_type(dtype):
    from sqlalchemy import Table, Column, INT, FLOAT, TEXT, BOOLEAN

    pytype = dtype.type

    if issubclass(pytype, np.floating):
        return FLOAT

    if issubclass(pytype, np.integer):
        #TODO: Refine integer size.
        return INT

    if issubclass(pytype, np.datetime64) or pytype is datetime:
        # Caution: np.datetime64 is also a subclass of np.number.
        return DATETIME

    if pytype is datetime.date:
        return DATE

    if issubclass(pytype, np.bool_):
        return BOOLEAN

    return TEXT

def _engine_create_table(frame, table_name, engine, keys=None, meta=None):
    from sqlalchemy import Table, Column
    if keys is None:
        keys = []
    if not meta:
        from sqlalchemy.schema import MetaData
        meta = MetaData(engine)
        meta.reflect(engine)

    safe_columns = [s.replace(' ', '_').strip() for s in frame.dtypes.index]  # may not be safe enough...
    column_types = map(_engine_lookup_type, frame.dtypes)

    columns = [(col_name, col_sqltype, col_name in keys)
                    for col_name, col_sqltype in zip(safe_columns, column_types)]
    columns = map(lambda (name, typ, pk): Column(name, typ, primary_key=pk), columns)

    table = Table(table_name, meta, *columns)

    table.create()

def _engine_get_table(table_name, engine, meta=None):
    if engine.has_table(table_name):
        if not meta:
            from sqlalchemy.schema import MetaData
            meta = MetaData(engine)
            meta.reflect(engine)
            return meta.tables[table_name]
    else:
        return None

def _engine_read_sql(sql, engine, params=None, index_col=None):

    if params is None:
        params = []

    try:
        result = engine.execute(sql, *params)
    except Exception as e:
        raise DatabaseError
    data = result.fetchall()
    columns = result.keys()

    df = DataFrame.from_records(data, columns=columns)
    if index_col is not None:
        df.set_index(index_col, inplace=True)
    return df

def _engine_read_table_name(table_name, engine, meta=None, index_col=None):
    table = _engine_get_table(table_name, engine=engine, meta=meta)

    if table is not None:
        sql_select = table.select()
        return _engine_read_sql(sql_select, engine=engine, index_col=index_col)
    else:
        raise ValueError("Table %s not found with %s." % table_name, engine)

def _engine_write_frame(frame, name, engine, if_exists='fail'):

    exists = engine.has_table(name)
    if exists:
        if if_exists == 'fail':
            raise ValueError("Table '%s' already exists." % name)
        elif if_exists == 'replace':
            _engine_drop_table(name)
    else:
        _engine_create_table(frame, name, engine=engine)

    _engine_write(frame, name, engine)
