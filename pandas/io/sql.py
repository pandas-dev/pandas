"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""
from datetime import datetime
import numpy as np
import traceback

from numpy import NaN

from pandas.core.datetools import format, to_datetime
from pandas.core.panel import pivot
from pandas.core.api import (DataFrame, DataMatrix, Series, Index,
                             isnull)

#-------------------------------------------------------------------------------
# Query formatting

_formatters = {
    datetime : lambda dt: "'%s'" % format(dt),
    str : lambda x: "'%s'" % x,
    np.str_ : lambda x: "'%s'" % x,
    unicode : lambda x: "'%s'" % x,
    float : lambda x: "%.8f" % x,
    int : lambda x: "%s" % x,
    type(None) : lambda x: "NULL",
    np.float64 : lambda x: "%.8f" % x,
    bool : lambda x: "'%s'" % x,
}

def format_query(sql, *args):
    """

    """
    processed_args = []
    for arg in args:
        if isinstance(arg, float) and isnull(arg):
            arg = None

        formatter = _formatters[type(arg)]
        processed_args.append(formatter(arg))

    return sql % tuple(processed_args)

#-------------------------------------------------------------------------------
# SQL Connection objects

class SQLConnection(object):
    """
    SQL Connection wrapper. Encapsulates parameters necessary to recreate
    connection upon database reset or failure.

    Parameters
    ----------
    host: string
    user: string
    password: string
    database: string
    """
    def __init__(self, host, user, password, database='.', trusted=False):
        self.host = host
        self.user = user
        self.password = password
        self.database = database
        self.trusted = trusted

        self._con = None
        self.connected = False

    @property
    def driver(self):
        raise Exception('No default driver')

    def connect(self):
        """

        """
        if self._con is not None:
            try:
                self._con.close()
            except Exception:
                pass

        if self.trusted:
            self._con = self.driver.connect(host=self.host,
                                            trusted=True,
                                            database=self.database)
        else:
            self._con = self.driver.connect(user=self.user,
                                            password=self.password,
                                            host=self.host,
                                            database=self.database)

    def execute(self, sql, retry=False):
        return execute(sql, self, retry=retry)

    def _get_con(self):
        """
        """
        if self._con is None:
            self.connect()
        return self._con

    def close(self):
        """
        """
        try:
            self._con.close()
        except Exception:
            pass

    def cursor(self, retry=True):
        """
        """
        try:
            return self._get_con().cursor()
        except Exception:
            if retry:
                self.connect()
                return self.cursor(retry=False)

    def commit(self):
        """
        """
        if self._con is None:
            raise Exception('Cannot commit with no connection created yet!')

        self._con.commit()

    def rollback(self):
        """
        """
        if self._con is None:
            raise Exception('Cannot commit with no connection created yet!')

        self._con.rollback()

    def getinfo(self, *args, **kwargs):
        """
        """
        return self._get_con().getinfo(*args, **kwargs)

class PymssqlConnection(SQLConnection):

    @property
    def driver(self):
        import pymssql
        return pymssql

class PyODBCConnection(SQLConnection):

    @property
    def driver(self):
        import pyodbc
        return pyodbc

    def connect(self):
        self._con = self.driver.connect(user=self.user, password=self.password,
                                        host=self.host, database=self.database,
                                        driver='{SQL Server}')

class SQLiteConnection(SQLConnection):

    def __init__(self, path):
        self.path = path
        self._con = None

    def connect(self):
        """

        """
        if self._con is not None:
            try:
                self._con.close()
            except Exception:
                pass

        self._con = self.driver.connect(self.path)

    @property
    def driver(self):
        import sqlite3
        return sqlite3

class ConnectionFactory(object):
    """
    SQL Connection Factory

    Parameters
    ----------
    host: string
        servername or path to database
    user: string
        Username (if needed)
    password: string (if needed)

    defaultDatabase: string
        Database to connect to by default
    """
    connectionClass = SQLConnection
    def __init__(self, host, user=None, password=None, defaultDatabase=None,
                 trusted=False):
        self.host = host
        self.user = user
        self.password = password
        self.trusted = trusted
        self.defaultDatabase = defaultDatabase
        self.connections = {}

    def __getitem__(self, database):
        return self.get_con(database)

    def is_user_authorized(self):
        return True

    def get_con(self, database=None, forceNew=False):
        if not self.is_user_authorized():
            raise Exception('Not authorized to use this connection')

        if database is None:
            database = self.defaultDatabase

        database = database.lower()

        if database in self.connections:
            con = self.connections[database]
        else:
            con = self.connectionClass(self.host,
                                       user=self.user,
                                       password=self.password,
                                       trusted=self.trusted,
                                       database=database)

            self.connections[database] = con

        if forceNew:
            con.connect()
        return con

class PymssqlFactory(ConnectionFactory):
    connectionClass = PymssqlConnection

# Note: with PyODBC, executing stored procedures may at times produce
# weird results. You may have to put "SET NOCOUNT ON;" at the
# beginning of the query

class PyODBCFactory(ConnectionFactory):
    connectionClass = PyODBCConnection

#-------------------------------------------------------------------------------
# Helper execution function

def execute(sql, con, retry=True, cur=None, params=()):
    """
    Execute the given SQL query using the provided connection object.

    Parameters
    ----------
    sql: string
        Query to be executed

    Returns
    -------
    Cursor object
    """
    providedCursor = False
    try:
        if cur is None:
            cur = con.cursor()
        else:
            providedCursor = True
            retry = False

        if len(params) == 0:
            cur.execute(sql)
        else:
            cur.execute(sql, params)
        return cur
    except ImportError:
        raise
    except Exception, e:
        excName = e.__class__.__name__
        if not providedCursor and excName in ('OperationalError', 'Error'):
            # connection wrapper
            if retry and isinstance(con, SQLConnection):
                con.connect()
                return execute(sql, con, retry=False)

        try:
            con.rollback()
        except Exception, e:
            pass

        print 'Error on sql %s' % sql
        raise

def _safe_fetch(cur):
    try:
        return cur.fetchall()
    except Exception, e:
        excName = e.__class__.__name__
        if excName == 'OperationalError':
            return []

def array_query(sql, con):
    """Returns results of query as a dict of numpy-arrays.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: DB connection object
    """
    cur = execute(sql, con)
    rows = _safe_fetch(cur)

    result = [np.array(x) for x in zip(*rows)]
    con.commit()

    return dict([(c[0], result[i] if len(result) > 0 else [])
                 for i, c in enumerate(cur.description)])

# def col_query(sql, con):
#     """Returns results of query as a dict of python lists.

#     Parameters
#     ----------
#     sql: string
#         SQL query to be executed
#     con: DB connection object, optional
#     """
#     cur = execute(sql, con)
#     rows = _safe_fetch(cur)

#     result = [list(x) for x in zip(*rows)]
#     con.commit()
#     if len(result) > 0:
#         return dict([(c[0], result[i]) for i, c in enumerate(cur.description)])
#     else:
#         return dict([(c[0], []) for c in cur.description])

def tquery(sql, con=None, cur=None, retry=True):
    """
    Returns list of tuples corresponding to each row in given sql
    query.

    If only one column selected, then plain list is returned.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: SQLConnection or DB API 2.0-compliant connection
    cur: DB API 2.0 cursor

    Provide a specific connection or a specific cursor if you are executing a
    lot of sequential statements and want to commit outside.
    """
    cur = execute(sql, con, cur=cur)
    result = _safe_fetch(cur)

    try:
        con.commit()
    except Exception, e:
        excName = e.__class__.__name__
        if excName == 'OperationalError':
            print 'Failed to commit, may need to restart interpreter'
        else:
            raise

        traceback.print_exc()
        if retry:
            return tquery(sql, con=con, retry=False)

    if result and len(result[0]) == 1:
        result = list(zip(*result)[0])
    elif result is None:
        result = []

    return result

def uquery(sql, con=None, cur=None, retry=True, params=()):
    """
    Does the same thing as tquery, but instead of returning results, it
    returns the number of rows affected.  Good for update queries.
    """
    cur = execute(sql, con, cur=cur, retry=retry, params=params)

    result = cur.rowcount
    try:
        con.commit()
    except Exception, e:
        excName = e.__class__.__name__
        if excName != 'OperationalError':
            raise

        traceback.print_exc()
        if retry:
            print 'Looks like your connection failed, reconnecting...'
            return uquery(sql, con, retry=False)
    return result

def frame_query(sql, con, indexField='Time', asDataMatrix=False):
    """
    Returns a DataFrame corresponding to the result set of the query
    string.

    Optionally provide an indexField parameter to use one of the
    columns as the index. Otherwise will be 0 to len(results) - 1.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: DB connection object, optional
    indexField: string, optional
        column name to use for the returned DataFrame object.
    """
    data = array_query(sql, con)
    if indexField is not None:
        try:
            idx = Index(data.pop(indexField))
        except KeyError:
            raise KeyError('indexField %s not found! %s' % (indexField, sql))
    else:
        idx = Index(np.arange(len(data.values()[0])))

    if asDataMatrix:
        return DataMatrix(data, index=idx)
    else:
        return DataFrame(data=data, index=idx)

def pivot_query(sql, rows, columns, values, con):
    """
    Returns DataFrame with columns corresponding to unique Item
    entries in the requested SQL query.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: SQLConnection
    """
    data = frame_query(sql, con)
    data = dict([(key.lower(), values) for key, values in data.iteritems()])

    pivoted = pivot(data[rows], data[columns], data[values])
    return pivoted
