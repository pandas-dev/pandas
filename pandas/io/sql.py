"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""
from datetime import datetime

import numpy as np
import traceback

from pandas.core.datetools import format as date_format
from pandas.core.api import DataFrame, isnull

#-------------------------------------------------------------------------------
# Helper execution function

def execute(sql, con, retry=True, cur=None, params=None):
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
    try:
        if cur is None:
            cur = con.cursor()

        if params is None:
            cur.execute(sql)
        else:
            cur.execute(sql, params)
        return cur
    except Exception:
        try:
            con.rollback()
        except Exception:  # pragma: no cover
            pass

        print 'Error on sql %s' % sql
        raise

def _safe_fetch(cur):
    try:
        return cur.fetchall()
    except Exception, e: # pragma: no cover
        excName = e.__class__.__name__
        if excName == 'OperationalError':
            return []

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

    if con is not None:
        try:
            con.commit()
        except Exception, e:
            excName = e.__class__.__name__
            if excName == 'OperationalError':  # pragma: no cover
                print 'Failed to commit, may need to restart interpreter'
            else:
                raise

            traceback.print_exc()
            if retry:
                return tquery(sql, con=con, retry=False)

    if result and len(result[0]) == 1:
        result = list(zip(*result)[0])
    elif result is None:  # pragma: no cover
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

def read_frame(sql, con, index_col=None):
    """
    Returns a DataFrame corresponding to the result set of the query
    string.

    Optionally provide an index_col parameter to use one of the
    columns as the index. Otherwise will be 0 to len(results) - 1.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: DB connection object, optional
    index_col: string, optional
        column name to use for the returned DataFrame object.
    """
    cur = execute(sql, con)
    rows = _safe_fetch(cur)
    con.commit()

    columns = [col_desc[0] for col_desc in cur.description]
    result = DataFrame.from_records(rows, columns=columns)

    if index_col is not None:
        result = result.set_index(index_col)

    return result

frame_query = read_frame

def write_frame(frame, name=None, con=None, flavor='sqlite'):
    """
    Write records stored in a DataFrame to SQLite. The index will currently be
    dropped
    """
    if flavor == 'sqlite':
        schema = get_sqlite_schema(frame, name)
    else:
        raise NotImplementedError

    con.execute(schema)

    wildcards = ','.join(['?'] * len(frame.columns))
    insert_sql = 'INSERT INTO %s VALUES (%s)' % (name, wildcards)
    data = [tuple(x) for x in frame.values]
    con.executemany(insert_sql, data)

def get_sqlite_schema(frame, name):
    template = """
CREATE TABLE %(name)s (
  %(columns)s
);"""

    column_types = []

    dtypes = frame.dtypes
    for k in dtypes.index:
        dt = dtypes[k]

        if issubclass(dt.type, (np.integer, np.bool_)):
            sqltype = 'INTEGER'
        elif issubclass(dt.type, np.floating):
            sqltype = 'REAL'
        else:
            sqltype = 'TEXT'

        column_types.append((k, sqltype))

    columns = ',\n  '.join('%s %s' % x for x in column_types)

    return template % {'name' : name, 'columns' : columns}




#-------------------------------------------------------------------------------
# Query formatting

_formatters = {
    datetime : lambda dt: "'%s'" % date_format(dt),
    str : lambda x: "'%s'" % x,
    np.str_ : lambda x: "'%s'" % x,
    unicode : lambda x: "'%s'" % x,
    float : lambda x: "%.8f" % x,
    int : lambda x: "%s" % x,
    type(None) : lambda x: "NULL",
    np.float64 : lambda x: "%.10f" % x,
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


