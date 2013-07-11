"""
Collection of query wrappers / abstractions to both facilitate data
retrieval and to reduce dependency on DB-specific API.
"""
from __future__ import print_function
from datetime import datetime, date
import warnings
from pandas.compat import range, lzip, map, zip, raise_with_traceback
import pandas.compat as compat
import numpy as np


from pandas.core.api import DataFrame
from pandas.core.base import PandasObject


class SQLAlchemyRequired(ImportError):
    pass


class LegacyMySQLConnection(Exception):
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


def execute(sql, con, cur=None, params=[], engine=None, flavor='sqlite'):
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
    flavor : string {sqlite, mysql} specifying the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2 connection.
    Returns
    -------
    Results Iterable
    """
    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    args = _convert_params(sql, params)
    return pandas_sql.execute(*args)


def tquery(sql, con, cur=None, params=[], engine=None, flavor='sqlite'):
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
    flavor : string {sqlite, mysql} specifying the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2 connection.
    """
    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    args = _convert_params(sql, params)
    return pandas_sql.tquery(*args)


def uquery(sql, con, cur=None, params=[], engine=None, flavor='sqlite'):
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
    flavor : string {sqlite, mysql} specifying the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2 connection.
    """
    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    args = _convert_params(sql, params)
    return pandas_sql.uquery(*args)


#------------------------------------------------------------------------------
# Read and write to DataFrames


def read_sql(sql, con, index_col=None, flavor='sqlite', coerce_float=True, params=[]):
    """
    Returns a DataFrame corresponding to the result set of the query
    string.

    Optionally provide an index_col parameter to use one of the
    columns as the index. Otherwise will be 0 to len(results) - 1.

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
    flavor : string {sqlite, mysql} specifying the flavor of SQL to use.
        Ignored when using SQLAlchemy engine. Required when using DBAPI2 connection.

    """

    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    return pandas_sql.read_sql(sql, index_col=index_col, params=params, coerce_float=coerce_float)


def to_sql(frame, name, con, flavor='sqlite', if_exists='fail'):
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
    flavor: {'sqlite', 'mysql', 'postgres'}, default 'sqlite', ignored when using engine
    if_exists: {'fail', 'replace', 'append'}, default 'fail'
        fail: If table exists, do nothing.
        replace: If table exists, drop it, recreate it, and insert data.
        append: If table exists, insert data. Create if does not exist.
    """
    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    pandas_sql.to_sql(frame, name, if_exists=if_exists)


# This is an awesome function
def read_table(table_name, con, meta=None, index_col=None, coerce_float=True):
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

    """
    pandas_sql = PandasSQLWithEngine(con, meta=meta)
    table = pandas_sql.get_table(table_name)

    if table is not None:
        sql_select = table.select()
        return pandas_sql.read_sql(sql_select, index_col=index_col, coerce_float=coerce_float)
    else:
        raise ValueError("Table %s not found with %s." % table_name, con)


def pandasSQL_builder(con, flavor=None, meta=None):
    """
    Convenience function to return the correct PandasSQL subclass based on the
    provided parameters
    """
    try:
        import sqlalchemy

        if isinstance(con, sqlalchemy.engine.Engine):
            return PandasSQLWithEngine(con, meta=meta)
        else:
<<<<<<< HEAD
            if_exists = 'fail'

    if if_exists not in ('fail', 'replace', 'append'):
        raise ValueError("'%s' is not valid for if_exists" % if_exists)

    exists = table_exists(name, con, flavor)
    if if_exists == 'fail' and exists:
        raise ValueError("Table '%s' already exists." % name)

    # creation/replacement dependent on the table existing and if_exist criteria
    create = None
    if exists:
        if if_exists == 'fail':
            raise ValueError("Table '%s' already exists." % name)
        elif if_exists == 'replace':
            cur = con.cursor()
            cur.execute("DROP TABLE %s;" % name)
            cur.close()
            create = get_schema(frame, name, flavor)
    else:
        create = get_schema(frame, name, flavor)

    if create is not None:
        cur = con.cursor()
        cur.execute(create)
=======
            warnings.warn("Not a valid SQLAlchemy engine, attempting to use as legacy DBAPI connection")
            if flavor is None:
                raise ValueError("""PandasSQL must be created with an SQLAlchemy engine
                    or a DBAPI2 connection and SQL flavour""")
            else:
                return PandasSQLWithCon(con, flavor)

    except ImportError:
        warnings.warn("SQLAlchemy not installed, using legacy mode")
        if flavor is None:
            raise SQLAlchemyRequired
        else:
            return PandasSQLWithCon(con, flavor)


class PandasSQL(PandasObject):
    """
    Subclasses Should define read_sql and to_sql
    """
    def read_sql(self, *args, **kwargs):
        raise ValueError("PandasSQL must be created with an engine,"
                         " connection or cursor.")

    def to_sql(self, *args, **kwargs):
        raise ValueError("PandasSQL must be created with an engine,"
                         " connection or cursor.")

    def _create_sql_schema(self, frame, name, keys):
        raise ValueError("PandasSQL must be created with an engine,"
                         " connection or cursor.")

    def _frame_from_data_and_columns(self, data, columns, index_col=None, coerce_float=True):
        df = DataFrame.from_records(data, columns=columns, coerce_float=coerce_float)
        if index_col is not None:
            df.set_index(index_col, inplace=True)
        return df

    def _safe_col_names(self, col_names):
        return [s.replace(' ', '_').strip() for s in col_names]  # may not be safe enough...


class PandasSQLWithEngine(PandasSQL):
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
        """Accepts same args as execute"""
        result = self.execute(*args, **kwargs)
        return result.fetchall()

    def uquery(self, *args, **kwargs):
        """Accepts same args as execute"""
        result = self.execute(*args, **kwargs)
        return result.rowcount

    def read_sql(self, sql, index_col=None, coerce_float=True, params=[]):
        args = _convert_params(sql, params)
        result = self.execute(*args)
        data = result.fetchall()
        columns = result.keys()

        return self._frame_from_data_and_columns(data, columns,
                                                 index_col=index_col,
                                                 coerce_float=coerce_float)

    def to_sql(self, frame, name, if_exists='fail'):
        if self.engine.has_table(name):
            if if_exists == 'fail':
                raise ValueError("Table '%s' already exists." % name)
            elif if_exists == 'replace':
                #TODO: this triggers a full refresh of metadata, could probably avoid this.
                self._drop_table(name)
                self._create_table(frame, name)
            elif if_exists == 'append':
                pass  # table exists and will automatically be appended to
        else:
            self._create_table(frame, name)
        self._write(frame, name)

    def _write(self, frame, table_name):
        table = self.get_table(table_name)
        ins = table.insert()
        # TODO: do this in one pass
        # TODO this should be done globally first (or work out how to pass np
        # dtypes to sql)

        def maybe_asscalar(i):
            try:
                return np.asscalar(i)
            except AttributeError:
                return i

        for t in frame.iterrows():
            self.engine.execute(ins, **dict((k, maybe_asscalar(v))
                                            for k, v in t[1].iteritems()))
            # TODO more efficient, I'm *sure* this was just working with tuples

    def has_table(self, name):
        return self.engine.has_table(name)

    def get_table(self, table_name):
        if self.engine.has_table(table_name):
            return self.meta.tables[table_name]
        else:
            return None

    def _drop_table(self, table_name):
        if self.engine.has_table(table_name):
            self.get_table(table_name).drop()
            self.meta.clear()
            self.meta.reflect()
            #print(table.exists())

    def _create_table(self, frame, table_name, keys=None):
        table = self._create_sqlalchemy_table(frame, table_name, keys)
        table.create()

    def _create_sql_schema(self, frame, table_name, keys=None):
        table = self._create_sqlalchemy_table(frame, table_name, keys)
        return str(table.compile())

    def _create_sqlalchemy_table(self, frame, table_name, keys=None):
        from sqlalchemy import Table, Column
        if keys is None:
            keys = []

        safe_columns = self._safe_col_names(frame.dtypes.index)
        column_types = map(self._lookup_type, frame.dtypes)

        columns = [(col_name, col_sqltype, col_name in keys)
                   for col_name, col_sqltype in zip(safe_columns, column_types)]

        columns = [Column(name, typ, primary_key=pk) for name, typ, pk in columns]

        return Table(table_name, self.meta, *columns)

    def _lookup_type(self, dtype):
        from sqlalchemy.types import Integer, Float, Text, Boolean, DateTime, Date

        pytype = dtype.type

        if issubclass(pytype, np.floating):
            return Float
        if issubclass(pytype, np.integer):
            # TODO: Refine integer size.
            return Integer
        if issubclass(pytype, np.datetime64) or pytype is datetime:
            # Caution: np.datetime64 is also a subclass of np.number.
            return DateTime
        if pytype is date:
            return Date
        if issubclass(pytype, np.bool_):
            return Boolean
        return Text


# ---- SQL without SQLAlchemy ---
# Flavour specific sql strings and handler class for access to DBs without SQLAlchemy installed

# SQL type convertions for each DB
_SQL_TYPES = {
    'text': {
        'mysql': 'VARCHAR (63)',
        'sqlite': 'TEXT',
        'postgres': 'text'
    },
    'float': {
        'mysql': 'FLOAT',
        'sqlite': 'REAL',
        'postgres': 'real'
    },
    'int': {
        'mysql': 'BIGINT',
        'sqlite': 'INTEGER',
        'postgres': 'integer'
    },
    'datetime': {
        'mysql': 'DATETIME',
        'sqlite': 'TIMESTAMP',
        'postgres': 'timestamp'
    },
    'date': {
        'mysql': 'DATE',
        'sqlite': 'TIMESTAMP',
        'postgres': 'date'
    },
    'bool': {
        'mysql': 'BOOLEAN',
        'sqlite': 'INTEGER',
        'postgres': 'boolean'
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
    },
    'postgres': {
        'br_l': '',
        'br_r': '',
        'wld': '?'
    }
}


class PandasSQLWithCon(PandasSQL):
    def __init__(self, con, flavor):
        self.con = con
        if flavor not in ['sqlite', 'mysql', 'postgres']:
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
        """
        Does the same thing as tquery, but instead of returning results, it
        returns the number of rows affected.  Good for update queries.
        """
        cur = self.execute(*args)
        return cur.rowcount

    def read_sql(self, sql, index_col=None, coerce_float=True, params=[], flavor='sqlite'):
        args = _convert_params(sql, params)
        cursor = self.execute(*args)
        columns = [col_desc[0] for col_desc in cursor.description]
        data = self._fetchall_as_list(cursor)
        cursor.close()

        return self._frame_from_data_and_columns(data, columns,
                                                 index_col=index_col,
                                                 coerce_float=coerce_float)

    def to_sql(self, frame, name, con=None, if_exists='fail'):
        """
        Write records stored in a DataFrame to a SQL database.

        Parameters
        ----------
        frame: DataFrame
        name: name of SQL table
        con: an open SQL database connection object
        flavor: {'sqlite', 'mysql', 'postgres'}, default 'sqlite'
        if_exists: {'fail', 'replace', 'append'}, default 'fail'
            fail: If table exists, do nothing.
            replace: If table exists, drop it, recreate it, and insert data.
            append: If table exists, insert data. Create if does not exist.
        """
        if self.has_table(name):
            if if_exists == 'fail':
                raise ValueError("Table '%s' already exists." % name)
            elif if_exists == 'replace':
                self._drop_table(name)
                self._create_table(frame, name)
            elif if_exists == "append":
                pass  # should just add...
        else:
            self._create_table(frame, name)

        self._write(frame, name)

    def _fetchall_as_list(self, cur):
        '''ensures result of fetchall is a list'''
        result = cur.fetchall()
        if not isinstance(result, list):
            result = list(result)
        return result

    def _write(self, frame, table_name):
        # Replace spaces in DataFrame column names with _.
        safe_names = self._safe_col_names(frame.columns)

        br_l = _SQL_SYMB[self.flavor]['br_l']  # left val quote char
        br_r = _SQL_SYMB[self.flavor]['br_r']  # right val quote char
        wld = _SQL_SYMB[self.flavor]['wld']  # wildcard char

        bracketed_names = [br_l + column + br_r for column in safe_names]
        col_names = ','.join(bracketed_names)
        wildcards = ','.join([wld] * len(safe_names))
        insert_query = 'INSERT INTO %s (%s) VALUES (%s)' % (
            table_name, col_names, wildcards)

        # pandas types are badly handled if there is only 1 col (Issue #3628)
        if len(frame.columns) != 1:
            data = [tuple(x) for x in frame.values]
        else:
            data = [tuple(x) for x in frame.values.tolist()]

        cur = self.con.cursor()
        cur.executemany(insert_query, data)
>>>>>>> 1259dca... ENH #4163 Use SQLAlchemy for DB abstraction
        cur.close()

    def _create_table(self, frame, name, keys=None):
        create_sql = self._create_sql_schema(frame, name, keys)
        self.execute(create_sql)

    def has_table(self, name):
        flavor_map = {
            'sqlite': ("SELECT name FROM sqlite_master "
                       "WHERE type='table' AND name='%s';") % name,
            'mysql': "SHOW TABLES LIKE '%s'" % name}
        query = flavor_map.get(self.flavor)
        if query is None:
            raise NotImplementedError
        return len(self.tquery(query)) > 0

    def _drop_table(self, name):
        # Previously this worried about connection tp cursor then closing...
        drop_sql = "DROP TABLE %s" % name
        self.execute(drop_sql)

    def _create_sql_schema(self, frame, table_name, keys=None):
        "Return a CREATE TABLE statement to suit the contents of a DataFrame."

        lookup_type = lambda dtype: self._get_sqltype(dtype.type)
        # Replace spaces in DataFrame column names with _.
        safe_columns = self._safe_col_names(frame.dtypes.index)

        column_types = lzip(safe_columns, map(lookup_type, frame.dtypes))

        br_l = _SQL_SYMB[self.flavor]['br_l']  # left val quote char
        br_r = _SQL_SYMB[self.flavor]['br_r']  # right val quote char
        col_template = br_l + '%s' + br_r + ' %s'
        columns = ',\n  '.join(col_template % x for x in column_types)

        keystr = ''
        if keys is not None:
            if isinstance(keys, compat.string_types):
                keys = (keys,)
            keystr = ', PRIMARY KEY (%s)' % ','.join(keys)
        template = """CREATE TABLE %(name)s (
                      %(columns)s
                      %(keystr)s
                      );"""
        create_statement = template % {'name': table_name, 'columns': columns,
                                       'keystr': keystr}
        return create_statement

    def _get_sqltype(self, pytype):
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

        return _SQL_TYPES[pytype_name][self.flavor]


# legacy names
def get_schema(frame, name, con=None, flavor='sqlite', engine=None):
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
    pandas_sql = pandasSQL_builder(con=con, flavor=flavor)
    return pandas_sql._create_sql_schema()



#TODO: add depreciation warnings
read_frame = read_sql
write_frame = to_sql

