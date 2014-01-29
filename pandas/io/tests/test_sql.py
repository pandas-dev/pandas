from __future__ import print_function
import sqlite3
import sys

import warnings

import nose

import numpy as np

from pandas.core.datetools import format as date_format
from pandas.core.api import DataFrame, isnull
from pandas.compat import StringIO, range, lrange
import pandas.compat as compat

import pandas.io.sql as sql
import pandas.util.testing as tm
from pandas import Series, Index, DataFrame
from datetime import datetime

_formatters = {
    datetime: lambda dt: "'%s'" % date_format(dt),
    str: lambda x: "'%s'" % x,
    np.str_: lambda x: "'%s'" % x,
    compat.text_type: lambda x: "'%s'" % x,
    compat.binary_type: lambda x: "'%s'" % x,
    float: lambda x: "%.8f" % x,
    int: lambda x: "%s" % x,
    type(None): lambda x: "NULL",
    np.float64: lambda x: "%.10f" % x,
    bool: lambda x: "'%s'" % x,
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

def _skip_if_no_MySQLdb():
    try:
        import MySQLdb
    except ImportError:
        raise nose.SkipTest('MySQLdb not installed, skipping')

class TestSQLite(tm.TestCase):

    def setUp(self):
        self.db = sqlite3.connect(':memory:')

    def test_basic(self):
        frame = tm.makeTimeDataFrame()
        self._check_roundtrip(frame)

    def test_write_row_by_row(self):
        frame = tm.makeTimeDataFrame()
        frame.ix[0, 0] = np.nan
        create_sql = sql.get_schema(frame, 'test', 'sqlite')
        cur = self.db.cursor()
        cur.execute(create_sql)

        cur = self.db.cursor()

        ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"
        for idx, row in frame.iterrows():
            fmt_sql = format_query(ins, *row)
            sql.tquery(fmt_sql, cur=cur)

        self.db.commit()

        result = sql.read_frame("select * from test", con=self.db)
        result.index = frame.index
        tm.assert_frame_equal(result, frame)

    def test_execute(self):
        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_schema(frame, 'test', 'sqlite')
        cur = self.db.cursor()
        cur.execute(create_sql)
        ins = "INSERT INTO test VALUES (?, ?, ?, ?)"

        row = frame.ix[0]
        sql.execute(ins, self.db, params=tuple(row))
        self.db.commit()

        result = sql.read_frame("select * from test", self.db)
        result.index = frame.index[:1]
        tm.assert_frame_equal(result, frame[:1])

    def test_schema(self):
        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_schema(frame, 'test', 'sqlite')
        lines = create_sql.splitlines()
        for l in lines:
            tokens = l.split(' ')
            if len(tokens) == 2 and tokens[0] == 'A':
                self.assert_(tokens[1] == 'DATETIME')

        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_schema(frame, 'test', 'sqlite', keys=['A', 'B'],)
        lines = create_sql.splitlines()
        self.assert_('PRIMARY KEY (A,B)' in create_sql)
        cur = self.db.cursor()
        cur.execute(create_sql)

    def test_execute_fail(self):
        create_sql = """
        CREATE TABLE test
        (
        a TEXT,
        b TEXT,
        c REAL,
        PRIMARY KEY (a, b)
        );
        """
        cur = self.db.cursor()
        cur.execute(create_sql)

        sql.execute('INSERT INTO test VALUES("foo", "bar", 1.234)', self.db)
        sql.execute('INSERT INTO test VALUES("foo", "baz", 2.567)', self.db)

        try:
            sys.stdout = StringIO()
            self.assertRaises(Exception, sql.execute,
                              'INSERT INTO test VALUES("foo", "bar", 7)',
                              self.db)
        finally:
            sys.stdout = sys.__stdout__

    def test_execute_closed_connection(self):
        create_sql = """
        CREATE TABLE test
        (
        a TEXT,
        b TEXT,
        c REAL,
        PRIMARY KEY (a, b)
        );
        """
        cur = self.db.cursor()
        cur.execute(create_sql)

        sql.execute('INSERT INTO test VALUES("foo", "bar", 1.234)', self.db)
        self.db.close()
        try:
            sys.stdout = StringIO()
            self.assertRaises(Exception, sql.tquery, "select * from test",
                              con=self.db)
        finally:
            sys.stdout = sys.__stdout__

    def test_na_roundtrip(self):
        pass

    def _check_roundtrip(self, frame):
        sql.write_frame(frame, name='test_table', con=self.db)
        result = sql.read_frame("select * from test_table", self.db)

        # HACK! Change this once indexes are handled properly.
        result.index = frame.index

        expected = frame
        tm.assert_frame_equal(result, expected)

        frame['txt'] = ['a'] * len(frame)
        frame2 = frame.copy()
        frame2['Idx'] = Index(lrange(len(frame2))) + 10
        sql.write_frame(frame2, name='test_table2', con=self.db)
        result = sql.read_frame("select * from test_table2", self.db,
                                index_col='Idx')
        expected = frame.copy()
        expected.index = Index(lrange(len(frame2))) + 10
        expected.index.name = 'Idx'
        print(expected.index.names)
        print(result.index.names)
        tm.assert_frame_equal(expected, result)

    def test_tquery(self):
        frame = tm.makeTimeDataFrame()
        sql.write_frame(frame, name='test_table', con=self.db)
        result = sql.tquery("select A from test_table", self.db)
        expected = frame.A
        result = Series(result, frame.index)
        tm.assert_series_equal(result, expected)

        try:
            sys.stdout = StringIO()
            self.assertRaises(sqlite3.OperationalError, sql.tquery,
                              'select * from blah', con=self.db)

            self.assertRaises(sqlite3.OperationalError, sql.tquery,
                              'select * from blah', con=self.db, retry=True)
        finally:
            sys.stdout = sys.__stdout__

    def test_uquery(self):
        frame = tm.makeTimeDataFrame()
        sql.write_frame(frame, name='test_table', con=self.db)
        stmt = 'INSERT INTO test_table VALUES(2.314, -123.1, 1.234, 2.3)'
        self.assertEqual(sql.uquery(stmt, con=self.db), 1)

        try:
            sys.stdout = StringIO()

            self.assertRaises(sqlite3.OperationalError, sql.tquery,
                              'insert into blah values (1)', con=self.db)

            self.assertRaises(sqlite3.OperationalError, sql.tquery,
                              'insert into blah values (1)', con=self.db,
                              retry=True)
        finally:
            sys.stdout = sys.__stdout__

    def test_keyword_as_column_names(self):
        '''
        '''
        df = DataFrame({'From':np.ones(5)})
        sql.write_frame(df, con = self.db, name = 'testkeywords')

    def test_onecolumn_of_integer(self):
        # GH 3628
        # a column_of_integers dataframe should transfer well to sql

        mono_df=DataFrame([1 , 2], columns=['c0'])
        sql.write_frame(mono_df, con = self.db, name = 'mono_df')
        # computing the sum via sql
        con_x=self.db
        the_sum=sum([my_c0[0] for  my_c0 in con_x.execute("select * from mono_df")])
        # it should not fail, and gives 3 ( Issue #3628 )
        self.assertEqual(the_sum , 3)

        result = sql.read_frame("select * from mono_df",con_x)
        tm.assert_frame_equal(result,mono_df)

    def test_if_exists(self):
        df_if_exists_1 = DataFrame({'col1': [1, 2], 'col2': ['A', 'B']})
        df_if_exists_2 = DataFrame({'col1': [3, 4, 5], 'col2': ['C', 'D', 'E']})
        table_name = 'table_if_exists'
        sql_select = "SELECT * FROM %s" % table_name

        def clean_up(test_table_to_drop):
            """
            Drops tables created from individual tests
            so no dependencies arise from sequential tests
            """
            if sql.table_exists(test_table_to_drop, self.db, flavor='sqlite'):
                cur = self.db.cursor()
                cur.execute("DROP TABLE %s" % test_table_to_drop)
                cur.close()

        # test if invalid value for if_exists raises appropriate error
        self.assertRaises(ValueError,
                          sql.write_frame,
                          frame=df_if_exists_1,
                          con=self.db,
                          name=table_name,
                          flavor='sqlite',
                          if_exists='notvalidvalue')
        clean_up(table_name)

        # test if_exists='fail'
        sql.write_frame(frame=df_if_exists_1, con=self.db, name=table_name,
                        flavor='sqlite', if_exists='fail')
        self.assertRaises(ValueError,
                          sql.write_frame,
                          frame=df_if_exists_1,
                          con=self.db,
                          name=table_name,
                          flavor='sqlite',
                          if_exists='fail')

        # test if_exists='replace'
        sql.write_frame(frame=df_if_exists_1, con=self.db, name=table_name,
                        flavor='sqlite', if_exists='replace')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(1, 'A'), (2, 'B')])
        sql.write_frame(frame=df_if_exists_2, con=self.db, name=table_name,
                        flavor='sqlite', if_exists='replace')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(3, 'C'), (4, 'D'), (5, 'E')])
        clean_up(table_name)
                        
        # test if_exists='append'
        sql.write_frame(frame=df_if_exists_1, con=self.db, name=table_name,
                        flavor='sqlite', if_exists='fail')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(1, 'A'), (2, 'B')])
        sql.write_frame(frame=df_if_exists_2, con=self.db, name=table_name,
                        flavor='sqlite', if_exists='append')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(1, 'A'), (2, 'B'), (3, 'C'), (4, 'D'), (5, 'E')])
        clean_up(table_name)


class TestMySQL(tm.TestCase):

    def setUp(self):
        _skip_if_no_MySQLdb()
        import MySQLdb
        try:
            # Try Travis defaults.
            # No real user should allow root access with a blank password.
            self.db = MySQLdb.connect(host='localhost', user='root', passwd='',
                                    db='pandas_nosetest')
        except:
            pass
        else:
            return
        try:
            self.db = MySQLdb.connect(read_default_group='pandas')
        except MySQLdb.ProgrammingError as e:
            raise nose.SkipTest(
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")
        except MySQLdb.Error as e:
            raise nose.SkipTest(
                "Cannot connect to database. "
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")

    def test_basic(self):
        _skip_if_no_MySQLdb()
        frame = tm.makeTimeDataFrame()
        self._check_roundtrip(frame)

    def test_write_row_by_row(self):
        _skip_if_no_MySQLdb()
        frame = tm.makeTimeDataFrame()
        frame.ix[0, 0] = np.nan
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = sql.get_schema(frame, 'test', 'mysql')
        cur = self.db.cursor()
        cur.execute(drop_sql)
        cur.execute(create_sql)
        ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"
        for idx, row in frame.iterrows():
            fmt_sql = format_query(ins, *row)
            sql.tquery(fmt_sql, cur=cur)

        self.db.commit()

        result = sql.read_frame("select * from test", con=self.db)
        result.index = frame.index
        tm.assert_frame_equal(result, frame)

    def test_execute(self):
        _skip_if_no_MySQLdb()
        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = sql.get_schema(frame, 'test', 'mysql')
        cur = self.db.cursor()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Unknown table.*")
            cur.execute(drop_sql)
        cur.execute(create_sql)
        ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"

        row = frame.ix[0]
        sql.execute(ins, self.db, params=tuple(row))
        self.db.commit()

        result = sql.read_frame("select * from test", self.db)
        result.index = frame.index[:1]
        tm.assert_frame_equal(result, frame[:1])

    def test_schema(self):
        _skip_if_no_MySQLdb()
        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_schema(frame, 'test', 'mysql')
        lines = create_sql.splitlines()
        for l in lines:
            tokens = l.split(' ')
            if len(tokens) == 2 and tokens[0] == 'A':
                self.assert_(tokens[1] == 'DATETIME')

        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = sql.get_schema(frame, 'test', 'mysql', keys=['A', 'B'],)
        lines = create_sql.splitlines()
        self.assert_('PRIMARY KEY (A,B)' in create_sql)
        cur = self.db.cursor()
        cur.execute(drop_sql)
        cur.execute(create_sql)

    def test_execute_fail(self):
        _skip_if_no_MySQLdb()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = """
        CREATE TABLE test
        (
        a TEXT,
        b TEXT,
        c REAL,
        PRIMARY KEY (a(5), b(5))
        );
        """
        cur = self.db.cursor()
        cur.execute(drop_sql)
        cur.execute(create_sql)

        sql.execute('INSERT INTO test VALUES("foo", "bar", 1.234)', self.db)
        sql.execute('INSERT INTO test VALUES("foo", "baz", 2.567)', self.db)

        try:
            sys.stdout = StringIO()
            self.assertRaises(Exception, sql.execute,
                              'INSERT INTO test VALUES("foo", "bar", 7)',
                              self.db)
        finally:
            sys.stdout = sys.__stdout__

    def test_execute_closed_connection(self):
        _skip_if_no_MySQLdb()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = """
        CREATE TABLE test
        (
        a TEXT,
        b TEXT,
        c REAL,
        PRIMARY KEY (a(5), b(5))
        );
        """
        cur = self.db.cursor()
        cur.execute(drop_sql)
        cur.execute(create_sql)

        sql.execute('INSERT INTO test VALUES("foo", "bar", 1.234)', self.db)
        self.db.close()
        try:
            sys.stdout = StringIO()
            self.assertRaises(Exception, sql.tquery, "select * from test",
                              con=self.db)
        finally:
            sys.stdout = sys.__stdout__

    def test_na_roundtrip(self):
        _skip_if_no_MySQLdb()
        pass

    def _check_roundtrip(self, frame):
        _skip_if_no_MySQLdb()
        drop_sql = "DROP TABLE IF EXISTS test_table"
        cur = self.db.cursor()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Unknown table.*")
            cur.execute(drop_sql)
        sql.write_frame(frame, name='test_table', con=self.db, flavor='mysql')
        result = sql.read_frame("select * from test_table", self.db)

        # HACK! Change this once indexes are handled properly.
        result.index = frame.index
        result.index.name = frame.index.name

        expected = frame
        tm.assert_frame_equal(result, expected)

        frame['txt'] = ['a'] * len(frame)
        frame2 = frame.copy()
        index = Index(lrange(len(frame2))) + 10
        frame2['Idx'] = index
        drop_sql = "DROP TABLE IF EXISTS test_table2"
        cur = self.db.cursor()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Unknown table.*")
            cur.execute(drop_sql)
        sql.write_frame(frame2, name='test_table2', con=self.db, flavor='mysql')
        result = sql.read_frame("select * from test_table2", self.db,
                                index_col='Idx')
        expected = frame.copy()

        # HACK! Change this once indexes are handled properly.
        expected.index = index
        expected.index.names = result.index.names
        tm.assert_frame_equal(expected, result)

    def test_tquery(self):
        try:
            import MySQLdb
        except ImportError:
            raise nose.SkipTest("no MySQLdb")
        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test_table"
        cur = self.db.cursor()
        cur.execute(drop_sql)
        sql.write_frame(frame, name='test_table', con=self.db, flavor='mysql')
        result = sql.tquery("select A from test_table", self.db)
        expected = frame.A
        result = Series(result, frame.index)
        tm.assert_series_equal(result, expected)

        try:
            sys.stdout = StringIO()
            self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
                              'select * from blah', con=self.db)

            self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
                              'select * from blah', con=self.db, retry=True)
        finally:
            sys.stdout = sys.__stdout__

    def test_uquery(self):
        try:
            import MySQLdb
        except ImportError:
            raise nose.SkipTest("no MySQLdb")
        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test_table"
        cur = self.db.cursor()
        cur.execute(drop_sql)
        sql.write_frame(frame, name='test_table', con=self.db, flavor='mysql')
        stmt = 'INSERT INTO test_table VALUES(2.314, -123.1, 1.234, 2.3)'
        self.assertEqual(sql.uquery(stmt, con=self.db), 1)

        try:
            sys.stdout = StringIO()

            self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
                              'insert into blah values (1)', con=self.db)

            self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
                              'insert into blah values (1)', con=self.db,
                              retry=True)
        finally:
            sys.stdout = sys.__stdout__

    def test_keyword_as_column_names(self):
        '''
        '''
        _skip_if_no_MySQLdb()
        df = DataFrame({'From':np.ones(5)})
        sql.write_frame(df, con = self.db, name = 'testkeywords',
                        if_exists='replace', flavor='mysql')

    def test_if_exists(self):
        _skip_if_no_MySQLdb()
        df_if_exists_1 = DataFrame({'col1': [1, 2], 'col2': ['A', 'B']})
        df_if_exists_2 = DataFrame({'col1': [3, 4, 5], 'col2': ['C', 'D', 'E']})
        table_name = 'table_if_exists'
        sql_select = "SELECT * FROM %s" % table_name

        def clean_up(test_table_to_drop):
            """
            Drops tables created from individual tests
            so no dependencies arise from sequential tests
            """
            if sql.table_exists(test_table_to_drop, self.db, flavor='mysql'):
                cur = self.db.cursor()
                cur.execute("DROP TABLE %s" % test_table_to_drop)
                cur.close()

        # test if invalid value for if_exists raises appropriate error
        self.assertRaises(ValueError,
                          sql.write_frame,
                          frame=df_if_exists_1,
                          con=self.db,
                          name=table_name,
                          flavor='mysql',
                          if_exists='notvalidvalue')
        clean_up(table_name)

        # test if_exists='fail'
        sql.write_frame(frame=df_if_exists_1, con=self.db, name=table_name,
                        flavor='mysql', if_exists='fail')
        self.assertRaises(ValueError,
                          sql.write_frame,
                          frame=df_if_exists_1,
                          con=self.db,
                          name=table_name,
                          flavor='mysql',
                          if_exists='fail')

        # test if_exists='replace'
        sql.write_frame(frame=df_if_exists_1, con=self.db, name=table_name,
                        flavor='mysql', if_exists='replace')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(1, 'A'), (2, 'B')])
        sql.write_frame(frame=df_if_exists_2, con=self.db, name=table_name,
                        flavor='mysql', if_exists='replace')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(3, 'C'), (4, 'D'), (5, 'E')])
        clean_up(table_name)
                        
        # test if_exists='append'
        sql.write_frame(frame=df_if_exists_1, con=self.db, name=table_name,
                        flavor='mysql', if_exists='fail')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(1, 'A'), (2, 'B')])
        sql.write_frame(frame=df_if_exists_2, con=self.db, name=table_name,
                        flavor='mysql', if_exists='append')
        self.assertEqual(sql.tquery(sql_select, con=self.db),
                         [(1, 'A'), (2, 'B'), (3, 'C'), (4, 'D'), (5, 'E')])
        clean_up(table_name)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
