from __future__ import with_statement
from pandas.util.py3compat import StringIO
import unittest
import sqlite3
import sys

import warnings

import nose

import numpy as np

from pandas.core.datetools import format as date_format

import pandas.io.sql as sql
import pandas.util.testing as tm
from pandas import Series, Index, DataFrame, isnull
from datetime import datetime

import sqlalchemy
import sqlite3  # try to import other db modules in their test classes


class TestSQLAlchemy(unittest.TestCase):

    def set_flavor_engine(self):
        # override for other db modules
        self.engine = sqlalchemy.create_engine('sqlite:///:memory:')

    def setUp(self):
        # this is overriden for other db modules        
        self.set_flavor_engine()

        # shared for all db modules
        self.meta = sqlalchemy.schema.MetaData(self.engine)
        self.drop_table('test')  # should already be done ?
        self.meta.reflect(self.engine) # not sure if this is different

        self.frame = tm.makeTimeDataFrame()

    def drop_table(self, table_name):
        drop_sql = "DROP TABLE IF EXISTS %s" % table_name
        self.engine.execute(drop_sql)

    def create_table(self, frame, table_name):
        create_sql = sql.get_schema(frame, table_name, self.engine.dialect.name)  # TODO let get_schema use engine
        self.engine.execute(create_sql)

    def get_table(self, table_name):
        self.meta.reflect(self.engine)
        return self.meta.tables[table_name]

    def tquery(self, fmt_sql, params=None):
        sql.tquery(fmt_sql, engine=self.engine, params=params)  # make tquery work with engine?

    def read_frame(self, fmt_sql):
        return sql.read_frame(fmt_sql, engine=self.engine)  # have read_frame use engine, pass select expression?



    def test_basic(self):
        self._check_roundtrip(self.frame)

    def test_write_row_by_row(self):
        self.frame.ix[0, 0] = np.nan
        self.create_table(self.frame, 'test')

        test_table = self.get_table('test')

        ins = test_table.insert() # INSERT INTO test VALUES (%s, %s, %s, %s)
        for idx, row in self.frame.iterrows():
            fmt_sql = tuple(row)
            self.tquery(ins, fmt_sql)

        select_test = str(test_table.select())  # SELECT * FROM test

        result = self.read_frame(select_test)

        result.index = self.frame.index
        tm.assert_frame_equal(result, self.frame)

    def test_execute(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        self.create_sql = sql.get_schema(self.frame, 'test', self.engine.dialect.name)
        self.engine.execute(create_sql)
        ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"

        row = frame.ix[0]
        sql.execute(ins, self.db, params=tuple(row))
        self.db.commit()

        result = sql.read_frame("select * from test", self.db)
        result.index = frame.index[:1]
        tm.assert_frame_equal(result, frame[:1])

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
        cur.execute(drop_sql)
        sql.write_frame(frame, name='test_table', con=self.db, flavor='mysql')
        result = sql.read_frame("select * from test_table", self.db)

        # HACK!
        result.index = frame.index

        expected = frame
        tm.assert_frame_equal(result, expected)

        frame['txt'] = ['a'] * len(frame)
        frame2 = frame.copy()
        frame2['Idx'] = Index(range(len(frame2))) + 10
        drop_sql = "DROP TABLE IF EXISTS test_table2"
        cur = self.db.cursor()
        cur.execute(drop_sql)
        sql.write_frame(frame2, name='test_table2', con=self.db, flavor='mysql')
        result = sql.read_frame("select * from test_table2", self.db,
                                index_col='Idx')
        expected = frame.copy()
        expected.index = Index(range(len(frame2))) + 10
        tm.assert_frame_equal(expected, result)

    def test_tquery(self):
        try:
            import MySQLdb
        except ImportError:
            raise nose.SkipTest
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
            raise nose.SkipTest
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


#########





'''

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
        frame2['Idx'] = Index(range(len(frame2))) + 10
        sql.write_frame(frame2, name='test_table2', con=self.db)
        result = sql.read_frame("select * from test_table2", self.db,
                                index_col='Idx')
        expected = frame.copy()
        expected.index = Index(range(len(frame2))) + 10
        expected.index.name = 'Idx'
        print expected.index.names
        print result.index.names
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
        df = DataFrame({'From':np.ones(5)})
        sql.write_frame(df, con = self.db, name = 'testkeywords')

    def test_onecolumn_of_integer(self):
        'GH 3628, a column_of_integers dataframe should transfer well to sql'
        mono_df=DataFrame([1 , 2], columns=['c0'])
        sql.write_frame(mono_df, con = self.db, name = 'mono_df')
        # computing the sum via sql
        con_x=self.db
        the_sum=sum([my_c0[0] for  my_c0 in con_x.execute("select * from mono_df")])
        # it should not fail, and gives 3 ( Issue #3628 )
        self.assertEqual(the_sum , 3)

        result = sql.read_frame("select * from mono_df",con_x)
        tm.assert_frame_equal(result,mono_df)
'''


class TestSQLA_pymysql(TestSQLAlchemy):
    def set_flavor_engine(self):
        # if can't import should skip all tests
        try:
            import pymysql
        except ImportError:
            # TODO skip all tests
            raise nose.SkipTest("pymysql was not installed")

        try:
            self.engine = sqlalchemy.create_engine("mysql+pymysql://root:@localhost/pandas_nosetest")
        except pymysql.Error, e:
            raise nose.SkipTest(
                "Cannot connect to database. "
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")
        except pymysql.ProgrammingError, e:
            raise nose.SkipTest(
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")


class TestSQLA_mysqldb(TestSQLAlchemy):
    def set_flavor_engine(self):
        # if can't import should skip all tests
        try:
            import MySQLdb
        except ImportError:
            raise nose.SkipTest("MySQLdb was not installed")

        try:
            self.engine = sqlalchemy.create_engine("mysql+MySQLdb://root:@localhost/pandas_nosetest")
        except mysql.Error, e:
            raise nose.SkipTest(
                "Cannot connect to database. "
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")
        except Mysqldb.ProgrammingError, e:
            raise nose.SkipTest(
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")