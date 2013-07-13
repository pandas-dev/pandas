from __future__ import print_function
import unittest
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
from pandas.io.sql import DatabaseError
import pandas.util.testing as tm
from pandas import Series, Index, DataFrame, isnull
from datetime import datetime

import sqlalchemy
import sqlite3  # try to import other db modules in their test classes

from sqlalchemy import Table, Column, INT, FLOAT, TEXT


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
        sql.PandasSQLWithEngine(engine=self.engine)._drop_table(table_name)

    def create_table(self, frame, table_name, keys=None):
        return sql.PandasSQLWithEngine(engine=self.engine)._create_table(frame, table_name, keys=None)

    def get_table(self, table_name):
        return sql.PandasSQLWithEngine(self.engine)._get_table(table_name)

    def tquery(self, fmt_sql, params=None, retry=False):
        sql.tquery(fmt_sql, engine=self.engine, params=params, retry=retry)

    def read_frame(self, fmt_sql=None):
        return sql.read_frame(fmt_sql, engine=self.engine)

    def _check_roundtrip(self, frame):
        self.drop_table('test')
        sql.write_frame(self.frame, 'test', engine=self.engine)
        result = sql._engine_read_table_name('test', engine=self.engine)

        # HACK!
        result.index = self.frame.index

        tm.assert_frame_equal(result, self.frame)

        self.frame['txt'] = ['a'] * len(self.frame)
        frame2 = self.frame.copy()
        frame2['Idx'] = Index(range(len(frame2))) + 10

        self.drop_table('test_table2')
        sql.write_frame(frame2, 'test_table2', engine=self.engine)
        result = sql._engine_read_table_name('test_table2', engine=self.engine, index_col='Idx')

        self.assertRaises(DatabaseError, self.tquery,
                          'insert into blah values (1)')

        self.assertRaises(DatabaseError, self.tquery,
                          'insert into blah values (1)',
                          retry=True)


    def test_basic(self):
        self._check_roundtrip(self.frame)

    # not sure what intention of this was?
    def test_na_roundtrip(self):
        pass

    def test_write_row_by_row(self):
        self.frame.ix[0, 0] = np.nan
        self.create_table(self.frame, 'test')

        test_table = self.get_table('test')

        ins = test_table.insert() # INSERT INTO test VALUES (%s, %s, %s, %s)
        for idx, row in self.frame.iterrows():
            values = tuple(row)
            sql.execute(ins.values(values), engine=self.engine)

        select_test = test_table.select()  # SELECT * FROM test

        result = self.read_frame(select_test)

        result.index = self.frame.index
        tm.assert_frame_equal(result, self.frame)

    def test_execute(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        self.create_table(self.frame, 'test')

        test_table = self.get_table('test')
        
        ins = test_table.insert() # INSERT INTO test VALUES (%s, %s, %s, %s)

        row = self.frame.ix[0]
        self.engine.execute(ins, **row)

        select_test = test_table.select()  # SELECT * FROM test
        result = self.read_frame(select_test)
        result.index = self.frame.index[:1]
        tm.assert_frame_equal(result, self.frame[:1])

    def test_execute_fail(self):
        """
        CREATE TABLE test
        (
        a TEXT,
        b TEXT,
        c REAL,
        PRIMARY KEY (a, b)
        );
        """
        from sqlalchemy import Table, Column, TEXT, REAL
        test_table = Table('test', self.meta,
                        Column('a', TEXT), Column('b', TEXT), Column('c', REAL))
        test_table.create()

        sql.execute('INSERT INTO test VALUES("foo", "bar", 1.234)', engine=self.engine)
        sql.execute('INSERT INTO test VALUES("foo", "baz", 2.567)', engine=self.engine)

        self.assertRaises(DatabaseError, sql.execute,
                              'INSERT INTO test VALUES("foo", "bar", 7)',
                              engine=self.engine)

    def test_tquery(self):
        self.drop_table('test_table')
        sql.write_frame(self.frame, 'test_table', engine=self.engine)
        result = sql.tquery("select A from test_table", engine=self.engine)
        expected = self.frame.A
        result = DataFrame(result, self.frame.index, columns=['A'])['A']
        tm.assert_series_equal(result, expected)

        self.assertRaises(DatabaseError, sql.tquery,
                           'select * from blah', engine=self.engine)

        self.assertRaises(DatabaseError, sql.tquery,
                           'select * from blah', engine=self.engine, retry=True)

    def test_uquery(self):
        self.drop_table('test_table')
        sql.write_frame(self.frame, 'test_table', engine=self.engine)

        ins = sql.PandasSQLWithEngine(self.engine)._get_table('test_table').insert()
        params = (2.314, -123.1, 1.234, 2.3)
        self.assertEqual(sql.uquery(ins, params=params, engine=self.engine), 1)

        self.assertRaises(DatabaseError, sql.uquery,
                            'insert into blah values (1)', engine=self.engine)

        self.assertRaises(DatabaseError, sql.tquery,
                            'insert into blah values (1)', engine=self.engine, retry=True)


    def test_onecolumn_of_integer(self):
        'GH 3628, a column_of_integers dataframe should transfer well to sql'
        mono_df = DataFrame([1 , 2], columns=['c0'])
        sql.write_frame(mono_df, 'mono_df', engine=self.engine)
        # computing the sum via sql
        select = sql.PandasSQLWithEngine(self.engine)._get_table('mono_df').select()
        the_sum = sum([my_c0[0] for  my_c0 in self.engine.execute(select)])
        # it should not fail, and gives 3 ( Issue #3628 )
        self.assertEqual(the_sum , 3)

        result = sql._engine_read_table_name('mono_df', engine=self.engine)
        tm.assert_frame_equal(result, mono_df)

    def test_keyword_as_column_names(self):
        df = DataFrame({'From':np.ones(5)})
        sql.write_frame(df, engine=self.engine, name='testkeywords',
                        flavor='mysql', if_exists='replace')


    # Not needed with engines, but add into con/cur tests later

    # def test_execute_closed_connection(self):
    #     create_sql = """
    #     CREATE TABLE test
    #     (
    #     a TEXT,
    #     b TEXT,
    #     c REAL,
    #     PRIMARY KEY (a, b)
    #     );
    #     """
    #     cur = self.db.cursor()
    #     cur.execute(create_sql)

    #     sql.execute('INSERT INTO test VALUES("foo", "bar", 1.234)', self.db)
    #     self.db.close()
    #     try:
    #         sys.stdout = StringIO()
    #         self.assertRaises(Exception, sql.tquery, "select * from test",
    #                           con=self.db)
    #     finally:
    #         sys.stdout = sys.__stdout__

    # def test_schema(self):
    #     create_sql = self.create_table(self.frame, 'test')[1]
    #     lines = create_sql.splitlines()
    #     for l in lines:
    #         tokens = l.split(' ')
    #         if len(tokens) == 2 and tokens[0] == 'A':
    #             self.assert_(tokens[1] == 'DATETIME')
    #     self.drop_table('test')
    #     create_sql = self.create_table(frame, 'test', keys=['A', 'B'])[1]
    #     self.assert_('PRIMARY KEY (A,B)' in create_sql)


class TestSQLA_pymysql(TestSQLAlchemy):
    def setUp(self):
        raise nose.SkipTest("MySQLdb was not installed")

    def set_flavor_engine(self):
        # if can't import should skip all tests
        try:
            import pymysql
        except ImportError:
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


class TestSQLA_MySQLdb(TestSQLAlchemy):
    def setUp(self):
        raise nose.SkipTest("MySQLdb was not installed")

    def set_flavor_engine(self):
        # if can't import should skip all tests
        try:
            import MySQLdb
        except ImportError:
            raise nose.SkipTest("MySQLdb was not installed")

        try:
            self.engine = sqlalchemy.create_engine("mysql+mysqldb://root:@localhost/pandas_nosetest")
        except MySQLdb.Error:
            raise nose.SkipTest(
                "Cannot connect to database. "
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")
        except MySQLdb.ProgrammingError:
            raise nose.SkipTest(
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")