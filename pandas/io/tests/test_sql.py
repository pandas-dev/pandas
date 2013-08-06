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
        sql._engine_drop_table(table_name, engine=self.engine)

    def create_table(self, frame, table_name, keys=None):
        return sql._engine_create_table(frame, table_name, keys=None, engine=self.engine)

    def get_table(self, table_name):
        return sql._engine_get_table(table_name, self.engine)

    def tquery(self, fmt_sql, params=None):
        sql.tquery(fmt_sql, engine=self.engine, params=params)

    def read_frame(self, fmt_sql=None):
        return sql.read_frame(fmt_sql, engine=self.engine) 

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
        self.engine.execute(ins, **dict(*row))

        select_test = test_table.select()  # SELECT * FROM test
        result = self.read_frame(select_test)
        result.index = self.frame.index[:1]
        tm.assert_frame_equal(result, self.frame[:1])

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

    def test_na_roundtrip(self):
        pass

    def _check_roundtrip(self, frame):
        self.drop_table('test')
        sql._engine_write_frame(self.frame, 'test', self.engine)
        result = sql._engine_read_table_name('test', engine=self.engine)

        # HACK!
        result.index = self.frame.index

        tm.assert_frame_equal(result, self.frame)

        self.frame['txt'] = ['a'] * len(self.frame)
        frame2 = self.frame.copy()
        frame2['Idx'] = Index(range(len(frame2))) + 10

        self.drop_table('test_table2')
        sql._engine_write_frame(frame2, 'test_table2', self.engine)
        result = sql._engine_read_table_name('test_table2', engine=self.engine, index_col='Idx')

        tm.assert_frame_equal(frame2, result)

    def test_tquery(self):
        self.drop_table('test_table')
        sql._engine_write_frame(self.frame, 'test_table', self.engine)
        result = sql.tquery("select A from test_table", engine=self.engine)
        expected = self.frame.A
        result = Series(result, self.frame.index)
        tm.assert_series_equal(result, expected)

        # try:
        #     sys.stdout = StringIO()
        #     self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
        #                       'select * from blah', con=self.db)

        #     self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
        #                       'select * from blah', con=self.db, retry=True)
        # finally:
        #     sys.stdout = sys.__stdout__

    def test_uquery(self):
        self.drop_table('test_table')
        sql._engine_write_frame(self.frame, 'test_table', self.engine)

        ins = sql._engine_get_table('test_table', self.engine).insert()
        params = (2.314, -123.1, 1.234, 2.3)
        self.assertEqual(sql.uquery(ins, params, engine=self.engine), 1)

        # try:
        #     sys.stdout = StringIO()

        #     self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
        #                       'insert into blah values (1)', con=self.db)

        #     self.assertRaises(MySQLdb.ProgrammingError, sql.tquery,
        #                       'insert into blah values (1)', con=self.db,
        #                       retry=True)
        # finally:
        #     sys.stdout = sys.__stdout__

    def test_keyword_as_column_names(self):
        df = DataFrame({'From':np.ones(5)})
        sql._engine_write_frame(df, 'testkeywords',
                        if_exists='replace', engine=self.engine)

    def test_onecolumn_of_integer(self):
        'GH 3628, a column_of_integers dataframe should transfer well to sql'
        mono_df = DataFrame([1 , 2], columns=['c0'])
        sql._engine_write_frame(mono_df, 'mono_df', self.engine)
        # computing the sum via sql
        select = sql._engine_get_table('mono_df', self.engine).select()
        the_sum = sum([my_c0[0] for  my_c0 in engine.execute(select)])
        # it should not fail, and gives 3 ( Issue #3628 )
        self.assertEqual(the_sum , 3)

        result = sql._engine_read_table_name('mono_df', engine=self.engine)
        tm.assert_frame_equal(result, mono_df)



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