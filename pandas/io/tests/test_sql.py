from __future__ import with_statement
from pandas.util.py3compat import StringIO
import unittest
import sqlite3
import sys

import warnings

import nose

import numpy as np

from pandas.core.datetools import format as date_format
from pandas.core.api import DataFrame, isnull

import pandas.io.sql as sql
import pandas.util.testing as tm
from pandas import Series, Index, DataFrame
from datetime import datetime

import sqlalchemy
import sqllite3  # try to import other db modules in their test classes

class TestSQLA_pymysql(TestSQLAlchemy):
    def set_flavor_engine(self):
        # if can't import should skip *all* tests
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
        # if can't import should skip *all* tests
        try:
            import MySQLdb
        except ImportError:
            # TODO skip all tests
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


class TestSQLAlchemy(unittest.TestCase):

    def set_flavor_engine(self):
        # override for other db modules
        self.engine = sqlalchemy.create_engine('sqlite:///:memory:')

    def setUp(self):
        # this is overriden for other db modules        
        set_flavor_engine(self)

        # shared for all db modules
        self.meta = sqlalchemy.schema.MetaData(self.engine)
        self.meta.reflect(engine) # not sure if this is different
        
        mytable = Table('test', meta,
            Column('col2', Integer),
            Column('col3', Integer),

            # table level CHECK constraint.  'name' is optional.
            CheckConstraint('col2 > col3 + 5', name='check1')
        )

        self.frame = tm.makeTimeDataFrame()


    def test_basic(self):
        self._check_roundtrip(self.frame)

    def test_write_row_by_row(self):
        self.frame.ix[0, 0] = np.nan
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        create_sql = sql.get_schema(frame, 'test', self.engine.dialect.name) # TODO let get_schema use engine
        #cur = self.db.cursor()
        #cur.execute(drop_sql)
        self.engine.execute(create_sql)

        ins = test
        test_table = self.meta.tables['test']
        ins = test_table.insert()
        # ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"
        for idx, row in self.frame.iterrows():
            fmt_sql = ins.values(row)  # does passing in a Series work here?
            sql.tquery(fmt_sql, cur=cur)  # make tquery work with engine?

        self.db.commit() # same for this

        select_test = str(test_table.select())

        result = sql.read_frame(select_test, con=self.db)  # have read_frame use engine, pass select expression?
        result.index = frame.index
        tm.assert_frame_equal(result, frame)

    def test_execute(self):
        _skip_if_no_MySQLdb()
        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = sql.get_schema(frame, 'test', 'mysql')
        cur = self.db.cursor()
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