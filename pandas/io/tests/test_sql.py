from __future__ import print_function
import unittest
import sqlite3
import csv
import os

import numpy as np

#from pandas.core.datetools import format as date_format
from pandas import DataFrame
from pandas.compat import range, lrange, iteritems


import pandas.io.sql as sql
import pandas.util.testing as tm

import sqlalchemy


class TestSQLAlchemy(unittest.TestCase):
    '''
    Test the sqlalchemy backend against an in-memory sqlite database.
    Assume that sqlalchemy takes case of the DB specifics
    '''

    def setUp(self):
        self.engine = sqlalchemy.create_engine('sqlite:///:memory:')
        self._load_iris_data(self.engine)

        self.test_frame_time = tm.makeTimeDataFrame()
        self._load_test1_data()

    def _load_iris_data(self, engine):
        self.dirpath = tm.get_data_path()
        iris_csv_file = os.path.join(self.dirpath, 'iris.csv')
        engine.execute("""CREATE TABLE iris (
                `SepalLength` REAL,
                `SepalWidth` REAL,
                `PetalLength` REAL,
                `PetalWidth` REAL,
                `Name` TEXT
            )""")

        with open(iris_csv_file, 'rU') as iris_csv:
            r = csv.reader(iris_csv)
            next(r)  # skip header row
            ins = """
                INSERT INTO iris
                VALUES(?, ?, ?, ?, ?)
                """
            for row in r:
                engine.execute(ins, *row)

    def _load_test1_data(self):
        test1_csv_file = os.path.join(self.dirpath, 'test1.csv')

        with open(test1_csv_file, 'rU') as test1_csv:
            dr = csv.DictReader(test1_csv)
            self.test_frame1 = DataFrame(list(dr))

    def _test_iris_loaded_frame(self, iris_frame):
        pytype = iris_frame.dtypes[0].type
        row = iris_frame.iloc[0]

        self.assertTrue(issubclass(pytype, np.floating), 'Loaded frame has incorrect type')
        tm.equalContents(row.values, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_read_sql(self):
        iris_frame = sql.read_sql("SELECT * FROM iris", con=self.engine)
        self._test_iris_loaded_frame(iris_frame)

    def test_read_table(self):
        iris_frame = sql.read_table("iris", con=self.engine)
        self._test_iris_loaded_frame(iris_frame)

    def test_to_sql(self):
        # Nuke table
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")

        sql.to_sql(self.test_frame1, 'test_frame1', con=self.engine)
        self.assertTrue(self.engine.has_table('test_frame1'), 'Table not written to DB')

        # Nuke table
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")

    def test_to_sql_fail(self):
        # Nuke table
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")

        sql.to_sql(self.test_frame1, 'test_frame1', con=self.engine, if_exists='fail')
        self.assertTrue(self.engine.has_table('test_frame1'), 'Table not written to DB')

        self.assertRaises(ValueError, sql.to_sql, self.test_frame1, 'test_frame1', con=self.engine, if_exists='fail')

        # Nuke table
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")

    def test_to_sql_replace(self):
        # Nuke table just in case
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.engine, if_exists='fail')
        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.engine, if_exists='replace')
        self.assertTrue(self.engine.has_table('test_frame1'), 'Table not written to DB')

        num_entries = len(self.test_frame1)

        result = self.engine.execute("SELECT count(*) AS count_1 FROM test_frame1").fetchone()
        num_rows = result[0]

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

        # Nuke table
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")

    def test_to_sql_append(self):
        # Nuke table just in case
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.engine, if_exists='fail')
        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.engine, if_exists='append')
        self.assertTrue(self.engine.has_table('test_frame1'), 'Table not written to DB')

        num_entries = 2*len(self.test_frame1)
        result = self.engine.execute("SELECT count(*) AS count_1 FROM test_frame1").fetchone()
        num_rows = result[0]

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

        # Nuke table
        self.engine.execute("DROP TABLE IF EXISTS test_frame1")

    def test_create_table(self):
        temp_engine = sqlalchemy.create_engine('sqlite:///:memory:')
        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithEngine(temp_engine)
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(temp_engine.has_table('temp_frame'), 'Table not written to DB')

    def test_drop_table(self):
        temp_engine = sqlalchemy.create_engine('sqlite:///:memory:')

        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithEngine(temp_engine)
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(temp_engine.has_table('temp_frame'), 'Table not written to DB')

        pandasSQL._drop_table('temp_frame')

        self.assertFalse(temp_engine.has_table('temp_frame'), 'Table not deleted from DB')

    def test_roundtrip(self):
        #temp_engine = sqlalchemy.create_engine('sqlite:///:memory:')

        sql.to_sql(self.test_frame1, 'test_frame_roundtrip', con=self.engine)
        result = sql.read_table('test_frame_roundtrip', con=self.engine)

        # HACK!
        result.index = self.test_frame1.index

        tm.assert_frame_equal(result, self.test_frame1)

    def test_execute_sql(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        iris_results = sql.execute("SELECT * FROM iris", con=self.engine)
        row = iris_results.fetchone()
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_tquery(self):
        iris_results = sql.tquery("SELECT * FROM iris", con=self.engine)
        row = iris_results[0]
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

# --- Test SQLITE fallback


class TestSQLite(unittest.TestCase):
    '''
    Test the sqlalchemy backend against an in-memory sqlite database.
    Assume that sqlalchemy takes case of the DB specifics
    '''

<<<<<<< HEAD
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
=======
    def setUp(self):
        self.conn = sqlite3.connect(':memory:')
        self.pandasSQL = sql.PandasSQLWithCon(self.conn, 'sqlite')

        self._load_iris_data(self.conn)

        self.test_frame_time = tm.makeTimeDataFrame()
        self._load_test1_data()

    def _load_iris_data(self, conn):
        self.dirpath = tm.get_data_path()
        iris_csv_file = os.path.join(self.dirpath, 'iris.csv')
        cur = conn.cursor()
        cur.execute("""CREATE TABLE iris (
                `SepalLength` REAL,
                `SepalWidth` REAL,
                `PetalLength` REAL,
                `PetalWidth` REAL,
                `Name` TEXT
            )""")

        with open(iris_csv_file, 'rU') as iris_csv:
            r = csv.reader(iris_csv)
            next(r)  # skip header row
            ins = """
                INSERT INTO iris
                VALUES(?, ?, ?, ?, ?)
                """
            for row in r:
                cur.execute(ins, row)
        conn.commit()

    def _load_test1_data(self):
        test1_csv_file = os.path.join(self.dirpath, 'test1.csv')

        with open(test1_csv_file, 'rU') as test1_csv:
            dr = csv.DictReader(test1_csv)
            self.test_frame1 = DataFrame(list(dr))

    def test_read_sql(self):
        iris_frame = sql.read_sql("SELECT * FROM iris", con=self.conn)
        pytype = iris_frame.dtypes[0].type
        row = iris_frame.iloc[0]

        self.assertTrue(issubclass(pytype, np.floating), 'Loaded frame has incorrect type')
        tm.equalContents(row.values, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_to_sql(self):
        # Nuke table
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()
        
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.conn, flavor='sqlite')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        # Nuke table
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()

    def test_to_sql_fail(self):
        # Nuke table
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.conn, if_exists='fail', flavor='sqlite')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        self.assertRaises(ValueError, sql.to_sql, self.test_frame1, 'test_frame1', con=self.conn, if_exists='fail')

        # Nuke table
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()

    def test_to_sql_replace(self):
        # Nuke table just in case
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.conn, if_exists='fail', flavor='sqlite')
        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.conn, if_exists='replace')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        num_entries = len(self.test_frame1)

        result = self.conn.execute("SELECT count(*) AS count_1 FROM test_frame1").fetchone()
        num_rows = result[0]

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")
>>>>>>> 1259dca... ENH #4163 Use SQLAlchemy for DB abstraction

        # Nuke table
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()

    def test_to_sql_append(self):
        # Nuke table just in case
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()

        sql.to_sql(self.test_frame1, 'test_frame1', con=self.conn, if_exists='fail', flavor='sqlite')

        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame1', con=self.conn, if_exists='append')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        num_entries = 2*len(self.test_frame1)
        result = self.conn.execute("SELECT count(*) AS count_1 FROM test_frame1").fetchone()
        num_rows = result[0]

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

        # Nuke table
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS test_frame1")
        self.conn.commit()

    def test_create_table(self):
        temp_conn = sqlite3.connect(':memory:')
        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithCon(temp_conn, 'sqlite')
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(pandasSQL.has_table('temp_frame'), 'Table not written to DB')

    def test_drop_table(self):
        temp_conn = sqlite3.connect(':memory:')

        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithCon(temp_conn, 'sqlite')
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(pandasSQL.has_table('temp_frame'), 'Table not written to DB')

        pandasSQL._drop_table('temp_frame')

        self.assertFalse(pandasSQL.has_table('temp_frame'), 'Table not deleted from DB')

    def test_roundtrip(self):

        sql.to_sql(self.test_frame1, 'test_frame_roundtrip', con=self.conn, flavor='sqlite')
        result = sql.read_sql('SELECT * FROM test_frame_roundtrip', con=self.conn, flavor='sqlite')

        # HACK!
        result.index = self.test_frame1.index

        tm.assert_frame_equal(result, self.test_frame1)

    def test_execute_sql(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        iris_results = sql.execute("SELECT * FROM iris", con=self.conn, flavor='sqlite')
        row = iris_results.fetchone()
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_tquery(self):
        iris_results = sql.tquery("SELECT * FROM iris", con=self.conn, flavor='sqlite')
        row = iris_results[0]
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])




"""
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
        except pymysql.Error as e:
            raise nose.SkipTest(
                "Cannot connect to database. "
                "Create a group of conn parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")
        except pymysql.ProgrammingError as e:
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
<<<<<<< HEAD
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
=======
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
"""
>>>>>>> 1259dca... ENH #4163 Use SQLAlchemy for DB abstraction
