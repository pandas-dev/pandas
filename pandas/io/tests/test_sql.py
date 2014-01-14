from __future__ import print_function
import unittest
import sqlite3
import csv
import os

import numpy as np

from pandas import DataFrame
from pandas.compat import range, lrange, iteritems
#from pandas.core.datetools import format as date_format

import pandas.io.sql as sql
import pandas.util.testing as tm


try:
    import sqlalchemy
    SQLALCHEMY_INSTALLED = True
except ImportError:
    SQLALCHEMY_INSTALLED = False


class PandasSQLTest(unittest.TestCase):
    """Base class with common private methods for
    SQLAlchemy and fallback case test suits"""

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
                engine.execute(ins, row)

    def _check_iris_loaded_frame(self, iris_frame):
        pytype = iris_frame.dtypes[0].type
        row = iris_frame.iloc[0]

        self.assertTrue(issubclass(pytype, np.floating), 'Loaded frame has incorrect type')
        tm.equalContents(row.values, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def _load_test1_data(self):
        test1_csv_file = os.path.join(self.dirpath, 'test1.csv')

        with open(test1_csv_file, 'rU') as test1_csv:
            dr = csv.DictReader(test1_csv)
            self.test_frame1 = DataFrame(list(dr))

    def _count_rows(self, table_name, con):
        result = con.execute("SELECT count(*) AS count_1 FROM %s" % table_name).fetchone()
        return result[0]

    def _read_sql_iris(self):
        iris_frame = self.pandasSQL.read_sql("SELECT * FROM iris")
        self._check_iris_loaded_frame(iris_frame)

    def _to_sql(self):
        # Nuke table
        self.drop_table('test_frame1', self.conn)

        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        # Nuke table
        self.drop_table('test_frame1', self.conn)

    def _to_sql_fail(self):
        # Nuke table
        self.drop_table('test_frame1', self.conn)

        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1', if_exists='fail')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        self.assertRaises(ValueError, self.pandasSQL.to_sql, self.test_frame1, 'test_frame1', if_exists='fail')

        # Nuke table
        self.drop_table('test_frame1', self.conn)

    def _to_sql_replace(self):
        # Nuke table just in case
        self.drop_table('test_frame1', self.conn)

        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1', if_exists='fail')
        # Add to table again
        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1', if_exists='replace')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        num_entries = len(self.test_frame1)
        num_rows = self._count_rows('test_frame1', self.conn)

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

        # Nuke table
        self.drop_table('test_frame1', self.conn)

    def _to_sql_append(self):
        # Nuke table just in case
        self.drop_table('test_frame1', self.conn)

        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1', if_exists='fail')

        # Add to table again
        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1', if_exists='append')
        self.assertTrue(self.pandasSQL.has_table('test_frame1'), 'Table not written to DB')

        num_entries = 2*len(self.test_frame1)
        num_rows = self._count_rows('test_frame1', self.conn)

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

        # Nuke table
        self.drop_table('test_frame1', self.conn)

    def _roundtrip(self):
        self.pandasSQL.to_sql(self.test_frame1, 'test_frame_roundtrip')
        result = self.pandasSQL.read_sql('SELECT * FROM test_frame_roundtrip')

        # HACK!
        result.index = self.test_frame1.index

        tm.assert_frame_equal(result, self.test_frame1)

    def _execute_sql(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        iris_results = self.pandasSQL.execute("SELECT * FROM iris")
        row = iris_results.fetchone()
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def _tquery(self):
        iris_results = self.pandasSQL.tquery("SELECT * FROM iris")
        row = iris_results[0]
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])


class TestSQLApi(PandasSQLTest):
    """Test the public API as it would be used
    directly, including legacy names

    Notes:
    flavor can always be passed even in SQLAlchemy mode,
    should be correctly ignored.

    we don't use drop_table because that isn't part of the public api

    """
    def connect(self):
        if SQLALCHEMY_INSTALLED:
            return sqlalchemy.create_engine('sqlite:///:memory:')
        else:
            return sqlite3.connect(':memory:')

    def setUp(self):
        self.conn = self.connect()
        self._load_iris_data(self.conn)
        self._load_test1_data()

    def test_read_sql_iris(self):
        iris_frame = sql.read_sql("SELECT * FROM iris", self.conn, flavor='sqlite')
        self._check_iris_loaded_frame(iris_frame)

    def test_legacy_read_frame(self):
        """Test legacy name read_frame"""
        iris_frame = sql.read_frame("SELECT * FROM iris", self.conn, flavor='sqlite')
        self._check_iris_loaded_frame(iris_frame)

    def test_to_sql(self):
        sql.to_sql(self.test_frame1, 'test_frame1', self.conn, flavor='sqlite')
        self.assertTrue(sql.has_table('test_frame1', self.conn, flavor='sqlite'), 'Table not written to DB')

    def test_to_sql_fail(self):
        sql.to_sql(self.test_frame1, 'test_frame2', self.conn, flavor='sqlite', if_exists='fail')
        self.assertTrue(sql.has_table('test_frame2', self.conn, flavor='sqlite'), 'Table not written to DB')

        self.assertRaises(ValueError, sql.to_sql, self.test_frame1, 'test_frame2', self.conn, flavor='sqlite', if_exists='fail')

    def test_to_sql_replace(self):
        sql.to_sql(self.test_frame1, 'test_frame3', self.conn, flavor='sqlite', if_exists='fail')
        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame3', self.conn, flavor='sqlite', if_exists='replace')
        self.assertTrue(sql.has_table('test_frame3', self.conn, flavor='sqlite'), 'Table not written to DB')

        num_entries = len(self.test_frame1)
        num_rows = self._count_rows('test_frame3', self.conn)

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

    def test_to_sql_append(self):
        sql.to_sql(self.test_frame1, 'test_frame4', self.conn, flavor='sqlite', if_exists='fail')

        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame4', self.conn, flavor='sqlite', if_exists='append')
        self.assertTrue(sql.has_table('test_frame4', self.conn, flavor='sqlite'), 'Table not written to DB')

        num_entries = 2*len(self.test_frame1)
        num_rows = self._count_rows('test_frame4', self.conn)

        self.assertEqual(num_rows, num_entries, "not the same number of rows as entries")

    def test_legacy_write_frame(self):
        """Test legacy write frame name.
        Assume that functionality is already tested above so just do quick check that it basically works"""
        sql.write_frame(self.test_frame1, 'test_frame_legacy', self.conn, flavor='sqlite')
        self.assertTrue(sql.has_table('test_frame_legacy', self.conn, flavor='sqlite'), 'Table not written to DB')

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


class TestSQLAlchemy(PandasSQLTest):
    '''
    Test the sqlalchemy backend against an in-memory sqlite database.
    Assume that sqlalchemy takes case of the DB specifics
    '''
    def connect(self):
        return sqlalchemy.create_engine('sqlite:///:memory:')

    def drop_table(self, table_name, conn):
        conn.execute("DROP TABLE IF EXISTS %s" % table_name)

    def setUp(self):
        # Skip this test if SQLAlchemy not available
        if not SQLALCHEMY_INSTALLED:
            raise unittest.SkipTest('SQLAlchemy not installed')

        self.conn = self.connect()
        self.pandasSQL = sql.PandasSQLWithEngine(self.conn)

        self._load_iris_data(self.conn)

        self._load_test1_data()

    def test_read_sql(self):
        self._read_sql_iris()

    def test_read_table(self):
        iris_frame = sql.read_table("iris", con=self.conn)
        self._check_iris_loaded_frame(iris_frame)

    def test_read_table_absent(self):
        self.assertRaises(ValueError, sql.read_table, "this_doesnt_exist", con=self.conn)

    def test_to_sql(self):
        self._to_sql()

    def test_to_sql_fail(self):
        self._to_sql_fail()

    def test_to_sql_replace(self):
        self._to_sql_replace()

    def test_to_sql_append(self):
        self._to_sql_append()

    def test_create_table(self):
        temp_conn = self.connect()
        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithEngine(temp_conn)
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(temp_conn.has_table('temp_frame'), 'Table not written to DB')

    def test_drop_table(self):
        temp_conn = self.connect()

        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithEngine(temp_conn)
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(temp_conn.has_table('temp_frame'), 'Table not written to DB')

        pandasSQL._drop_table('temp_frame')

        self.assertFalse(temp_conn.has_table('temp_frame'), 'Table not deleted from DB')

    def test_roundtrip(self):
        self._roundtrip()

    def test_execute_sql(self):
        self._execute_sql()

    def test_tquery(self):
        self._tquery()


# --- Test SQLITE fallback


class TestSQLite(PandasSQLTest):
    '''
    Test the sqlalchemy backend against an in-memory sqlite database.
    Assume that sqlalchemy takes case of the DB specifics
    '''
    def connect(self):
        return sqlite3.connect(':memory:')

    def drop_table(self, table_name, conn):
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS %s" % table_name)
        self.conn.commit()

    def setUp(self):
        self.conn = self.connect()
        self.pandasSQL = sql.PandasSQLWithCon(self.conn, 'sqlite')

        self._load_iris_data(self.conn)

        self._load_test1_data()

    def test_invalid_flavor(self):
        self.assertRaises(NotImplementedError, sql.PandasSQLWithCon, self.conn, 'oracle')

    def test_read_sql(self):
        self._read_sql_iris()

    def test_to_sql(self):
        self._to_sql()

    def test_to_sql_fail(self):
        self._to_sql_fail()

    def test_to_sql_replace(self):
        self._to_sql_replace()

    def test_to_sql_append(self):
        self._to_sql_append()

    def test_create_table(self):
        temp_conn = self.connect()
        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithCon(temp_conn, 'sqlite')
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(pandasSQL.has_table('temp_frame'), 'Table not written to DB')

    def test_drop_table(self):
        temp_conn = self.connect()

        temp_frame = DataFrame({'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLWithCon(temp_conn, 'sqlite')
        pandasSQL._create_table(temp_frame, 'temp_frame')

        self.assertTrue(pandasSQL.has_table('temp_frame'), 'Table not written to DB')

        pandasSQL._drop_table('temp_frame')

        self.assertFalse(pandasSQL.has_table('temp_frame'), 'Table not deleted from DB')

    def test_roundtrip(self):
        self._roundtrip()

    def test_execute_sql(self):
        self._execute_sql()

    def test_tquery(self):
        self._tquery()


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
