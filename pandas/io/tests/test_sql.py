from __future__ import print_function
import unittest
import sqlite3
import csv
import os

import nose
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

SQL_STRINGS = {
    'create_iris': {
        'sqlite': """CREATE TABLE iris (
                `SepalLength` REAL,
                `SepalWidth` REAL,
                `PetalLength` REAL,
                `PetalWidth` REAL,
                `Name` TEXT
            )""",
        'mysql': """CREATE TABLE iris (
                `SepalLength` DOUBLE,
                `SepalWidth` DOUBLE,
                `PetalLength` DOUBLE,
                `PetalWidth` DOUBLE,
                `Name` VARCHAR(200)
            )"""
    },
    'insert_iris': {
        'sqlite': """INSERT INTO iris VALUES(?, ?, ?, ?, ?)""",
        'mysql': """INSERT INTO iris VALUES(%s, %s, %s, %s, "%s");"""
    },
    'create_test_types': {
        'sqlite': """CREATE TABLE types_test_data (
                    `TextCol` TEXT,
                    `DateCol` TEXT,
                    `IntDateCol` INTEGER,
                    `FloatCol` REAL,
                    `IntCol` INTEGER,
                    `BoolCol` INTEGER,
                    `IntColWithNull` INTEGER,
                    `BoolColWithNull` INTEGER
                )""",
        'mysql': """CREATE TABLE types_test_data (
                    `TextCol` TEXT,
                    `DateCol` DATETIME,
                    `IntDateCol` INTEGER,
                    `FloatCol` DOUBLE,
                    `IntCol` INTEGER,
                    `BoolCol` BOOLEAN,
                    `IntColWithNull` INTEGER,
                    `BoolColWithNull` BOOLEAN
                )"""
    },
    'insert_test_types': {
        'sqlite': """
                INSERT INTO types_test_data
                VALUES(?, ?, ?, ?, ?, ?, ?, ?)
                """,
        'mysql': """
                INSERT INTO types_test_data
                VALUES("%s", %s, %s, %s, %s, %s, %s, %s)
                """
    }
}


class PandasSQLTest(unittest.TestCase):

    """Base class with common private methods for
    SQLAlchemy and fallback cases.
    """

    def drop_table(self, table_name):
        self._get_exec().execute("DROP TABLE IF EXISTS %s" % table_name)

    def _get_exec(self):
        if hasattr(self.conn, 'execute'):
            return self.conn
        else:
            return self.conn.cursor()

    def _load_iris_data(self):
        iris_csv_file = os.path.join(tm.get_data_path(), 'iris.csv')

        self.drop_table('iris')
        self._get_exec().execute(SQL_STRINGS['create_iris'][self.flavor])

        with open(iris_csv_file, 'rU') as iris_csv:
            r = csv.reader(iris_csv)
            next(r)  # skip header row
            ins = SQL_STRINGS['insert_iris'][self.flavor]

            for row in r:
                self._get_exec().execute(ins, row)

    def _check_iris_loaded_frame(self, iris_frame):
        pytype = iris_frame.dtypes[0].type
        row = iris_frame.iloc[0]

        self.assertTrue(
            issubclass(pytype, np.floating), 'Loaded frame has incorrect type')
        tm.equalContents(row.values, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def _load_test1_data(self):
        columns = ['index', 'A', 'B', 'C', 'D']
        data = [(
            '2000-01-03 00:00:00', 0.980268513777, 3.68573087906, -0.364216805298, -1.15973806169),
            ('2000-01-04 00:00:00', 1.04791624281, -
             0.0412318367011, -0.16181208307, 0.212549316967),
            ('2000-01-05 00:00:00', 0.498580885705,
             0.731167677815, -0.537677223318, 1.34627041952),
            ('2000-01-06 00:00:00', 1.12020151869, 1.56762092543, 0.00364077397681, 0.67525259227)]

        self.test_frame1 = DataFrame(data, columns=columns)

    def _load_raw_sql(self):
        self.drop_table('types_test_data')
        self._get_exec().execute(SQL_STRINGS['create_test_types'][self.flavor])
        ins = SQL_STRINGS['insert_test_types'][self.flavor]

        data = [(
            'first', '2000-01-03 00:00:00', 535852800, 10.10, 1, False, 1, False),
            ('first', '2000-01-04 00:00:00', 1356998400, 10.10, 1, False, None, None)]
        for d in data:
            self._get_exec().execute(ins, d)

    def _count_rows(self, table_name):
        result = self._get_exec().execute(
            "SELECT count(*) AS count_1 FROM %s" % table_name).fetchone()
        return result[0]

    def _read_sql_iris(self):
        iris_frame = self.pandasSQL.read_sql("SELECT * FROM iris")
        self._check_iris_loaded_frame(iris_frame)

    def _to_sql(self):
        self.drop_table('test_frame1')

        self.pandasSQL.to_sql(self.test_frame1, 'test_frame1')
        self.assertTrue(self.pandasSQL.has_table(
            'test_frame1'), 'Table not written to DB')

        # Nuke table
        self.drop_table('test_frame1')

    def _to_sql_fail(self):
        self.drop_table('test_frame1')

        self.pandasSQL.to_sql(
            self.test_frame1, 'test_frame1', if_exists='fail')
        self.assertTrue(self.pandasSQL.has_table(
            'test_frame1'), 'Table not written to DB')

        self.assertRaises(ValueError, self.pandasSQL.to_sql,
                          self.test_frame1, 'test_frame1', if_exists='fail')

        self.drop_table('test_frame1')

    def _to_sql_replace(self):
        self.drop_table('test_frame1')

        self.pandasSQL.to_sql(
            self.test_frame1, 'test_frame1', if_exists='fail')
        # Add to table again
        self.pandasSQL.to_sql(
            self.test_frame1, 'test_frame1', if_exists='replace')
        self.assertTrue(self.pandasSQL.has_table(
            'test_frame1'), 'Table not written to DB')

        num_entries = len(self.test_frame1)
        num_rows = self._count_rows('test_frame1')

        self.assertEqual(
            num_rows, num_entries, "not the same number of rows as entries")

        self.drop_table('test_frame1')

    def _to_sql_append(self):
        # Nuke table just in case
        self.drop_table('test_frame1')

        self.pandasSQL.to_sql(
            self.test_frame1, 'test_frame1', if_exists='fail')

        # Add to table again
        self.pandasSQL.to_sql(
            self.test_frame1, 'test_frame1', if_exists='append')
        self.assertTrue(self.pandasSQL.has_table(
            'test_frame1'), 'Table not written to DB')

        num_entries = 2 * len(self.test_frame1)
        num_rows = self._count_rows('test_frame1')

        self.assertEqual(
            num_rows, num_entries, "not the same number of rows as entries")

        self.drop_table('test_frame1')

    def _roundtrip(self):
        self.drop_table('test_frame_roundtrip')
        self.pandasSQL.to_sql(self.test_frame1, 'test_frame_roundtrip')
        result = self.pandasSQL.read_sql('SELECT * FROM test_frame_roundtrip')

        result.set_index('pandas_index', inplace=True)
        #result.index.astype(int)

        result.index.name = None

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
    flavor = 'sqlite'

    def connect(self):
        if SQLALCHEMY_INSTALLED:
            return sqlalchemy.create_engine('sqlite:///:memory:')
        else:
            return sqlite3.connect(':memory:')

    def setUp(self):
        self.conn = self.connect()
        self._load_iris_data()
        self._load_test1_data()
        self._load_raw_sql()

    def test_read_sql_iris(self):
        iris_frame = sql.read_sql(
            "SELECT * FROM iris", self.conn, flavor='sqlite')
        self._check_iris_loaded_frame(iris_frame)

    def test_legacy_read_frame(self):
        """Test legacy name read_frame"""
        iris_frame = sql.read_frame(
            "SELECT * FROM iris", self.conn, flavor='sqlite')
        self._check_iris_loaded_frame(iris_frame)

    def test_to_sql(self):
        sql.to_sql(self.test_frame1, 'test_frame1', self.conn, flavor='sqlite')
        self.assertTrue(
            sql.has_table('test_frame1', self.conn, flavor='sqlite'), 'Table not written to DB')

    def test_to_sql_fail(self):
        sql.to_sql(self.test_frame1, 'test_frame2',
                   self.conn, flavor='sqlite', if_exists='fail')
        self.assertTrue(
            sql.has_table('test_frame2', self.conn, flavor='sqlite'), 'Table not written to DB')

        self.assertRaises(ValueError, sql.to_sql, self.test_frame1,
                          'test_frame2', self.conn, flavor='sqlite', if_exists='fail')

    def test_to_sql_replace(self):
        sql.to_sql(self.test_frame1, 'test_frame3',
                   self.conn, flavor='sqlite', if_exists='fail')
        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame3',
                   self.conn, flavor='sqlite', if_exists='replace')
        self.assertTrue(
            sql.has_table('test_frame3', self.conn, flavor='sqlite'), 'Table not written to DB')

        num_entries = len(self.test_frame1)
        num_rows = self._count_rows('test_frame3')

        self.assertEqual(
            num_rows, num_entries, "not the same number of rows as entries")

    def test_to_sql_append(self):
        sql.to_sql(self.test_frame1, 'test_frame4',
                   self.conn, flavor='sqlite', if_exists='fail')

        # Add to table again
        sql.to_sql(self.test_frame1, 'test_frame4',
                   self.conn, flavor='sqlite', if_exists='append')
        self.assertTrue(
            sql.has_table('test_frame4', self.conn, flavor='sqlite'), 'Table not written to DB')

        num_entries = 2 * len(self.test_frame1)
        num_rows = self._count_rows('test_frame4')

        self.assertEqual(
            num_rows, num_entries, "not the same number of rows as entries")

    def test_legacy_write_frame(self):
        """Test legacy write frame name.
        Assume that functionality is already tested above so just do quick check that it basically works"""
        sql.write_frame(
            self.test_frame1, 'test_frame_legacy', self.conn, flavor='sqlite')
        self.assertTrue(
            sql.has_table('test_frame_legacy', self.conn, flavor='sqlite'), 'Table not written to DB')

    def test_roundtrip(self):
        sql.to_sql(self.test_frame1, 'test_frame_roundtrip',
                   con=self.conn, flavor='sqlite')
        result = sql.read_sql(
            'SELECT * FROM test_frame_roundtrip', con=self.conn, flavor='sqlite')

        # HACK!
        result.index = self.test_frame1.index
        result.set_index('pandas_index', inplace=True)
        result.index.astype(int)
        result.index.name = None
        tm.assert_frame_equal(result, self.test_frame1)

    def test_execute_sql(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        iris_results = sql.execute(
            "SELECT * FROM iris", con=self.conn, flavor='sqlite')
        row = iris_results.fetchone()
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_tquery(self):
        iris_results = sql.tquery(
            "SELECT * FROM iris", con=self.conn, flavor='sqlite')
        row = iris_results[0]
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_date_parsing(self):
        """ Test date parsing in read_sql """
        # No Parsing
        df = sql.read_sql(
            "SELECT * FROM types_test_data", self.conn, flavor='sqlite')
        self.assertFalse(
            issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

        df = sql.read_sql("SELECT * FROM types_test_data",
                          self.conn, flavor='sqlite', parse_dates=['DateCol'])
        self.assertTrue(
            issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

        df = sql.read_sql("SELECT * FROM types_test_data", self.conn,
                          flavor='sqlite', parse_dates={'DateCol': '%Y-%m-%d %H:%M:%S'})
        self.assertTrue(
            issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

        df = sql.read_sql("SELECT * FROM types_test_data",
                          self.conn, flavor='sqlite', parse_dates=['IntDateCol'])
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_sql("SELECT * FROM types_test_data",
                          self.conn, flavor='sqlite', parse_dates={'IntDateCol': 's'})
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")


class TestSQLAlchemy(PandasSQLTest):

    '''
    Test the sqlalchemy backend against an in-memory sqlite database.
    Assume that sqlalchemy takes case of the DB specifics
    '''
    flavor = 'sqlite'

    def connect(self):
        return sqlalchemy.create_engine('sqlite:///:memory:')

    def setUp(self):
        # Skip this test if SQLAlchemy not available
        if not SQLALCHEMY_INSTALLED:
            raise nose.SkipTest('SQLAlchemy not installed')

        self.conn = self.connect()
        self.pandasSQL = sql.PandasSQLAlchemy(self.conn)

        self._load_iris_data()
        self._load_raw_sql()

        self._load_test1_data()

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
        temp_frame = DataFrame(
            {'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLAlchemy(temp_conn)
        pandasSQL.to_sql(temp_frame, 'temp_frame')

        self.assertTrue(
            temp_conn.has_table('temp_frame'), 'Table not written to DB')

    def test_drop_table(self):
        temp_conn = self.connect()

        temp_frame = DataFrame(
            {'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        pandasSQL = sql.PandasSQLAlchemy(temp_conn)
        pandasSQL.to_sql(temp_frame, 'temp_frame')

        self.assertTrue(
            temp_conn.has_table('temp_frame'), 'Table not written to DB')

        pandasSQL.drop_table('temp_frame')

        self.assertFalse(
            temp_conn.has_table('temp_frame'), 'Table not deleted from DB')

    def test_roundtrip(self):
        self._roundtrip()

    def test_execute_sql(self):
        self._execute_sql()

    def test_read_table(self):
        iris_frame = sql.read_table("iris", con=self.conn)
        self._check_iris_loaded_frame(iris_frame)

    def test_read_table_columns(self):
        iris_frame = sql.read_table(
            "iris", con=self.conn, columns=['SepalLength', 'SepalLength'])
        tm.equalContents(
            iris_frame.columns.values, ['SepalLength', 'SepalLength'])

    def test_read_table_absent(self):
        self.assertRaises(
            ValueError, sql.read_table, "this_doesnt_exist", con=self.conn)

    def test_default_type_convertion(self):
        """ Test default type conversion"""
        df = sql.read_table("types_test_data", self.conn)
        self.assertTrue(
            issubclass(df.FloatCol.dtype.type, np.floating), "FloatCol loaded with incorrect type")
        self.assertTrue(
            issubclass(df.IntCol.dtype.type, np.integer), "IntCol loaded with incorrect type")
        self.assertTrue(
            issubclass(df.BoolCol.dtype.type, np.integer), "BoolCol loaded with incorrect type")

        # Int column with NA values stays as float
        self.assertTrue(issubclass(df.IntColWithNull.dtype.type, np.floating),
                        "IntColWithNull loaded with incorrect type")
        # Non-native Bool column with NA values stays as float
        self.assertTrue(
            issubclass(df.BoolColWithNull.dtype.type, np.floating), "BoolCol loaded with incorrect type")

    def test_default_date_load(self):
        df = sql.read_table("types_test_data", self.conn)

        # IMPORTANT - sqlite has no native date type, so shouldn't parse, but
        # MySQL SHOULD be converted.
        self.assertFalse(
            issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

    def test_date_parsing(self):
        """ Test date parsing """
        # No Parsing
        df = sql.read_table("types_test_data", self.conn)

        df = sql.read_table(
            "types_test_data", self.conn, parse_dates=['DateCol'])
        self.assertTrue(
            issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

        df = sql.read_table(
            "types_test_data", self.conn, parse_dates={'DateCol': '%Y-%m-%d %H:%M:%S'})
        self.assertTrue(
            issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

        df = sql.read_table("types_test_data", self.conn, parse_dates={
                            'DateCol': {'format': '%Y-%m-%d %H:%M:%S'}})
        self.assertTrue(issubclass(df.DateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_table(
            "types_test_data", self.conn, parse_dates=['IntDateCol'])
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_table(
            "types_test_data", self.conn, parse_dates={'IntDateCol': 's'})
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_table(
            "types_test_data", self.conn, parse_dates={'IntDateCol': {'unit': 's'}})
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")


# --- Test SQLITE fallback
class TestSQLite(PandasSQLTest):

    '''
    Test the sqlalchemy backend against an in-memory sqlite database.
    Assume that sqlalchemy takes case of the DB specifics
    '''
    flavor = 'sqlite'

    def connect(self):
        return sqlite3.connect(':memory:')

    def drop_table(self, table_name):
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS %s" % table_name)
        self.conn.commit()

    def setUp(self):
        self.conn = self.connect()
        self.pandasSQL = sql.PandasSQLLegacy(self.conn, 'sqlite')

        self._load_iris_data()

        self._load_test1_data()

    def test_invalid_flavor(self):
        self.assertRaises(
            NotImplementedError, sql.PandasSQLLegacy, self.conn, 'oracle')

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

    def test_create_and_drop_table(self):
        temp_frame = DataFrame(
            {'one': [1., 2., 3., 4.], 'two': [4., 3., 2., 1.]})

        self.pandasSQL.to_sql(temp_frame, 'drop_test_frame')

        self.assertTrue(self.pandasSQL.has_table(
            'drop_test_frame'), 'Table not written to DB')

        self.pandasSQL.drop_table('drop_test_frame')

        self.assertFalse(self.pandasSQL.has_table(
            'drop_test_frame'), 'Table not deleted from DB')

    def test_roundtrip(self):
        self._roundtrip()

    def test_execute_sql(self):
        self._execute_sql()

    def test_tquery(self):
        self._tquery()


class TestMySQL(TestSQLite):
    flavor = 'mysql'

    def drop_table(self, table_name):
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS %s" % table_name)
        self.conn.commit()

    def _count_rows(self, table_name):
        cur = self._get_exec()
        cur.execute(
            "SELECT count(*) AS count_1 FROM %s" % table_name)
        rows = cur.fetchall()
        return rows[0][0]

    def connect(self):
        return self.driver.connect(host='127.0.0.1', user='root', passwd='', db='pandas_nosetest')

    def setUp(self):
        try:
            import pymysql
            self.driver = pymysql

        except ImportError:
            raise nose.SkipTest

        self.conn = self.connect()
        self.pandasSQL = sql.PandasSQLLegacy(self.conn, 'mysql')

        self._load_iris_data()
        self._load_test1_data()

    def tearDown(self):
        c = self.conn.cursor()
        c.execute('SHOW TABLES')
        for table in c.fetchall():
            c.execute('DROP TABLE %s' % table[0])
        self.conn.commit()
        self.conn.close()


class TestMySQLAlchemy(TestSQLAlchemy):
        flavor = 'mysql'

        def connect(self):
            return sqlalchemy.create_engine(
                'mysql+{driver}://root@localhost/pandas_nosetest'.format(driver=self.driver))

        def setUp(self):
            if not SQLALCHEMY_INSTALLED:
                raise nose.SkipTest('SQLAlchemy not installed')

            try:
                import pymysql
                self.driver = 'pymysql'

            except ImportError:
                raise nose.SkipTest

            self.conn = self.connect()
            self.pandasSQL = sql.PandasSQLAlchemy(self.conn)

            self._load_iris_data()
            self._load_raw_sql()

            self._load_test1_data()

        def tearDown(self):
            c = self.conn.execute('SHOW TABLES')
            for table in c.fetchall():
                self.conn.execute('DROP TABLE %s' % table[0])

        def test_default_date_load(self):
            df = sql.read_table("types_test_data", self.conn)

            # IMPORTANT - sqlite has no native date type, so shouldn't parse,
            # but MySQL SHOULD be converted.
            self.assertTrue(
                issubclass(df.DateCol.dtype.type, np.datetime64), "DateCol loaded with incorrect type")

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
