"""SQL io tests

The SQL tests are broken down in different classes:

- `PandasSQLTest`: base class with common methods for all test classes
- Tests for the public API (only tests with sqlite3)
    - `_TestSQLApi` base class
    - `TestSQLApi`: test the public API with sqlalchemy engine
    - `TesySQLLegacyApi`: test the public API with DBAPI connection
- Tests for the different SQL flavors (flavor specific type conversions)
    - Tests for the sqlalchemy mode: `_TestSQLAlchemy` is the base class with
      common methods, the different tested flavors (sqlite3, MySQL, PostgreSQL)
      derive from the base class
    - Tests for the legacy mode (`TestSQLiteLegacy` and `TestMySQLLegacy`)

"""

from __future__ import print_function
import unittest
import sqlite3
import csv
import os
import sys

import nose
import warnings
import numpy as np

from datetime import datetime

from pandas import DataFrame, Series, Index, MultiIndex, isnull
from pandas import date_range, to_datetime, to_timedelta
import pandas.compat as compat
from pandas.compat import StringIO, range, lrange, string_types
from pandas.core.datetools import format as date_format

import pandas.io.sql as sql
import pandas.util.testing as tm
from pandas import _np_version_under1p7


try:
    import sqlalchemy
    SQLALCHEMY_INSTALLED = True
except ImportError:
    SQLALCHEMY_INSTALLED = False

SQL_STRINGS = {
    'create_iris': {
        'sqlite': """CREATE TABLE iris (
                "SepalLength" REAL,
                "SepalWidth" REAL,
                "PetalLength" REAL,
                "PetalWidth" REAL,
                "Name" TEXT
            )""",
        'mysql': """CREATE TABLE iris (
                `SepalLength` DOUBLE,
                `SepalWidth` DOUBLE,
                `PetalLength` DOUBLE,
                `PetalWidth` DOUBLE,
                `Name` VARCHAR(200)
            )""",
        'postgresql': """CREATE TABLE iris (
                "SepalLength" DOUBLE PRECISION,
                "SepalWidth" DOUBLE PRECISION,
                "PetalLength" DOUBLE PRECISION,
                "PetalWidth" DOUBLE PRECISION,
                "Name" VARCHAR(200)
            )"""
    },
    'insert_iris': {
        'sqlite': """INSERT INTO iris VALUES(?, ?, ?, ?, ?)""",
        'mysql': """INSERT INTO iris VALUES(%s, %s, %s, %s, "%s");""",
        'postgresql': """INSERT INTO iris VALUES(%s, %s, %s, %s, %s);"""
    },
    'create_test_types': {
        'sqlite': """CREATE TABLE types_test_data (
                    "TextCol" TEXT,
                    "DateCol" TEXT,
                    "IntDateCol" INTEGER,
                    "FloatCol" REAL,
                    "IntCol" INTEGER,
                    "BoolCol" INTEGER,
                    "IntColWithNull" INTEGER,
                    "BoolColWithNull" INTEGER
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
                )""",
        'postgresql': """CREATE TABLE types_test_data (
                    "TextCol" TEXT,
                    "DateCol" TIMESTAMP,
                    "IntDateCol" INTEGER,
                    "FloatCol" DOUBLE PRECISION,
                    "IntCol" INTEGER,
                    "BoolCol" BOOLEAN,
                    "IntColWithNull" INTEGER,
                    "BoolColWithNull" BOOLEAN
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
                """,
        'postgresql': """
                INSERT INTO types_test_data
                VALUES(%s, %s, %s, %s, %s, %s, %s, %s)
                """
    },
    'read_parameters': {
        'sqlite': "SELECT * FROM iris WHERE Name=? AND SepalLength=?",
        'mysql': 'SELECT * FROM iris WHERE `Name`="%s" AND `SepalLength`=%s',
        'postgresql': 'SELECT * FROM iris WHERE "Name"=%s AND "SepalLength"=%s'
    },
    'read_named_parameters': {
        'sqlite': """
                SELECT * FROM iris WHERE Name=:name AND SepalLength=:length
                """,
        'mysql': """
                SELECT * FROM iris WHERE
                `Name`="%(name)s" AND `SepalLength`=%(length)s
                """,
        'postgresql': """
                SELECT * FROM iris WHERE
                "Name"=%(name)s AND "SepalLength"=%(length)s
                """
    }
}


class PandasSQLTest(unittest.TestCase):
    """
    Base class with common private methods for SQLAlchemy and fallback cases.

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

    def _read_sql_iris_parameter(self):
        query = SQL_STRINGS['read_parameters'][self.flavor]
        params = ['Iris-setosa', 5.1]
        iris_frame = self.pandasSQL.read_sql(query, params=params)
        self._check_iris_loaded_frame(iris_frame)

    def _read_sql_iris_named_parameter(self):
        query = SQL_STRINGS['read_named_parameters'][self.flavor]
        params = {'name': 'Iris-setosa', 'length': 5.1}
        iris_frame = self.pandasSQL.read_sql(query, params=params)
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

        result.set_index('level_0', inplace=True)
        # result.index.astype(int)

        result.index.name = None

        tm.assert_frame_equal(result, self.test_frame1)

    def _execute_sql(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        iris_results = self.pandasSQL.execute("SELECT * FROM iris")
        row = iris_results.fetchone()
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])


#------------------------------------------------------------------------------
#--- Testing the public API

class _TestSQLApi(PandasSQLTest):

    """
    Base class to test the public API.

    From this two classes are derived to run these tests for both the
    sqlalchemy mode (`TestSQLApi`) and the legacy mode (`TestSQLLegacyApi`).
    These tests are run with sqlite3. Specific tests for the different
    sql flavours are included in `_TestSQLAlchemy`.

    Notes:
    flavor can always be passed even in SQLAlchemy mode,
    should be correctly ignored.

    we don't use drop_table because that isn't part of the public api

    """
    flavor = 'sqlite'

    def setUp(self):
        self.conn = self.connect()
        self._load_iris_data()
        self._load_test1_data()
        self._load_raw_sql()

    def test_read_sql_iris(self):
        iris_frame = sql.read_sql_query(
            "SELECT * FROM iris", self.conn)
        self._check_iris_loaded_frame(iris_frame)

    def test_legacy_read_frame(self):
        with tm.assert_produces_warning(FutureWarning):
            iris_frame = sql.read_frame(
                "SELECT * FROM iris", self.conn)
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
            sql.has_table('test_frame3', self.conn, flavor='sqlite'),
            'Table not written to DB')

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
            sql.has_table('test_frame4', self.conn, flavor='sqlite'),
            'Table not written to DB')

        num_entries = 2 * len(self.test_frame1)
        num_rows = self._count_rows('test_frame4')

        self.assertEqual(
            num_rows, num_entries, "not the same number of rows as entries")

    def test_to_sql_series(self):
        s = Series(np.arange(5, dtype='int64'), name='series')
        sql.to_sql(s, "test_series", self.conn, flavor='sqlite', index=False)
        s2 = sql.read_sql_query("SELECT * FROM test_series", self.conn)
        tm.assert_frame_equal(s.to_frame(), s2)

    def test_to_sql_panel(self):
        panel = tm.makePanel()
        self.assertRaises(NotImplementedError, sql.to_sql, panel,
                          'test_panel', self.conn, flavor='sqlite')

    def test_legacy_write_frame(self):
        # Assume that functionality is already tested above so just do
        # quick check that it basically works
        with tm.assert_produces_warning(FutureWarning):
            sql.write_frame(self.test_frame1, 'test_frame_legacy', self.conn,
                            flavor='sqlite')

        self.assertTrue(
            sql.has_table('test_frame_legacy', self.conn, flavor='sqlite'),
            'Table not written to DB')

    def test_roundtrip(self):
        sql.to_sql(self.test_frame1, 'test_frame_roundtrip',
                   con=self.conn, flavor='sqlite')
        result = sql.read_sql_query(
            'SELECT * FROM test_frame_roundtrip',
            con=self.conn)

        # HACK!
        result.index = self.test_frame1.index
        result.set_index('level_0', inplace=True)
        result.index.astype(int)
        result.index.name = None
        tm.assert_frame_equal(result, self.test_frame1)

    def test_execute_sql(self):
        # drop_sql = "DROP TABLE IF EXISTS test"  # should already be done
        iris_results = sql.execute("SELECT * FROM iris", con=self.conn)
        row = iris_results.fetchone()
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_date_parsing(self):
        # Test date parsing in read_sq
        # No Parsing
        df = sql.read_sql_query("SELECT * FROM types_test_data", self.conn)
        self.assertFalse(
            issubclass(df.DateCol.dtype.type, np.datetime64),
            "DateCol loaded with incorrect type")

        df = sql.read_sql_query("SELECT * FROM types_test_data", self.conn,
                                parse_dates=['DateCol'])
        self.assertTrue(
            issubclass(df.DateCol.dtype.type, np.datetime64),
            "DateCol loaded with incorrect type")

        df = sql.read_sql_query("SELECT * FROM types_test_data", self.conn,
                                parse_dates={'DateCol': '%Y-%m-%d %H:%M:%S'})
        self.assertTrue(
            issubclass(df.DateCol.dtype.type, np.datetime64),
            "DateCol loaded with incorrect type")

        df = sql.read_sql_query("SELECT * FROM types_test_data", self.conn,
                                parse_dates=['IntDateCol'])

        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_sql_query("SELECT * FROM types_test_data", self.conn,
                                parse_dates={'IntDateCol': 's'})

        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

    def test_date_and_index(self):
        # Test case where same column appears in parse_date and index_col

        df = sql.read_sql_query("SELECT * FROM types_test_data", self.conn,
                                index_col='DateCol',
                                parse_dates=['DateCol', 'IntDateCol'])

        self.assertTrue(issubclass(df.index.dtype.type, np.datetime64),
                        "DateCol loaded with incorrect type")

        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

    def test_timedelta(self):
        # see #6921
        if _np_version_under1p7:
            raise nose.SkipTest("test only valid in numpy >= 1.7")

        df = to_timedelta(Series(['00:00:01', '00:00:03'], name='foo')).to_frame()
        with tm.assert_produces_warning(UserWarning):
            df.to_sql('test_timedelta', self.conn)
        result = sql.read_sql_query('SELECT * FROM test_timedelta', self.conn)
        tm.assert_series_equal(result['foo'], df['foo'].astype('int64'))

    def test_to_sql_index_label(self):
        temp_frame = DataFrame({'col1': range(4)})

        # no index name, defaults to 'index'
        sql.to_sql(temp_frame, 'test_index_label', self.conn)
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[0], 'index')

        # specifying index_label
        sql.to_sql(temp_frame, 'test_index_label', self.conn,
                   if_exists='replace', index_label='other_label')
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[0], 'other_label',
                         "Specified index_label not written to database")

        # using the index name
        temp_frame.index.name = 'index_name'
        sql.to_sql(temp_frame, 'test_index_label', self.conn,
                   if_exists='replace')
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[0], 'index_name',
                         "Index name not written to database")

        # has index name, but specifying index_label
        sql.to_sql(temp_frame, 'test_index_label', self.conn,
                   if_exists='replace', index_label='other_label')
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[0], 'other_label',
                         "Specified index_label not written to database")

    def test_to_sql_index_label_multiindex(self):
        temp_frame = DataFrame({'col1': range(4)},
            index=MultiIndex.from_product([('A0', 'A1'), ('B0', 'B1')]))

        # no index name, defaults to 'level_0' and 'level_1'
        sql.to_sql(temp_frame, 'test_index_label', self.conn)
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[0], 'level_0')
        self.assertEqual(frame.columns[1], 'level_1')

        # specifying index_label
        sql.to_sql(temp_frame, 'test_index_label', self.conn,
                   if_exists='replace', index_label=['A', 'B'])
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[:2].tolist(), ['A', 'B'],
                         "Specified index_labels not written to database")

        # using the index name
        temp_frame.index.names = ['A', 'B']
        sql.to_sql(temp_frame, 'test_index_label', self.conn,
                   if_exists='replace')
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[:2].tolist(), ['A', 'B'],
                         "Index names not written to database")

        # has index name, but specifying index_label
        sql.to_sql(temp_frame, 'test_index_label', self.conn,
                   if_exists='replace', index_label=['C', 'D'])
        frame = sql.read_sql_query('SELECT * FROM test_index_label', self.conn)
        self.assertEqual(frame.columns[:2].tolist(), ['C', 'D'],
                         "Specified index_labels not written to database")

        # wrong length of index_label
        self.assertRaises(ValueError, sql.to_sql, temp_frame,
                          'test_index_label', self.conn, if_exists='replace',
                          index_label='C')

    def test_integer_col_names(self):
        df = DataFrame([[1, 2], [3, 4]], columns=[0, 1])
        sql.to_sql(df, "test_frame_integer_col_names", self.conn,
                   if_exists='replace')

    def test_get_schema(self):
        create_sql = sql.get_schema(self.test_frame1, 'test', 'sqlite',
                                    con=self.conn)
        self.assertTrue('CREATE' in create_sql)


class TestSQLApi(_TestSQLApi):
    """
    Test the public API as it would be used directly

    Tests for `read_sql_table` are included here, as this is specific for the
    sqlalchemy mode.

    """
    flavor = 'sqlite'

    def connect(self):
        if SQLALCHEMY_INSTALLED:
            return sqlalchemy.create_engine('sqlite:///:memory:')
        else:
            raise nose.SkipTest('SQLAlchemy not installed')

    def test_read_table_columns(self):
        # test columns argument in read_table
        sql.to_sql(self.test_frame1, 'test_frame', self.conn)

        cols = ['A', 'B']
        result = sql.read_sql_table('test_frame', self.conn, columns=cols)
        self.assertEqual(result.columns.tolist(), cols,
                         "Columns not correctly selected")

    def test_read_table_index_col(self):
        # test columns argument in read_table
        sql.to_sql(self.test_frame1, 'test_frame', self.conn)

        result = sql.read_sql_table('test_frame', self.conn, index_col="index")
        self.assertEqual(result.index.names, ["index"],
                         "index_col not correctly set")

        result = sql.read_sql_table('test_frame', self.conn, index_col=["A", "B"])
        self.assertEqual(result.index.names, ["A", "B"],
                         "index_col not correctly set")

        result = sql.read_sql_table('test_frame', self.conn, index_col=["A", "B"],
                                columns=["C", "D"])
        self.assertEqual(result.index.names, ["A", "B"],
                         "index_col not correctly set")
        self.assertEqual(result.columns.tolist(), ["C", "D"],
                         "columns not set correctly whith index_col")

    def test_read_sql_delegate(self):
        iris_frame1 = sql.read_sql_query(
            "SELECT * FROM iris", self.conn)
        iris_frame2 = sql.read_sql(
            "SELECT * FROM iris", self.conn)
        tm.assert_frame_equal(iris_frame1, iris_frame2,
                              "read_sql and read_sql_query have not the same"
                              " result with a query")

        iris_frame1 = sql.read_sql_table('iris', self.conn)
        iris_frame2 = sql.read_sql('iris', self.conn)
        tm.assert_frame_equal(iris_frame1, iris_frame2)


class TestSQLLegacyApi(_TestSQLApi):
    """
    Test the public legacy API

    """
    flavor = 'sqlite'

    def connect(self, database=":memory:"):
        return sqlite3.connect(database)

    def _load_test2_data(self):
        columns = ['index', 'A', 'B']
        data = [(
            '2000-01-03 00:00:00', 2 ** 31 - 1, -1.987670),
            ('2000-01-04 00:00:00', -29, -0.0412318367011),
            ('2000-01-05 00:00:00', 20000, 0.731167677815),
            ('2000-01-06 00:00:00', -290867, 1.56762092543)]

        self.test_frame2 = DataFrame(data, columns=columns)

    def test_sql_open_close(self):
        # Test if the IO in the database still work if the connection closed
        # between the writing and reading (as in many real situations).

        self._load_test2_data()

        with tm.ensure_clean() as name:

            conn = self.connect(name)
            sql.to_sql(self.test_frame2, "test_frame2_legacy", conn,
                       flavor="sqlite", index=False)
            conn.close()

            conn = self.connect(name)
            result = sql.read_sql_query("SELECT * FROM test_frame2_legacy;",
                                        conn)
            conn.close()

        tm.assert_frame_equal(self.test_frame2, result)

    def test_read_sql_delegate(self):
        iris_frame1 = sql.read_sql_query("SELECT * FROM iris", self.conn)
        iris_frame2 = sql.read_sql("SELECT * FROM iris", self.conn)
        tm.assert_frame_equal(iris_frame1, iris_frame2,
                              "read_sql and read_sql_query have not the same"
                              " result with a query")

        self.assertRaises(ValueError, sql.read_sql, 'iris', self.conn)

    def test_safe_names_warning(self):
        # GH 6798
        df = DataFrame([[1, 2], [3, 4]], columns=['a', 'b '])  # has a space
        # warns on create table with spaces in names
        with tm.assert_produces_warning():
            sql.to_sql(df, "test_frame3_legacy", self.conn,
                       flavor="sqlite", index=False)

    def test_get_schema2(self):
        # without providing a connection object (available for backwards comp)
        create_sql = sql.get_schema(self.test_frame1, 'test', 'sqlite')
        self.assertTrue('CREATE' in create_sql)

    def test_tquery(self):
        with tm.assert_produces_warning(FutureWarning):
            iris_results = sql.tquery("SELECT * FROM iris", con=self.conn)
        row = iris_results[0]
        tm.equalContents(row, [5.1, 3.5, 1.4, 0.2, 'Iris-setosa'])

    def test_uquery(self):
        with tm.assert_produces_warning(FutureWarning):
            rows = sql.uquery("SELECT * FROM iris LIMIT 1", con=self.conn)
        self.assertEqual(rows, -1)


#------------------------------------------------------------------------------
#--- Database flavor specific tests


class _TestSQLAlchemy(PandasSQLTest):
    """
    Base class for testing the sqlalchemy backend.

    Subclasses for specific database types are created below. Tests that
    deviate for each flavor are overwritten there.

    """
    flavor = None

    def setUp(self):
        self.setup_import()
        self.setup_driver()
        self.setup_connect()

        self._load_iris_data()
        self._load_raw_sql()
        self._load_test1_data()

    def setup_import(self):
        # Skip this test if SQLAlchemy not available
        if not SQLALCHEMY_INSTALLED:
            raise nose.SkipTest('SQLAlchemy not installed')

    def setup_driver(self):
        raise NotImplementedError()

    def connect(self):
        raise NotImplementedError()

    def setup_connect(self):
        try:
            self.conn = self.connect()
            self.pandasSQL = sql.PandasSQLAlchemy(self.conn)
        except sqlalchemy.exc.OperationalError:
            raise nose.SkipTest("Can't connect to {0} server".format(self.flavor))

    def tearDown(self):
        raise NotImplementedError()

    def test_aread_sql(self):
        self._read_sql_iris()

    def test_read_sql_parameter(self):
        self._read_sql_iris_parameter()

    def test_read_sql_named_parameter(self):
        self._read_sql_iris_named_parameter()

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
        iris_frame = sql.read_sql_table("iris", con=self.conn)
        self._check_iris_loaded_frame(iris_frame)

    def test_read_table_columns(self):
        iris_frame = sql.read_sql_table(
            "iris", con=self.conn, columns=['SepalLength', 'SepalLength'])
        tm.equalContents(
            iris_frame.columns.values, ['SepalLength', 'SepalLength'])

    def test_read_table_absent(self):
        self.assertRaises(
            ValueError, sql.read_sql_table, "this_doesnt_exist", con=self.conn)

    def test_default_type_conversion(self):
        df = sql.read_sql_table("types_test_data", self.conn)

        self.assertTrue(issubclass(df.FloatCol.dtype.type, np.floating),
                        "FloatCol loaded with incorrect type")
        self.assertTrue(issubclass(df.IntCol.dtype.type, np.integer),
                        "IntCol loaded with incorrect type")
        self.assertTrue(issubclass(df.BoolCol.dtype.type, np.bool_),
                        "BoolCol loaded with incorrect type")

        # Int column with NA values stays as float
        self.assertTrue(issubclass(df.IntColWithNull.dtype.type, np.floating),
                        "IntColWithNull loaded with incorrect type")
        # Bool column with NA values becomes object
        self.assertTrue(issubclass(df.BoolColWithNull.dtype.type, np.object),
                        "BoolColWithNull loaded with incorrect type")

    def test_default_date_load(self):
        df = sql.read_sql_table("types_test_data", self.conn)

        # IMPORTANT - sqlite has no native date type, so shouldn't parse, but
        # MySQL SHOULD be converted.
        self.assertTrue(issubclass(df.DateCol.dtype.type, np.datetime64),
                        "DateCol loaded with incorrect type")

    def test_date_parsing(self):
        # No Parsing
        df = sql.read_sql_table("types_test_data", self.conn)

        df = sql.read_sql_table("types_test_data", self.conn,
                                parse_dates=['DateCol'])
        self.assertTrue(issubclass(df.DateCol.dtype.type, np.datetime64),
                        "DateCol loaded with incorrect type")

        df = sql.read_sql_table("types_test_data", self.conn,
                                parse_dates={'DateCol': '%Y-%m-%d %H:%M:%S'})
        self.assertTrue(issubclass(df.DateCol.dtype.type, np.datetime64),
                        "DateCol loaded with incorrect type")

        df = sql.read_sql_table("types_test_data", self.conn, parse_dates={
                            'DateCol': {'format': '%Y-%m-%d %H:%M:%S'}})
        self.assertTrue(issubclass(df.DateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_sql_table(
            "types_test_data", self.conn, parse_dates=['IntDateCol'])
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_sql_table(
            "types_test_data", self.conn, parse_dates={'IntDateCol': 's'})
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

        df = sql.read_sql_table(
            "types_test_data", self.conn, parse_dates={'IntDateCol': {'unit': 's'}})
        self.assertTrue(issubclass(df.IntDateCol.dtype.type, np.datetime64),
                        "IntDateCol loaded with incorrect type")

    def test_datetime(self):
        if self.driver == 'pymysql':
             raise nose.SkipTest('writing datetime not working with pymysql')

        df = DataFrame({'A': date_range('2013-01-01 09:00:00', periods=3),
                        'B': np.arange(3.0)})
        df.to_sql('test_datetime', self.conn)

        # with read_table -> type information from schema used
        result = sql.read_sql_table('test_datetime', self.conn)
        result = result.drop('index', axis=1)
        tm.assert_frame_equal(result, df)

        # with read_sql -> no type information -> sqlite has no native
        result = sql.read_sql_query('SELECT * FROM test_datetime', self.conn)
        result = result.drop('index', axis=1)
        if self.flavor == 'sqlite':
            self.assertTrue(isinstance(result.loc[0, 'A'], string_types))
            result['A'] = to_datetime(result['A'])
            tm.assert_frame_equal(result, df)
        else:
            tm.assert_frame_equal(result, df)

    def test_datetime_NaT(self):
        # status:
        # - postgresql: gives error on inserting "0001-255-255T00:00:00"
        # - sqlite3: works, but reading it with query returns '-001--1--1 -1:-1:-1.-00001'

        if self.driver == 'pymysql':
            raise nose.SkipTest('writing datetime not working with pymysql')
        if self.driver == 'psycopg2':
            raise nose.SkipTest('writing datetime NaT not working with psycopg2')
        if self.flavor == 'sqlite':
            raise nose.SkipTest('reading datetime NaT not working with sqlite')

        df = DataFrame({'A': date_range('2013-01-01 09:00:00', periods=3),
                        'B': np.arange(3.0)})
        df.loc[1, 'A'] = np.nan
        df.to_sql('test_datetime', self.conn, index=False)

        # with read_table -> type information from schema used
        result = sql.read_sql_table('test_datetime', self.conn)
        tm.assert_frame_equal(result, df)

        # with read_sql -> no type information -> sqlite has no native
        result = sql.read_sql_query('SELECT * FROM test_datetime', self.conn)
        if self.flavor == 'sqlite':
            self.assertTrue(isinstance(result.loc[0, 'A'], string_types))
            result['A'] = to_datetime(result['A'], coerce=True)
            tm.assert_frame_equal(result, df)
        else:
            tm.assert_frame_equal(result, df)

    def test_mixed_dtype_insert(self):
        # see GH6509
        s1 = Series(2**25 + 1,dtype=np.int32)
        s2 = Series(0.0,dtype=np.float32)
        df = DataFrame({'s1': s1, 's2': s2})

        # write and read again
        df.to_sql("test_read_write", self.conn, index=False)
        df2 = sql.read_sql_table("test_read_write", self.conn)

        tm.assert_frame_equal(df, df2, check_dtype=False, check_exact=True)

    def test_nan_numeric(self):
        if self.driver == 'pymysql':
            raise nose.SkipTest('writing NaNs not working with pymysql')

        # NaNs in numeric float column
        df = DataFrame({'A':[0, 1, 2], 'B':[0.2, np.nan, 5.6]})
        df.to_sql('test_nan', self.conn, index=False)

        # with read_table
        result = sql.read_sql_table('test_nan', self.conn)
        tm.assert_frame_equal(result, df)

        # with read_sql
        result = sql.read_sql_query('SELECT * FROM test_nan', self.conn)
        tm.assert_frame_equal(result, df)

    def test_nan_fullcolumn(self):
        if self.driver == 'pymysql':
            raise nose.SkipTest('writing NaNs not working with pymysql')

        # full NaN column (numeric float column)
        df = DataFrame({'A':[0, 1, 2], 'B':[np.nan, np.nan, np.nan]})
        df.to_sql('test_nan', self.conn, index=False)

        if self.flavor == 'sqlite':
            df['B'] = df['B'].astype('object')
            df['B'] = None

        # with read_table
        result = sql.read_sql_table('test_nan', self.conn)
        tm.assert_frame_equal(result, df)

        # with read_sql
        result = sql.read_sql_query('SELECT * FROM test_nan', self.conn)
        tm.assert_frame_equal(result, df)

    def test_nan_string(self):
        if self.driver == 'pymysql':
             raise nose.SkipTest('writing NaNs not working with pymysql')

        # NaNs in string column
        df = DataFrame({'A':[0, 1, 2], 'B':['a', 'b', np.nan]})
        df.to_sql('test_nan', self.conn, index=False)

        if self.flavor == 'sqlite':
            df.loc[2, 'B'] = None
        elif self.flavor == 'postgresql':
            df = df.fillna('NaN')

        # with read_table
        result = sql.read_sql_table('test_nan', self.conn)
        tm.assert_frame_equal(result, df)

        # with read_sql
        result = sql.read_sql_query('SELECT * FROM test_nan', self.conn)
        tm.assert_frame_equal(result, df)


class TestSQLiteAlchemy(_TestSQLAlchemy):
    """
    Test the sqlalchemy backend against an in-memory sqlite database.

    """
    flavor = 'sqlite'

    def connect(self):
        return sqlalchemy.create_engine('sqlite:///:memory:')

    def setup_driver(self):
        # sqlite3 is built-in
        self.driver = None

    def tearDown(self):
        # in memory so tables should not be removed explicitly
        pass

    def test_default_type_conversion(self):
        df = sql.read_sql_table("types_test_data", self.conn)

        self.assertTrue(issubclass(df.FloatCol.dtype.type, np.floating),
                        "FloatCol loaded with incorrect type")
        self.assertTrue(issubclass(df.IntCol.dtype.type, np.integer),
                        "IntCol loaded with incorrect type")
        # sqlite has no boolean type, so integer type is returned
        self.assertTrue(issubclass(df.BoolCol.dtype.type, np.integer),
                        "BoolCol loaded with incorrect type")

        # Int column with NA values stays as float
        self.assertTrue(issubclass(df.IntColWithNull.dtype.type, np.floating),
                        "IntColWithNull loaded with incorrect type")
        # Non-native Bool column with NA values stays as float
        self.assertTrue(issubclass(df.BoolColWithNull.dtype.type, np.floating),
                        "BoolColWithNull loaded with incorrect type")

    def test_default_date_load(self):
        df = sql.read_sql_table("types_test_data", self.conn)

        # IMPORTANT - sqlite has no native date type, so shouldn't parse, but
        self.assertFalse(issubclass(df.DateCol.dtype.type, np.datetime64),
                         "DateCol loaded with incorrect type")


class TestMySQLAlchemy(_TestSQLAlchemy):
    """
    Test the sqlalchemy backend against an MySQL database.

    """
    flavor = 'mysql'

    def connect(self):
        return sqlalchemy.create_engine(
            'mysql+{driver}://root@localhost/pandas_nosetest'.format(driver=self.driver))

    def setup_driver(self):
        try:
            import pymysql
            self.driver = 'pymysql'
        except ImportError:
            raise nose.SkipTest('pymysql not installed')

    def tearDown(self):
        c = self.conn.execute('SHOW TABLES')
        for table in c.fetchall():
            self.conn.execute('DROP TABLE %s' % table[0])

    def test_default_type_conversion(self):
        df = sql.read_sql_table("types_test_data", self.conn)

        self.assertTrue(issubclass(df.FloatCol.dtype.type, np.floating),
                        "FloatCol loaded with incorrect type")
        self.assertTrue(issubclass(df.IntCol.dtype.type, np.integer),
                        "IntCol loaded with incorrect type")
        # MySQL has no real BOOL type (it's an alias for TINYINT)
        self.assertTrue(issubclass(df.BoolCol.dtype.type, np.integer),
                        "BoolCol loaded with incorrect type")

        # Int column with NA values stays as float
        self.assertTrue(issubclass(df.IntColWithNull.dtype.type, np.floating),
                        "IntColWithNull loaded with incorrect type")
        # Bool column with NA = int column with NA values => becomes float
        self.assertTrue(issubclass(df.BoolColWithNull.dtype.type, np.floating),
                        "BoolColWithNull loaded with incorrect type")


class TestPostgreSQLAlchemy(_TestSQLAlchemy):
    """
    Test the sqlalchemy backend against an PostgreSQL database.

    """
    flavor = 'postgresql'

    def connect(self):
        return sqlalchemy.create_engine(
            'postgresql+{driver}://postgres@localhost/pandas_nosetest'.format(driver=self.driver))

    def setup_driver(self):
        try:
            import psycopg2
            self.driver = 'psycopg2'
        except ImportError:
            raise nose.SkipTest('psycopg2 not installed')

    def tearDown(self):
        c = self.conn.execute(
            "SELECT table_name FROM information_schema.tables"
            " WHERE table_schema = 'public'")
        for table in c.fetchall():
            self.conn.execute("DROP TABLE %s" % table[0])


#------------------------------------------------------------------------------
#--- Test Sqlite / MySQL fallback

class TestSQLiteLegacy(PandasSQLTest):
    """
    Test the legacy mode against an in-memory sqlite database.

    """
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

    def test_read_sql_parameter(self):
        self._read_sql_iris_parameter()

    def test_read_sql_named_parameter(self):
        self._read_sql_iris_named_parameter()

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

        self.assertTrue(self.pandasSQL.has_table('drop_test_frame'),
                        'Table not written to DB')

        self.pandasSQL.drop_table('drop_test_frame')

        self.assertFalse(self.pandasSQL.has_table('drop_test_frame'),
                         'Table not deleted from DB')

    def test_roundtrip(self):
        self._roundtrip()

    def test_execute_sql(self):
        self._execute_sql()


class TestMySQLLegacy(TestSQLiteLegacy):
    """
    Test the legacy mode against a MySQL database.

    """
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
            raise nose.SkipTest('pymysql not installed')

        try:
            self.conn = self.connect()
        except self.driver.err.OperationalError:
            raise nose.SkipTest("Can't connect to MySQL server")

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

    def test_a_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            sql.to_sql(self.test_frame1, 'test_frame1', self.conn,
                       flavor='mysql')
        self.assertTrue(
            sql.has_table('test_frame1', self.conn, flavor='mysql'),
            'Table not written to DB')


#------------------------------------------------------------------------------
#--- Old tests from 0.13.1 (before refactor using sqlalchemy)


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

def _skip_if_no_pymysql():
    try:
        import pymysql
    except ImportError:
        raise nose.SkipTest('pymysql not installed, skipping')


class TestXSQLite(tm.TestCase):

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
                self.assertTrue(tokens[1] == 'DATETIME')

        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_schema(frame, 'test', 'sqlite', keys=['A', 'B'],)
        lines = create_sql.splitlines()
        self.assertTrue('PRIMARY KEY (A,B)' in create_sql)
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
            self.assertRaises(sql.DatabaseError, sql.tquery,
                              'select * from blah', con=self.db)

            self.assertRaises(sql.DatabaseError, sql.tquery,
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

            self.assertRaises(sql.DatabaseError, sql.tquery,
                              'insert into blah values (1)', con=self.db)

            self.assertRaises(sql.DatabaseError, sql.tquery,
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


class TestXMySQL(tm.TestCase):

    def setUp(self):
        _skip_if_no_pymysql()
        import pymysql
        try:
            # Try Travis defaults.
            # No real user should allow root access with a blank password.
            self.db = pymysql.connect(host='localhost', user='root', passwd='',
                                    db='pandas_nosetest')
        except:
            pass
        else:
            return
        try:
            self.db = pymysql.connect(read_default_group='pandas')
        except pymysql.ProgrammingError as e:
            raise nose.SkipTest(
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")
        except pymysql.Error as e:
            raise nose.SkipTest(
                "Cannot connect to database. "
                "Create a group of connection parameters under the heading "
                "[pandas] in your system's mysql default file, "
                "typically located at ~/.my.cnf or /etc/.my.cnf. ")

    def test_basic(self):
        _skip_if_no_pymysql()
        frame = tm.makeTimeDataFrame()
        self._check_roundtrip(frame)

    def test_write_row_by_row(self):

        _skip_if_no_pymysql()
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
        _skip_if_no_pymysql()
        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = sql.get_schema(frame, 'test', 'mysql')
        cur = self.db.cursor()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Unknown table.*")
            cur.execute(drop_sql)
        cur.execute(create_sql)
        ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"

        row = frame.ix[0].values.tolist()
        sql.execute(ins, self.db, params=tuple(row))
        self.db.commit()

        result = sql.read_frame("select * from test", self.db)
        result.index = frame.index[:1]
        tm.assert_frame_equal(result, frame[:1])

    def test_schema(self):
        _skip_if_no_pymysql()
        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_schema(frame, 'test', 'mysql')
        lines = create_sql.splitlines()
        for l in lines:
            tokens = l.split(' ')
            if len(tokens) == 2 and tokens[0] == 'A':
                self.assertTrue(tokens[1] == 'DATETIME')

        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test"
        create_sql = sql.get_schema(frame, 'test', 'mysql', keys=['A', 'B'],)
        lines = create_sql.splitlines()
        self.assertTrue('PRIMARY KEY (A,B)' in create_sql)
        cur = self.db.cursor()
        cur.execute(drop_sql)
        cur.execute(create_sql)

    def test_execute_fail(self):
        _skip_if_no_pymysql()
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
        _skip_if_no_pymysql()
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
        _skip_if_no_pymysql()
        pass

    def _check_roundtrip(self, frame):
        _skip_if_no_pymysql()
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
            import pymysql
        except ImportError:
            raise nose.SkipTest("no pymysql")
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
            self.assertRaises(sql.DatabaseError, sql.tquery,
                              'select * from blah', con=self.db)

            self.assertRaises(sql.DatabaseError, sql.tquery,
                              'select * from blah', con=self.db, retry=True)
        finally:
            sys.stdout = sys.__stdout__

    def test_uquery(self):
        try:
            import pymysql
        except ImportError:
            raise nose.SkipTest("no pymysql")
        frame = tm.makeTimeDataFrame()
        drop_sql = "DROP TABLE IF EXISTS test_table"
        cur = self.db.cursor()
        cur.execute(drop_sql)
        sql.write_frame(frame, name='test_table', con=self.db, flavor='mysql')
        stmt = 'INSERT INTO test_table VALUES(2.314, -123.1, 1.234, 2.3)'
        self.assertEqual(sql.uquery(stmt, con=self.db), 1)

        try:
            sys.stdout = StringIO()

            self.assertRaises(sql.DatabaseError, sql.tquery,
                              'insert into blah values (1)', con=self.db)

            self.assertRaises(sql.DatabaseError, sql.tquery,
                              'insert into blah values (1)', con=self.db,
                              retry=True)
        finally:
            sys.stdout = sys.__stdout__

    def test_keyword_as_column_names(self):
        '''
        '''
        _skip_if_no_pymysql()
        df = DataFrame({'From':np.ones(5)})
        sql.write_frame(df, con = self.db, name = 'testkeywords',
                        if_exists='replace', flavor='mysql')

    def test_if_exists(self):
        _skip_if_no_pymysql()
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
