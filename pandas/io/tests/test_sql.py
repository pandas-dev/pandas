from cStringIO import StringIO
import unittest
import sqlite3
import sys

import pandas.io.sql as sql
import pandas.util.testing as tm

class TestSQLite(unittest.TestCase):

    def setUp(self):
        self.db = sqlite3.connect(':memory:')

    def test_basic(self):
        frame = tm.makeTimeDataFrame()
        self._check_roundtrip(frame)

    def test_write_row_by_row(self):
        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_sqlite_schema(frame, 'test')
        self.db.execute(create_sql)

        cur = self.db.cursor()

        ins = "INSERT INTO test VALUES (%s, %s, %s, %s)"
        for idx, row in frame.iterrows():
            fmt_sql = sql.format_query(ins, *row)
            sql.tquery(fmt_sql, cur=cur)

        self.db.commit()

        result = sql.read_frame("select * from test", con=self.db)
        result.index = frame.index
        tm.assert_frame_equal(result, frame)

    def test_execute(self):
        frame = tm.makeTimeDataFrame()
        create_sql = sql.get_sqlite_schema(frame, 'test')
        self.db.execute(create_sql)
        ins = "INSERT INTO test VALUES (?, ?, ?, ?)"

        row = frame.ix[0]
        sql.execute(ins, self.db, params=tuple(row))
        self.db.commit()

        result = sql.read_frame("select * from test", self.db)
        result.index = frame.index[:1]
        tm.assert_frame_equal(result, frame[:1])



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
        self.db.execute(create_sql)

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
        self.db.execute(create_sql)

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

        # HACK!
        result.index = frame.index

        expected = frame
        tm.assert_frame_equal(result, expected)

if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

