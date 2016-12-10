import sqlalchemy
from .pandas_vb_common import *
import sqlite3
from sqlalchemy import create_engine


#-------------------------------------------------------------------------------
# to_sql

class WriteSQL(object):
    goal_time = 0.2

    def setup(self):
        self.engine = create_engine('sqlite:///:memory:')
        self.con = sqlite3.connect(':memory:')
        self.index = tm.makeStringIndex(10000)
        self.df = DataFrame({'float1': randn(10000), 'float2': randn(10000), 'string1': (['foo'] * 10000), 'bool1': ([True] * 10000), 'int1': np.random.randint(0, 100000, size=10000), }, index=self.index)

    def time_fallback(self):
        self.df.to_sql('test1', self.con, if_exists='replace')

    def time_sqlalchemy(self):
        self.df.to_sql('test1', self.engine, if_exists='replace')


#-------------------------------------------------------------------------------
# read_sql

class ReadSQL(object):
    goal_time = 0.2

    def setup(self):
        self.engine = create_engine('sqlite:///:memory:')
        self.con = sqlite3.connect(':memory:')
        self.index = tm.makeStringIndex(10000)
        self.df = DataFrame({'float1': randn(10000), 'float2': randn(10000), 'string1': (['foo'] * 10000), 'bool1': ([True] * 10000), 'int1': np.random.randint(0, 100000, size=10000), }, index=self.index)
        self.df.to_sql('test2', self.engine, if_exists='replace')
        self.df.to_sql('test2', self.con, if_exists='replace')

    def time_read_query_fallback(self):
        read_sql_query('SELECT * FROM test2', self.con)

    def time_read_query_sqlalchemy(self):
        read_sql_query('SELECT * FROM test2', self.engine)

    def time_read_table_sqlalchemy(self):
        read_sql_table('test2', self.engine)


#-------------------------------------------------------------------------------
# type specific write

class WriteSQLTypes(object):
    goal_time = 0.2

    def setup(self):
        self.engine = create_engine('sqlite:///:memory:')
        self.con = sqlite3.connect(':memory:')
        self.df = DataFrame({'float': randn(10000), 'string': (['foo'] * 10000), 'bool': ([True] * 10000), 'datetime': date_range('2000-01-01', periods=10000, freq='s'), })
        self.df.loc[1000:3000, 'float'] = np.nan

    def time_string_fallback(self):
        self.df[['string']].to_sql('test_string', self.con, if_exists='replace')

    def time_string_sqlalchemy(self):
        self.df[['string']].to_sql('test_string', self.engine, if_exists='replace')

    def time_float_fallback(self):
        self.df[['float']].to_sql('test_float', self.con, if_exists='replace')

    def time_float_sqlalchemy(self):
        self.df[['float']].to_sql('test_float', self.engine, if_exists='replace')

    def time_datetime_sqlalchemy(self):
        self.df[['datetime']].to_sql('test_datetime', self.engine, if_exists='replace')


#-------------------------------------------------------------------------------
# type specific read

class ReadSQLTypes(object):
    goal_time = 0.2

    def setup(self):
        self.engine = create_engine('sqlite:///:memory:')
        self.con = sqlite3.connect(':memory:')
        self.df = DataFrame({'float': randn(10000), 'datetime': date_range('2000-01-01', periods=10000, freq='s'), })
        self.df['datetime_string'] = self.df['datetime'].map(str)
        self.df.to_sql('test_type', self.engine, if_exists='replace')
        self.df[['float', 'datetime_string']].to_sql('test_type', self.con, if_exists='replace')

    def time_datetime_read_and_parse_sqlalchemy(self):
        read_sql_table('test_type', self.engine, columns=['datetime_string'], parse_dates=['datetime_string'])

    def time_datetime_read_as_native_sqlalchemy(self):
        read_sql_table('test_type', self.engine, columns=['datetime'])

    def time_float_read_query_fallback(self):
        read_sql_query('SELECT float FROM test_type', self.con)

    def time_float_read_query_sqlalchemy(self):
        read_sql_query('SELECT float FROM test_type', self.engine)

    def time_float_read_table_sqlalchemy(self):
        read_sql_table('test_type', self.engine, columns=['float'])
