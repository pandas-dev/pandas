import numpy as np
from collections import defaultdict
import gc
import time
from pandas import DataFrame
from pandas.util.testing import rands
from pandas.compat import range, zip
import random

N = 10000

indices = np.array([rands(10) for _ in range(N)], dtype='O')
indices2 = np.array([rands(10) for _ in range(N)], dtype='O')
key = np.tile(indices[:8000], 10)
key2 = np.tile(indices2[:8000], 10)

left = DataFrame({'key': key, 'key2': key2,
                  'value': np.random.randn(80000)})
right = DataFrame({'key': indices[2000:], 'key2': indices2[2000:],
                   'value2': np.random.randn(8000)})

# right2 = right.append(right, ignore_index=True)
# right = right2

# random.shuffle(key2)
# indices2 = indices.copy()
# random.shuffle(indices2)

# Prepare Database
import sqlite3
create_sql_indexes = True

conn = sqlite3.connect(':memory:')
conn.execute(
    'create table left( key varchar(10), key2 varchar(10), value int);')
conn.execute(
    'create table right( key varchar(10), key2 varchar(10), value2 int);')
conn.executemany('insert into left values (?, ?, ?)',
                 zip(key, key2, left['value']))
conn.executemany('insert into right values (?, ?, ?)',
                 zip(right['key'], right['key2'], right['value2']))

# Create Indices
if create_sql_indexes:
    conn.execute('create index left_ix on left(key, key2)')
    conn.execute('create index right_ix on right(key, key2)')


join_methods = ['inner', 'left outer', 'left']  # others not supported
sql_results = DataFrame(index=join_methods, columns=[False])
niter = 5
for sort in [False]:
    for join_method in join_methods:
        sql = """CREATE TABLE test as select *
        from left
           %s join right
             on left.key=right.key
               and left.key2 = right.key2;""" % join_method
        sql = """select *
        from left
           %s join right
             on left.key=right.key
               and left.key2 = right.key2;""" % join_method

        if sort:
            sql = '%s order by key, key2' % sql
        f = lambda: list(conn.execute(sql))  # list fetches results
        g = lambda: conn.execute(sql)  # list fetches results
        gc.disable()
        start = time.time()
        # for _ in range(niter):
        g()
        elapsed = (time.time() - start) / niter
        gc.enable()

        cur = conn.execute("DROP TABLE test")
        conn.commit()

        sql_results[sort][join_method] = elapsed
        sql_results.columns = ['sqlite3']  # ['dont_sort', 'sort']
        sql_results.index = ['inner', 'outer', 'left']

        sql = """select *
        from left
           inner join right
             on left.key=right.key
               and left.key2 = right.key2;"""
