from pandas import *
from pandas.util.testing import rands
import random

N = 10000
ngroups = 10

def get_test_data(ngroups=100, n=N):
    unique_groups = range(ngroups)
    arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
                         dtype=object)

    random.shuffle(arr)
    return arr

# aggregate multiple columns
# df = DataFrame({'key1' : get_test_data(ngroups=ngroups),
#                 'key2' : get_test_data(ngroups=ngroups),
#                 'data1' : np.random.randn(N),
#                 'data2' : np.random.randn(N)})

# df2 = DataFrame({'key1'  : get_test_data(ngroups=ngroups, n=N//10),
#                  'key2'  : get_test_data(ngroups=ngroups//2, n=N//10),
#                  'value' : np.random.randn(N // 10)})
# result = merge.merge(df, df2, on='key2')

from collections import defaultdict
import gc
import time
from pandas.util.testing import rands
N = 10000
indices = np.array([rands(10) for _ in xrange(N)], dtype='O')

key = np.tile(indices, 10)
key2 = key.copy()
random.shuffle(key2)
indices2 = indices.copy()
random.shuffle(indices2)
left = DataFrame({'key' : key, 'key2':key2,
                  'value' : np.random.randn(100000)})
right = DataFrame({'key': indices, 'key2':indices2,
                   'value2' : np.random.randn(10000)})
join_methods = ['inner', 'outer', 'left', 'right']
results = DataFrame(index=join_methods, columns=[False, True])
niter = 10
for sort in [False, True]:
    for join_method in join_methods:
        f = lambda: merge(left, right, how=join_method, sort=sort)
        gc.disable()
        start = time.time()
        for _ in xrange(niter):
            f()
        elapsed = (time.time() - start) / niter
        gc.enable()
        results[sort][join_method] = elapsed
results.columns = ['dont_sort', 'sort']


# R results
from StringIO import StringIO
r_results = read_table(StringIO("""dont_sort   sort
inner    0.2297 0.2286
outer    1.1811 1.2843
left     0.6706 0.7766
right    0.2995 0.3371
"""), sep='\s+')

sort_results = DataFrame.from_items([('pandas', results['sort']),
                                     ('R', r_results['sort'])])
sort_results['Ratio'] = sort_results['R'] / sort_results['pandas']


nosort_results = DataFrame.from_items([('pandas', results['dont_sort']),
                                       ('R', r_results['dont_sort'])])
nosort_results['Ratio'] = sort_results['R'] / sort_results['pandas']

