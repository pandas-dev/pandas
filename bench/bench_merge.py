import random
import gc
import time
from pandas import *
from pandas.compat import range, lrange, StringIO
from pandas.util.testing import rands

N = 10000
ngroups = 10


def get_test_data(ngroups=100, n=N):
    unique_groups = lrange(ngroups)
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

N = 10000

indices = np.array([rands(10) for _ in range(N)], dtype='O')
indices2 = np.array([rands(10) for _ in range(N)], dtype='O')
key = np.tile(indices[:8000], 10)
key2 = np.tile(indices2[:8000], 10)

left = DataFrame({'key': key, 'key2': key2,
                  'value': np.random.randn(80000)})
right = DataFrame({'key': indices[2000:], 'key2': indices2[2000:],
                   'value2': np.random.randn(8000)})

right2 = right.append(right, ignore_index=True)


join_methods = ['inner', 'outer', 'left', 'right']
results = DataFrame(index=join_methods, columns=[False, True])
niter = 10
for sort in [False, True]:
    for join_method in join_methods:
        f = lambda: merge(left, right, how=join_method, sort=sort)
        gc.disable()
        start = time.time()
        for _ in range(niter):
            f()
        elapsed = (time.time() - start) / niter
        gc.enable()
        results[sort][join_method] = elapsed
# results.columns = ['pandas']
results.columns = ['dont_sort', 'sort']


# R results
# many to one
r_results = read_table(StringIO("""      base::merge   plyr data.table
inner      0.2475 0.1183     0.1100
outer      0.4213 0.1916     0.2090
left       0.2998 0.1188     0.0572
right      0.3102 0.0536     0.0376
"""), sep='\s+')

presults = results[['dont_sort']].rename(columns={'dont_sort': 'pandas'})
all_results = presults.join(r_results)

all_results = all_results.div(all_results['pandas'], axis=0)

all_results = all_results.ix[:, ['pandas', 'data.table', 'plyr',
                                 'base::merge']]

sort_results = DataFrame.from_items([('pandas', results['sort']),
                                     ('R', r_results['base::merge'])])
sort_results['Ratio'] = sort_results['R'] / sort_results['pandas']


nosort_results = DataFrame.from_items([('pandas', results['dont_sort']),
                                       ('R', r_results['base::merge'])])
nosort_results['Ratio'] = nosort_results['R'] / nosort_results['pandas']

# many to many

# many to one
r_results = read_table(StringIO("""base::merge   plyr data.table
inner      0.4610 0.1276     0.1269
outer      0.9195 0.1881     0.2725
left       0.6559 0.1257     0.0678
right      0.6425 0.0522     0.0428
"""), sep='\s+')

all_results = presults.join(r_results)
all_results = all_results.div(all_results['pandas'], axis=0)
all_results = all_results.ix[:, ['pandas', 'data.table', 'plyr',
                                 'base::merge']]
