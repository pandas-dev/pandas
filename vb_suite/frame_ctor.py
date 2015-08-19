from vbench.benchmark import Benchmark
from datetime import datetime
try:
    import pandas.tseries.offsets as offsets
except:
    import pandas.core.datetools as offsets

common_setup = """from pandas_vb_common import *
try:
    from pandas.tseries.offsets import *
except:
    from pandas.core.datetools import *
"""

#----------------------------------------------------------------------
# Creation from nested dict

setup = common_setup + """
N, K = 5000, 50
index = tm.makeStringIndex(N)
columns = tm.makeStringIndex(K)
frame = DataFrame(np.random.randn(N, K), index=index, columns=columns)

try:
    data = frame.to_dict()
except:
    data = frame.toDict()

some_dict = data.values()[0]
dict_list = [dict(zip(columns, row)) for row in frame.values]
"""

frame_ctor_nested_dict = Benchmark("DataFrame(data)", setup)

# From JSON-like stuff
frame_ctor_list_of_dict = Benchmark("DataFrame(dict_list)", setup,
                                    start_date=datetime(2011, 12, 20))

series_ctor_from_dict = Benchmark("Series(some_dict)", setup)

# nested dict, integer indexes, regression described in #621
setup = common_setup + """
data = dict((i,dict((j,float(j)) for j in xrange(100))) for i in xrange(2000))
"""
frame_ctor_nested_dict_int64 = Benchmark("DataFrame(data)", setup)

# dynamically generate benchmarks for every offset
#
# get_period_count & get_index_for_offset are there because blindly taking each
# offset times 1000 can easily go out of Timestamp bounds and raise errors.
dynamic_benchmarks = {}
n_steps = [1, 2]
offset_kwargs = {'WeekOfMonth': {'weekday': 1, 'week': 1},
                 'LastWeekOfMonth': {'weekday': 1, 'week': 1},
                 'FY5253': {'startingMonth': 1, 'weekday': 1},
                 'FY5253Quarter': {'qtr_with_extra_week': 1, 'startingMonth': 1, 'weekday': 1}}

offset_extra_cases = {'FY5253': {'variation': ['nearest', 'last']},
                      'FY5253Quarter': {'variation': ['nearest', 'last']}}

for offset in offsets.__all__:
    for n in n_steps:
        kwargs = {}
        if offset in offset_kwargs:
            kwargs = offset_kwargs[offset]

        if offset in offset_extra_cases:
            extras = offset_extra_cases[offset]
        else:
            extras = {'': ['']}

        for extra_arg in extras:
            for extra in extras[extra_arg]:
                if extra:
                    kwargs[extra_arg] = extra
                setup = common_setup + """

def get_period_count(start_date, off):
    ten_offsets_in_days = ((start_date + off * 10) - start_date).days
    if ten_offsets_in_days == 0:
        return 1000
    else:
        return min(9 * ((Timestamp.max - start_date).days //
                        ten_offsets_in_days),
                   1000)

def get_index_for_offset(off):
    start_date = Timestamp('1/1/1900')
    return date_range(start_date,
                      periods=min(1000, get_period_count(start_date, off)),
                      freq=off)

idx = get_index_for_offset({}({}, **{}))
df = DataFrame(np.random.randn(len(idx),10), index=idx)
d = dict([ (col,df[col]) for col in df.columns ])
""".format(offset, n, kwargs)
                key = 'frame_ctor_dtindex_{}x{}'.format(offset, n)
                if extra:
                    key += '__{}_{}'.format(extra_arg, extra)
                dynamic_benchmarks[key] = Benchmark("DataFrame(d)", setup, name=key)

# Have to stuff them in globals() so vbench detects them
globals().update(dynamic_benchmarks)

# from a mi-series
setup = common_setup + """
mi = MultiIndex.from_tuples([(x,y) for x in range(100) for y in range(100)])
s = Series(randn(10000), index=mi)
"""
frame_from_series = Benchmark("DataFrame(s)", setup)

#----------------------------------------------------------------------
# get_numeric_data

setup = common_setup + """
df = DataFrame(randn(10000, 25))
df['foo'] = 'bar'
df['bar'] = 'baz'
df = df.consolidate()
"""

frame_get_numeric_data = Benchmark('df._get_numeric_data()', setup,
                                   start_date=datetime(2011, 11, 1))
