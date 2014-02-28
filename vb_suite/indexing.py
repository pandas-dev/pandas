from vbench.benchmark import Benchmark
from datetime import datetime

SECTION = 'Indexing and scalar value access'

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# Series.__getitem__, get_value, __getitem__(slice)

setup = common_setup + """
tm.N = 1000
ts = tm.makeTimeSeries()
dt = ts.index[500]
"""
statement = "ts[dt]"

bm_getitem = Benchmark(statement, setup, ncalls=100000,
                       name='series_getitem_scalar')

setup = common_setup + """
index = [tm.rands(10) for _ in xrange(1000)]
s = Series(np.random.rand(1000), index=index)
idx = index[100]
"""
statement = "s.get_value(idx)"
bm_df_getitem3 = Benchmark(statement, setup,
                           name='series_get_value',
                           start_date=datetime(2011, 11, 12))


setup = common_setup + """
index = tm.makeStringIndex(1000000)
s = Series(np.random.rand(1000000), index=index)
"""
series_getitem_pos_slice = Benchmark("s[:800000]", setup,
                                     name="series_getitem_pos_slice")


setup = common_setup + """
index = tm.makeStringIndex(1000000)
s = Series(np.random.rand(1000000), index=index)
lbl = s.index[800000]
"""
series_getitem_label_slice = Benchmark("s[:lbl]", setup,
                                       name="series_getitem_label_slice")


#----------------------------------------------------------------------
# DataFrame __getitem__

setup = common_setup + """
index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = DataFrame(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem = Benchmark(statement, setup,
                          name='dataframe_getitem_scalar')

setup = common_setup + """
try:
    klass = DataMatrix
except:
    klass = DataFrame

index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = klass(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem2 = Benchmark(statement, setup,
                           name='datamatrix_getitem_scalar')


#----------------------------------------------------------------------
# ix get scalar

setup = common_setup + """
index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = DataFrame(np.random.randn(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""

indexing_frame_get_value_ix = Benchmark("df.ix[idx,col]", setup,
                                        name='indexing_frame_get_value_ix',
                                        start_date=datetime(2011, 11, 12))

indexing_frame_get_value = Benchmark("df.get_value(idx,col)", setup,
                                     name='indexing_frame_get_value',
                                     start_date=datetime(2011, 11, 12))

setup = common_setup + """
mi = MultiIndex.from_tuples([(x,y) for x in range(1000) for y in range(1000)])
s =  Series(np.random.randn(1000000), index=mi)
"""

series_xs_mi_ix = Benchmark("s.ix[999]", setup,
                            name='series_xs_mi_ix',
                            start_date=datetime(2013, 1, 1))

setup = common_setup + """
mi = MultiIndex.from_tuples([(x,y) for x in range(1000) for y in range(1000)])
s =  Series(np.random.randn(1000000), index=mi)
df = DataFrame(s)
"""

frame_xs_mi_ix = Benchmark("df.ix[999]", setup,
                           name='frame_xs_mi_ix',
                           start_date=datetime(2013, 1, 1))

#----------------------------------------------------------------------
# Boolean DataFrame row selection

setup = common_setup + """
df  = DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
indexer = df['B'] > 0
obj_indexer = indexer.astype('O')
"""
indexing_dataframe_boolean_rows = \
    Benchmark("df[indexer]", setup, name='indexing_dataframe_boolean_rows')

indexing_dataframe_boolean_rows_object = \
    Benchmark("df[obj_indexer]", setup,
              name='indexing_dataframe_boolean_rows_object')

setup = common_setup + """
df  = DataFrame(np.random.randn(50000, 100))
df2 = DataFrame(np.random.randn(50000, 100))
"""
indexing_dataframe_boolean = \
    Benchmark("df > df2", setup, name='indexing_dataframe_boolean',
              start_date=datetime(2012, 1, 1))

setup = common_setup + """
import pandas.computation.expressions as expr
df  = DataFrame(np.random.randn(50000, 100))
df2 = DataFrame(np.random.randn(50000, 100))
expr.set_numexpr_threads(1)
"""

indexing_dataframe_boolean_st = \
    Benchmark("df > df2", setup, name='indexing_dataframe_boolean_st',cleanup="expr.set_numexpr_threads()",
              start_date=datetime(2013, 2, 26))


setup = common_setup + """
import pandas.computation.expressions as expr
df  = DataFrame(np.random.randn(50000, 100))
df2 = DataFrame(np.random.randn(50000, 100))
expr.set_use_numexpr(False)
"""

indexing_dataframe_boolean_no_ne = \
    Benchmark("df > df2", setup, name='indexing_dataframe_boolean_no_ne',cleanup="expr.set_use_numexpr(True)",
              start_date=datetime(2013, 2, 26))
#----------------------------------------------------------------------
# MultiIndex sortlevel

setup = common_setup + """
a = np.repeat(np.arange(100), 1000)
b = np.tile(np.arange(1000), 100)
midx = MultiIndex.from_arrays([a, b])
midx = midx.take(np.random.permutation(np.arange(100000)))
"""
sort_level_zero = Benchmark("midx.sortlevel(0)", setup,
                            start_date=datetime(2012, 1, 1))
sort_level_one = Benchmark("midx.sortlevel(1)", setup,
                           start_date=datetime(2012, 1, 1))

#----------------------------------------------------------------------
# Panel subset selection

setup = common_setup + """
p = Panel(np.random.randn(100, 100, 100))
inds = range(0, 100, 10)
"""

indexing_panel_subset = Benchmark('p.ix[inds, inds, inds]', setup,
                                  start_date=datetime(2012, 1, 1))

#----------------------------------------------------------------------
# Iloc

setup = common_setup + """
df = DataFrame({'A' : [0.1] * 3000, 'B' : [1] * 3000})
idx = np.array(range(30)) * 99
df2 = DataFrame({'A' : [0.1] * 1000, 'B' : [1] * 1000})
df2 = concat([df2, 2*df2, 3*df2])
"""

frame_iloc_dups = Benchmark('df2.iloc[idx]', setup,
                            start_date=datetime(2013, 1, 1))

frame_loc_dups = Benchmark('df2.loc[idx]', setup,
                            start_date=datetime(2013, 1, 1))

setup = common_setup + """
df = DataFrame(dict( A = [ 'foo'] * 1000000))
"""

frame_iloc_big = Benchmark('df.iloc[:100,0]', setup,
                            start_date=datetime(2013, 1, 1))


#----------------------------------------------------------------------
# Thorough checks of all containers and all indexing types

import pandas.util.testing as tm

setup_template = common_setup + """
import sys

try:
    make_index = tm.makeCustomIndexWithCache
except AttributeError:
    MAX_ENTRIES = 1000000
    _indices = {}

    def makeCustomIndexWithCache(nentries, **kwargs):
        assert nentries < MAX_ENTRIES

        key = tuple(kwargs.items())
        try:
            full_idx = _indices[key]
        except KeyError:
            full_idx = _indices[key] = tm.makeCustomIndex(nentries=MAX_ENTRIES,
                                                          **kwargs)
        return full_idx[:nentries]

    make_index = tm.makeCustomIndexWithCache = makeCustomIndexWithCache

obj = %(class_name)s(%(ctor_args)s)

pos = -1
axis = obj._get_axis(%(axis)r)
label = axis[pos]
arr_pos = np.arange(int(len(axis) / 2))
arr_label = axis[arr_pos].values
mask = tm.np.arange(len(axis)) %% 3 == 0
series_mask = Series(mask)
"""

def generate_index_benchmarks(klass, idx_type, shape):
    if not isinstance(shape, tuple):
        shape = (shape,)
    ndim = len(shape)

    if not isinstance(idx_type, tuple):
        idx_types = tuple([idx_type] * ndim)
    else:
        assert len(idx_type) == ndim
        idx_types = idx_type

    axes = klass._AXIS_ORDERS
    ctor_args = ',\n    '.join([
        '%s=make_index(idx_type=%r, nentries=%s, nlevels=1)' % v
        for v in zip(axes, idx_types, shape)])

    def get_benchmark_name(indexer, axis):
        shape_type_str = 'x'.join([str(s) + str(t)
                                   for s, t in zip(shape, idx_types)])

        components = ['indexing_', klass.__name__.lower(), indexer,
                      shape_type_str]
        if axis is not None:
            components.append("ax%s" % axis)

        return '_'.join(components)

    def make_suffix(attrname, indexer_str, axis):
        if axis is not None:
            indexers = [':,'] * ndim
            indexers[axis] = indexer_str + ','
            indexer_str = ''.join(indexers)
        return '%s[%s]' % (attrname, indexer_str)

    benchmarked_axes = set([None, 0, ndim - 1])

    result = {}
    for axis in benchmarked_axes:
        for params in [
            {'indexer': 'basic_pos',
             'suffix': make_suffix('.iloc', 'pos', axis)},
            {'indexer': 'basic_label',
             'suffix': make_suffix('.loc', 'label', axis)},

            {'indexer': 'slice_pos',
             'suffix': make_suffix('.iloc', ':pos', axis)},
            {'indexer': 'slice_label',
             'suffix': make_suffix('.loc', ':label', axis)},

            {'indexer': 'arr_pos',
             'suffix': make_suffix('.iloc', 'arr_pos', axis)},
            {'indexer': 'arr_label',
             'suffix': make_suffix('.loc', 'arr_label', axis)},

            {'indexer': 'iloc_mask',
             'suffix': make_suffix('.iloc', 'mask', axis)},
            {'indexer': 'loc_mask',
             'suffix': make_suffix('.loc', 'mask', axis)}, ]:

            b = Benchmark('obj%s' % params['suffix'],
                          setup_template % {
                              'class_name': klass.__name__,
                              'ctor_args': ctor_args, 'axis': axis or 0},
                          name=get_benchmark_name(params['indexer'], axis))
            result[b.name] = b

    return result

globals().update(generate_index_benchmarks(tm.Series, 's', 100000))
globals().update(generate_index_benchmarks(tm.DataFrame, 's', (10, 100000)))
globals().update(generate_index_benchmarks(tm.DataFrame, 's', (100000, 10)))
globals().update(generate_index_benchmarks(tm.Panel, 's', (100000, 10, 10)))
globals().update(generate_index_benchmarks(tm.Panel, 's', (10, 10, 100000)))
