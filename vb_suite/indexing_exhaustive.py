#----------------------------------------------------------------------
# Thorough checks of all containers and all indexing types

from vbench.benchmark import Benchmark

SECTION = 'Exhaustive check of indexing and scalar value access'

common_setup = """from pandas_vb_common import *
"""


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
