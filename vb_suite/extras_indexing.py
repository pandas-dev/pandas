#----------------------------------------------------------------------
# Thorough checks of all containers and all indexing types

from vbench.benchmark import Benchmark

SECTION = 'Exhaustive check of indexing and scalar value access'

common_setup = """from pandas_vb_common import *
"""


import pandas.util.testing as tm

MAX_ENTRIES = 100000

# FIXME: makeCustomIndexWithCache reimplements (sort of) tm.makeCustomIndex,
# because the latter doesn't offer customization of date/period index
# frequencies and integer index offset.

setup_template = common_setup + """
import sys
import pandas as pd

try:
    make_index = tm.makeCustomIndexWithCache
except AttributeError:
    MAX_ENTRIES = %(MAX_ENTRIES)s
    _indices = {}

    def makeCustomIndexWithCache(nentries, idx_type):
        assert nentries <= MAX_ENTRIES

        key = idx_type
        try:
            full_idx = _indices[key]
        except KeyError:
            if idx_type == 'mi':
                full_idx = tm.makeCustomIndex(nentries=MAX_ENTRIES, nlevels=2)
            elif idx_type == 'dt':
                full_idx = pd.date_range('2000-01-01', periods=MAX_ENTRIES, freq='T')
            elif idx_type == 'p':
                full_idx = pd.period_range('2000-01-01', periods=MAX_ENTRIES, freq='T')
            elif idx_type == 's':
                full_idx = tm.makeStringIndex(k=MAX_ENTRIES)
            elif idx_type == 'u':
                full_idx = tm.makeUnicodeIndex(k=MAX_ENTRIES)
            elif idx_type == 'i':
                full_idx = pd.Index(np.arange(MAX_ENTRIES) + MAX_ENTRIES)
            elif idx_type == 'f':
                full_idx = tm.makeFloatIndex(MAX_ENTRIES)
            else:
                raise ValueError('Wrong idx type: %%s' %% idx_type)

            _indices[key] = full_idx

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

# generate_index_benchmarks(
#                     klass, long_axis=axis, idx_type=idx_type, is_dup=is_dup)


def generate_index_benchmarks(klass, idx_type, long_axis):
    ndim = klass().ndim

    shape = [10] * ndim
    shape[long_axis] = MAX_ENTRIES
    shape = tuple(shape)

    types = ['i'] * ndim
    types[long_axis] = idx_type
    types = tuple(types)

    axes = klass._AXIS_ORDERS
    ctor_args = ',\n    '.join([
        '%s=make_index(nentries=%r, idx_type=%r)' % v
        for v in zip(axes, shape, types)])

    def get_benchmark_name(indexer, axis):
        shape_type_str = 'x'.join([str(s) + str(t)
                                   for s, t in zip(shape, types)])

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
                              'ctor_args': ctor_args, 'axis': axis or 0,
                              'MAX_ENTRIES': MAX_ENTRIES},
                          name=get_benchmark_name(params['indexer'], axis))
            result[b.name] = b

    return result

# Benchmarks are generated as follows: given a container type, generate an
# instance of it with one of the axes long enough to produce statistically
# significant timing values and try different kinds of indexing on it.
#
# Generated benchmark set involves a cartesian product of
# - container types
# - designated "long" axis (minor or major one)
# - "long" axis type (string, integer, datetime, period, multiindex)
# - indexer type (positional, slice, fancy, etc.)
# - indexer axis (indexing is not limited to "long" axis)
# - label/positional indexer
#
# FIXME: add multiindex indexers?
# FIXME: add non-unique axes?
# FIXME: add non-unique non-monotonic axes?
for klass in (tm.Series, tm.DataFrame, tm.Panel):
    for axis in set([0, klass().ndim - 1]):
        for idx_type in ('s', 'i', 'dt', 'p', 'mi'):
                bms = generate_index_benchmarks(
                    klass, long_axis=axis, idx_type=idx_type)
                globals().update(bms)
