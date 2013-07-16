# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import itertools

import numpy as np

from pandas.core.series import Series
from pandas.core.frame import DataFrame

from pandas.core.categorical import Categorical
from pandas.core.common import (notnull, _ensure_platform_int, _maybe_promote,
                                isnull)
from pandas.core.groupby import (get_group_index, _compress_group_index,
                                 decons_group_index)
import pandas.core.common as com
import pandas.algos as algos

from pandas.core.index import MultiIndex


class ReshapeError(Exception):
    pass


class _Unstacker(object):
    """
    Helper class to unstack data / pivot with multi-level index

    Parameters
    ----------
    level : int or str, default last level
        Level to "unstack". Accepts a name for the level.

    Examples
    --------
    >>> import pandas as pd
    >>> index = pd.MultiIndex.from_tuples([('one', 'a'), ('one', 'b'),
    ...                                    ('two', 'a'), ('two', 'b')])
    >>> s = pd.Series(np.arange(1.0, 5.0), index=index)
    >>> s
    one  a   1
         b   2
    two  a   3
         b   4
    dtype: float64

    >>> s.unstack(level=-1)
         a   b
    one  1  2
    two  3  4

    >>> s.unstack(level=0)
       one  two
    a  1   2
    b  3   4

    Returns
    -------
    unstacked : DataFrame
    """
    def __init__(self, values, index, level=-1, value_columns=None):
        if values.ndim == 1:
            values = values[:, np.newaxis]
        self.values = values
        self.value_columns = value_columns

        if value_columns is None and values.shape[1] != 1:  # pragma: no cover
            raise ValueError('must pass column labels for multi-column data')

        self.index = index
        self.level = self.index._get_level_number(level)

        levels = index.levels
        labels = index.labels
        def _make_index(lev,lab):
            i = lev.__class__(_make_index_array_level(lev.values,lab))
            i.name = lev.name
            return i

        self.new_index_levels = list([ _make_index(lev,lab) for lev,lab in zip(levels,labels) ])
        self.new_index_names = list(index.names)

        self.removed_name = self.new_index_names.pop(self.level)
        self.removed_level = self.new_index_levels.pop(self.level)

        self._make_sorted_values_labels()
        self._make_selectors()

    def _make_sorted_values_labels(self):
        v = self.level

        labs = self.index.labels
        levs = self.index.levels
        to_sort = labs[:v] + labs[v + 1:] + [labs[v]]
        sizes = [len(x) for x in levs[:v] + levs[v + 1:] + [levs[v]]]

        comp_index, obs_ids = get_compressed_ids(to_sort, sizes)

        # group_index = get_group_index(to_sort, sizes)
        # comp_index, obs_ids = _compress_group_index(group_index)

        ngroups = len(obs_ids)

        indexer = algos.groupsort_indexer(comp_index, ngroups)[0]
        indexer = _ensure_platform_int(indexer)

        self.sorted_values = com.take_nd(self.values, indexer, axis=0)
        self.sorted_labels = [l.take(indexer) for l in to_sort]

    def _make_selectors(self):
        new_levels = self.new_index_levels

        # make the mask
        remaining_labels = self.sorted_labels[:-1]
        level_sizes = [len(x) for x in new_levels]

        comp_index, obs_ids = get_compressed_ids(remaining_labels, level_sizes)
        ngroups = len(obs_ids)

        comp_index = _ensure_platform_int(comp_index)
        stride = self.index.levshape[self.level]
        self.full_shape = ngroups, stride

        selector = self.sorted_labels[-1] + stride * comp_index
        mask = np.zeros(np.prod(self.full_shape), dtype=bool)
        mask.put(selector, True)

        if mask.sum() < len(self.index):
            raise ReshapeError('Index contains duplicate entries, '
                               'cannot reshape')

        self.group_index = comp_index
        self.mask = mask
        self.unique_groups = obs_ids
        self.compressor = comp_index.searchsorted(np.arange(ngroups))

    def get_result(self):
        # TODO: find a better way than this masking business

        values, value_mask = self.get_new_values()
        columns = self.get_new_columns()
        index = self.get_new_index()

        # filter out missing levels
        if values.shape[1] > 0:
            col_inds, obs_ids = _compress_group_index(self.sorted_labels[-1])
            # rare case, level values not observed
            if len(obs_ids) < self.full_shape[1]:
                inds = (value_mask.sum(0) > 0).nonzero()[0]
                values = com.take_nd(values, inds, axis=1)
                columns = columns[inds]

        # we might have a missing index
        if len(index) != values.shape[0]:
            mask = isnull(index)
            if mask.any():
                l = np.arange(len(index))
                values, orig_values = np.empty((len(index),values.shape[1])), values
                values.fill(np.nan)
                values_indexer = com._ensure_int64(l[~mask])
                for i, j in enumerate(values_indexer):
                    values[j] = orig_values[i]
            else:
                index = index.take(self.unique_groups)

        return DataFrame(values, index=index, columns=columns)

    def get_new_values(self):
        values = self.values

        # place the values
        length, width = self.full_shape
        stride = values.shape[1]
        result_width = width * stride
        result_shape = (length, result_width)

        # if our mask is all True, then we can use our existing dtype
        if self.mask.all():
            dtype = values.dtype
            new_values = np.empty(result_shape, dtype=dtype)
        else:
            dtype, fill_value = _maybe_promote(values.dtype)
            new_values = np.empty(result_shape, dtype=dtype)
            new_values.fill(fill_value)

        new_mask = np.zeros(result_shape, dtype=bool)

        # is there a simpler / faster way of doing this?
        for i in xrange(values.shape[1]):
            chunk = new_values[:, i * width: (i + 1) * width]
            mask_chunk = new_mask[:, i * width: (i + 1) * width]

            chunk.flat[self.mask] = self.sorted_values[:, i]
            mask_chunk.flat[self.mask] = True

        return new_values, new_mask

    def get_new_columns(self):
        if self.value_columns is None:
            return self.removed_level

        stride = len(self.removed_level)
        width = len(self.value_columns)
        propagator = np.repeat(np.arange(width), stride)
        if isinstance(self.value_columns, MultiIndex):
            new_levels = self.value_columns.levels + [self.removed_level]
            new_names = self.value_columns.names + [self.removed_name]

            new_labels = [lab.take(propagator)
                          for lab in self.value_columns.labels]
            new_labels.append(np.tile(np.arange(stride), width))
        else:
            new_levels = [self.value_columns, self.removed_level]
            new_names = [self.value_columns.name, self.removed_name]

            new_labels = []

            new_labels.append(propagator)
            new_labels.append(np.tile(np.arange(stride), width))

        return MultiIndex(levels=new_levels, labels=new_labels,
                          names=new_names)

    def get_new_index(self):
        result_labels = []
        for cur in self.sorted_labels[:-1]:
            labels = cur.take(self.compressor)
            labels = _make_index_array_level(labels,cur)
            result_labels.append(labels)

        # construct the new index
        if len(self.new_index_levels) == 1:
            new_index = self.new_index_levels[0]
            new_index.name = self.new_index_names[0]
        else:
            new_index = MultiIndex(levels=self.new_index_levels,
                                   labels=result_labels,
                                   names=self.new_index_names)

        return new_index


def _make_index_array_level(lev,lab):
    """ create the combined index array, preserving nans, return an array """
    mask = lab == -1
    if not mask.any():
        return lev

    l = np.arange(len(lab))
    mask_labels  = np.empty(len(mask[mask]),dtype=object)
    mask_labels.fill(np.nan)
    mask_indexer = com._ensure_int64(l[mask])

    labels = lev
    labels_indexer = com._ensure_int64(l[~mask])

    new_labels = np.empty(tuple([len(lab)]),dtype=object)
    new_labels[labels_indexer] = labels
    new_labels[mask_indexer]   = mask_labels

    return new_labels

def _unstack_multiple(data, clocs):
    if len(clocs) == 0:
        return data

    # NOTE: This doesn't deal with hierarchical columns yet

    index = data.index

    clocs = [index._get_level_number(i) for i in clocs]

    rlocs = [i for i in range(index.nlevels) if i not in clocs]

    clevels = [index.levels[i] for i in clocs]
    clabels = [index.labels[i] for i in clocs]
    cnames = [index.names[i] for i in clocs]
    rlevels = [index.levels[i] for i in rlocs]
    rlabels = [index.labels[i] for i in rlocs]
    rnames = [index.names[i] for i in rlocs]

    shape = [len(x) for x in clevels]
    group_index = get_group_index(clabels, shape)

    comp_ids, obs_ids = _compress_group_index(group_index, sort=False)
    recons_labels = decons_group_index(obs_ids, shape)

    dummy_index = MultiIndex(levels=rlevels + [obs_ids],
                             labels=rlabels + [comp_ids],
                             names=rnames + ['__placeholder__'])

    if isinstance(data, Series):
        dummy = Series(data.values, index=dummy_index)
        unstacked = dummy.unstack('__placeholder__')
        new_levels = clevels
        new_names = cnames
        new_labels = recons_labels
    else:
        if isinstance(data.columns, MultiIndex):
            result = data
            for i in range(len(clocs)):
                val = clocs[i]
                result = result.unstack(val)
                clocs = [val if i > val else val - 1 for val in clocs]

            return result

        dummy = DataFrame(data.values, index=dummy_index,
                          columns=data.columns)

        unstacked = dummy.unstack('__placeholder__')
        if isinstance(unstacked, Series):
            unstcols = unstacked.index
        else:
            unstcols = unstacked.columns
        new_levels = [unstcols.levels[0]] + clevels
        new_names = [data.columns.name] + cnames

        new_labels = [unstcols.labels[0]]
        for rec in recons_labels:
            new_labels.append(rec.take(unstcols.labels[-1]))

    new_columns = MultiIndex(levels=new_levels, labels=new_labels,
                             names=new_names)

    if isinstance(unstacked, Series):
        unstacked.index = new_columns
    else:
        unstacked.columns = new_columns

    return unstacked


def pivot(self, index=None, columns=None, values=None):
    """
    See DataFrame.pivot
    """
    if values is None:
        indexed = self.set_index([index, columns])
        return indexed.unstack(columns)
    else:
        indexed = Series(self[values].values,
                         index=[self[index], self[columns]])
        return indexed.unstack(columns)


def pivot_simple(index, columns, values):
    """
    Produce 'pivot' table based on 3 columns of this DataFrame.
    Uses unique values from index / columns and fills with values.

    Parameters
    ----------
    index : ndarray
        Labels to use to make new frame's index
    columns : ndarray
        Labels to use to make new frame's columns
    values : ndarray
        Values to use for populating new frame's values

    Note
    ----
    Obviously, all 3 of the input arguments must have the same length

    Returns
    -------
    DataFrame
    """
    if (len(index) != len(columns)) or (len(columns) != len(values)):
        raise AssertionError('Length of index, columns, and values must be the'
                             ' same')

    if len(index) == 0:
        return DataFrame(index=[])

    hindex = MultiIndex.from_arrays([index, columns])
    series = Series(values.ravel(), index=hindex)
    series = series.sortlevel(0)
    return series.unstack()


def _slow_pivot(index, columns, values):
    """
    Produce 'pivot' table based on 3 columns of this DataFrame.
    Uses unique values from index / columns and fills with values.

    Parameters
    ----------
    index : string or object
        Column name to use to make new frame's index
    columns : string or object
        Column name to use to make new frame's columns
    values : string or object
        Column name to use for populating new frame's values

    Could benefit from some Cython here.
    """
    tree = {}
    for i, (idx, col) in enumerate(itertools.izip(index, columns)):
        if col not in tree:
            tree[col] = {}
        branch = tree[col]
        branch[idx] = values[i]

    return DataFrame(tree)


def unstack(obj, level):
    if isinstance(level, (tuple, list)):
        return _unstack_multiple(obj, level)

    if isinstance(obj, DataFrame):
        if isinstance(obj.index, MultiIndex):
            return _unstack_frame(obj, level)
        else:
            return obj.T.stack(dropna=False)
    else:
        unstacker = _Unstacker(obj.values, obj.index, level=level)
        return unstacker.get_result()


def _unstack_frame(obj, level):
    from pandas.core.internals import BlockManager, make_block

    if obj._is_mixed_type:
        unstacker = _Unstacker(np.empty(obj.shape, dtype=bool),  # dummy
                               obj.index, level=level,
                               value_columns=obj.columns)
        new_columns = unstacker.get_new_columns()
        new_index = unstacker.get_new_index()
        new_axes = [new_columns, new_index]

        new_blocks = []
        mask_blocks = []
        for blk in obj._data.blocks:
            bunstacker = _Unstacker(blk.values.T, obj.index, level=level,
                                    value_columns=blk.items)
            new_items = bunstacker.get_new_columns()
            new_values, mask = bunstacker.get_new_values()

            mblk = make_block(mask.T, new_items, new_columns)
            mask_blocks.append(mblk)

            newb = make_block(new_values.T, new_items, new_columns)
            new_blocks.append(newb)

        result = DataFrame(BlockManager(new_blocks, new_axes))
        mask_frame = DataFrame(BlockManager(mask_blocks, new_axes))
        return result.ix[:, mask_frame.sum(0) > 0]
    else:
        unstacker = _Unstacker(obj.values, obj.index, level=level,
                               value_columns=obj.columns)
        return unstacker.get_result()


def get_compressed_ids(labels, sizes):
    # no overflow
    if com._long_prod(sizes) < 2 ** 63:
        group_index = get_group_index(labels, sizes)
        comp_index, obs_ids = _compress_group_index(group_index)
    else:
        n = len(labels[0])
        mask = np.zeros(n, dtype=bool)
        for v in labels:
            mask |= v < 0

        while com._long_prod(sizes) >= 2 ** 63:
            i = len(sizes)
            while com._long_prod(sizes[:i]) >= 2 ** 63:
                i -= 1

            rem_index, rem_ids = get_compressed_ids(labels[:i],
                                                    sizes[:i])
            sizes = [len(rem_ids)] + sizes[i:]
            labels = [rem_index] + labels[i:]

        return get_compressed_ids(labels, sizes)

    return comp_index, obs_ids


def stack(frame, level=-1, dropna=True):
    """
    Convert DataFrame to Series with multi-level Index. Columns become the
    second level of the resulting hierarchical index

    Returns
    -------
    stacked : Series
    """
    N, K = frame.shape
    if isinstance(level, int) and level < 0:
        level += frame.columns.nlevels

    level = frame.columns._get_level_number(level)

    if isinstance(frame.columns, MultiIndex):
        return _stack_multi_columns(frame, level=level, dropna=dropna)
    elif isinstance(frame.index, MultiIndex):
        new_levels = list(frame.index.levels)
        new_levels.append(frame.columns)

        new_labels = [lab.repeat(K) for lab in frame.index.labels]
        new_labels.append(np.tile(np.arange(K), N).ravel())

        new_names = list(frame.index.names)
        new_names.append(frame.columns.name)
        new_index = MultiIndex(levels=new_levels, labels=new_labels,
                               names=new_names)
    else:
        ilabels = np.arange(N).repeat(K)
        clabels = np.tile(np.arange(K), N).ravel()
        new_index = MultiIndex(levels=[frame.index, frame.columns],
                               labels=[ilabels, clabels],
                               names=[frame.index.name, frame.columns.name])

    new_values = frame.values.ravel()
    if dropna:
        mask = notnull(new_values)
        new_values = new_values[mask]
        new_index = new_index[mask]
    return Series(new_values, index=new_index)


def _stack_multi_columns(frame, level=-1, dropna=True):
    this = frame.copy()

    # this makes life much simpler
    if level != frame.columns.nlevels - 1:
        # roll levels to put selected level at end
        roll_columns = this.columns
        for i in range(level, frame.columns.nlevels - 1):
            roll_columns = roll_columns.swaplevel(i, i + 1)
        this.columns = roll_columns

    if not this.columns.is_lexsorted():
        this = this.sortlevel(0, axis=1)

    # tuple list excluding level for grouping columns
    if len(frame.columns.levels) > 2:
        tuples = zip(*[lev.values.take(lab)
                       for lev, lab in zip(this.columns.levels[:-1],
                                           this.columns.labels[:-1])])
        unique_groups = [key for key, _ in itertools.groupby(tuples)]
        new_names = this.columns.names[:-1]
        new_columns = MultiIndex.from_tuples(unique_groups, names=new_names)
    else:
        new_columns = unique_groups = this.columns.levels[0]

    # time to ravel the values
    new_data = {}
    level_vals = this.columns.levels[-1]
    levsize = len(level_vals)
    drop_cols = []
    for key in unique_groups:
        loc = this.columns.get_loc(key)
        slice_len = loc.stop - loc.start
        # can make more efficient?

        if slice_len == 0:
            drop_cols.append(key)
            continue
        elif slice_len != levsize:
            chunk = this.ix[:, this.columns[loc]]
            chunk.columns = level_vals.take(chunk.columns.labels[-1])
            value_slice = chunk.reindex(columns=level_vals).values
        else:
            if frame._is_mixed_type:
                value_slice = this.ix[:, this.columns[loc]].values
            else:
                value_slice = this.values[:, loc]

        new_data[key] = value_slice.ravel()

    if len(drop_cols) > 0:
        new_columns = new_columns - drop_cols

    N = len(this)

    if isinstance(this.index, MultiIndex):
        new_levels = list(this.index.levels)
        new_names = list(this.index.names)
        new_labels = [lab.repeat(levsize) for lab in this.index.labels]
    else:
        new_levels = [this.index]
        new_labels = [np.arange(N).repeat(levsize)]
        new_names = [this.index.name]  # something better?

    new_levels.append(frame.columns.levels[level])
    new_labels.append(np.tile(np.arange(levsize), N))
    new_names.append(frame.columns.names[level])

    new_index = MultiIndex(levels=new_levels, labels=new_labels,
                           names=new_names)

    result = DataFrame(new_data, index=new_index, columns=new_columns)

    # more efficient way to go about this? can do the whole masking biz but
    # will only save a small amount of time...
    if dropna:
        result = result.dropna(axis=0, how='all')

    return result


def melt(frame, id_vars=None, value_vars=None,
         var_name=None, value_name='value', col_level=None):
    """
    "Unpivots" a DataFrame from wide format to long format, optionally leaving
    id variables set

    Parameters
    ----------
    frame : DataFrame
    id_vars : tuple, list, or ndarray
    value_vars : tuple, list, or ndarray
    var_name : scalar, if None uses frame.column.name or 'variable'
    value_name : scalar, default 'value'
    col_level : scalar, if columns are a MultiIndex then use this level to melt

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'A': {0: 'a', 1: 'b', 2: 'c'},
    ...                    'B': {0: 1, 1: 3, 2: 5},
    ...                    'C': {0: 2, 1: 4, 2: 6}})

    >>> df
       A  B  C
    0  a  1  2
    1  b  3  4
    2  c  5  6

    >>> melt(df, id_vars=['A'], value_vars=['B'])
       A variable  value
    0  a        B      1
    1  b        B      3
    2  c        B      5

    >>> melt(df, id_vars=['A'], value_vars=['B'],
    ... var_name='myVarname', value_name='myValname')
       A myVarname  myValname
    0  a         B          1
    1  b         B          3
    2  c         B          5

    >>> df.columns = [list('ABC'), list('DEF')]

    >>> melt(df, col_level=0, id_vars=['A'], value_vars=['B'])
       A variable  value
    0  a        B      1
    1  b        B      3
    2  c        B      5

    >>> melt(df, id_vars=[('A', 'D')], value_vars=[('B', 'E')])
      (A, D) variable_0 variable_1  value
    0      a          B          E      1
    1      b          B          E      3
    2      c          B          E      5

    """
    # TODO: what about the existing index?
    if id_vars is not None:
        if not isinstance(id_vars, (tuple, list, np.ndarray)):
            id_vars = [id_vars]
        else:
            id_vars = list(id_vars)
    else:
        id_vars = []

    if value_vars is not None:
        if not isinstance(value_vars, (tuple, list, np.ndarray)):
            value_vars = [value_vars]
        frame = frame.ix[:, id_vars + value_vars]
    else:
        frame = frame.copy()

    if col_level is not None:  # allow list or other?
        frame.columns = frame.columns.get_level_values(col_level) #  frame is a copy

    if var_name is None:
        if isinstance(frame.columns, MultiIndex):
            if len(frame.columns.names) == len(set(frame.columns.names)):
                var_name = frame.columns.names
            else:
                var_name = ['variable_%s' % i for i in
                            xrange(len(frame.columns.names))]
        else:
            var_name = [frame.columns.name if frame.columns.name is not None
                        else 'variable']
    if isinstance(var_name, basestring):
        var_name = [var_name]

    N, K = frame.shape
    K -= len(id_vars)

    mdata = {}
    for col in id_vars:
        mdata[col] = np.tile(frame.pop(col).values, K)

    mcolumns = id_vars + var_name + [value_name]

    mdata[value_name] = frame.values.ravel('F')
    for i, col in enumerate(var_name):
        # asanyarray will keep the columns as an Index
        mdata[col] = np.asanyarray(frame.columns.get_level_values(i)).repeat(N)

    return DataFrame(mdata, columns=mcolumns)


def lreshape(data, groups, dropna=True, label=None):
    """
    Reshape long-format data to wide. Generalized inverse of DataFrame.pivot

    Parameters
    ----------
    data : DataFrame
    groups : dict
        {new_name : list_of_columns}
    dropna : boolean, default True

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.DataFrame({'hr1': [514, 573], 'hr2': [545, 526],
    ...                      'team': ['Red Sox', 'Yankees'],
    ...                      'year1': [2007, 2008], 'year2': [2008, 2008]})
    >>> data
       hr1  hr2     team  year1  year2
    0  514  545  Red Sox   2007   2008
    1  573  526  Yankees   2007   2008

    >>> pd.lreshape(data, {'year': ['year1', 'year2'], 'hr': ['hr1', 'hr2']})
          team   hr  year
    0  Red Sox  514  2007
    1  Yankees  573  2007
    2  Red Sox  545  2008
    3  Yankees  526  2008

    Returns
    -------
    reshaped : DataFrame
    """
    if isinstance(groups, dict):
        keys = groups.keys()
        values = groups.values()
    else:
        keys, values = zip(*groups)

    all_cols = list(set.union(*[set(x) for x in values]))
    id_cols = list(data.columns.diff(all_cols))

    K = len(values[0])

    for seq in values:
        if len(seq) != K:
            raise ValueError('All column lists must be same length')

    mdata = {}
    pivot_cols = []

    for target, names in zip(keys, values):
        mdata[target] = com._concat_compat([data[col].values for col in names])
        pivot_cols.append(target)

    for col in id_cols:
        mdata[col] = np.tile(data[col].values, K)

    if dropna:
        mask = np.ones(len(mdata[pivot_cols[0]]), dtype=bool)
        for c in pivot_cols:
            mask &= notnull(mdata[c])
        if not mask.all():
            mdata = dict((k, v[mask]) for k, v in mdata.iteritems())

    return DataFrame(mdata, columns=id_cols + pivot_cols)


def convert_dummies(data, cat_variables, prefix_sep='_'):
    """
    Compute DataFrame with specified columns converted to dummy variables (0 /
    1). Result columns will be prefixed with the column name, then the level
    name, e.g. 'A_foo' for column A and level foo

    Parameters
    ----------
    data : DataFrame
    cat_variables : list-like
        Must be column names in the DataFrame
    prefix_sep : string, default '_'
        String to use to separate column name from dummy level

    Returns
    -------
    dummies : DataFrame
    """
    result = data.drop(cat_variables, axis=1)
    for variable in cat_variables:
        dummies = get_dummies(data[variable], prefix=variable,
                              prefix_sep=prefix_sep)
        result = result.join(dummies)
    return result


def get_dummies(data, prefix=None, prefix_sep='_'):
    """
    Convert categorical variable into dummy/indicator variables

    Parameters
    ----------
    data : array-like or Series
    prefix : string, default None
        String to append DataFrame column names
    prefix_sep : string, default '_'
        If appending prefix, separator/delimiter to use

    Returns
    -------
    dummies : DataFrame
    """
    cat = Categorical.from_array(np.asarray(data))
    dummy_mat = np.eye(len(cat.levels)).take(cat.labels, axis=0)

    if prefix is not None:
        dummy_cols = ['%s%s%s' % (prefix, prefix_sep, str(v))
                      for v in cat.levels]
    else:
        dummy_cols = cat.levels

    if isinstance(data, Series):
        index = data.index
    else:
        index = None

    return DataFrame(dummy_mat, index=index, columns=dummy_cols)


def make_axis_dummies(frame, axis='minor', transform=None):
    """
    Construct 1-0 dummy variables corresponding to designated axis
    labels

    Parameters
    ----------
    frame : DataFrame
    axis : {'major', 'minor'}, default 'minor'
    transform : function, default None
        Function to apply to axis labels first. For example, to
        get "day of week" dummies in a time series regression
        you might call::

            make_axis_dummies(panel, axis='major',
                              transform=lambda d: d.weekday())
    Returns
    -------
    dummies : DataFrame
        Column names taken from chosen axis
    """
    numbers = {
        'major': 0,
        'minor': 1
    }
    num = numbers.get(axis, axis)

    items = frame.index.levels[num]
    labels = frame.index.labels[num]
    if transform is not None:
        mapped_items = items.map(transform)
        cat = Categorical.from_array(mapped_items.take(labels))
        labels = cat.labels
        items = cat.levels

    values = np.eye(len(items), dtype=float)
    values = values.take(labels, axis=0)

    return DataFrame(values, columns=items, index=frame.index)


def block2d_to_blocknd(values, items, shape, labels, ref_items=None):
    """ pivot to the labels shape """
    from pandas.core.internals import make_block
    panel_shape = (len(items),) + shape

    # TODO: lexsort depth needs to be 2!!

    # Create observation selection vector using major and minor
    # labels, for converting to panel format.
    selector = factor_indexer(shape[1:], labels)
    mask = np.zeros(np.prod(shape), dtype=bool)
    mask.put(selector, True)

    if mask.all():
        pvalues = np.empty(panel_shape, dtype=values.dtype)
    else:
        dtype, fill_value = _maybe_promote(values.dtype)
        pvalues = np.empty(panel_shape, dtype=dtype)
        pvalues.fill(fill_value)

    values = values
    for i in xrange(len(items)):
        pvalues[i].flat[mask] = values[:, i]

    if ref_items is None:
        ref_items = items

    return make_block(pvalues, items, ref_items)


def factor_indexer(shape, labels):
    """ given a tuple of shape and a list of Categorical labels, return the expanded label indexer """
    mult = np.array(shape)[::-1].cumprod()[::-1]
    return com._ensure_platform_int(np.sum(np.array(labels).T * np.append(mult, [1]), axis=1).T)
