# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import itertools

import numpy as np

from pandas.core.series import Series
from pandas.core.frame import DataFrame
from pandas.core.panel import Panel

from pandas.core.common import notnull
from pandas.core.groupby import get_group_index
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
    >>> s
    one  a   1.
    one  b   2.
    two  a   3.
    two  b   4.

    >>> s.unstack(level=-1)
         a   b
    one  1.  2.
    two  3.  4.

    >>> s.unstack(level=0)
       one  two
    a  1.   2.
    b  3.   4.

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

        self.new_index_levels = list(index.levels)
        self.new_index_names = list(index.names)

        self.removed_name = self.new_index_names.pop(self.level)
        self.removed_level = self.new_index_levels.pop(self.level)

        v = self.level
        lshape = self.index.levshape
        self.full_shape = np.prod(lshape[:v] + lshape[v+1:]), lshape[v]

        self._make_sorted_values_labels()
        self._make_selectors()

    def _make_sorted_values_labels(self):
        v = self.level

        labs = self.index.labels
        to_sort = labs[:v] + labs[v+1:] + [labs[v]]
        indexer = np.lexsort(to_sort[::-1])

        self.sorted_values = self.values.take(indexer, axis=0)
        self.sorted_labels = [l.take(indexer) for l in to_sort]

    def _make_selectors(self):
        new_levels = self.new_index_levels

        # make the mask
        group_index = get_group_index(self.sorted_labels,
                                      [len(x) for x in new_levels])

        group_mask = np.zeros(self.full_shape[0], dtype=bool)
        group_mask.put(group_index, True)

        stride = self.index.levshape[self.level]
        selector = self.sorted_labels[-1] + stride * group_index
        mask = np.zeros(np.prod(self.full_shape), dtype=bool)
        mask.put(selector, True)

        # compress labels
        unique_groups = np.arange(self.full_shape[0])[group_mask]
        compressor = group_index.searchsorted(unique_groups)

        if mask.sum() < len(self.index):
            raise ReshapeError('Index contains duplicate entries, '
                               'cannot reshape')

        self.group_mask = group_mask
        self.group_index = group_index
        self.mask = mask
        self.unique_groups = unique_groups
        self.compressor = compressor

    def get_result(self):
        # TODO: find a better way than this masking business

        values, value_mask = self.get_new_values()
        columns = self.get_new_columns()
        index = self.get_new_index()

        # filter out missing levels
        if values.shape[1] > 0:
            mask = value_mask.sum(0) > 0
            values = values[:, mask]
            columns = columns[mask]

        return DataFrame(values, index=index, columns=columns)

    def get_new_values(self):
        return self._reshape_values(self.values)

    def _reshape_values(self, values):
        values = self.values
        # place the values
        length, width = self.full_shape
        stride = values.shape[1]
        result_width = width * stride

        new_values = np.empty((length, result_width), dtype=values.dtype)
        new_mask = np.zeros((length, result_width), dtype=bool)

        if issubclass(values.dtype.type, np.integer):
            new_values = new_values.astype(float)

        new_values.fill(np.nan)

        # is there a simpler / faster way of doing this?
        for i in xrange(self.values.shape[1]):
            chunk = new_values[:, i * width : (i + 1) * width]
            mask_chunk = new_mask[:, i * width : (i + 1) * width]

            chunk.flat[self.mask] = self.sorted_values[:, i]
            mask_chunk.flat[self.mask] = True

        new_values = new_values.take(self.unique_groups, axis=0)
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
            result_labels.append(cur.take(self.compressor))

        # construct the new index
        if len(self.new_index_levels) == 1:
            new_index = self.new_index_levels[0].take(self.unique_groups)
            new_index.name = self.new_index_names[0]
        else:
            new_index = MultiIndex(levels=self.new_index_levels,
                                   labels=result_labels,
                                   names=self.new_index_names)

        return new_index

def pivot(self, index=None, columns=None, values=None):
    """
    See DataFrame.pivot
    """
    index_vals = self[index]
    column_vals = self[columns]
    mindex = MultiIndex.from_arrays([index_vals, column_vals],
                                    names=[index, columns])

    if values is None:
        items = self.columns - [index, columns]
        mat = self.reindex(columns=items).values
    else:
        items = [values]
        mat = np.atleast_2d(self[values].values).T

    stacked = DataFrame(mat, index=mindex, columns=items)

    if not mindex.is_lexsorted():
        stacked = stacked.sortlevel(level=0)

    unstacked = stacked.unstack()
    if values is not None:
        unstacked.columns = unstacked.columns.droplevel(0)
    return unstacked

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
    assert(len(index) == len(columns) == len(values))

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
        unstacker = _Unstacker(np.empty(obj.shape, dtype=bool), # dummy
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
        return _stack_multi_columns(frame, level=level, dropna=True)
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
    for key in unique_groups:
        loc = this.columns.get_loc(key)

        # can make more efficient?
        if loc.stop - loc.start != levsize:
            chunk = this.ix[:, this.columns[loc]]
            chunk.columns = level_vals.take(chunk.columns.labels[-1])
            value_slice = chunk.reindex(columns=level_vals).values
        else:
            if frame._is_mixed_type:
                value_slice = this.ix[:, this.columns[loc]].values
            else:
                value_slice = this.values[:, loc]

        new_data[key] = value_slice.ravel()

    N = len(this)

    if isinstance(this.index, MultiIndex):
        new_levels = list(this.index.levels)
        new_names = list(this.index.names)
        new_labels = [lab.repeat(levsize) for lab in this.index.labels]
    else:
        new_levels = [this.index]
        new_labels = [np.arange(N).repeat(levsize)]
        new_names = [this.index.name] # something better?

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


def melt(frame, id_vars=None, value_vars=None):
    """
    "Unpivots" a DataFrame from wide format to long format, optionally leaving
    id variables set

    Parameters
    ----------
    frame : DataFrame
    id_vars :
    value_vars :

    Examples
    --------
    >>> df
    A B C
    a 1 2
    b 3 4
    c 5 6

    >>> melt(df, id_vars=['A'])
    A variable value
    a B        1
    b B        3
    c B        5
    a C        2
    b C        4
    c C        6
    """
    # TODO: what about the existing index?

    N, K = frame.shape

    mdata = {}

    if id_vars is not None:
        id_vars = list(id_vars)
        frame = frame.copy()
        K -= len(id_vars)
        for col in id_vars:
            mdata[col] = np.tile(frame.pop(col).values, K)
    else:
        id_vars = []

    mcolumns = id_vars + ['variable', 'value']

    mdata['value'] = frame.values.ravel('F')
    mdata['variable'] = np.asarray(frame.columns).repeat(N)
    return DataFrame(mdata, columns=mcolumns)

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
        dummies = make_column_dummies(data, variable, prefix=True,
                                      prefix_sep=prefix_sep)
        result = result.join(dummies)
    return result

def make_column_dummies(data, column, prefix=False, prefix_sep='_'):
    from pandas import Factor
    factor = Factor(data[column].values)
    dummy_mat = np.eye(len(factor.levels)).take(factor.labels, axis=0)

    if prefix:
        dummy_cols = ['%s%s%s' % (column, prefix_sep, str(v))
                      for v in factor.levels]
    else:
        dummy_cols = factor.levels
    dummies = DataFrame(dummy_mat, index=data.index, columns=dummy_cols)
    return dummies

def make_axis_dummies(frame, axis='minor', transform=None):
    """
    Construct 1-0 dummy variables corresponding to designated axis
    labels

    Parameters
    ----------
    axis : {'major', 'minor'}, default 'minor'
    transform : function, default None
        Function to apply to axis labels first. For example, to
        get "day of week" dummies in a time series regression you might
        call:
            make_axis_dummies(panel, axis='major',
                              transform=lambda d: d.weekday())
    Returns
    -------
    dummies : DataFrame
        Column names taken from chosen axis
    """
    from pandas import Factor

    numbers = {
        'major' : 0,
        'minor' : 1
    }
    num = numbers.get(axis, axis)

    items = frame.index.levels[num]
    labels = frame.index.labels[num]
    if transform is not None:
        mapped_items = items.map(transform)
        factor = Factor(mapped_items.take(labels))
        labels = factor.labels
        items = factor.levels

    values = np.eye(len(items), dtype=float)
    values = values.take(labels, axis=0)

    return DataFrame(values, columns=items, index=frame.index)

def block2d_to_block3d(values, items, shape, major_labels, minor_labels,
                       ref_items=None):
    """
    Developer method for pivoting DataFrame -> Panel. Used in HDFStore and
    DataFrame.to_panel
    """
    from pandas.core.internals import make_block
    panel_shape = (len(items),) + shape

    # TODO: lexsort depth needs to be 2!!

    # Create observation selection vector using major and minor
    # labels, for converting to panel format.
    selector = minor_labels + shape[1] * major_labels
    mask = np.zeros(np.prod(shape), dtype=bool)
    mask.put(selector, True)

    pvalues = np.empty(panel_shape, dtype=values.dtype)
    if not issubclass(pvalues.dtype.type, np.integer):
        pvalues.fill(np.nan)

    values = values
    for i in xrange(len(items)):
        pvalues[i].flat[mask] = values[:, i]

    if ref_items is None:
        ref_items = items

    return make_block(pvalues, items, ref_items)
