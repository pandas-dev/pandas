import collections
from textwrap import dedent
from functools import partial

import numpy as np

from pandas.compat import map
import pandas.core.common as com
from pandas.core.series import Series
from pandas.core.generic import _shared_docs
from pandas.core.groupby.groupby import GroupBy
from pandas.core.groupby import base
from pandas.util._decorators import Appender
from pandas.core.index import MultiIndex, CategoricalIndex
from pandas.core.arrays.categorical import Categorical
from pandas.core.frame import DataFrame
from pandas.core.groupby.groupby import NDFrameGroupBy
from pandas.core.groupby.series import SeriesGroupBy
from pandas._libs.lib import count_level_2d
from pandas.plotting._core import boxplot_frame_groupby
from pandas.core.internals import BlockManager, make_block


class DataFrameGroupBy(NDFrameGroupBy):

    _apply_whitelist = base.dataframe_apply_whitelist

    #
    # Make class defs of attributes on DataFrameGroupBy whitelist.
    for _def_str in base.whitelist_method_generator(
            GroupBy, DataFrame, _apply_whitelist):
        exec(_def_str)

    _block_agg_axis = 1

    _agg_doc = dedent("""
    Examples
    --------

    >>> df = pd.DataFrame({'A': [1, 1, 2, 2],
    ...                    'B': [1, 2, 3, 4],
    ...                    'C': np.random.randn(4)})

    >>> df
       A  B         C
    0  1  1  0.362838
    1  1  2  0.227877
    2  2  3  1.267767
    3  2  4 -0.562860

    The aggregation is for each column.

    >>> df.groupby('A').agg('min')
       B         C
    A
    1  1  0.227877
    2  3 -0.562860

    Multiple aggregations

    >>> df.groupby('A').agg(['min', 'max'])
        B             C
      min max       min       max
    A
    1   1   2  0.227877  0.362838
    2   3   4 -0.562860  1.267767

    Select a column for aggregation

    >>> df.groupby('A').B.agg(['min', 'max'])
       min  max
    A
    1    1    2
    2    3    4

    Different aggregations per column

    >>> df.groupby('A').agg({'B': ['min', 'max'], 'C': 'sum'})
        B             C
      min max       sum
    A
    1   1   2  0.590716
    2   3   4  0.704907

    See also
    --------
    pandas.DataFrame.groupby.apply
    pandas.DataFrame.groupby.transform
    pandas.DataFrame.aggregate

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        klass='DataFrame',
        versionadded='',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):
        return super(DataFrameGroupBy, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    def _gotitem(self, key, ndim, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """

        if ndim == 2:
            if subset is None:
                subset = self.obj
            return DataFrameGroupBy(subset, self.grouper, selection=key,
                                    grouper=self.grouper,
                                    exclusions=self.exclusions,
                                    as_index=self.as_index)
        elif ndim == 1:
            if subset is None:
                subset = self.obj[key]
            return SeriesGroupBy(subset, selection=key,
                                 grouper=self.grouper)

        raise AssertionError("invalid ndim for _gotitem")

    def _wrap_generic_output(self, result, obj):
        result_index = self.grouper.levels[0]

        if self.axis == 0:
            return DataFrame(result, index=obj.columns,
                             columns=result_index).T
        else:
            return DataFrame(result, index=obj.index,
                             columns=result_index)

    def _get_data_to_aggregate(self):
        obj = self._obj_with_exclusions
        if self.axis == 1:
            return obj.T._data, 1
        else:
            return obj._data, 1

    def _insert_inaxis_grouper_inplace(self, result):
        # zip in reverse so we can always insert at loc 0
        izip = zip(* map(reversed, (
            self.grouper.names,
            self.grouper.get_group_levels(),
            [grp.in_axis for grp in self.grouper.groupings])))

        for name, lev, in_axis in izip:
            if in_axis:
                result.insert(0, name, lev)

    def _wrap_aggregated_output(self, output, names=None):
        agg_axis = 0 if self.axis == 1 else 1
        agg_labels = self._obj_with_exclusions._get_axis(agg_axis)

        output_keys = self._decide_output_index(output, agg_labels)

        if not self.as_index:
            result = DataFrame(output, columns=output_keys)
            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
        else:
            index = self.grouper.result_index
            result = DataFrame(output, index=index, columns=output_keys)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _wrap_transformed_output(self, output, names=None):
        return DataFrame(output, index=self.obj.index)

    def _wrap_agged_blocks(self, items, blocks):
        if not self.as_index:
            index = np.arange(blocks[0].values.shape[-1])
            mgr = BlockManager(blocks, [items, index])
            result = DataFrame(mgr)

            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
        else:
            index = self.grouper.result_index
            mgr = BlockManager(blocks, [items, index])
            result = DataFrame(mgr)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _reindex_output(self, result):
        """
        If we have categorical groupers, then we want to make sure that
        we have a fully reindex-output to the levels. These may have not
        participated in the groupings (e.g. may have all been
        nan groups);

        This can re-expand the output space
        """

        # we need to re-expand the output space to accomodate all values
        # whether observed or not in the cartesian product of our groupes
        groupings = self.grouper.groupings
        if groupings is None:
            return result
        elif len(groupings) == 1:
            return result

        # if we only care about the observed values
        # we are done
        elif self.observed:
            return result

        # reindexing only applies to a Categorical grouper
        elif not any(isinstance(ping.grouper, (Categorical, CategoricalIndex))
                     for ping in groupings):
            return result

        levels_list = [ping.group_index for ping in groupings]
        index, _ = MultiIndex.from_product(
            levels_list, names=self.grouper.names).sortlevel()

        if self.as_index:
            d = {self.obj._get_axis_name(self.axis): index, 'copy': False}
            return result.reindex(**d)

        # GH 13204
        # Here, the categorical in-axis groupers, which need to be fully
        # expanded, are columns in `result`. An idea is to do:
        # result = result.set_index(self.grouper.names)
        #                .reindex(index).reset_index()
        # but special care has to be taken because of possible not-in-axis
        # groupers.
        # So, we manually select and drop the in-axis grouper columns,
        # reindex `result`, and then reset the in-axis grouper columns.

        # Select in-axis groupers
        in_axis_grps = [(i, ping.name) for (i, ping)
                        in enumerate(groupings) if ping.in_axis]
        g_nums, g_names = zip(*in_axis_grps)

        result = result.drop(labels=list(g_names), axis=1)

        # Set a temp index and reindex (possibly expanding)
        result = result.set_index(self.grouper.result_index
                                  ).reindex(index, copy=False)

        # Reset in-axis grouper columns
        # (using level numbers `g_nums` because level names may not be unique)
        result = result.reset_index(level=g_nums)

        return result.reset_index(drop=True)

    def _iterate_column_groupbys(self):
        for i, colname in enumerate(self._selected_obj.columns):
            yield colname, SeriesGroupBy(self._selected_obj.iloc[:, i],
                                         selection=colname,
                                         grouper=self.grouper,
                                         exclusions=self.exclusions)

    def _apply_to_column_groupbys(self, func):
        from pandas.core.reshape.concat import concat
        return concat(
            (func(col_groupby) for _, col_groupby
             in self._iterate_column_groupbys()),
            keys=self._selected_obj.columns, axis=1)

    def _fill(self, direction, limit=None):
        """Overridden method to join grouped columns in output"""
        res = super(DataFrameGroupBy, self)._fill(direction, limit=limit)
        output = collections.OrderedDict(
            (grp.name, grp.grouper) for grp in self.grouper.groupings)

        from pandas import concat
        return concat((self._wrap_transformed_output(output), res), axis=1)

    def count(self):
        """ Compute count of group, excluding missing values """
        from pandas.core.dtypes.missing import _isna_ndarraylike as _isna

        data, _ = self._get_data_to_aggregate()
        ids, _, ngroups = self.grouper.group_info
        mask = ids != -1

        val = ((mask & ~_isna(np.atleast_2d(blk.get_values())))
               for blk in data.blocks)
        loc = (blk.mgr_locs for blk in data.blocks)

        counter = partial(count_level_2d, labels=ids, max_bin=ngroups, axis=1)
        blk = map(make_block, map(counter, val), loc)

        return self._wrap_agged_blocks(data.items, list(blk))

    def nunique(self, dropna=True):
        """
        Return DataFrame with number of distinct observations per group for
        each column.

        .. versionadded:: 0.20.0

        Parameters
        ----------
        dropna : boolean, default True
            Don't include NaN in the counts.

        Returns
        -------
        nunique: DataFrame

        Examples
        --------
        >>> df = pd.DataFrame({'id': ['spam', 'egg', 'egg', 'spam',
        ...                           'ham', 'ham'],
        ...                    'value1': [1, 5, 5, 2, 5, 5],
        ...                    'value2': list('abbaxy')})
        >>> df
             id  value1 value2
        0  spam       1      a
        1   egg       5      b
        2   egg       5      b
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y

        >>> df.groupby('id').nunique()
            id  value1  value2
        id
        egg    1       1       1
        ham    1       1       2
        spam   1       2       1

        # check for rows with the same id but conflicting values
        >>> df.groupby('id').filter(lambda g: (g.nunique() > 1).any())
             id  value1 value2
        0  spam       1      a
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y
        """

        obj = self._selected_obj

        def groupby_series(obj, col=None):
            return SeriesGroupBy(obj,
                                 selection=col,
                                 grouper=self.grouper).nunique(dropna=dropna)

        if isinstance(obj, Series):
            results = groupby_series(obj)
        else:
            from pandas.core.reshape.concat import concat
            results = [groupby_series(obj[col], col) for col in obj.columns]
            results = concat(results, axis=1)

        if not self.as_index:
            results.index = com._default_index(len(results))
        return results

    boxplot = boxplot_frame_groupby
