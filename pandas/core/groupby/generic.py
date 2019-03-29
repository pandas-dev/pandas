"""
Define the SeriesGroupBy and DataFrameGroupBy
classes that hold the groupby interfaces (and some implementations).

These are user facing as the result of the ``df.groupby(...)`` operations,
which here returns a DataFrameGroupBy object.
"""

import collections
import copy
from functools import partial
from textwrap import dedent
import warnings

import numpy as np

from pandas._libs import Timestamp, lib
import pandas.compat as compat
from pandas.compat import lzip
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, Substitution

from pandas.core.dtypes.cast import maybe_downcast_to_dtype
from pandas.core.dtypes.common import (
    ensure_int64, ensure_platform_int, is_bool, is_datetimelike,
    is_integer_dtype, is_interval_dtype, is_numeric_dtype, is_scalar)
from pandas.core.dtypes.missing import isna, notna

import pandas.core.algorithms as algorithms
from pandas.core.arrays import Categorical
from pandas.core.base import DataError, SpecificationError
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame, _shared_docs
from pandas.core.groupby import base
from pandas.core.groupby.groupby import (
    GroupBy, _apply_docs, _transform_template)
from pandas.core.index import CategoricalIndex, Index, MultiIndex
import pandas.core.indexes.base as ibase
from pandas.core.internals import BlockManager, make_block
from pandas.core.series import Series
from pandas.core.sparse.frame import SparseDataFrame

from pandas.plotting._core import boxplot_frame_groupby


class NDFrameGroupBy(GroupBy):

    def _iterate_slices(self):
        if self.axis == 0:
            # kludge
            if self._selection is None:
                slice_axis = self.obj.columns
            else:
                slice_axis = self._selection_list
            slicer = lambda x: self.obj[x]
        else:
            slice_axis = self.obj.index
            slicer = self.obj.xs

        for val in slice_axis:
            if val in self.exclusions:
                continue
            yield val, slicer(val)

    def _cython_agg_general(self, how, alt=None, numeric_only=True,
                            min_count=-1):
        new_items, new_blocks = self._cython_agg_blocks(
            how, alt=alt, numeric_only=numeric_only, min_count=min_count)
        return self._wrap_agged_blocks(new_items, new_blocks)

    def _wrap_agged_blocks(self, items, blocks):
        obj = self._obj_with_exclusions

        new_axes = list(obj._data.axes)

        # more kludge
        if self.axis == 0:
            new_axes[0], new_axes[1] = new_axes[1], self.grouper.result_index
        else:
            new_axes[self.axis] = self.grouper.result_index

        # Make sure block manager integrity check passes.
        assert new_axes[0].equals(items)
        new_axes[0] = items

        mgr = BlockManager(blocks, new_axes)

        new_obj = type(obj)(mgr)

        return self._post_process_cython_aggregate(new_obj)

    _block_agg_axis = 0

    def _cython_agg_blocks(self, how, alt=None, numeric_only=True,
                           min_count=-1):
        # TODO: the actual managing of mgr_locs is a PITA
        # here, it should happen via BlockManager.combine

        data, agg_axis = self._get_data_to_aggregate()

        if numeric_only:
            data = data.get_numeric_data(copy=False)

        new_blocks = []
        new_items = []
        deleted_items = []
        for block in data.blocks:

            locs = block.mgr_locs.as_array
            try:
                result, _ = self.grouper.aggregate(
                    block.values, how, axis=agg_axis, min_count=min_count)
            except NotImplementedError:
                # generally if we have numeric_only=False
                # and non-applicable functions
                # try to python agg

                if alt is None:
                    # we cannot perform the operation
                    # in an alternate way, exclude the block
                    deleted_items.append(locs)
                    continue

                # call our grouper again with only this block
                from pandas.core.groupby.groupby import groupby

                obj = self.obj[data.items[locs]]
                s = groupby(obj, self.grouper)
                result = s.aggregate(lambda x: alt(x, axis=self.axis))

            finally:

                # see if we can cast the block back to the original dtype
                result = block._try_coerce_and_cast_result(result)
                newb = block.make_block(result)

            new_items.append(locs)
            new_blocks.append(newb)

        if len(new_blocks) == 0:
            raise DataError('No numeric types to aggregate')

        # reset the locs in the blocks to correspond to our
        # current ordering
        indexer = np.concatenate(new_items)
        new_items = data.items.take(np.sort(indexer))

        if len(deleted_items):

            # we need to adjust the indexer to account for the
            # items we have removed
            # really should be done in internals :<

            deleted = np.concatenate(deleted_items)
            ai = np.arange(len(data))
            mask = np.zeros(len(data))
            mask[deleted] = 1
            indexer = (ai - mask.cumsum())[indexer]

        offset = 0
        for b in new_blocks:
            loc = len(b.mgr_locs)
            b.mgr_locs = indexer[offset:(offset + loc)]
            offset += loc

        return new_items, new_blocks

    def _get_data_to_aggregate(self):
        obj = self._obj_with_exclusions
        if self.axis == 0:
            return obj.swapaxes(0, 1)._data, 1
        else:
            return obj._data, self.axis

    def _post_process_cython_aggregate(self, obj):
        # undoing kludge from below
        if self.axis == 0:
            obj = obj.swapaxes(0, 1)
        return obj

    def aggregate(self, arg, *args, **kwargs):

        _level = kwargs.pop('_level', None)
        result, how = self._aggregate(arg, _level=_level, *args, **kwargs)
        if how is None:
            return result

        if result is None:

            # grouper specific aggregations
            if self.grouper.nkeys > 1:
                return self._python_agg_general(arg, *args, **kwargs)
            else:

                # try to treat as if we are passing a list
                try:
                    assert not args and not kwargs
                    result = self._aggregate_multiple_funcs(
                        [arg], _level=_level, _axis=self.axis)

                    result.columns = Index(
                        result.columns.levels[0],
                        name=self._selected_obj.columns.name)

                    if isinstance(self.obj, SparseDataFrame):
                        # Backwards compat for groupby.agg() with sparse
                        # values. concat no longer converts DataFrame[Sparse]
                        # to SparseDataFrame, so we do it here.
                        result = SparseDataFrame(result._data)
                except Exception:
                    result = self._aggregate_generic(arg, *args, **kwargs)

        if not self.as_index:
            self._insert_inaxis_grouper_inplace(result)
            result.index = np.arange(len(result))

        return result._convert(datetime=True)

    agg = aggregate

    def _aggregate_generic(self, func, *args, **kwargs):
        if self.grouper.nkeys != 1:
            raise AssertionError('Number of keys must be 1')

        axis = self.axis
        obj = self._obj_with_exclusions

        result = collections.OrderedDict()
        if axis != obj._info_axis_number:
            try:
                for name, data in self:
                    result[name] = self._try_cast(func(data, *args, **kwargs),
                                                  data)
            except Exception:
                return self._aggregate_item_by_item(func, *args, **kwargs)
        else:
            for name in self.indices:
                try:
                    data = self.get_group(name, obj=obj)
                    result[name] = self._try_cast(func(data, *args, **kwargs),
                                                  data)
                except Exception:
                    wrapper = lambda x: func(x, *args, **kwargs)
                    result[name] = data.apply(wrapper, axis=axis)

        return self._wrap_generic_output(result, obj)

    def _wrap_aggregated_output(self, output, names=None):
        raise AbstractMethodError(self)

    def _aggregate_item_by_item(self, func, *args, **kwargs):
        # only for axis==0

        obj = self._obj_with_exclusions
        result = collections.OrderedDict()
        cannot_agg = []
        errors = None
        for item in obj:
            try:
                data = obj[item]
                colg = SeriesGroupBy(data, selection=item,
                                     grouper=self.grouper)
                result[item] = self._try_cast(
                    colg.aggregate(func, *args, **kwargs), data)
            except ValueError:
                cannot_agg.append(item)
                continue
            except TypeError as e:
                cannot_agg.append(item)
                errors = e
                continue

        result_columns = obj.columns
        if cannot_agg:
            result_columns = result_columns.drop(cannot_agg)

            # GH6337
            if not len(result_columns) and errors is not None:
                raise errors

        return DataFrame(result, columns=result_columns)

    def _decide_output_index(self, output, labels):
        if len(output) == len(labels):
            output_keys = labels
        else:
            output_keys = sorted(output)
            try:
                output_keys.sort()
            except Exception:  # pragma: no cover
                pass

            if isinstance(labels, MultiIndex):
                output_keys = MultiIndex.from_tuples(output_keys,
                                                     names=labels.names)

        return output_keys

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        from pandas.core.index import _all_indexes_same
        from pandas.core.tools.numeric import to_numeric

        if len(keys) == 0:
            return DataFrame(index=keys)

        key_names = self.grouper.names

        # GH12824.
        def first_not_none(values):
            try:
                return next(com._not_none(*values))
            except StopIteration:
                return None

        v = first_not_none(values)

        if v is None:
            # GH9684. If all values are None, then this will throw an error.
            # We'd prefer it return an empty dataframe.
            return DataFrame()
        elif isinstance(v, DataFrame):
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        elif self.grouper.groupings is not None:
            if len(self.grouper.groupings) > 1:
                key_index = self.grouper.result_index

            else:
                ping = self.grouper.groupings[0]
                if len(keys) == ping.ngroups:
                    key_index = ping.group_index
                    key_index.name = key_names[0]

                    key_lookup = Index(keys)
                    indexer = key_lookup.get_indexer(key_index)

                    # reorder the values
                    values = [values[i] for i in indexer]
                else:

                    key_index = Index(keys, name=key_names[0])

                # don't use the key indexer
                if not self.as_index:
                    key_index = None

            # make Nones an empty object
            v = first_not_none(values)
            if v is None:
                return DataFrame()
            elif isinstance(v, NDFrame):
                values = [
                    x if x is not None else
                    v._constructor(**v._construct_axes_dict())
                    for x in values
                ]

            v = values[0]

            if isinstance(v, (np.ndarray, Index, Series)):
                if isinstance(v, Series):
                    applied_index = self._selected_obj._get_axis(self.axis)
                    all_indexed_same = _all_indexes_same([
                        x.index for x in values
                    ])
                    singular_series = (len(values) == 1 and
                                       applied_index.nlevels == 1)

                    # GH3596
                    # provide a reduction (Frame -> Series) if groups are
                    # unique
                    if self.squeeze:

                        # assign the name to this series
                        if singular_series:
                            values[0].name = keys[0]

                            # GH2893
                            # we have series in the values array, we want to
                            # produce a series:
                            # if any of the sub-series are not indexed the same
                            # OR we don't have a multi-index and we have only a
                            # single values
                            return self._concat_objects(
                                keys, values, not_indexed_same=not_indexed_same
                            )

                        # still a series
                        # path added as of GH 5545
                        elif all_indexed_same:
                            from pandas.core.reshape.concat import concat
                            return concat(values)

                    if not all_indexed_same:
                        # GH 8467
                        return self._concat_objects(
                            keys, values, not_indexed_same=True,
                        )

                try:
                    if self.axis == 0:
                        # GH6124 if the list of Series have a consistent name,
                        # then propagate that name to the result.
                        index = v.index.copy()
                        if index.name is None:
                            # Only propagate the series name to the result
                            # if all series have a consistent name.  If the
                            # series do not have a consistent name, do
                            # nothing.
                            names = {v.name for v in values}
                            if len(names) == 1:
                                index.name = list(names)[0]

                        # normally use vstack as its faster than concat
                        # and if we have mi-columns
                        if (isinstance(v.index, MultiIndex) or
                                key_index is None or
                                isinstance(key_index, MultiIndex)):
                            stacked_values = np.vstack([
                                np.asarray(v) for v in values
                            ])
                            result = DataFrame(stacked_values, index=key_index,
                                               columns=index)
                        else:
                            # GH5788 instead of stacking; concat gets the
                            # dtypes correct
                            from pandas.core.reshape.concat import concat
                            result = concat(values, keys=key_index,
                                            names=key_index.names,
                                            axis=self.axis).unstack()
                            result.columns = index
                    else:
                        stacked_values = np.vstack([np.asarray(v)
                                                    for v in values])
                        result = DataFrame(stacked_values.T, index=v.index,
                                           columns=key_index)

                except (ValueError, AttributeError):
                    # GH1738: values is list of arrays of unequal lengths fall
                    # through to the outer else caluse
                    return Series(values, index=key_index,
                                  name=self._selection_name)

                # if we have date/time like in the original, then coerce dates
                # as we are stacking can easily have object dtypes here
                so = self._selected_obj
                if (so.ndim == 2 and so.dtypes.apply(is_datetimelike).any()):
                    result = result.apply(
                        lambda x: to_numeric(x, errors='ignore'))
                    date_cols = self._selected_obj.select_dtypes(
                        include=['datetime', 'timedelta']).columns
                    date_cols = date_cols.intersection(result.columns)
                    result[date_cols] = (result[date_cols]
                                         ._convert(datetime=True,
                                                   coerce=True))
                else:
                    result = result._convert(datetime=True)

                return self._reindex_output(result)

            # values are not series or array-like but scalars
            else:
                # only coerce dates if we find at least 1 datetime
                coerce = any(isinstance(x, Timestamp) for x in values)
                # self._selection_name not passed through to Series as the
                # result should not take the name of original selection
                # of columns
                return (Series(values, index=key_index)
                        ._convert(datetime=True,
                                  coerce=coerce))

        else:
            # Handle cases like BinGrouper
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)

    def _transform_general(self, func, *args, **kwargs):
        from pandas.core.reshape.concat import concat

        applied = []
        obj = self._obj_with_exclusions
        gen = self.grouper.get_iterator(obj, axis=self.axis)
        fast_path, slow_path = self._define_paths(func, *args, **kwargs)

        path = None
        for name, group in gen:
            object.__setattr__(group, 'name', name)

            if path is None:
                # Try slow path and fast path.
                try:
                    path, res = self._choose_path(fast_path, slow_path, group)
                except TypeError:
                    return self._transform_item_by_item(obj, fast_path)
                except ValueError:
                    msg = 'transform must return a scalar value for each group'
                    raise ValueError(msg)
            else:
                res = path(group)

            if isinstance(res, Series):

                # we need to broadcast across the
                # other dimension; this will preserve dtypes
                # GH14457
                if not np.prod(group.shape):
                    continue
                elif res.index.is_(obj.index):
                    r = concat([res] * len(group.columns), axis=1)
                    r.columns = group.columns
                    r.index = group.index
                else:
                    r = DataFrame(
                        np.concatenate([res.values] * len(group.index)
                                       ).reshape(group.shape),
                        columns=group.columns, index=group.index)

                applied.append(r)
            else:
                applied.append(res)

        concat_index = obj.columns if self.axis == 0 else obj.index
        concatenated = concat(applied, join_axes=[concat_index],
                              axis=self.axis, verify_integrity=False)
        return self._set_result_index_ordered(concatenated)

    @Substitution(klass='DataFrame', selected='')
    @Appender(_transform_template)
    def transform(self, func, *args, **kwargs):

        # optimized transforms
        func = self._is_cython_func(func) or func
        if isinstance(func, str):
            if func in base.cython_transforms:
                # cythonized transform
                return getattr(self, func)(*args, **kwargs)
            else:
                # cythonized aggregation and merge
                result = getattr(self, func)(*args, **kwargs)
        else:
            return self._transform_general(func, *args, **kwargs)

        # a reduction transform
        if not isinstance(result, DataFrame):
            return self._transform_general(func, *args, **kwargs)

        obj = self._obj_with_exclusions

        # nuiscance columns
        if not result.columns.equals(obj.columns):
            return self._transform_general(func, *args, **kwargs)

        return self._transform_fast(result, obj, func)

    def _transform_fast(self, result, obj, func_nm):
        """
        Fast transform path for aggregations
        """
        # if there were groups with no observations (Categorical only?)
        # try casting data to original dtype
        cast = self._transform_should_cast(func_nm)

        # for each col, reshape to to size of original frame
        # by take operation
        ids, _, ngroup = self.grouper.group_info
        output = []
        for i, _ in enumerate(result.columns):
            res = algorithms.take_1d(result.iloc[:, i].values, ids)
            if cast:
                res = self._try_cast(res, obj.iloc[:, i])
            output.append(res)

        return DataFrame._from_arrays(output, columns=result.columns,
                                      index=obj.index)

    def _define_paths(self, func, *args, **kwargs):
        if isinstance(func, str):
            fast_path = lambda group: getattr(group, func)(*args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: getattr(x, func)(*args, **kwargs), axis=self.axis)
        else:
            fast_path = lambda group: func(group, *args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: func(x, *args, **kwargs), axis=self.axis)
        return fast_path, slow_path

    def _choose_path(self, fast_path, slow_path, group):
        path = slow_path
        res = slow_path(group)

        # if we make it here, test if we can use the fast path
        try:
            res_fast = fast_path(group)

            # verify fast path does not change columns (and names), otherwise
            # its results cannot be joined with those of the slow path
            if res_fast.columns != group.columns:
                return path, res
            # verify numerical equality with the slow path
            if res.shape == res_fast.shape:
                res_r = res.values.ravel()
                res_fast_r = res_fast.values.ravel()
                mask = notna(res_r)
                if (res_r[mask] == res_fast_r[mask]).all():
                    path = fast_path
        except Exception:
            pass
        return path, res

    def _transform_item_by_item(self, obj, wrapper):
        # iterate through columns
        output = {}
        inds = []
        for i, col in enumerate(obj):
            try:
                output[col] = self[col].transform(wrapper)
                inds.append(i)
            except Exception:
                pass

        if len(output) == 0:  # pragma: no cover
            raise TypeError('Transform function invalid for data types')

        columns = obj.columns
        if len(output) < len(obj.columns):
            columns = columns.take(inds)

        return DataFrame(output, index=obj.index, columns=columns)

    def filter(self, func, dropna=True, *args, **kwargs):  # noqa
        """
        Return a copy of a DataFrame excluding elements from groups that
        do not satisfy the boolean criterion specified by func.

        Parameters
        ----------
        f : function
            Function to apply to each subframe. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            if False, groups that evaluate False are filled with NaNs.

        Returns
        -------
        filtered : DataFrame

        Notes
        -----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Examples
        --------
        >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
        ...                           'foo', 'bar'],
        ...                    'B' : [1, 2, 3, 4, 5, 6],
        ...                    'C' : [2.0, 5., 8., 1., 2., 9.]})
        >>> grouped = df.groupby('A')
        >>> grouped.filter(lambda x: x['B'].mean() > 3.)
             A  B    C
        1  bar  2  5.0
        3  bar  4  1.0
        5  bar  6  9.0
        """

        indices = []

        obj = self._selected_obj
        gen = self.grouper.get_iterator(obj, axis=self.axis)

        for name, group in gen:
            object.__setattr__(group, 'name', name)

            res = func(group, *args, **kwargs)

            try:
                res = res.squeeze()
            except AttributeError:  # allow e.g., scalars and frames to pass
                pass

            # interpret the result of the filter
            if is_bool(res) or (is_scalar(res) and isna(res)):
                if res and notna(res):
                    indices.append(self._get_index(name))
            else:
                # non scalars aren't allowed
                raise TypeError("filter function returned a %s, "
                                "but expected a scalar bool" %
                                type(res).__name__)

        return self._apply_filter(indices, dropna)


class SeriesGroupBy(GroupBy):
    #
    # Make class defs of attributes on SeriesGroupBy whitelist

    _apply_whitelist = base.series_apply_whitelist
    for _def_str in base.whitelist_method_generator(
            GroupBy, Series, _apply_whitelist):
        exec(_def_str)

    @property
    def _selection_name(self):
        """
        since we are a series, we by definition only have
        a single name, but may be the result of a selection or
        the name of our object
        """
        if self._selection is None:
            return self.obj.name
        else:
            return self._selection

    _agg_see_also_doc = dedent("""
    See Also
    --------
    pandas.Series.groupby.apply
    pandas.Series.groupby.transform
    pandas.Series.aggregate
    """)

    _agg_examples_doc = dedent("""
    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4])

    >>> s
    0    1
    1    2
    2    3
    3    4
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).min()
    1    1
    2    3
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).agg('min')
    1    1
    2    3
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).agg(['min', 'max'])
       min  max
    1    1    2
    2    3    4
    """)

    @Appender(_apply_docs['template']
              .format(input='series',
                      examples=_apply_docs['series_examples']))
    def apply(self, func, *args, **kwargs):
        return super(SeriesGroupBy, self).apply(func, *args, **kwargs)

    @Substitution(see_also=_agg_see_also_doc,
                  examples=_agg_examples_doc,
                  versionadded='',
                  klass='Series',
                  axis='')
    @Appender(_shared_docs['aggregate'])
    def aggregate(self, func_or_funcs, *args, **kwargs):
        _level = kwargs.pop('_level', None)
        if isinstance(func_or_funcs, str):
            return getattr(self, func_or_funcs)(*args, **kwargs)

        if isinstance(func_or_funcs, compat.Iterable):
            # Catch instances of lists / tuples
            # but not the class list / tuple itself.
            ret = self._aggregate_multiple_funcs(func_or_funcs,
                                                 (_level or 0) + 1)
        else:
            cyfunc = self._is_cython_func(func_or_funcs)
            if cyfunc and not args and not kwargs:
                return getattr(self, cyfunc)()

            if self.grouper.nkeys > 1:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)

            try:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)
            except Exception:
                result = self._aggregate_named(func_or_funcs, *args, **kwargs)

            index = Index(sorted(result), name=self.grouper.names[0])
            ret = Series(result, index=index)

        if not self.as_index:  # pragma: no cover
            print('Warning, ignoring as_index=True')

        # _level handled at higher
        if not _level and isinstance(ret, dict):
            from pandas import concat
            ret = concat(ret, axis=1)
        return ret

    agg = aggregate

    def _aggregate_multiple_funcs(self, arg, _level):
        if isinstance(arg, dict):

            # show the deprecation, but only if we
            # have not shown a higher level one
            # GH 15931
            if isinstance(self._selected_obj, Series) and _level <= 1:
                warnings.warn(
                    ("using a dict on a Series for aggregation\n"
                     "is deprecated and will be removed in a future "
                     "version"),
                    FutureWarning, stacklevel=3)

            columns = list(arg.keys())
            arg = list(arg.items())
        elif any(isinstance(x, (tuple, list)) for x in arg):
            arg = [(x, x) if not isinstance(x, (tuple, list)) else x
                   for x in arg]

            # indicated column order
            columns = lzip(*arg)[0]
        else:
            # list of functions / function names
            columns = []
            for f in arg:
                if isinstance(f, str):
                    columns.append(f)
                else:
                    # protect against callables without names
                    columns.append(com.get_callable_name(f))
            arg = lzip(columns, arg)

        results = collections.OrderedDict()
        for name, func in arg:
            obj = self
            if name in results:
                raise SpecificationError(
                    'Function names must be unique, found multiple named '
                    '{}'.format(name))

            # reset the cache so that we
            # only include the named selection
            if name in self._selected_obj:
                obj = copy.copy(obj)
                obj._reset_cache()
                obj._selection = name
            results[name] = obj.aggregate(func)

        if any(isinstance(x, DataFrame) for x in compat.itervalues(results)):
            # let higher level handle
            if _level:
                return results

        return DataFrame(results, columns=columns)

    def _wrap_output(self, output, index, names=None):
        """ common agg/transform wrapping logic """
        output = output[self._selection_name]

        if names is not None:
            return DataFrame(output, index=index, columns=names)
        else:
            name = self._selection_name
            if name is None:
                name = self._selected_obj.name
            return Series(output, index=index, name=name)

    def _wrap_aggregated_output(self, output, names=None):
        return self._wrap_output(output=output,
                                 index=self.grouper.result_index,
                                 names=names)

    def _wrap_transformed_output(self, output, names=None):
        return self._wrap_output(output=output,
                                 index=self.obj.index,
                                 names=names)

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            # GH #6265
            return Series([], name=self._selection_name, index=keys)

        def _get_index():
            if self.grouper.nkeys > 1:
                index = MultiIndex.from_tuples(keys, names=self.grouper.names)
            else:
                index = Index(keys, name=self.grouper.names[0])
            return index

        if isinstance(values[0], dict):
            # GH #823
            index = _get_index()
            result = DataFrame(values, index=index).stack()
            result.name = self._selection_name
            return result

        if isinstance(values[0], (Series, dict)):
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        elif isinstance(values[0], DataFrame):
            # possible that Series -> DataFrame by applied function
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        else:
            # GH #6265
            return Series(values, index=_get_index(),
                          name=self._selection_name)

    def _aggregate_named(self, func, *args, **kwargs):
        result = collections.OrderedDict()

        for name, group in self:
            group.name = name
            output = func(group, *args, **kwargs)
            if isinstance(output, (Series, Index, np.ndarray)):
                raise Exception('Must produce aggregated value')
            result[name] = self._try_cast(output, group)

        return result

    @Substitution(klass='Series', selected='A.')
    @Appender(_transform_template)
    def transform(self, func, *args, **kwargs):
        func = self._is_cython_func(func) or func

        # if string function
        if isinstance(func, str):
            if func in base.cython_transforms:
                # cythonized transform
                return getattr(self, func)(*args, **kwargs)
            else:
                # cythonized aggregation and merge
                return self._transform_fast(
                    lambda: getattr(self, func)(*args, **kwargs), func)

        # reg transform
        klass = self._selected_obj.__class__
        results = []
        wrapper = lambda x: func(x, *args, **kwargs)
        for name, group in self:
            object.__setattr__(group, 'name', name)
            res = wrapper(group)

            if hasattr(res, 'values'):
                res = res.values

            indexer = self._get_index(name)
            s = klass(res, indexer)
            results.append(s)

        from pandas.core.reshape.concat import concat
        result = concat(results).sort_index()

        # we will only try to coerce the result type if
        # we have a numeric dtype, as these are *always* udfs
        # the cython take a different path (and casting)
        dtype = self._selected_obj.dtype
        if is_numeric_dtype(dtype):
            result = maybe_downcast_to_dtype(result, dtype)

        result.name = self._selected_obj.name
        result.index = self._selected_obj.index
        return result

    def _transform_fast(self, func, func_nm):
        """
        fast version of transform, only applicable to
        builtin/cythonizable functions
        """
        if isinstance(func, str):
            func = getattr(self, func)

        ids, _, ngroup = self.grouper.group_info
        cast = self._transform_should_cast(func_nm)
        out = algorithms.take_1d(func()._values, ids)
        if cast:
            out = self._try_cast(out, self.obj)
        return Series(out, index=self.obj.index, name=self.obj.name)

    def filter(self, func, dropna=True, *args, **kwargs):  # noqa
        """
        Return a copy of a Series excluding elements from groups that
        do not satisfy the boolean criterion specified by func.

        Parameters
        ----------
        func : function
            To apply to each group. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            if False, groups that evaluate False are filled with NaNs.

        Examples
        --------
        >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
        ...                           'foo', 'bar'],
        ...                    'B' : [1, 2, 3, 4, 5, 6],
        ...                    'C' : [2.0, 5., 8., 1., 2., 9.]})
        >>> grouped = df.groupby('A')
        >>> df.groupby('A').B.filter(lambda x: x.mean() > 3.)
        1    2
        3    4
        5    6
        Name: B, dtype: int64

        Returns
        -------
        filtered : Series
        """
        if isinstance(func, str):
            wrapper = lambda x: getattr(x, func)(*args, **kwargs)
        else:
            wrapper = lambda x: func(x, *args, **kwargs)

        # Interpret np.nan as False.
        def true_and_notna(x, *args, **kwargs):
            b = wrapper(x, *args, **kwargs)
            return b and notna(b)

        try:
            indices = [self._get_index(name) for name, group in self
                       if true_and_notna(group)]
        except ValueError:
            raise TypeError("the filter must return a boolean result")
        except TypeError:
            raise TypeError("the filter must return a boolean result")

        filtered = self._apply_filter(indices, dropna)
        return filtered

    def nunique(self, dropna=True):
        """
        Return number of unique elements in the group.
        """
        ids, _, _ = self.grouper.group_info

        val = self.obj.get_values()

        try:
            sorter = np.lexsort((val, ids))
        except TypeError:  # catches object dtypes
            msg = 'val.dtype must be object, got {}'.format(val.dtype)
            assert val.dtype == object, msg
            val, _ = algorithms.factorize(val, sort=False)
            sorter = np.lexsort((val, ids))
            _isna = lambda a: a == -1
        else:
            _isna = isna

        ids, val = ids[sorter], val[sorter]

        # group boundaries are where group ids change
        # unique observations are where sorted values change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]
        inc = np.r_[1, val[1:] != val[:-1]]

        # 1st item of each group is a new unique observation
        mask = _isna(val)
        if dropna:
            inc[idx] = 1
            inc[mask] = 0
        else:
            inc[mask & np.r_[False, mask[:-1]]] = 0
            inc[idx] = 1

        out = np.add.reduceat(inc, idx).astype('int64', copy=False)
        if len(ids):
            # NaN/NaT group exists if the head of ids is -1,
            # so remove it from res and exclude its index from idx
            if ids[0] == -1:
                res = out[1:]
                idx = idx[np.flatnonzero(idx)]
            else:
                res = out
        else:
            res = out[1:]
        ri = self.grouper.result_index

        # we might have duplications among the bins
        if len(res) != len(ri):
            res, out = np.zeros(len(ri), dtype=out.dtype), res
            res[ids[idx]] = out

        return Series(res,
                      index=ri,
                      name=self._selection_name)

    @Appender(Series.describe.__doc__)
    def describe(self, **kwargs):
        result = self.apply(lambda x: x.describe(**kwargs))
        if self.axis == 1:
            return result.T
        return result.unstack()

    def value_counts(self, normalize=False, sort=True, ascending=False,
                     bins=None, dropna=True):

        from pandas.core.reshape.tile import cut
        from pandas.core.reshape.merge import _get_join_indexers

        if bins is not None and not np.iterable(bins):
            # scalar bins cannot be done at top level
            # in a backward compatible way
            return self.apply(Series.value_counts,
                              normalize=normalize,
                              sort=sort,
                              ascending=ascending,
                              bins=bins)

        ids, _, _ = self.grouper.group_info
        val = self.obj.get_values()

        # groupby removes null keys from groupings
        mask = ids != -1
        ids, val = ids[mask], val[mask]

        if bins is None:
            lab, lev = algorithms.factorize(val, sort=True)
            llab = lambda lab, inc: lab[inc]
        else:

            # lab is a Categorical with categories an IntervalIndex
            lab = cut(Series(val), bins, include_lowest=True)
            lev = lab.cat.categories
            lab = lev.take(lab.cat.codes)
            llab = lambda lab, inc: lab[inc]._multiindex.codes[-1]

        if is_interval_dtype(lab):
            # TODO: should we do this inside II?
            sorter = np.lexsort((lab.left, lab.right, ids))
        else:
            sorter = np.lexsort((lab, ids))

        ids, lab = ids[sorter], lab[sorter]

        # group boundaries are where group ids change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]

        # new values are where sorted labels change
        lchanges = llab(lab, slice(1, None)) != llab(lab, slice(None, -1))
        inc = np.r_[True, lchanges]
        inc[idx] = True  # group boundaries are also new values
        out = np.diff(np.nonzero(np.r_[inc, True])[0])  # value counts

        # num. of times each group should be repeated
        rep = partial(np.repeat, repeats=np.add.reduceat(inc, idx))

        # multi-index components
        labels = list(map(rep, self.grouper.recons_labels)) + [llab(lab, inc)]
        levels = [ping.group_index for ping in self.grouper.groupings] + [lev]
        names = self.grouper.names + [self._selection_name]

        if dropna:
            mask = labels[-1] != -1
            if mask.all():
                dropna = False
            else:
                out, labels = out[mask], [label[mask] for label in labels]

        if normalize:
            out = out.astype('float')
            d = np.diff(np.r_[idx, len(ids)])
            if dropna:
                m = ids[lab == -1]
                np.add.at(d, m, -1)
                acc = rep(d)[mask]
            else:
                acc = rep(d)
            out /= acc

        if sort and bins is None:
            cat = ids[inc][mask] if dropna else ids[inc]
            sorter = np.lexsort((out if ascending else -out, cat))
            out, labels[-1] = out[sorter], labels[-1][sorter]

        if bins is None:
            mi = MultiIndex(levels=levels, codes=labels, names=names,
                            verify_integrity=False)

            if is_integer_dtype(out):
                out = ensure_int64(out)
            return Series(out, index=mi, name=self._selection_name)

        # for compat. with libgroupby.value_counts need to ensure every
        # bin is present at every index level, null filled with zeros
        diff = np.zeros(len(out), dtype='bool')
        for lab in labels[:-1]:
            diff |= np.r_[True, lab[1:] != lab[:-1]]

        ncat, nbin = diff.sum(), len(levels[-1])

        left = [np.repeat(np.arange(ncat), nbin),
                np.tile(np.arange(nbin), ncat)]

        right = [diff.cumsum() - 1, labels[-1]]

        _, idx = _get_join_indexers(left, right, sort=False, how='left')
        out = np.where(idx != -1, out[idx], 0)

        if sort:
            sorter = np.lexsort((out if ascending else -out, left[0]))
            out, left[-1] = out[sorter], left[-1][sorter]

        # build the multi-index w/ full levels
        codes = list(map(lambda lab: np.repeat(lab[diff], nbin), labels[:-1]))
        codes.append(left[-1])

        mi = MultiIndex(levels=levels, codes=codes, names=names,
                        verify_integrity=False)

        if is_integer_dtype(out):
            out = ensure_int64(out)
        return Series(out, index=mi, name=self._selection_name)

    def count(self):
        """ Compute count of group, excluding missing values """
        ids, _, ngroups = self.grouper.group_info
        val = self.obj.get_values()

        mask = (ids != -1) & ~isna(val)
        ids = ensure_platform_int(ids)
        minlength = ngroups or 0
        out = np.bincount(ids[mask], minlength=minlength)

        return Series(out,
                      index=self.grouper.result_index,
                      name=self._selection_name,
                      dtype='int64')

    def _apply_to_column_groupbys(self, func):
        """ return a pass thru """
        return func(self)

    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None):
        """Calcuate pct_change of each value to previous entry in group"""
        # TODO: Remove this conditional when #23918 is fixed
        if freq:
            return self.apply(lambda x: x.pct_change(periods=periods,
                                                     fill_method=fill_method,
                                                     limit=limit, freq=freq))
        filled = getattr(self, fill_method)(limit=limit)
        fill_grp = filled.groupby(self.grouper.labels)
        shifted = fill_grp.shift(periods=periods, freq=freq)

        return (filled / shifted) - 1


class DataFrameGroupBy(NDFrameGroupBy):

    _apply_whitelist = base.dataframe_apply_whitelist

    #
    # Make class defs of attributes on DataFrameGroupBy whitelist.
    for _def_str in base.whitelist_method_generator(
            GroupBy, DataFrame, _apply_whitelist):
        exec(_def_str)

    _block_agg_axis = 1

    _agg_see_also_doc = dedent("""
    See Also
    --------
    pandas.DataFrame.groupby.apply
    pandas.DataFrame.groupby.transform
    pandas.DataFrame.aggregate
    """)

    _agg_examples_doc = dedent("""
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
    """)

    @Substitution(see_also=_agg_see_also_doc,
                  examples=_agg_examples_doc,
                  versionadded='',
                  klass='DataFrame',
                  axis='')
    @Appender(_shared_docs['aggregate'])
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
                                    as_index=self.as_index,
                                    observed=self.observed)
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
        in_axis_grps = ((i, ping.name) for (i, ping)
                        in enumerate(groupings) if ping.in_axis)
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

        counter = partial(
            lib.count_level_2d, labels=ids, max_bin=ngroups, axis=1)
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
            results.columns.names = obj.columns.names

        if not self.as_index:
            results.index = ibase.default_index(len(results))
        return results

    boxplot = boxplot_frame_groupby
