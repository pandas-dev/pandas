import collections
import warnings
import copy
from textwrap import dedent
from functools import partial

import numpy as np

from pandas.compat import lzip, map
from pandas import compat
import pandas.core.common as com
from pandas.core.series import Series
from pandas.core.generic import _shared_docs
from pandas.core.groupby.groupby import (
    GroupBy, _apply_docs, _transform_template)
from pandas.core.groupby import base
from pandas.util._decorators import Substitution, Appender
from pandas.core.dtypes.common import (
    is_numeric_dtype,
    is_integer_dtype,
    is_interval_dtype,
    _ensure_platform_int,
    _ensure_int64)
from pandas.core.dtypes.missing import isna, notna
from pandas.core.index import Index, MultiIndex
import pandas.core.algorithms as algorithms
from pandas.core.frame import DataFrame
from pandas.core.dtypes.cast import maybe_downcast_to_dtype
from pandas.core.base import SpecificationError


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

    _agg_doc = dedent("""
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

    See also
    --------
    pandas.Series.groupby.apply
    pandas.Series.groupby.transform
    pandas.Series.aggregate

    """)

    @Appender(_apply_docs['template']
              .format(input='series',
                      examples=_apply_docs['series_examples']))
    def apply(self, func, *args, **kwargs):
        return super(SeriesGroupBy, self).apply(func, *args, **kwargs)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        klass='Series',
        versionadded='',
        axis=''))
    def aggregate(self, func_or_funcs, *args, **kwargs):
        _level = kwargs.pop('_level', None)
        if isinstance(func_or_funcs, compat.string_types):
            return getattr(self, func_or_funcs)(*args, **kwargs)

        if isinstance(func_or_funcs, collections.Iterable):
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
                if isinstance(f, compat.string_types):
                    columns.append(f)
                else:
                    # protect against callables without names
                    columns.append(com._get_callable_name(f))
            arg = lzip(columns, arg)

        results = {}
        for name, func in arg:
            obj = self
            if name in results:
                raise SpecificationError('Function names must be unique, '
                                         'found multiple named %s' % name)

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
        result = {}

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
        if isinstance(func, compat.string_types):
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
        if isinstance(func, compat.string_types):
            func = getattr(self, func)

        ids, _, ngroup = self.grouper.group_info
        cast = self._transform_should_cast(func_nm)
        out = algorithms.take_1d(func().values, ids)
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
        if isinstance(func, compat.string_types):
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
        """ Returns number of unique elements in the group """
        ids, _, _ = self.grouper.group_info

        val = self.obj.get_values()

        try:
            sorter = np.lexsort((val, ids))
        except TypeError:  # catches object dtypes
            assert val.dtype == object, \
                'val.dtype must be object, got %s' % val.dtype
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
            llab = lambda lab, inc: lab[inc]._multiindex.labels[-1]

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
            mi = MultiIndex(levels=levels, labels=labels, names=names,
                            verify_integrity=False)

            if is_integer_dtype(out):
                out = _ensure_int64(out)
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
        labels = list(map(lambda lab: np.repeat(lab[diff], nbin), labels[:-1]))
        labels.append(left[-1])

        mi = MultiIndex(levels=levels, labels=labels, names=names,
                        verify_integrity=False)

        if is_integer_dtype(out):
            out = _ensure_int64(out)
        return Series(out, index=mi, name=self._selection_name)

    def count(self):
        """ Compute count of group, excluding missing values """
        ids, _, ngroups = self.grouper.group_info
        val = self.obj.get_values()

        mask = (ids != -1) & ~isna(val)
        ids = _ensure_platform_int(ids)
        out = np.bincount(ids[mask], minlength=ngroups or 0)

        return Series(out,
                      index=self.grouper.result_index,
                      name=self._selection_name,
                      dtype='int64')

    def _apply_to_column_groupbys(self, func):
        """ return a pass thru """
        return func(self)

    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None):
        """Calculate percent change of each value to previous entry in group"""
        filled = getattr(self, fill_method)(limit=limit)
        shifted = filled.shift(periods=periods, freq=freq)

        return (filled / shifted) - 1
