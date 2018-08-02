"""
Provide basic components for groupby. These defintiions
hold the whitelist of methods that are exposed on the
SeriesGroupBy and the DataFrameGroupBy objects.
"""

import types
from pandas.util._decorators import make_signature
from pandas.core.dtypes.common import is_scalar, is_list_like


class GroupByMixin(object):
    """ provide the groupby facilities to the mixed object """

    @staticmethod
    def _dispatch(name, *args, **kwargs):
        """ dispatch to apply """

        def outer(self, *args, **kwargs):
            def f(x):
                x = self._shallow_copy(x, groupby=self._groupby)
                return getattr(x, name)(*args, **kwargs)
            return self._groupby.apply(f)
        outer.__name__ = name
        return outer

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
        # create a new object to prevent aliasing
        if subset is None:
            subset = self.obj

        # we need to make a shallow copy of ourselves
        # with the same groupby
        kwargs = {attr: getattr(self, attr) for attr in self._attributes}
        self = self.__class__(subset,
                              groupby=self._groupby[key],
                              parent=self,
                              **kwargs)
        self._reset_cache()
        if subset.ndim == 2:
            if is_scalar(key) and key in subset or is_list_like(key):
                self._selection = key
        return self


# special case to prevent duplicate plots when catching exceptions when
# forwarding methods from NDFrames
plotting_methods = frozenset(['plot', 'boxplot', 'hist'])

common_apply_whitelist = frozenset([
    'last', 'first',
    'head', 'tail', 'median',
    'mean', 'sum', 'min', 'max',
    'cumcount', 'ngroup',
    'resample',
    'rank', 'quantile',
    'fillna',
    'mad',
    'any', 'all',
    'take',
    'idxmax', 'idxmin',
    'shift', 'tshift',
    'ffill', 'bfill',
    'pct_change', 'skew',
    'corr', 'cov', 'diff',
]) | plotting_methods

series_apply_whitelist = ((common_apply_whitelist |
                           {'nlargest', 'nsmallest',
                            'is_monotonic_increasing',
                            'is_monotonic_decreasing'}) -
                          {'boxplot'}) | frozenset(['dtype', 'unique'])

dataframe_apply_whitelist = ((common_apply_whitelist |
                              frozenset(['dtypes', 'corrwith'])) -
                             {'boxplot'})

cython_transforms = frozenset(['cumprod', 'cumsum', 'shift',
                               'cummin', 'cummax'])

cython_cast_blacklist = frozenset(['rank', 'count', 'size'])


def whitelist_method_generator(base, klass, whitelist):
    """
    Yields all GroupBy member defs for DataFrame/Series names in whitelist.

    Parameters
    ----------
    base : class
        base class
    klass : class
        class where members are defined.
        Should be Series or DataFrame
    whitelist : list
        list of names of klass methods to be constructed

    Returns
    -------
    The generator yields a sequence of strings, each suitable for exec'ing,
    that define implementations of the named methods for DataFrameGroupBy
    or SeriesGroupBy.

    Since we don't want to override methods explicitly defined in the
    base class, any such name is skipped.
    """

    method_wrapper_template = \
        """def %(name)s(%(sig)s) :
    \"""
    %(doc)s
    \"""
    f = %(self)s.__getattr__('%(name)s')
    return f(%(args)s)"""
    property_wrapper_template = \
        """@property
def %(name)s(self) :
    \"""
    %(doc)s
    \"""
    return self.__getattr__('%(name)s')"""

    for name in whitelist:
        # don't override anything that was explicitly defined
        # in the base class
        if hasattr(base, name):
            continue
        # ugly, but we need the name string itself in the method.
        f = getattr(klass, name)
        doc = f.__doc__
        doc = doc if type(doc) == str else ''
        if isinstance(f, types.MethodType):
            wrapper_template = method_wrapper_template
            decl, args = make_signature(f)
            # pass args by name to f because otherwise
            # GroupBy._make_wrapper won't know whether
            # we passed in an axis parameter.
            args_by_name = ['{0}={0}'.format(arg) for arg in args[1:]]
            params = {'name': name,
                      'doc': doc,
                      'sig': ','.join(decl),
                      'self': args[0],
                      'args': ','.join(args_by_name)}
        else:
            wrapper_template = property_wrapper_template
            params = {'name': name, 'doc': doc}
        yield wrapper_template % params
