import inspect
import numpy as np
from pandas import Series

# pylint: disable=C0111
# pylint: disable=W0163


# boilerplate -----------------------------------------------------------------

def _get_methods(obj, names_only=True):
    members = inspect.getmembers(obj)
    methods = []
    for member in members:
        if inspect.ismethod(member[1]):
            if names_only:
                methods.append(member[0])
            else:
                methods.append(member)
    return methods


def smoke(name):
    def _smoke(factory):
        getattr(factory(), name)()
    return _smoke


def series_smoker(smoker):
    """
    Decorator to register a function as a Series method smoker.

    Parameters
    ----------
    smoker: Name of a smoker or a smoker function.
        In the later case, the smoker function name should be of the form
        `smoke_<smoker_name>`. Smoker name will be derived from the function
        name.

    Example
    -------
    @series_smoker('something')
    def method_that_smokes(factory):
        pass

    @series_smoker
    def smoke_something(factory):
        pass
    """
    def register(smoker):
        if name in smokers:
            other = smokers[name]
            if isinstance(other, list):
                other.append(smoker)
            else:
                smokers[name] = [other, smoker]
        else:
            smokers[name] = smoker
        return smoker

    if hasattr(smoker, '__call__'):
        if smoker.__name__.startswith('smoke_'):
            name = smoker.__name__[6:]
        else:
            raise ValueError(("Can not derive smoker name from %s."
                % smoker.__name__))
        return register(smoker)
    else:
        name = smoker
        return register

# Series factories ------------------------------------------------------------

# A series factory is an argless callable object returning a Series object.
# Every smoker will get every factory as input. See also `test_smoke_series`.
# TODO write some interesting factories, for now just lambda stuff.

factories = [lambda: Series(),
             lambda: Series(['foo', 'bar', 'baz']),
             lambda: Series(np.random.randint(-1e9, 1e9, 51)),
             lambda: Series(np.random.randn(10))]

# Series smokers --------------------------------------------------------------

# smokers = dict of (Series method name, smoker function)
# Every method of Series should be in the smokers dict. If not test
# :func:`test_all_smoked` will fail. The decorator :func:`series_smoker` can be
# used to register a smoker.  If there is no smoker for a given Series method,
# register it with None as smoker function.

smokers = {'__repr__': smoke('__repr__'),
           '_agg_by_level': None,
           '_binop': None,
           '_check_bool_indexer': None,
           '_get_axis': None,
           '_get_axis_name': None,
           '_get_axis_number': None,
           '_get_repr': None,
           '_get_values': None,
           '_get_values_tuple': None,
           '_get_with': None,
           '_reindex_indexer': None,
           '_repr_footer': None,
           '_set_labels': None,
           '_set_values': None,
           '_set_with': None,
           '_tidy_repr': None,
           'value_counts': smoke('value_counts')}

@series_smoker('__add__')
def smoke_dunder_add(factory):
    pass
    # check if data is numeric and do
    #for factory_other in factories:
    #    factory() + factory_other()

@series_smoker('__and__')
def smoke_dunder_and(factory):
    pass

@series_smoker('__array_finalize__')
def smoke_dunder_array_finalize(factory):
    pass

@series_smoker('__contains__')
def smoke_dunder_contains(factory):
    s = factory()
    for key in [1, 1.0, 'Abcd', u'\u03B1', True, None, ('foo', -3.14)]:
        key in s
    for key in s.values:
        key in s

@series_smoker('__div__')
def smoke_dunder_div(factory):
    pass

@series_smoker('__eq__')
def smoke_dunder_eq(factory):
    for factory_other in factories:
        factory() == factory_other()

@series_smoker('__floordiv__')
def smoke_dunder_floordiv(factory):
    pass

@series_smoker('__ge__')
def smoke_dunder_ge(factory):
    pass

@series_smoker('__getitem__')
def smoke_dunder_getitem(factory):
    pass

@series_smoker('__getslice__')
def smoke_dunder_getslice(factory):
    pass

@series_smoker('__gt__')
def smoke_dunder_gt(factory):
    pass

@series_smoker('__hash__')
def smoke_dunder_hash(factory):
    pass

@series_smoker('__iadd__')
def smoke_dunder_iadd(factory):
    pass

@series_smoker('__idiv__')
def smoke_dunder_idiv(factory):
    pass

@series_smoker('__ifloordiv__')
def smoke_dunder_ifloordiv(factory):
    pass

@series_smoker('__imul__')
def smoke_dunder_imul(factory):
    pass

@series_smoker('__init__')
def smoke_dunder_init(factory):
    pass

@series_smoker('__ipow__')
def smoke_dunder_ipow(factory):
    pass

@series_smoker('__isub__')
def smoke_dunder_isub(factory):
    pass

@series_smoker('__iter__')
def smoke_dunder_iter(factory):
    pass

@series_smoker('__itruediv__')
def smoke_dunder_itruediv(factory):
    pass

@series_smoker('__le__')
def smoke_dunder_le(factory):
    pass

@series_smoker('__lt__')
def smoke_dunder_lt(factory):
    pass

@series_smoker('__mul__')
def smoke_dunder_mul(factory):
    pass

@series_smoker('__ne__')
def smoke_dunder_ne(factory):
    pass

@series_smoker('__or__')
def smoke_dunder_or(factory):
    pass

@series_smoker('__pow__')
def smoke_dunder_pow(factory):
    pass

@series_smoker('__radd__')
def smoke_dunder_radd(factory):
    pass

@series_smoker('__rdiv__')
def smoke_dunder_rdiv(factory):
    pass

@series_smoker('__reduce__')
def smoke_dunder_reduce(factory):
    pass

@series_smoker('__rfloordiv__')
def smoke_dunder_rfloordiv(factory):
    pass

@series_smoker('__rmul__')
def smoke_dunder_rmul(factory):
    pass

@series_smoker('__rpow__')
def smoke_dunder_rpow(factory):
    pass

@series_smoker('__rsub__')
def smoke_dunder_rsub(factory):
    pass

@series_smoker('__rtruediv__')
def smoke_dunder_rtruediv(factory):
    pass

@series_smoker('__setitem__')
def smoke_dunder_setitem(factory):
    pass

@series_smoker('__setslice__')
def smoke_dunder_setslice(factory):
    pass

@series_smoker('__setstate__')
def smoke_dunder_setstate(factory):
    pass

@series_smoker('__str__')
def smoke_dunder_str(factory):
    pass

@series_smoker('__sub__')
def smoke_dunder_sub(factory):
    pass

@series_smoker('__truediv__')
def smoke_dunder_truediv(factory):
    pass

@series_smoker('__xor__')
def smoke_dunder_xor(factory):
    pass

@series_smoker
def smoke_abs(factory):
    pass

@series_smoker
def smoke_add(factory):
    pass

@series_smoker
def smoke_align(factory):
    pass

@series_smoker
def smoke_all(factory):
    pass

@series_smoker
def smoke_any(factory):
    pass

@series_smoker
def smoke_append(factory):
    pass

@series_smoker
def smoke_apply(factory):
    pass

@series_smoker
def smoke_argsort(factory):
    pass

@series_smoker
def smoke_asfreq(factory):
    pass

@series_smoker
def smoke_asof(factory):
    pass

@series_smoker
def smoke_astype(factory):
    pass

@series_smoker
def smoke_autocorr(factory):
    pass

@series_smoker
def smoke_between(factory):
    pass

@series_smoker
def smoke_clip(factory):
    pass

@series_smoker
def smoke_clip_lower(factory):
    pass

@series_smoker
def smoke_clip_upper(factory):
    pass

@series_smoker
def smoke_combine(factory):
    pass

@series_smoker
def smoke_combine_first(factory):
    pass

@series_smoker
def smoke_copy(factory):
    pass

@series_smoker
def smoke_corr(factory):
    pass

@series_smoker
def smoke_count(factory):
    pass

@series_smoker
def smoke_cov(factory):
    pass

@series_smoker
def smoke_cummax(factory):
    pass

@series_smoker
def smoke_cummin(factory):
    pass

@series_smoker
def smoke_cumprod(factory):
    pass

@series_smoker
def smoke_cumsum(factory):
    pass

@series_smoker
def smoke_describe(factory):
    pass

@series_smoker
def smoke_diff(factory):
    pass

@series_smoker
def smoke_div(factory):
    pass

@series_smoker
def smoke_drop(factory):
    pass

@series_smoker
def smoke_dropna(factory):
    pass

@series_smoker
def smoke_fillna(factory):
    pass

@series_smoker
def smoke_first(factory):
    pass

@series_smoker
def smoke_first_valid_index(factory):
    pass

@series_smoker
def smoke_from_csv(factory):
    pass

@series_smoker
def smoke_get(factory):
    pass

@series_smoker
def smoke_get_value(factory):
    pass

@series_smoker
def smoke_groupby(factory):
    pass

@series_smoker
def smoke_head(factory):
    pass

@series_smoker
def smoke_hist(factory):
    pass

@series_smoker
def smoke_idxmax(factory):
    pass

@series_smoker
def smoke_idxmin(factory):
    pass

@series_smoker
def smoke_iget(factory):
    pass

@series_smoker
def smoke_iget_value(factory):
    pass

@series_smoker
def smoke_interpolate(factory):
    pass

@series_smoker
def smoke_irow(factory):
    pass

@series_smoker
def smoke_isin(factory):
    pass

@series_smoker
def smoke_isnull(factory):
    pass

@series_smoker
def smoke_iteritems(factory):
    pass

@series_smoker
def smoke_iterkv(factory):
    pass

@series_smoker
def smoke_keys(factory):
    pass

@series_smoker
def smoke_kurt(factory):
    pass

@series_smoker
def smoke_last(factory):
    pass

@series_smoker
def smoke_last_valid_index(factory):
    pass

@series_smoker
def smoke_load(factory):
    pass

@series_smoker
def smoke_mad(factory):
    pass

@series_smoker
def smoke_map(factory):
    pass

@series_smoker
def smoke_max(factory):
    pass

@series_smoker
def smoke_mean(factory):
    pass

@series_smoker
def smoke_median(factory):
    pass

@series_smoker
def smoke_min(factory):
    pass

@series_smoker
def smoke_mul(factory):
    pass

@series_smoker
def smoke_notnull(factory):
    pass

@series_smoker
def smoke_nunique(factory):
    pass

@series_smoker
def smoke_order(factory):
    pass

@series_smoker
def smoke_pct_change(factory):
    pass

@series_smoker
def smoke_plot(factory):
    pass

@series_smoker
def smoke_prod(factory):
    pass

@series_smoker
def smoke_ptp(factory):
    pass

@series_smoker
def smoke_quantile(factory):
    pass

@series_smoker
def smoke_rank(factory):
    pass

@series_smoker
def smoke_reindex(factory):
    pass

@series_smoker
def smoke_reindex_like(factory):
    pass

@series_smoker
def smoke_rename(factory):
    pass

@series_smoker
def smoke_reorder_levels(factory):
    pass

@series_smoker
def smoke_repeat(factory):
    pass

@series_smoker
def smoke_replace(factory):
    pass

@series_smoker
def smoke_resample(factory):
    pass

@series_smoker
def smoke_reset_index(factory):
    pass

@series_smoker
def smoke_reshape(factory):
    pass

@series_smoker
def smoke_round(factory):
    pass

@series_smoker
def smoke_save(factory):
    pass

@series_smoker
def smoke_select(factory):
    pass

@series_smoker
def smoke_set_value(factory):
    pass

@series_smoker
def smoke_shift(factory):
    pass

@series_smoker
def smoke_skew(factory):
    pass

@series_smoker
def smoke_sort(factory):
    pass

@series_smoker
def smoke_sort_index(factory):
    pass

@series_smoker
def smoke_sortlevel(factory):
    pass

@series_smoker
def smoke_std(factory):
    pass

@series_smoker
def smoke_sub(factory):
    pass

@series_smoker
def smoke_sum(factory):
    pass

@series_smoker
def smoke_swaplevel(factory):
    pass

@series_smoker
def smoke_tail(factory):
    pass

@series_smoker
def smoke_take(factory):
    pass

@series_smoker
def smoke_to_csv(factory):
    pass

@series_smoker
def smoke_to_dict(factory):
    pass

@series_smoker
def smoke_to_sparse(factory):
    pass

@series_smoker
def smoke_to_string(factory):
    pass

@series_smoker
def smoke_truncate(factory):
    pass

@series_smoker
def smoke_tshift(factory):
    pass

@series_smoker
def smoke_tz_convert(factory):
    pass

@series_smoker
def smoke_tz_localize(factory):
    pass

@series_smoker
def smoke_unique(factory):
    pass

@series_smoker
def smoke_unstack(factory):
    pass

@series_smoker
def smoke_update(factory):
    pass

@series_smoker
def smoke_valid(factory):
    pass

@series_smoker
def smoke_var(factory):
    pass


# tests -----------------------------------------------------------------------

def test_smoke_series():
    def _light(name, factory_id):
        smokers[name](factories[factory_id])
    for name, smoker in smokers.iteritems():
        if smoker is not None:
            for factory_id in range(len(factories)):
                yield _light, name, factory_id

def test_all_smoke():
    """
    Test if every Series method has a smoker.
    """
    methods = _get_methods(Series)
    unsmoked = []
    for method in methods:
        if not method in smokers:
            unsmoked.append(method)
    # FIXME fail if unsmoked is not empty
    print unsmoked
