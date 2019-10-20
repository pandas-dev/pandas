from collections import OrderedDict
import functools
from typing import Any, Sequence

from pandas.compat import PY36

from pandas.core.dtypes.common import is_dict_like, is_list_like

import pandas.core.common as com
from pandas.core.index import Index


def _is_multi_agg_with_relabel(**kwargs):
    """
    Check whether kwargs passed to .agg look like multi-agg with relabeling.

    Parameters
    ----------
    **kwargs : dict

    Returns
    -------
    bool

    Examples
    --------
    >>> _is_multi_agg_with_relabel(a='max')
    False
    >>> _is_multi_agg_with_relabel(a_max=('a', 'max'),
    ...                            a_min=('a', 'min'))
    True
    >>> _is_multi_agg_with_relabel()
    False
    """
    return all(isinstance(v, tuple) and len(v) == 2 for v in kwargs.values()) and kwargs


def _normalize_keyword_aggregation(kwargs):
    """
    Normalize user-provided "named aggregation" kwargs.

    Transforms from the new ``Dict[str, NamedAgg]`` style kwargs
    to the old OrderedDict[str, List[scalar]]].

    Parameters
    ----------
    kwargs : dict

    Returns
    -------
    aggspec : dict
        The transformed kwargs.
    columns : List[str]
        The user-provided keys.
    col_idx_order : List[int]
        List of columns indices.

    Examples
    --------
    >>> _normalize_keyword_aggregation({'output': ('input', 'sum')})
    (OrderedDict([('input', ['sum'])]), ('output',), [('input', 'sum')])
    """
    if not PY36:
        kwargs = OrderedDict(sorted(kwargs.items()))

    # Normalize the aggregation functions as Dict[column, List[func]],
    # process normally, then fixup the names.
    # TODO(Py35): When we drop python 3.5, change this to
    # defaultdict(list)
    # TODO: aggspec type: typing.OrderedDict[str, List[AggScalar]]
    # May be hitting https://github.com/python/mypy/issues/5958
    # saying it doesn't have an attribute __name__
    aggspec = OrderedDict()
    order = []
    columns, pairs = list(zip(*kwargs.items()))

    for name, (column, aggfunc) in zip(columns, pairs):
        if column in aggspec:
            aggspec[column].append(aggfunc)
        else:
            aggspec[column] = [aggfunc]
        order.append((column, com.get_callable_name(aggfunc) or aggfunc))

    # uniquify aggfunc name if duplicated in order list
    uniquified_order = _make_unique(order)

    # GH 25719, due to aggspec will change the order of assigned columns in aggregation
    # uniquified_aggspec will store uniquified order list and will compare it with order
    # based on index
    aggspec_order = [
        (column, com.get_callable_name(aggfunc) or aggfunc)
        for column, aggfuncs in aggspec.items()
        for aggfunc in aggfuncs
    ]
    uniquified_aggspec = _make_unique(aggspec_order)

    # get the new indice of columns by comparison
    col_idx_order = Index(uniquified_aggspec).get_indexer(uniquified_order)
    return aggspec, columns, col_idx_order


def _managle_lambda_list(aggfuncs: Sequence[Any]) -> Sequence[Any]:
    """
    Possibly mangle a list of aggfuncs.

    Parameters
    ----------
    aggfuncs : Sequence

    Returns
    -------
    mangled: list-like
        A new AggSpec sequence, where lambdas have been converted
        to have unique names.

    Notes
    -----
    If just one aggfunc is passed, the name will not be mangled.
    """
    if len(aggfuncs) <= 1:
        # don't mangle for .agg([lambda x: .])
        return aggfuncs
    i = 0
    mangled_aggfuncs = []
    for aggfunc in aggfuncs:
        if com.get_callable_name(aggfunc) == "<lambda>":
            aggfunc = functools.partial(aggfunc)
            aggfunc.__name__ = "<lambda_{}>".format(i)
            i += 1
        mangled_aggfuncs.append(aggfunc)

    return mangled_aggfuncs


def _maybe_mangle_lambdas(agg_spec: Any) -> Any:
    """
    Make new lambdas with unique names.

    Parameters
    ----------
    agg_spec : Any
        An argument to GroupBy.agg.
        Non-dict-like `agg_spec` are pass through as is.
        For dict-like `agg_spec` a new spec is returned
        with name-mangled lambdas.

    Returns
    -------
    mangled : Any
        Same type as the input.

    Examples
    --------
    >>> _maybe_mangle_lambdas('sum')
    'sum'

    >>> _maybe_mangle_lambdas([lambda: 1, lambda: 2])  # doctest: +SKIP
    [<function __main__.<lambda_0>,
     <function pandas...._make_lambda.<locals>.f(*args, **kwargs)>]
    """
    is_dict = is_dict_like(agg_spec)
    if not (is_dict or is_list_like(agg_spec)):
        return agg_spec
    mangled_aggspec = type(agg_spec)()  # dict or OrderdDict

    if is_dict:
        for key, aggfuncs in agg_spec.items():
            if is_list_like(aggfuncs) and not is_dict_like(aggfuncs):
                mangled_aggfuncs = _managle_lambda_list(aggfuncs)
            else:
                mangled_aggfuncs = aggfuncs

            mangled_aggspec[key] = mangled_aggfuncs
    else:
        mangled_aggspec = _managle_lambda_list(agg_spec)

    return mangled_aggspec


def _make_unique(seq):
    """Uniquify aggfunc name of the pairs in the order list

    Examples:
    --------
    >>> _make_unique([('a', '<lambda>'), ('a', '<lambda>'), ('b', '<lambda>')])
    [('a', '<lambda>_0'), ('a', '<lambda>_1'), ('b', '<lambda>')]
    """
    return [
        (pair[0], "_".join([pair[1], str(seq[:i].count(pair))]))
        if seq.count(pair) > 1
        else pair
        for i, pair in enumerate(seq)
    ]
