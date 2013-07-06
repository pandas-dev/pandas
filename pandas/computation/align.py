from functools import partial, wraps
from itertools import izip

import numpy as np

import pandas as pd
import pandas.core.common as com
from pandas.computation.ops import is_const
from pandas.computation.common import flatten


def _align_core_single_unary_op(term):
    if isinstance(term.value, np.ndarray) and not com.is_series(term.value):
        typ = partial(np.asanyarray, dtype=term.value.dtype)
    else:
        typ = type(term.value)
    ret = typ,

    if not hasattr(term.value, 'axes'):
        ret += None,
    else:
        ret += _zip_axes_from_type(typ, term.value.axes),
    return ret


def _zip_axes_from_type(typ, new_axes):
    axes = {}
    for ax_ind, ax_name in typ._AXIS_NAMES.iteritems():
        axes[ax_name] = new_axes[ax_ind]
    return axes


def _maybe_promote_shape(values, naxes):
    # test to see if we have an array else leave since must be a number
    if not isinstance(values, np.ndarray):
        return values

    ndims = values.ndim
    if ndims > naxes:
        raise AssertionError('cannot have more dims than axes, '
                             '{0} > {1}'.format(ndims, naxes))
    if ndims == naxes:
        return values

    ndim = set(xrange(ndims))
    nax = set(xrange(naxes))

    axes_slice = [slice(None)] * naxes

    # symmetric difference of numaxes and ndims
    slices = nax - ndim

    if ndims == naxes:
        if slices:
            raise AssertionError('slices should be empty if ndims == naxes '
                                 '{0}'.format(slices))
    else:
        if not slices:
            raise AssertionError('slices should NOT be empty if ndim != naxes '
                                 '{0}'.format(slices))

    for sl in slices:
        axes_slice[sl] = np.newaxis

    return values[tuple(axes_slice)]


def _any_pandas_objects(terms):
    """Check a sequence of terms for instances of PandasObject."""
    return any(com.is_pd_obj(term.value) for term in terms)


def _filter_special_cases(f):
    @wraps(f)
    def wrapper(terms):
        # single unary operand
        if len(terms) == 1:
            return _align_core_single_unary_op(terms[0])

        term_values = (term.value for term in terms)
        # only scalars
        if all(isinstance(term.value, pd.Index) or term.isscalar for term in
               terms):
            return np.result_type(*term_values), None

        # single element ndarrays
        all_has_size = all(hasattr(term.value, 'size') for term in terms)
        if all_has_size and all(term.value.size == 1 for term in terms):
            return np.result_type(*term_values), None

        # no pandas so just punt to the evaluator
        if not _any_pandas_objects(terms):
            return np.result_type(*term_values), None

        return f(terms)
    return wrapper


@_filter_special_cases
def _align_core(terms):
    term_index = [i for i, term in enumerate(terms) if hasattr(term.value,
                                                               'axes')]
    term_dims = [terms[i].value.ndim for i in term_index]
    ndims = pd.Series(dict(zip(term_index, term_dims)))

    # initial axes are the axes of the largest-axis'd term
    biggest = terms[ndims.idxmax()].value
    typ = biggest._constructor
    axes = biggest.axes
    naxes = len(axes)

    for term in (terms[i] for i in term_index):
        for axis, items in enumerate(term.value.axes):
            if com.is_series(term.value) and naxes > 1:
                ax, itm = naxes - 1, term.value.index
            else:
                ax, itm = axis, items
            axes[ax] = axes[ax].join(itm, how='outer')

    for i, ndim in ndims.iteritems():
        for axis, items in izip(xrange(ndim), axes):
            ti = terms[i].value

            if hasattr(ti, 'reindex_axis'):
                transpose = com.is_series(ti) and naxes > 1

                if transpose:
                    f = partial(ti.reindex, index=axes[naxes - 1], copy=False)
                else:
                    f = partial(ti.reindex_axis, items, axis=axis, copy=False)

                if pd.lib.is_bool_array(ti.values):
                    r = f(fill_value=True)
                else:
                    r = f()

                terms[i].update(r)

        res = _maybe_promote_shape(terms[i].value.T if transpose else
                                   terms[i].value, naxes)
        res = res.T if transpose else res

        try:
            v = res.values
        except AttributeError:
            v = res
        terms[i].update(v)

    return typ, _zip_axes_from_type(typ, axes)


def _filter_terms(flat):
    # numeric literals
    literals = set(filter(is_const, flat))

    # these are strings which are variable names
    names = set(flat) - literals

    # literals are not names and names are not literals, so intersection should
    # be empty
    if literals & names:
        raise ValueError('literals cannot be names and names cannot be '
                         'literals')
    return names, literals


def _align(terms):
    """Align a set of terms"""
    # flatten the parse tree (a nested list, really)
    terms = list(flatten(terms))

    # if all resolved variables are numeric scalars
    if all(term.isscalar for term in terms):
        return np.result_type(*(term.value for term in terms)).type, None

    # perform the main alignment
    typ, axes = _align_core(terms)
    return typ, axes


def _reconstruct_object(typ, obj, axes, dtype):
    """Reconstruct an object given its type, raw value, and possibly empty
    (None) axes.

    Parameters
    ----------
    typ : object
        A type
    obj : object
        The value to use in the type constructor
    axes : dict
        The axes to use to construct the resulting pandas object

    Returns
    -------
    reconst : typ
        An object of type ``typ`` with the value `obj` and possible axes
        `axes`.
    """
    #import ipdb; ipdb.set_trace()
    try:
        typ = typ.type
    except AttributeError:
        pass

    if (not isinstance(typ, partial) and
        issubclass(typ, pd.core.generic.PandasObject)):
        return typ(obj, dtype=dtype, **axes)

    ret_value = typ(obj).astype(dtype)

    try:
        ret = ret_value.item()
    except ValueError:
        ret = ret_value
    return ret
