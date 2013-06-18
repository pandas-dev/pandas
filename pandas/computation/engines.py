import abc
import functools
from functools import partial
from itertools import izip

import numpy as np

import pandas as pd
import pandas.core.common as com
from pandas.computation.ops import _resolve_name, _update_names
from pandas.computation.common import flatten


def _align_core_single_unary_op(term):
    if isinstance(term, np.ndarray) and not com.is_series(term):
        typ = np.asanyarray
    else:
        typ = type(term)
    ret = typ, [term]

    if not hasattr(term, 'axes'):
        ret += None,
    else:
        ret += _zip_axes_from_type(typ, term.axes),
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

    # symmetric difference
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
    return any(com.is_pd_obj(term) for term in terms)


def _filter_special_cases(f):
    @functools.wraps(f)
    def wrapper(terms):
        # need to ensure that terms is not an iterator
        terms = list(terms)

        ## special cases

        # single unary operand
        if len(terms) == 1:
            return _align_core_single_unary_op(terms[0])

        # only scalars
        elif all(np.isscalar(term) for term in terms):
            return np.result_type(*terms), terms, None

        # single element ndarrays
        all_has_size = all(hasattr(term, 'size') for term in terms)
        if (all_has_size and all(term.size == 1 for term in terms)):
            return np.result_type(*terms), terms, None

        # no pandas so just punt to the evaluator
        if not _any_pandas_objects(terms):
            return np.result_type(*terms), terms, None

        return f(terms)
    return wrapper


@_filter_special_cases
def _align_core(terms):
    term_index = [i for i, term in enumerate(terms) if hasattr(term, 'axes')]
    term_dims = [terms[i].ndim for i in term_index]
    ndims = pd.Series(dict(zip(term_index, term_dims)))

    # initial axes are the axes of the largest-axis'd term
    biggest = terms[ndims.idxmax()]
    typ = biggest._constructor
    axes = biggest.axes
    naxes = len(axes)

    for i in term_index:
        for axis, items in enumerate(terms[i].axes):
            if com.is_series(terms[i]) and naxes > 1:
                axes[naxes - 1] = axes[naxes - 1].join(terms[i].index,
                                                       how='outer')
            else:
                axes[axis] = axes[axis].join(items, how='outer')

    for i, ndim in ndims.iteritems():
        for axis, items in izip(xrange(ndim), axes):
            ti = terms[i]  # needed here because we modify it in the inner loop

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

                terms[i] = r

        res = _maybe_promote_shape(terms[i].T if transpose else terms[i],
                                   naxes)
        res = res.T if transpose else res

        try:
            terms[i] = res.values
        except AttributeError:
            terms[i] = res

    return typ, terms, _zip_axes_from_type(typ, axes)


def _filter_terms(flat):
    # numeric literals
    literals = filter(lambda string: not com.is_string(string), flat)
    literals_set = set(literals)

    # these are strings which are variable names
    names = filter(com.is_string, flat)
    names_set = set(names)

    # literals are not names and names are not literals, by definition
    if literals_set & names_set:
        raise ValueError('literals cannot be names and names cannot be '
                         'literals')
    return names, literals


def _align(terms, env):
    # flatten the parse tree (a nested list)
    flat = list(flatten(terms))

    # separate names and literals
    names, literals = _filter_terms(flat)

    if not names:  # only literals so just promote to a common type
        return np.result_type(*literals).type, None

    # get the variables out
    resolve_in_env = partial(_resolve_name, env)
    resolved = map(resolve_in_env, names)

    # if all resolved variables are numeric scalars
    if all(np.isscalar(rsv) for rsv in resolved):
        return np.result_type(*resolved).type, None

    # perform the main alignment
    typ, resolved, axes = _align_core(resolved)

    # put the aligned arrays back in the table
    _update_names(env, dict(izip(names, resolved)))

    # we need this to reconstruct things after evaluation since we CANNOT
    # depend on the array interface
    return typ, axes


def _reconstruct_object(typ, obj, axes):
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
    try:
        # handle numpy dtypes
        typ = typ.type
    except AttributeError:
        pass

    if typ != np.asanyarray and issubclass(typ, pd.core.generic.PandasObject):
        return typ(obj, **axes)

    ret_value = typ(obj)

    try:
        return ret_value.item()
    except (AttributeError, ValueError):
        return ret_value


class AbstractEngine(object):
    """"""
    __metaclass__ = abc.ABCMeta

    has_neg_frac = False

    def __init__(self, expr):
        self.expr = expr
        self.aligned_axes = None
        self.result_type = None

    @abc.abstractmethod
    def convert(self):
        """Convert an expression for evaluation."""
        pass

    def evaluate(self, env):
        if not self._is_aligned:
            self.result_type, self.aligned_axes = _align(self.expr.terms, env)

        res = self._evaluate(env)
        return _reconstruct_object(self.result_type, res, self.aligned_axes)

    @property
    def _is_aligned(self):
        return self.aligned_axes is not None and self.result_type is not None

    @abc.abstractmethod
    def _evaluate(self, env):
        """Return an evaluated expression."""
        pass


class NumExprEngine(AbstractEngine):
    """NumExpr engine class"""
    has_neg_frac = True

    def __init__(self, expr):
        super(NumExprEngine, self).__init__(expr)

    def convert(self):
        """Return a string"""
        return '%s' % self.expr

    def _evaluate(self, env):
        import numexpr as ne

        try:
            return ne.evaluate(self.convert(), local_dict=env.locals,
                               global_dict=env.globals,
                               truediv=self.expr.truediv)
        except KeyError as e:
            raise NameError('{0!r} is not defined'.format(e.message))


class PythonEngine(AbstractEngine):
    """Use NumPy even if numexpr is installed"""
    has_neg_frac = False

    def __init__(self, expr):
        super(PythonEngine, self).__init__(expr)

    def convert(self):
        pass

    def evaluate(self, env):
        return self.expr(env)

    def _evaluate(self, env):
        pass


_engines = {'numexpr': NumExprEngine, 'python': PythonEngine}
