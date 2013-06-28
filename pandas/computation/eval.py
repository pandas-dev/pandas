#!/usr/bin/env python

import sys
import numbers

import numpy as np

from pandas.computation.expr import Expr, Scope
from pandas.computation.engines import _engines


def eval(expr, engine='numexpr', truediv=True, local_dict=None,
         global_dict=None):
    """Evaluate a Python expression as a string.

    Parameters
    ----------
    expr : string or Expr object
    engine : string, optional, default 'numexpr'
        The engine to use to evaluate the passed expression
    truediv : bool, optional, default True
    local_dict : dict or None, optional, default None
    global_dict : dict or None, optional, default None

    Returns
    -------
    obj : ndarray, scalar, DataFrame, Series, or Panel
    """
    # make sure we're passed a valid engine
    if not engine in _engines:
        raise KeyError('Invalid engine {0} passed, valid engines are'
                       ' {1}'.format(_engines.keys()))

    eng = _engines[engine]

    if isinstance(expr, basestring):
        frame = sys._getframe(1)

        # get the globals and locals
        gbl, lcl = (global_dict or frame.f_globals,
                    local_dict or frame.f_locals)

        try:
            # shallow copy the scope so we don't overwrite anything
            env = Scope(gbl.copy(), lcl.copy())
        finally:
            del frame
        parsed_expr = Expr(expr, engine, env, truediv)
    elif isinstance(expr, Expr):
        parsed_expr = expr
    else:
        raise TypeError("eval only accepts strings and Expr objects, you "
                        "passed a {0!r}".format(expr.__class__.__name__))


    # construct the engine and evaluate
    ret = eng(parsed_expr).evaluate()

    # sanity check for a number
    if np.isscalar(ret):
        if not isinstance(ret, (np.number, numbers.Number, np.bool_, bool)):
            raise TypeError('scalar result must be numeric or bool, type is '
                            '{0!r}'.format(ret.__class__.__name__))
    return ret
