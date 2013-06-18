#!/usr/bin/env python

import sys
import numbers

import numpy as np

from pandas.computation.expr import Expr, Scope
from pandas.computation.engines import _engines


def eval(expr, engine='numexpr', truediv=True, local_dict=None,
         global_dict=None):
    # make sure we're passed a valid engine
    if not engine in _engines:
        raise KeyError('Invalid engine {0} passed, valid engines are'
                       ' {1}'.format(_engines.keys()))

    # 1 up in the call stack for locals/globals; see the documentation for the
    # inspect module for why you must decrease the refcount of frame at all
    # costs
    frame = sys._getframe(1)

    try:
        # get the globals and locals
        gbl, lcl = global_dict or frame.f_globals, local_dict or frame.f_locals

        # shallow copy the scope so we don't overwrite everything
        env = Scope(gbl.copy(), lcl.copy())

        # parse the expression
        parsed_expr = Expr(expr, engine, env, truediv)

        # choose the engine
        eng = _engines[engine]

        # construct the engine and evaluate
        ret = eng(parsed_expr).evaluate(env)
    finally:
        del frame

    # sanity check for a number
    if np.isscalar(ret):
        if not isinstance(ret, (np.number, numbers.Number, np.bool_, bool)):
            raise TypeError('scalar result must be numeric or bool, type is '
                            '{0!r}'.format(ret.__class__.__name__))
    return ret
