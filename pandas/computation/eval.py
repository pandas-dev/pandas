#!/usr/bin/env python

import sys
import numbers

import numpy as np

from pandas.computation.expr import Expr, Scope
from pandas.computation.engines import _engines


def eval(expr, engine='numexpr', truediv=True, local_dict=None,
         global_dict=None):
    """Evaluate a Python expression as a string using various backends.

    The following arithmetic operations are supported: +, -, *, /, **, %, //
    (python engine only) along with the following boolean operations: | (or), &
    (and), and ~ (not). All Pandas objects are supported and behave as they
    would with in-Python evaluation.

    Parameters
    ----------
    expr : string or Expr object
        The expression to evaluate. This can be either a string or an ``Expr``
        object.
    engine : string, optional, default 'numexpr', {'python', 'numexpr', 'pytables'}
        The engine used to evaluate the expression. Supported engines are

        - 'numexpr': This default engine evaluates pandas objects using numexpr
                     for large speed ups in complex expressions with large
                     frames.
        - 'python': Performs operations as if you had eval'd in top level
                    python
        - 'pytables': Engine used for evaluating expressions for selection of
                      objects from PyTables HDF5 tables.

    truediv : bool, optional, default True
        Whether to use true division, like in Python >= 3
    local_dict : dict or None, optional, default None
        A dictionary of local variables, taken from locals() by default.
    global_dict : dict or None, optional, default None
        A dictionary of global variables, taken from globals() by default.

    Returns
    -------
    obj : ndarray, scalar, DataFrame, Series, or Panel

    Notes
    -----
    The benefits of using ``eval`` are that very large frames that are terms in
    long expressions are sped up, sometimes by as much as 10x.
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
