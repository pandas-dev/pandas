#!/usr/bin/env python

import numbers

import numpy as np

import six

from pandas.computation.expr import Expr, Scope
from pandas.computation.engines import _engines


def eval(expr, engine='numexpr', truediv=True, local_dict=None,
         global_dict=None, resolvers=None):
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
    engine : string, optional, default 'numexpr', {'python', 'numexpr' }
        The engine used to evaluate the expression. Supported engines are

        - 'numexpr': This default engine evaluates pandas objects using numexpr
                     for large speed ups in complex expressions with large
                     frames.
        - 'python': Performs operations as if you had eval'd in top level
                    python

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
    * The benefits of using ``eval`` are that very large frames that are terms in
      long expressions are sped up, sometimes by as much as 10x.

    See :ref:`Enhancing performance <enhancingperf.eval>` for more details.
    """
    # make sure we're passed a valid engine
    if not engine in _engines:
        raise KeyError('Invalid engine {0} passed, valid engines are'
                       ' {1}'.format(_engines.keys()))

    eng = _engines[engine]

    if isinstance(expr, six.string_types):
        # need to go 2 up in the call stack from the constructor
        env = Scope(global_dict, local_dict, frame_level=2,
                    resolvers=resolvers)
        parsed_expr = Expr(expr, engine, env, truediv)
    elif isinstance(expr, Expr):
        parsed_expr = expr
    else:
        raise TypeError("eval only accepts strings and Expr objects, you "
                        "passed a {0!r}".format(expr.__class__.__name__))


    # construct the engine and evaluate
    ret = eng(parsed_expr).evaluate()

    # sanity check for a number
    # TODO: eventually take out
    if np.isscalar(ret):
        if not isinstance(ret, (np.number, np.bool_, numbers.Number)):
            raise TypeError('scalar result must be numeric or bool, passed '
                            'type is {0!r}'.format(ret.__class__.__name__))
    return ret
