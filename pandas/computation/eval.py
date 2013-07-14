#!/usr/bin/env python

import numbers

import numpy as np

from pandas.computation.expr import Expr, Scope, _parsers
from pandas.computation.engines import _engines


def _check_engine(engine):
    if engine not in _engines:
        raise KeyError('Invalid engine {0!r} passed, valid engines are'
                       ' {1}'.format(engine, _engines.keys()))
    if engine == 'numexpr':
        try:
            import numexpr
        except ImportError:
            raise ImportError("'numexpr' not found. Cannot use "
                              "engine='numexpr' if 'numexpr' is not installed")


def _check_parser(parser):
    if parser not in _parsers:
        raise KeyError('Invalid parser {0!r} passed, valid parsers are'
                       ' {1}'.format(parser, _parsers.keys()))



def eval(expr, parser='pandas', engine='numexpr', truediv=True,
         local_dict=None, global_dict=None, resolvers=None):
    """Evaluate a Python expression as a string using various backends.

    The following arithmetic operations are supported: +, -, *, /, **, %, //
    (python engine only) along with the following boolean operations: | (or), &
    (and), and ~ (not). Series and DataFrame objects are supported and behave
    as they would with in-Python evaluation.

    Parameters
    ----------
    expr : string or Expr object
        The expression to evaluate. This can be either a string or an ``Expr``
        object.
    parser : str, optional, default 'pandas', {'pandas', 'python'}
        The parser to use to construct the syntax tree from the expression. The
        default of 'pandas' parses code slightly different than standard
        Python. See the :ref:`enhancing performance <enhancingperf.eval>`
        documentation for more details.
    engine : string, optional, default 'numexpr', {'python', 'numexpr'}
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
    obj : ndarray, scalar, DataFrame, Series

    Notes
    -----
    See :ref:`Enhancing performance <enhancingperf.eval>` for more details.
    """
    # make sure we're passed a valid engine
    _check_engine(engine)
    _check_parser(parser)

    eng = _engines[engine]

    if isinstance(expr, basestring):
        # need to go 2 up in the call stack from the constructor
        env = Scope(global_dict, local_dict, frame_level=2,
                    resolvers=resolvers)
        parsed_expr = Expr(expr, engine, parser, env, truediv)
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
