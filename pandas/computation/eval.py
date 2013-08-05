#!/usr/bin/env python

import numbers
import numpy as np

from pandas.compat import string_types
from pandas.computation.expr import Expr, _parsers, _ensure_scope
from pandas.computation.engines import _engines


def _check_engine(engine):
    """make sure a valid engine is passed"""
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
    """make sure a valid parser is passed"""
    if parser not in _parsers:
        raise KeyError('Invalid parser {0!r} passed, valid parsers are'
                       ' {1}'.format(parser, _parsers.keys()))


def eval(expr, parser='pandas', engine='numexpr', truediv=True,
         local_dict=None, global_dict=None, resolvers=None, level=2):
    """Evaluate a Python expression as a string using various backends.

    The following arithmetic operations are supported: ``+``, ``-``, ``*``,
    ``/``, ``**``, ``%``, ``//`` (python engine only) along with the following
    boolean operations: ``|`` (or), ``&`` (and), and ``~`` (not).
    Additionally, the ``'pandas'`` parser allows the use of :keyword:`and`,
    :keyword:`or`, and :keyword:`not` with the same semantics as the
    corresponding bitwise operators.  :class:`~pandas.Series` and
    :class:`~pandas.DataFrame` objects are supported and behave as they would
    with plain ol' Python evaluation.

    Parameters
    ----------
    expr : string
        The expression to evaluate.
    parser : string, default 'pandas', {'pandas', 'python'}
        The parser to use to construct the syntax tree from the expression. The
        default of 'pandas' parses code slightly different than standard
        Python. See the :ref:`enhancing performance <enhancingperf.eval>`
        documentation for more details.
    engine : string, default 'numexpr', {'python', 'numexpr'}

        The engine used to evaluate the expression. Supported engines are

        - ``'numexpr'``: This default engine evaluates pandas objects using
                         numexpr for large speed ups in complex expressions
                         with large frames.
        - ``'python'``: Performs operations as if you had ``eval``'d in top
                        level python. This engine is generally not that useful.

    truediv : bool, default True
        Whether to use true division, like in Python >= 3
    local_dict : dict or None, default None
        A dictionary of local variables, taken from locals() by default.
    global_dict : dict or None, default None
        A dictionary of global variables, taken from globals() by default.
    resolvers : dict of dict-like or None, default None
        A dictionary of dict-like object (specifically they must implement the
        ``get`` method) that you can use to inject an additional collection of
        namespaces to use for variable lookup. This is used in the
        :meth:`~pandas.DataFrame.query` method to inject the
        :attr:`~pandas.DataFrame.index` and :attr:`~pandas.DataFrame.columns`
        variables that refer to their respective :class:`~pandas.DataFrame`
        instance attributes.
    level : int, default 2
        The number of prior stack frames to traverse and add to the current
        scope.

    Returns
    -------
    ret : ndarray, numeric scalar, :class:`~pandas.DataFrame`, :class:`~pandas.Series`

    Notes
    -----
    The ``dtype`` of any objects involved in an arithmetic ``%`` operation are
    recursively cast to ``float64``.

    See the :ref:`enhancing performance <enhancingperf.eval>` documentation for
    more details.

    See Also
    --------
    pandas.DataFrame.query
    """
    # make sure we're passed a valid engine and parser
    _check_engine(engine)
    _check_parser(parser)

    env = _ensure_scope(global_dict=global_dict, local_dict=local_dict,
                        resolvers=resolvers, level=level)

    if isinstance(expr, string_types):
        parsed_expr = Expr(expr, engine=engine, parser=parser, env=env,
                           truediv=truediv)
    else:
        raise TypeError("eval only accepts strings, you passed an object of "
                        "type {0!r}".format(expr.__class__.__name__))

    # construct the engine and evaluate
    eng = _engines[engine]
    eng_inst = eng(parsed_expr)
    ret = eng_inst.evaluate()

    # sanity check for a number if it's a scalar result
    # TODO: eventually take out
    if np.isscalar(ret):
        if not isinstance(ret, (np.number, np.bool_, numbers.Number)):
            raise TypeError('scalar result must be numeric or bool, return'
                            ' type is {0!r}'.format(ret.__class__.__name__))
    return ret
