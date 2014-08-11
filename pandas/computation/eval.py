#!/usr/bin/env python

"""Top level ``eval`` module.
"""

import tokenize
from pandas.core import common as com
from pandas.computation.expr import Expr, _parsers, tokenize_string
from pandas.computation.scope import _ensure_scope
from pandas.compat import DeepChainMap, builtins
from pandas.computation.engines import _engines
from distutils.version import LooseVersion


def _check_engine(engine):
    """Make sure a valid engine is passed.

    Parameters
    ----------
    engine : str

    Raises
    ------
    KeyError
      * If an invalid engine is passed
    ImportError
      * If numexpr was requested but doesn't exist
    """
    if engine not in _engines:
        raise KeyError('Invalid engine {0!r} passed, valid engines are'
                       ' {1}'.format(engine, list(_engines.keys())))

    # TODO: validate this in a more general way (thinking of future engines
    # that won't necessarily be import-able)
    # Could potentially be done on engine instantiation
    if engine == 'numexpr':
        try:
            import numexpr
        except ImportError:
            raise ImportError("'numexpr' not found. Cannot use "
                              "engine='numexpr' for query/eval "
                              "if 'numexpr' is not installed")
        else:
            ne_version = numexpr.__version__
            if ne_version < LooseVersion('2.1'):
                raise ImportError("'numexpr' version is %s, "
                                  "must be >= 2.1" % ne_version)


def _check_parser(parser):
    """Make sure a valid parser is passed.

    Parameters
    ----------
    parser : str

    Raises
    ------
    KeyError
      * If an invalid parser is passed
    """
    if parser not in _parsers:
        raise KeyError('Invalid parser {0!r} passed, valid parsers are'
                       ' {1}'.format(parser, _parsers.keys()))


def _check_resolvers(resolvers):
    if resolvers is not None:
        for resolver in resolvers:
            if not hasattr(resolver, '__getitem__'):
                name = type(resolver).__name__
                raise TypeError('Resolver of type %r does not implement '
                                'the __getitem__ method' % name)


def _check_expression(expr):
    """Make sure an expression is not an empty string

    Parameters
    ----------
    expr : object
        An object that can be converted to a string

    Raises
    ------
    ValueError
      * If expr is an empty string
    """
    if not expr:
        raise ValueError("expr cannot be an empty string")


def _convert_expression(expr):
    """Convert an object to an expression.

    Thus function converts an object to an expression (a unicode string) and
    checks to make sure it isn't empty after conversion. This is used to
    convert operators to their string representation for recursive calls to
    :func:`~pandas.eval`.

    Parameters
    ----------
    expr : object
        The object to be converted to a string.

    Returns
    -------
    s : unicode
        The string representation of an object.

    Raises
    ------
    ValueError
      * If the expression is empty.
    """
    s = com.pprint_thing(expr)
    _check_expression(s)
    return s


def _check_for_locals(expr, stack_level, parser):
    at_top_of_stack = stack_level == 0
    not_pandas_parser = parser != 'pandas'

    if not_pandas_parser:
        msg = "The '@' prefix is only supported by the pandas parser"
    elif at_top_of_stack:
        msg = ("The '@' prefix is not allowed in "
               "top-level eval calls, \nplease refer to "
               "your variables by name without the '@' "
               "prefix")

    if at_top_of_stack or not_pandas_parser:
        for toknum, tokval in tokenize_string(expr):
            if toknum == tokenize.OP and tokval == '@':
                raise SyntaxError(msg)


def eval(expr, parser='pandas', engine='numexpr', truediv=True,
         local_dict=None, global_dict=None, resolvers=(), level=0,
         target=None):
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
    expr : str or unicode
        The expression to evaluate. This string cannot contain any Python
        `statements
        <http://docs.python.org/2/reference/simple_stmts.html#simple-statements>`__,
        only Python `expressions
        <http://docs.python.org/2/reference/simple_stmts.html#expression-statements>`__.
    parser : string, default 'pandas', {'pandas', 'python'}
        The parser to use to construct the syntax tree from the expression. The
        default of ``'pandas'`` parses code slightly different than standard
        Python. Alternatively, you can parse an expression using the
        ``'python'`` parser to retain strict Python semantics.  See the
        :ref:`enhancing performance <enhancingperf.eval>` documentation for
        more details.
    engine : string, default 'numexpr', {'python', 'numexpr'}

        The engine used to evaluate the expression. Supported engines are

        - ``'numexpr'``: This default engine evaluates pandas objects using
                         numexpr for large speed ups in complex expressions
                         with large frames.
        - ``'python'``: Performs operations as if you had ``eval``'d in top
                        level python. This engine is generally not that useful.

        More backends may be available in the future.

    truediv : bool, optional
        Whether to use true division, like in Python >= 3
    local_dict : dict or None, optional
        A dictionary of local variables, taken from locals() by default.
    global_dict : dict or None, optional
        A dictionary of global variables, taken from globals() by default.
    resolvers : list of dict-like or None, optional
        A list of objects implementing the ``__getitem__`` special method that
        you can use to inject an additional collection of namespaces to use for
        variable lookup. For example, this is used in the
        :meth:`~pandas.DataFrame.query` method to inject the
        :attr:`~pandas.DataFrame.index` and :attr:`~pandas.DataFrame.columns`
        variables that refer to their respective :class:`~pandas.DataFrame`
        instance attributes.
    level : int, optional
        The number of prior stack frames to traverse and add to the current
        scope. Most users will **not** need to change this parameter.
    target : a target object for assignment, optional, default is None
        essentially this is a passed in resolver

    Returns
    -------
    ndarray, numeric scalar, DataFrame, Series

    Notes
    -----
    The ``dtype`` of any objects involved in an arithmetic ``%`` operation are
    recursively cast to ``float64``.

    See the :ref:`enhancing performance <enhancingperf.eval>` documentation for
    more details.

    See Also
    --------
    pandas.DataFrame.query
    pandas.DataFrame.eval
    """
    expr = _convert_expression(expr)
    _check_engine(engine)
    _check_parser(parser)
    _check_resolvers(resolvers)
    _check_for_locals(expr, level, parser)

    # get our (possibly passed-in) scope
    level += 1
    env = _ensure_scope(level, global_dict=global_dict,
                        local_dict=local_dict, resolvers=resolvers,
                        target=target)

    parsed_expr = Expr(expr, engine=engine, parser=parser, env=env,
                       truediv=truediv)

    # construct the engine and evaluate the parsed expression
    eng = _engines[engine]
    eng_inst = eng(parsed_expr)
    ret = eng_inst.evaluate()

    # assign if needed
    if env.target is not None and parsed_expr.assigner is not None:
        env.target[parsed_expr.assigner] = ret
        return None

    return ret
