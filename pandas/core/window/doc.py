"""Any shareable docstring components for rolling/expanding/ewm"""

from __future__ import annotations

from textwrap import dedent

from pandas.core.shared_docs import _shared_docs

_shared_docs = dict(**_shared_docs)


def create_section_header(header: str) -> str:
    """Create numpydoc section header"""
    return f"{header}\n{'-' * len(header)}\n"


template_header = "\nCalculate the {window_method} {aggregation_description}.\n\n"

template_returns = dedent(
    """
    Series or DataFrame
        Return type is the same as the original object with ``np.float64`` dtype.\n
    """
).replace("\n", "", 1)

template_see_also = dedent(
    """
    Series.{window_method} : Calling {window_method} with Series data.
    DataFrame.{window_method} : Calling {window_method} with DataFrames.
    Series.{agg_method} : Aggregating {agg_method} for Series.
    DataFrame.{agg_method} : Aggregating {agg_method} for DataFrame.\n
    """
).replace("\n", "", 1)

kwargs_numeric_only = dedent(
    """
    numeric_only : bool, default False
        Include only float, int, boolean columns.

        .. versionadded:: 1.5.0\n
    """
).replace("\n", "", 1)

kwargs_scipy = dedent(
    """
    **kwargs
        Keyword arguments to configure the ``SciPy`` weighted window type.\n
    """
).replace("\n", "", 1)

window_apply_parameters = dedent(
    """
    func : function
        Must produce a single value from an ndarray input if ``raw=True``
        or a single value from a Series if ``raw=False``. Can also accept a
        Numba JIT function with ``engine='numba'`` specified.

    raw : bool, default False
        * ``False`` : passes each row or column as a Series to the
          function.
        * ``True`` : the passed function will receive ndarray
          objects instead.
          If you are just applying a NumPy reduction function this will
          achieve much better performance.

    engine : str, default None
        * ``'cython'`` : Runs rolling apply through C-extensions from cython.
        * ``'numba'`` : Runs rolling apply through JIT compiled code from numba.
          Only available when ``raw`` is set to ``True``.
        * ``None`` : Defaults to ``'cython'`` or globally setting ``compute.use_numba``

    engine_kwargs : dict, default None
        * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
        * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
          and ``parallel`` dictionary keys. The values must either be ``True`` or
          ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
          ``{{'nopython': True, 'nogil': False, 'parallel': False}}`` and will be
          applied to both the ``func`` and the ``apply`` rolling aggregation.

    args : tuple, default None
        Positional arguments to be passed into func.

    kwargs : dict, default None
        Keyword arguments to be passed into func.\n
    """
).replace("\n", "", 1)

template_pipe = """
Apply a ``func`` with arguments to this %(klass)s object and return its result.

Use `.pipe` when you want to improve readability by chaining together
functions that expect Series, DataFrames, GroupBy, Rolling, Expanding or Resampler
objects.
Instead of writing

>>> h = lambda x, arg2, arg3: x + 1 - arg2 * arg3
>>> g = lambda x, arg1: x * 5 / arg1
>>> f = lambda x: x ** 4
>>> df = pd.DataFrame({'A': [1, 2, 3, 4]}, index=pd.date_range('2012-08-02', periods=4))
>>> h(g(f(df.rolling('2D')), arg1=1), arg2=2, arg3=3)  # doctest: +SKIP

You can write

>>> (df.rolling('2D')
...    .pipe(f)
...    .pipe(g, arg1=1)
...    .pipe(h, arg2=2, arg3=3))  # doctest: +SKIP

which is much more readable.

Parameters
----------
func : callable or tuple of (callable, str)
    Function to apply to this %(klass)s object or, alternatively,
    a `(callable, data_keyword)` tuple where `data_keyword` is a
    string indicating the keyword of `callable` that expects the
    %(klass)s object.
*args : iterable, optional
       Positional arguments passed into `func`.
**kwargs : dict, optional
         A dictionary of keyword arguments passed into `func`.

Returns
-------
%(klass)s
    The original object with the function `func` applied.

See Also
--------
Series.pipe : Apply a function with arguments to a series.
DataFrame.pipe: Apply a function with arguments to a dataframe.
apply : Apply function to each group instead of to the
    full %(klass)s object.

Notes
-----
See more `here
<https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#piping-function-calls>`_

Examples
--------
%(examples)s
"""

numba_notes = (
    "See :ref:`window.numba_engine` and :ref:`enhancingperf.numba` for "
    "extended documentation and performance considerations for the Numba engine.\n\n"
)


def window_agg_numba_parameters(version: str = "1.3") -> str:
    return (
        dedent(
            """
    engine : str, default None
        * ``'cython'`` : Runs the operation through C-extensions from cython.
        * ``'numba'`` : Runs the operation through JIT compiled code from numba.
        * ``None`` : Defaults to ``'cython'`` or globally setting ``compute.use_numba``

          .. versionadded:: {version}.0

    engine_kwargs : dict, default None
        * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
        * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
          and ``parallel`` dictionary keys. The values must either be ``True`` or
          ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
          ``{{'nopython': True, 'nogil': False, 'parallel': False}}``

          .. versionadded:: {version}.0\n
    """
        )
        .replace("\n", "", 1)
        .replace("{version}", version)
    )
