from __future__ import annotations

_shared_docs: dict[str, str] = {}

_shared_docs[
    "aggregate"
] = """
Aggregate using one or more operations over the specified axis.

Parameters
----------
func : function, str, list or dict
    Function to use for aggregating the data. If a function, must either
    work when passed a {klass} or when passed to {klass}.apply.

    Accepted combinations are:

    - function
    - string function name
    - list of functions and/or function names, e.g. ``[np.sum, 'mean']``
    - dict of axis labels -> functions, function names or list of such.
{axis}
*args
    Positional arguments to pass to `func`.
**kwargs
    Keyword arguments to pass to `func`.

Returns
-------
scalar, Series or DataFrame

    The return can be:

    * scalar : when Series.agg is called with single function
    * Series : when DataFrame.agg is called with a single function
    * DataFrame : when DataFrame.agg is called with several functions

    Return scalar, Series or DataFrame.
{see_also}
Notes
-----
The aggregation operations are always performed over an axis, either the
index (default) or the column axis. This behavior is different from
`numpy` aggregation functions (`mean`, `median`, `prod`, `sum`, `std`,
`var`), where the default is to compute the aggregation of the flattened
array, e.g., ``numpy.mean(arr_2d)`` as opposed to
``numpy.mean(arr_2d, axis=0)``.

`agg` is an alias for `aggregate`. Use the alias.

Functions that mutate the passed object can produce unexpected
behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
for more details.

A passed user-defined-function will be passed a Series for evaluation.
{examples}"""

_shared_docs[
    "compare"
] = """
Compare to another {klass} and show the differences.

.. versionadded:: 1.1.0

Parameters
----------
other : {klass}
    Object to compare with.

align_axis : {{0 or 'index', 1 or 'columns'}}, default 1
    Determine which axis to align the comparison on.

    * 0, or 'index' : Resulting differences are stacked vertically
        with rows drawn alternately from self and other.
    * 1, or 'columns' : Resulting differences are aligned horizontally
        with columns drawn alternately from self and other.

keep_shape : bool, default False
    If true, all rows and columns are kept.
    Otherwise, only the ones with different values are kept.

keep_equal : bool, default False
    If true, the result keeps values that are equal.
    Otherwise, equal values are shown as NaNs.

result_names : tuple, default ('self', 'other')
    Set the dataframes names in the comparison.

    .. versionadded:: 1.5.0
"""

_shared_docs[
    "groupby"
] = """
Group %(klass)s using a mapper or by a Series of columns.

A groupby operation involves some combination of splitting the
object, applying a function, and combining the results. This can be
used to group large amounts of data and compute operations on these
groups.

Parameters
----------
by : mapping, function, label, pd.Grouper or list of such
    Used to determine the groups for the groupby.
    If ``by`` is a function, it's called on each value of the object's
    index. If a dict or Series is passed, the Series or dict VALUES
    will be used to determine the groups (the Series' values are first
    aligned; see ``.align()`` method). If a list or ndarray of length
    equal to the selected axis is passed (see the `groupby user guide
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#splitting-an-object-into-groups>`_),
    the values are used as-is to determine the groups. A label or list
    of labels may be passed to group by the columns in ``self``.
    Notice that a tuple is interpreted as a (single) key.
axis : {0 or 'index', 1 or 'columns'}, default 0
    Split along rows (0) or columns (1). For `Series` this parameter
    is unused and defaults to 0.
level : int, level name, or sequence of such, default None
    If the axis is a MultiIndex (hierarchical), group by a particular
    level or levels. Do not specify both ``by`` and ``level``.
as_index : bool, default True
    Return object with group labels as the
    index. Only relevant for DataFrame input. as_index=False is
    effectively "SQL-style" grouped output. This argument has no effect
    on filtrations (see the `filtrations in the user guide
    <https://pandas.pydata.org/docs/dev/user_guide/groupby.html#filtration>`_),
    such as ``head()``, ``tail()``, ``nth()`` and in transformations
    (see the `transformations in the user guide
    <https://pandas.pydata.org/docs/dev/user_guide/groupby.html#transformation>`_).
sort : bool, default True
    Sort group keys. Get better performance by turning this off.
    Note this does not influence the order of observations within each
    group. Groupby preserves the order of rows within each group.
    This argument has no effect on filtrations (see the `filtrations in the user guide
    <https://pandas.pydata.org/docs/dev/user_guide/groupby.html#filtration>`_),
    such as ``head()``, ``tail()``, ``nth()`` and in transformations
    (see the `transformations in the user guide
    <https://pandas.pydata.org/docs/dev/user_guide/groupby.html#transformation>`_).

    .. versionchanged:: 2.0.0

        Specifying ``sort=False`` with an ordered categorical grouper will no
        longer sort the values.

group_keys : bool, default True
    When calling apply and the ``by`` argument produces a like-indexed
    (i.e. :ref:`a transform <groupby.transform>`) result, add group keys to
    index to identify pieces. By default group keys are not included
    when the result's index (and column) labels match the inputs, and
    are included otherwise.

    .. versionchanged:: 1.5.0

       Warns that ``group_keys`` will no longer be ignored when the
       result from ``apply`` is a like-indexed Series or DataFrame.
       Specify ``group_keys`` explicitly to include the group keys or
       not.

    .. versionchanged:: 2.0.0

       ``group_keys`` now defaults to ``True``.

observed : bool, default False
    This only applies if any of the groupers are Categoricals.
    If True: only show observed values for categorical groupers.
    If False: show all values for categorical groupers.

    .. deprecated:: 2.1.0

        The default value will change to True in a future version of pandas.

dropna : bool, default True
    If True, and if group keys contain NA values, NA values together
    with row/column will be dropped.
    If False, NA values will also be treated as the key in groups.

    .. versionadded:: 1.1.0

Returns
-------
%(klass)sGroupBy
    Returns a groupby object that contains information about the groups.

See Also
--------
resample : Convenience method for frequency conversion and resampling
    of time series.

Notes
-----
See the `user guide
<https://pandas.pydata.org/pandas-docs/stable/groupby.html>`__ for more
detailed usage and examples, including splitting an object into groups,
iterating through groups, selecting a group, aggregation, and more.
"""

_shared_docs[
    "melt"
] = """
Unpivot a DataFrame from wide to long format, optionally leaving identifiers set.

This function is useful to massage a DataFrame into a format where one
or more columns are identifier variables (`id_vars`), while all other
columns, considered measured variables (`value_vars`), are "unpivoted" to
the row axis, leaving just two non-identifier columns, 'variable' and
'value'.

Parameters
----------
id_vars : tuple, list, or ndarray, optional
    Column(s) to use as identifier variables.
value_vars : tuple, list, or ndarray, optional
    Column(s) to unpivot. If not specified, uses all columns that
    are not set as `id_vars`.
var_name : scalar
    Name to use for the 'variable' column. If None it uses
    ``frame.columns.name`` or 'variable'.
value_name : scalar, default 'value'
    Name to use for the 'value' column.
col_level : int or str, optional
    If columns are a MultiIndex then use this level to melt.
ignore_index : bool, default True
    If True, original index is ignored. If False, the original index is retained.
    Index labels will be repeated as necessary.

    .. versionadded:: 1.1.0

Returns
-------
DataFrame
    Unpivoted DataFrame.

See Also
--------
%(other)s : Identical method.
pivot_table : Create a spreadsheet-style pivot table as a DataFrame.
DataFrame.pivot : Return reshaped DataFrame organized
    by given index / column values.
DataFrame.explode : Explode a DataFrame from list-like
        columns to long format.

Notes
-----
Reference :ref:`the user guide <reshaping.melt>` for more examples.

Examples
--------
>>> df = pd.DataFrame({'A': {0: 'a', 1: 'b', 2: 'c'},
...                    'B': {0: 1, 1: 3, 2: 5},
...                    'C': {0: 2, 1: 4, 2: 6}})
>>> df
   A  B  C
0  a  1  2
1  b  3  4
2  c  5  6

>>> %(caller)sid_vars=['A'], value_vars=['B'])
   A variable  value
0  a        B      1
1  b        B      3
2  c        B      5

>>> %(caller)sid_vars=['A'], value_vars=['B', 'C'])
   A variable  value
0  a        B      1
1  b        B      3
2  c        B      5
3  a        C      2
4  b        C      4
5  c        C      6

The names of 'variable' and 'value' columns can be customized:

>>> %(caller)sid_vars=['A'], value_vars=['B'],
...         var_name='myVarname', value_name='myValname')
   A myVarname  myValname
0  a         B          1
1  b         B          3
2  c         B          5

Original index values can be kept around:

>>> %(caller)sid_vars=['A'], value_vars=['B', 'C'], ignore_index=False)
   A variable  value
0  a        B      1
1  b        B      3
2  c        B      5
0  a        C      2
1  b        C      4
2  c        C      6

If you have multi-index columns:

>>> df.columns = [list('ABC'), list('DEF')]
>>> df
   A  B  C
   D  E  F
0  a  1  2
1  b  3  4
2  c  5  6

>>> %(caller)scol_level=0, id_vars=['A'], value_vars=['B'])
   A variable  value
0  a        B      1
1  b        B      3
2  c        B      5

>>> %(caller)sid_vars=[('A', 'D')], value_vars=[('B', 'E')])
  (A, D) variable_0 variable_1  value
0      a          B          E      1
1      b          B          E      3
2      c          B          E      5
"""

_shared_docs[
    "transform"
] = """
Call ``func`` on self producing a {klass} with the same axis shape as self.

Parameters
----------
func : function, str, list-like or dict-like
    Function to use for transforming the data. If a function, must either
    work when passed a {klass} or when passed to {klass}.apply. If func
    is both list-like and dict-like, dict-like behavior takes precedence.

    Accepted combinations are:

    - function
    - string function name
    - list-like of functions and/or function names, e.g. ``[np.exp, 'sqrt']``
    - dict-like of axis labels -> functions, function names or list-like of such.
{axis}
*args
    Positional arguments to pass to `func`.
**kwargs
    Keyword arguments to pass to `func`.

Returns
-------
{klass}
    A {klass} that must have the same length as self.

Raises
------
ValueError : If the returned {klass} has a different length than self.

See Also
--------
{klass}.agg : Only perform aggregating type operations.
{klass}.apply : Invoke function on a {klass}.

Notes
-----
Functions that mutate the passed object can produce unexpected
behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
for more details.

Examples
--------
>>> df = pd.DataFrame({{'A': range(3), 'B': range(1, 4)}})
>>> df
   A  B
0  0  1
1  1  2
2  2  3
>>> df.transform(lambda x: x + 1)
   A  B
0  1  2
1  2  3
2  3  4

Even though the resulting {klass} must have the same length as the
input {klass}, it is possible to provide several input functions:

>>> s = pd.Series(range(3))
>>> s
0    0
1    1
2    2
dtype: int64
>>> s.transform([np.sqrt, np.exp])
       sqrt        exp
0  0.000000   1.000000
1  1.000000   2.718282
2  1.414214   7.389056

You can call transform on a GroupBy object:

>>> df = pd.DataFrame({{
...     "Date": [
...         "2015-05-08", "2015-05-07", "2015-05-06", "2015-05-05",
...         "2015-05-08", "2015-05-07", "2015-05-06", "2015-05-05"],
...     "Data": [5, 8, 6, 1, 50, 100, 60, 120],
... }})
>>> df
         Date  Data
0  2015-05-08     5
1  2015-05-07     8
2  2015-05-06     6
3  2015-05-05     1
4  2015-05-08    50
5  2015-05-07   100
6  2015-05-06    60
7  2015-05-05   120
>>> df.groupby('Date')['Data'].transform('sum')
0     55
1    108
2     66
3    121
4     55
5    108
6     66
7    121
Name: Data, dtype: int64

>>> df = pd.DataFrame({{
...     "c": [1, 1, 1, 2, 2, 2, 2],
...     "type": ["m", "n", "o", "m", "m", "n", "n"]
... }})
>>> df
   c type
0  1    m
1  1    n
2  1    o
3  2    m
4  2    m
5  2    n
6  2    n
>>> df['size'] = df.groupby('c')['type'].transform(len)
>>> df
   c type size
0  1    m    3
1  1    n    3
2  1    o    3
3  2    m    4
4  2    m    4
5  2    n    4
6  2    n    4
"""

_shared_docs[
    "storage_options"
] = """storage_options : dict, optional
    Extra options that make sense for a particular storage connection, e.g.
    host, port, username, password, etc. For HTTP(S) URLs the key-value pairs
    are forwarded to ``urllib.request.Request`` as header options. For other
    URLs (e.g. starting with "s3://", and "gcs://") the key-value pairs are
    forwarded to ``fsspec.open``. Please see ``fsspec`` and ``urllib`` for more
    details, and for more examples on storage options refer `here
    <https://pandas.pydata.org/docs/user_guide/io.html?
    highlight=storage_options#reading-writing-remote-files>`_."""

_shared_docs[
    "compression_options"
] = """compression : str or dict, default 'infer'
    For on-the-fly compression of the output data. If 'infer' and '%s' is
    path-like, then detect compression from the following extensions: '.gz',
    '.bz2', '.zip', '.xz', '.zst', '.tar', '.tar.gz', '.tar.xz' or '.tar.bz2'
    (otherwise no compression).
    Set to ``None`` for no compression.
    Can also be a dict with key ``'method'`` set
    to one of {``'zip'``, ``'gzip'``, ``'bz2'``, ``'zstd'``, ``'tar'``} and other
    key-value pairs are forwarded to
    ``zipfile.ZipFile``, ``gzip.GzipFile``,
    ``bz2.BZ2File``, ``zstandard.ZstdCompressor`` or
    ``tarfile.TarFile``, respectively.
    As an example, the following could be passed for faster compression and to create
    a reproducible gzip archive:
    ``compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}``.

    .. versionadded:: 1.5.0
        Added support for `.tar` files."""

_shared_docs[
    "decompression_options"
] = """compression : str or dict, default 'infer'
    For on-the-fly decompression of on-disk data. If 'infer' and '%s' is
    path-like, then detect compression from the following extensions: '.gz',
    '.bz2', '.zip', '.xz', '.zst', '.tar', '.tar.gz', '.tar.xz' or '.tar.bz2'
    (otherwise no compression).
    If using 'zip' or 'tar', the ZIP file must contain only one data file to be read in.
    Set to ``None`` for no decompression.
    Can also be a dict with key ``'method'`` set
    to one of {``'zip'``, ``'gzip'``, ``'bz2'``, ``'zstd'``, ``'tar'``} and other
    key-value pairs are forwarded to
    ``zipfile.ZipFile``, ``gzip.GzipFile``,
    ``bz2.BZ2File``, ``zstandard.ZstdDecompressor`` or
    ``tarfile.TarFile``, respectively.
    As an example, the following could be passed for Zstandard decompression using a
    custom compression dictionary:
    ``compression={'method': 'zstd', 'dict_data': my_compression_dict}``.

    .. versionadded:: 1.5.0
        Added support for `.tar` files."""

_shared_docs[
    "replace"
] = """
    Replace values given in `to_replace` with `value`.

    Values of the {klass} are replaced with other values dynamically.
    {replace_iloc}

    Parameters
    ----------
    to_replace : str, regex, list, dict, Series, int, float, or None
        How to find the values that will be replaced.

        * numeric, str or regex:

            - numeric: numeric values equal to `to_replace` will be
              replaced with `value`
            - str: string exactly matching `to_replace` will be replaced
              with `value`
            - regex: regexs matching `to_replace` will be replaced with
              `value`

        * list of str, regex, or numeric:

            - First, if `to_replace` and `value` are both lists, they
              **must** be the same length.
            - Second, if ``regex=True`` then all of the strings in **both**
              lists will be interpreted as regexs otherwise they will match
              directly. This doesn't matter much for `value` since there
              are only a few possible substitution regexes you can use.
            - str, regex and numeric rules apply as above.

        * dict:

            - Dicts can be used to specify different replacement values
              for different existing values. For example,
              ``{{'a': 'b', 'y': 'z'}}`` replaces the value 'a' with 'b' and
              'y' with 'z'. To use a dict in this way, the optional `value`
              parameter should not be given.
            - For a DataFrame a dict can specify that different values
              should be replaced in different columns. For example,
              ``{{'a': 1, 'b': 'z'}}`` looks for the value 1 in column 'a'
              and the value 'z' in column 'b' and replaces these values
              with whatever is specified in `value`. The `value` parameter
              should not be ``None`` in this case. You can treat this as a
              special case of passing two lists except that you are
              specifying the column to search in.
            - For a DataFrame nested dictionaries, e.g.,
              ``{{'a': {{'b': np.nan}}}}``, are read as follows: look in column
              'a' for the value 'b' and replace it with NaN. The optional `value`
              parameter should not be specified to use a nested dict in this
              way. You can nest regular expressions as well. Note that
              column names (the top-level dictionary keys in a nested
              dictionary) **cannot** be regular expressions.

        * None:

            - This means that the `regex` argument must be a string,
              compiled regular expression, or list, dict, ndarray or
              Series of such elements. If `value` is also ``None`` then
              this **must** be a nested dictionary or Series.

        See the examples section for examples of each of these.
    value : scalar, dict, list, str, regex, default None
        Value to replace any values matching `to_replace` with.
        For a DataFrame a dict of values can be used to specify which
        value to use for each column (columns not in the dict will not be
        filled). Regular expressions, strings and lists or dicts of such
        objects are also allowed.
    {inplace}
    limit : int, default None
        Maximum size gap to forward or backward fill.
    regex : bool or same types as `to_replace`, default False
        Whether to interpret `to_replace` and/or `value` as regular
        expressions. If this is ``True`` then `to_replace` *must* be a
        string. Alternatively, this could be a regular expression or a
        list, dict, or array of regular expressions in which case
        `to_replace` must be ``None``.
    method : {{'pad', 'ffill', 'bfill'}}
        The method to use when for replacement, when `to_replace` is a
        scalar, list or tuple and `value` is ``None``.

    Returns
    -------
    {klass}
        Object after replacement.

    Raises
    ------
    AssertionError
        * If `regex` is not a ``bool`` and `to_replace` is not
          ``None``.

    TypeError
        * If `to_replace` is not a scalar, array-like, ``dict``, or ``None``
        * If `to_replace` is a ``dict`` and `value` is not a ``list``,
          ``dict``, ``ndarray``, or ``Series``
        * If `to_replace` is ``None`` and `regex` is not compilable
          into a regular expression or is a list, dict, ndarray, or
          Series.
        * When replacing multiple ``bool`` or ``datetime64`` objects and
          the arguments to `to_replace` does not match the type of the
          value being replaced

    ValueError
        * If a ``list`` or an ``ndarray`` is passed to `to_replace` and
          `value` but they are not the same length.

    See Also
    --------
    {klass}.fillna : Fill NA values.
    {klass}.where : Replace values based on boolean condition.
    DataFrame.applymap: Apply a function to a Dataframe elementwise.
    Series.map: Map values of Series according to an input mapping or function.
    Series.str.replace : Simple string replacement.

    Notes
    -----
    * Regex substitution is performed under the hood with ``re.sub``. The
      rules for substitution for ``re.sub`` are the same.
    * Regular expressions will only substitute on strings, meaning you
      cannot provide, for example, a regular expression matching floating
      point numbers and expect the columns in your frame that have a
      numeric dtype to be matched. However, if those floating point
      numbers *are* strings, then you can do this.
    * This method has *a lot* of options. You are encouraged to experiment
      and play with this method to gain intuition about how it works.
    * When dict is used as the `to_replace` value, it is like
      key(s) in the dict are the to_replace part and
      value(s) in the dict are the value parameter.

    Examples
    --------

    **Scalar `to_replace` and `value`**

    >>> s = pd.Series([1, 2, 3, 4, 5])
    >>> s.replace(1, 5)
    0    5
    1    2
    2    3
    3    4
    4    5
    dtype: int64

    >>> df = pd.DataFrame({{'A': [0, 1, 2, 3, 4],
    ...                    'B': [5, 6, 7, 8, 9],
    ...                    'C': ['a', 'b', 'c', 'd', 'e']}})
    >>> df.replace(0, 5)
        A  B  C
    0  5  5  a
    1  1  6  b
    2  2  7  c
    3  3  8  d
    4  4  9  e

    **List-like `to_replace`**

    >>> df.replace([0, 1, 2, 3], 4)
        A  B  C
    0  4  5  a
    1  4  6  b
    2  4  7  c
    3  4  8  d
    4  4  9  e

    >>> df.replace([0, 1, 2, 3], [4, 3, 2, 1])
        A  B  C
    0  4  5  a
    1  3  6  b
    2  2  7  c
    3  1  8  d
    4  4  9  e

    >>> s.replace([1, 2], method='bfill')
    0    3
    1    3
    2    3
    3    4
    4    5
    dtype: int64

    **dict-like `to_replace`**

    >>> df.replace({{0: 10, 1: 100}})
            A  B  C
    0   10  5  a
    1  100  6  b
    2    2  7  c
    3    3  8  d
    4    4  9  e

    >>> df.replace({{'A': 0, 'B': 5}}, 100)
            A    B  C
    0  100  100  a
    1    1    6  b
    2    2    7  c
    3    3    8  d
    4    4    9  e

    >>> df.replace({{'A': {{0: 100, 4: 400}}}})
            A  B  C
    0  100  5  a
    1    1  6  b
    2    2  7  c
    3    3  8  d
    4  400  9  e

    **Regular expression `to_replace`**

    >>> df = pd.DataFrame({{'A': ['bat', 'foo', 'bait'],
    ...                    'B': ['abc', 'bar', 'xyz']}})
    >>> df.replace(to_replace=r'^ba.$', value='new', regex=True)
            A    B
    0   new  abc
    1   foo  new
    2  bait  xyz

    >>> df.replace({{'A': r'^ba.$'}}, {{'A': 'new'}}, regex=True)
            A    B
    0   new  abc
    1   foo  bar
    2  bait  xyz

    >>> df.replace(regex=r'^ba.$', value='new')
            A    B
    0   new  abc
    1   foo  new
    2  bait  xyz

    >>> df.replace(regex={{r'^ba.$': 'new', 'foo': 'xyz'}})
            A    B
    0   new  abc
    1   xyz  new
    2  bait  xyz

    >>> df.replace(regex=[r'^ba.$', 'foo'], value='new')
            A    B
    0   new  abc
    1   new  new
    2  bait  xyz

    Compare the behavior of ``s.replace({{'a': None}})`` and
    ``s.replace('a', None)`` to understand the peculiarities
    of the `to_replace` parameter:

    >>> s = pd.Series([10, 'a', 'a', 'b', 'a'])

    When one uses a dict as the `to_replace` value, it is like the
    value(s) in the dict are equal to the `value` parameter.
    ``s.replace({{'a': None}})`` is equivalent to
    ``s.replace(to_replace={{'a': None}}, value=None, method=None)``:

    >>> s.replace({{'a': None}})
    0      10
    1    None
    2    None
    3       b
    4    None
    dtype: object

    When ``value`` is not explicitly passed and `to_replace` is a scalar, list
    or tuple, `replace` uses the method parameter (default 'pad') to do the
    replacement. So this is why the 'a' values are being replaced by 10
    in rows 1 and 2 and 'b' in row 4 in this case.

    >>> s.replace('a')
    0    10
    1    10
    2    10
    3     b
    4     b
    dtype: object

    On the other hand, if ``None`` is explicitly passed for ``value``, it will
    be respected:

    >>> s.replace('a', None)
    0      10
    1    None
    2    None
    3       b
    4    None
    dtype: object

        .. versionchanged:: 1.4.0
            Previously the explicit ``None`` was silently ignored.
"""

_shared_docs[
    "idxmin"
] = """
    Return index of first occurrence of minimum over requested axis.

    NA/null values are excluded.

    Parameters
    ----------
    axis : {{0 or 'index', 1 or 'columns'}}, default 0
        The axis to use. 0 or 'index' for row-wise, 1 or 'columns' for column-wise.
    skipna : bool, default True
        Exclude NA/null values. If an entire row/column is NA, the result
        will be NA.
    numeric_only : bool, default {numeric_only_default}
        Include only `float`, `int` or `boolean` data.

        .. versionadded:: 1.5.0

    Returns
    -------
    Series
        Indexes of minima along the specified axis.

    Raises
    ------
    ValueError
        * If the row/column is empty

    See Also
    --------
    Series.idxmin : Return index of the minimum element.

    Notes
    -----
    This method is the DataFrame version of ``ndarray.argmin``.

    Examples
    --------
    Consider a dataset containing food consumption in Argentina.

    >>> df = pd.DataFrame({{'consumption': [10.51, 103.11, 55.48],
    ...                     'co2_emissions': [37.2, 19.66, 1712]}},
    ...                   index=['Pork', 'Wheat Products', 'Beef'])

    >>> df
                    consumption  co2_emissions
    Pork                  10.51         37.20
    Wheat Products       103.11         19.66
    Beef                  55.48       1712.00

    By default, it returns the index for the minimum value in each column.

    >>> df.idxmin()
    consumption                Pork
    co2_emissions    Wheat Products
    dtype: object

    To return the index for the minimum value in each row, use ``axis="columns"``.

    >>> df.idxmin(axis="columns")
    Pork                consumption
    Wheat Products    co2_emissions
    Beef                consumption
    dtype: object
"""

_shared_docs[
    "idxmax"
] = """
    Return index of first occurrence of maximum over requested axis.

    NA/null values are excluded.

    Parameters
    ----------
    axis : {{0 or 'index', 1 or 'columns'}}, default 0
        The axis to use. 0 or 'index' for row-wise, 1 or 'columns' for column-wise.
    skipna : bool, default True
        Exclude NA/null values. If an entire row/column is NA, the result
        will be NA.
    numeric_only : bool, default {numeric_only_default}
        Include only `float`, `int` or `boolean` data.

        .. versionadded:: 1.5.0

    Returns
    -------
    Series
        Indexes of maxima along the specified axis.

    Raises
    ------
    ValueError
        * If the row/column is empty

    See Also
    --------
    Series.idxmax : Return index of the maximum element.

    Notes
    -----
    This method is the DataFrame version of ``ndarray.argmax``.

    Examples
    --------
    Consider a dataset containing food consumption in Argentina.

    >>> df = pd.DataFrame({{'consumption': [10.51, 103.11, 55.48],
    ...                     'co2_emissions': [37.2, 19.66, 1712]}},
    ...                   index=['Pork', 'Wheat Products', 'Beef'])

    >>> df
                    consumption  co2_emissions
    Pork                  10.51         37.20
    Wheat Products       103.11         19.66
    Beef                  55.48       1712.00

    By default, it returns the index for the maximum value in each column.

    >>> df.idxmax()
    consumption     Wheat Products
    co2_emissions             Beef
    dtype: object

    To return the index for the maximum value in each row, use ``axis="columns"``.

    >>> df.idxmax(axis="columns")
    Pork              co2_emissions
    Wheat Products     consumption
    Beef              co2_emissions
    dtype: object
"""

window_doc = """
Provide rolling window calculations.

Parameters
----------
window : int, timedelta, str, offset, or BaseIndexer subclass
    Size of the moving window.

    If an integer, the fixed number of observations used for
    each window.

    If a timedelta, str, or offset, the time period of each window. Each
    window will be a variable sized based on the observations included in
    the time-period. This is only valid for datetimelike indexes.
    To learn more about the offsets & frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.

    If a BaseIndexer subclass, the window boundaries
    based on the defined ``get_window_bounds`` method. Additional rolling
    keyword arguments, namely ``min_periods``, ``center``, ``closed`` and
    ``step`` will be passed to ``get_window_bounds``.

min_periods : int, default None
    Minimum number of observations in window required to have a value;
    otherwise, result is ``np.nan``.

    For a window that is specified by an offset, ``min_periods`` will default to 1.

    For a window that is specified by an integer, ``min_periods`` will default
    to the size of the window.

center : bool, default False
    If False, set the window labels as the right edge of the window index.

    If True, set the window labels as the center of the window index.

win_type : str, default None
    If ``None``, all points are evenly weighted.

    If a string, it must be a valid `scipy.signal window function
    <https://docs.scipy.org/doc/scipy/reference/signal.windows.html#module-scipy.signal.windows>`__.

    Certain Scipy window types require additional parameters to be passed
    in the aggregation function. The additional parameters must match
    the keywords specified in the Scipy window type method signature.

on : str, optional
    For a DataFrame, a column label or Index level on which
    to calculate the rolling window, rather than the DataFrame's index.

    Provided integer column is ignored and excluded from result since
    an integer index is not used to calculate the rolling window.

axis : int or str, default 0
    If ``0`` or ``'index'``, roll across the rows.

    If ``1`` or ``'columns'``, roll across the columns.

    For `Series` this parameter is unused and defaults to 0.

closed : str, default None
    If ``'right'``, the first point in the window is excluded from calculations.

    If ``'left'``, the last point in the window is excluded from calculations.

    If ``'both'``, the no points in the window are excluded from calculations.

    If ``'neither'``, the first and last points in the window are excluded
    from calculations.

    Default ``None`` (``'right'``).

    .. versionchanged:: 1.2.0

        The closed parameter with fixed windows is now supported.

step : int, default None

    .. versionadded:: 1.5.0

    Evaluate the window at every ``step`` result, equivalent to slicing as
    ``[::step]``. ``window`` must be an integer. Using a step argument other
    than None or 1 will produce a result with a different shape than the input.

method : str {'single', 'table'}, default 'single'

    .. versionadded:: 1.3.0

    Execute the rolling operation per single column or row (``'single'``)
    or over the entire object (``'table'``).

    This argument is only implemented when specifying ``engine='numba'``
    in the method call.

Returns
-------
``Window`` subclass if a ``win_type`` is passed

``Rolling`` subclass if ``win_type`` is not passed

See Also
--------
expanding : Provides expanding transformations.
ewm : Provides exponential weighted functions.

Notes
-----
See :ref:`Windowing Operations <window.generic>` for further usage details
and examples.

Examples
--------
>>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]})
>>> df
     B
0  0.0
1  1.0
2  2.0
3  NaN
4  4.0

**window**

Rolling sum with a window length of 2 observations.

>>> df.rolling(2).sum()
     B
0  NaN
1  1.0
2  3.0
3  NaN
4  NaN

Rolling sum with a window span of 2 seconds.

>>> df_time = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]},
...                        index = [pd.Timestamp('20130101 09:00:00'),
...                                 pd.Timestamp('20130101 09:00:02'),
...                                 pd.Timestamp('20130101 09:00:03'),
...                                 pd.Timestamp('20130101 09:00:05'),
...                                 pd.Timestamp('20130101 09:00:06')])

>>> df_time
                       B
2013-01-01 09:00:00  0.0
2013-01-01 09:00:02  1.0
2013-01-01 09:00:03  2.0
2013-01-01 09:00:05  NaN
2013-01-01 09:00:06  4.0

>>> df_time.rolling('2s').sum()
                       B
2013-01-01 09:00:00  0.0
2013-01-01 09:00:02  1.0
2013-01-01 09:00:03  3.0
2013-01-01 09:00:05  NaN
2013-01-01 09:00:06  4.0

Rolling sum with forward looking windows with 2 observations.

>>> indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=2)
>>> df.rolling(window=indexer, min_periods=1).sum()
     B
0  1.0
1  3.0
2  2.0
3  4.0
4  4.0

**min_periods**

Rolling sum with a window length of 2 observations, but only needs a minimum of 1
observation to calculate a value.

>>> df.rolling(2, min_periods=1).sum()
     B
0  0.0
1  1.0
2  3.0
3  2.0
4  4.0

**center**

Rolling sum with the result assigned to the center of the window index.

>>> df.rolling(3, min_periods=1, center=True).sum()
     B
0  1.0
1  3.0
2  3.0
3  6.0
4  4.0

>>> df.rolling(3, min_periods=1, center=False).sum()
     B
0  0.0
1  1.0
2  3.0
3  3.0
4  6.0

**step**

Rolling sum with a window length of 2 observations, minimum of 1 observation to
calculate a value, and a step of 2.

>>> df.rolling(2, min_periods=1, step=2).sum()
     B
0  0.0
2  3.0
4  4.0

**win_type**

Rolling sum with a window length of 2, using the Scipy ``'gaussian'``
window type. ``std`` is required in the aggregation function.

>>> df.rolling(2, win_type='gaussian').sum(std=3)
          B
0       NaN
1  0.986207
2  2.958621
3       NaN
4       NaN

**on**

Rolling sum with a window length of 2 days.

>>> df = pd.DataFrame({
...     'A': [pd.to_datetime('2020-01-01'),
...           pd.to_datetime('2020-01-01'),
...           pd.to_datetime('2020-01-02'),],
...     'B': [1, 2, 3], },
...     index=pd.date_range('2020', periods=3))

>>> df
                    A  B
2020-01-01 2020-01-01  1
2020-01-02 2020-01-01  2
2020-01-03 2020-01-02  3

>>> df.rolling('2D', on='A').sum()
                    A    B
2020-01-01 2020-01-01  1.0
2020-01-02 2020-01-01  3.0
2020-01-03 2020-01-02  6.0
"""


expanding_doc = """
Provide expanding window calculations.

Parameters
----------
min_periods : int, default 1
    Minimum number of observations in window required to have a value;
    otherwise, result is ``np.nan``.

axis : int or str, default 0
    If ``0`` or ``'index'``, roll across the rows.

    If ``1`` or ``'columns'``, roll across the columns.

    For `Series` this parameter is unused and defaults to 0.

method : str {'single', 'table'}, default 'single'
    Execute the rolling operation per single column or row (``'single'``)
    or over the entire object (``'table'``).

    This argument is only implemented when specifying ``engine='numba'``
    in the method call.

    .. versionadded:: 1.3.0

Returns
-------
``Expanding`` subclass

See Also
--------
rolling : Provides rolling window calculations.
ewm : Provides exponential weighted functions.

Notes
-----
See :ref:`Windowing Operations <window.expanding>` for further usage details
and examples.

Examples
--------
>>> df = pd.DataFrame({"B": [0, 1, 2, np.nan, 4]})
>>> df
     B
0  0.0
1  1.0
2  2.0
3  NaN
4  4.0

**min_periods**

Expanding sum with 1 vs 3 observations needed to calculate a value.

>>> df.expanding(1).sum()
     B
0  0.0
1  1.0
2  3.0
3  3.0
4  7.0
>>> df.expanding(3).sum()
     B
0  NaN
1  NaN
2  3.0
3  3.0
4  7.0
"""

exponential_moving_window_doc = r"""
Provide exponentially weighted (EW) calculations.

Exactly one of ``com``, ``span``, ``halflife``, or ``alpha`` must be
provided if ``times`` is not provided. If ``times`` is provided,
``halflife`` and one of ``com``, ``span`` or ``alpha`` may be provided.

Parameters
----------
com : float, optional
    Specify decay in terms of center of mass

    :math:`\alpha = 1 / (1 + com)`, for :math:`com \geq 0`.

span : float, optional
    Specify decay in terms of span

    :math:`\alpha = 2 / (span + 1)`, for :math:`span \geq 1`.

halflife : float, str, timedelta, optional
    Specify decay in terms of half-life

    :math:`\alpha = 1 - \exp\left(-\ln(2) / halflife\right)`, for
    :math:`halflife > 0`.

    If ``times`` is specified, a timedelta convertible unit over which an
    observation decays to half its value. Only applicable to ``mean()``,
    and halflife value will not apply to the other functions.

    .. versionadded:: 1.1.0

alpha : float, optional
    Specify smoothing factor :math:`\alpha` directly

    :math:`0 < \alpha \leq 1`.

min_periods : int, default 0
    Minimum number of observations in window required to have a value;
    otherwise, result is ``np.nan``.

adjust : bool, default True
    Divide by decaying adjustment factor in beginning periods to account
    for imbalance in relative weightings (viewing EWMA as a moving average).

    - When ``adjust=True`` (default), the EW function is calculated using weights
      :math:`w_i = (1 - \alpha)^i`. For example, the EW moving average of the series
      [:math:`x_0, x_1, ..., x_t`] would be:

    .. math::
        y_t = \frac{x_t + (1 - \alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ... + (1 -
        \alpha)^t x_0}{1 + (1 - \alpha) + (1 - \alpha)^2 + ... + (1 - \alpha)^t}

    - When ``adjust=False``, the exponentially weighted function is calculated
      recursively:

    .. math::
        \begin{split}
            y_0 &= x_0\\
            y_t &= (1 - \alpha) y_{t-1} + \alpha x_t,
        \end{split}
ignore_na : bool, default False
    Ignore missing values when calculating weights.

    - When ``ignore_na=False`` (default), weights are based on absolute positions.
      For example, the weights of :math:`x_0` and :math:`x_2` used in calculating
      the final weighted average of [:math:`x_0`, None, :math:`x_2`] are
      :math:`(1-\alpha)^2` and :math:`1` if ``adjust=True``, and
      :math:`(1-\alpha)^2` and :math:`\alpha` if ``adjust=False``.

    - When ``ignore_na=True``, weights are based
      on relative positions. For example, the weights of :math:`x_0` and :math:`x_2`
      used in calculating the final weighted average of
      [:math:`x_0`, None, :math:`x_2`] are :math:`1-\alpha` and :math:`1` if
      ``adjust=True``, and :math:`1-\alpha` and :math:`\alpha` if ``adjust=False``.

axis : {0, 1}, default 0
    If ``0`` or ``'index'``, calculate across the rows.

    If ``1`` or ``'columns'``, calculate across the columns.

    For `Series` this parameter is unused and defaults to 0.

times : np.ndarray, Series, default None

    .. versionadded:: 1.1.0

    Only applicable to ``mean()``.

    Times corresponding to the observations. Must be monotonically increasing and
    ``datetime64[ns]`` dtype.

    If 1-D array like, a sequence with the same shape as the observations.

method : str {'single', 'table'}, default 'single'
    .. versionadded:: 1.4.0

    Execute the rolling operation per single column or row (``'single'``)
    or over the entire object (``'table'``).

    This argument is only implemented when specifying ``engine='numba'``
    in the method call.

    Only applicable to ``mean()``

Returns
-------
``ExponentialMovingWindow`` subclass

See Also
--------
rolling : Provides rolling window calculations.
expanding : Provides expanding transformations.

Notes
-----
See :ref:`Windowing Operations <window.exponentially_weighted>`
for further usage details and examples.

Examples
--------
>>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]})
>>> df
     B
0  0.0
1  1.0
2  2.0
3  NaN
4  4.0

>>> df.ewm(com=0.5).mean()
          B
0  0.000000
1  0.750000
2  1.615385
3  1.615385
4  3.670213
>>> df.ewm(alpha=2 / 3).mean()
          B
0  0.000000
1  0.750000
2  1.615385
3  1.615385
4  3.670213

**adjust**

>>> df.ewm(com=0.5, adjust=True).mean()
          B
0  0.000000
1  0.750000
2  1.615385
3  1.615385
4  3.670213
>>> df.ewm(com=0.5, adjust=False).mean()
          B
0  0.000000
1  0.666667
2  1.555556
3  1.555556
4  3.650794

**ignore_na**

>>> df.ewm(com=0.5, ignore_na=True).mean()
          B
0  0.000000
1  0.750000
2  1.615385
3  1.615385
4  3.225000
>>> df.ewm(com=0.5, ignore_na=False).mean()
          B
0  0.000000
1  0.750000
2  1.615385
3  1.615385
4  3.670213

**times**

Exponentially weighted mean with weights calculated with a timedelta ``halflife``
relative to ``times``.

>>> times = ['2020-01-01', '2020-01-03', '2020-01-10', '2020-01-15', '2020-01-17']
>>> df.ewm(halflife='4 days', times=pd.DatetimeIndex(times)).mean()
          B
0  0.000000
1  0.585786
2  1.523889
3  1.523889
4  3.233686
"""
