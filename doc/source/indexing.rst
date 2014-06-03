.. _indexing:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   np.random.seed(123456)
   from pandas import *
   options.display.max_rows=15
   import pandas as pd
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import range, zip

***************************
Indexing and Selecting Data
***************************

The axis labeling information in pandas objects serves many purposes:

  - Identifies data (i.e. provides *metadata*) using known indicators,
    important for analysis, visualization, and interactive console display
  - Enables automatic and explicit data alignment
  - Allows intuitive getting and setting of subsets of the data set

In this section, we will focus on the final point: namely, how to slice, dice,
and generally get and set subsets of pandas objects. The primary focus will be
on Series and DataFrame as they have received more development attention in
this area. Expect more work to be invested higher-dimensional data structures
(including ``Panel``) in the future, especially in label-based advanced
indexing.

.. note::

   The Python and NumPy indexing operators ``[]`` and attribute operator ``.``
   provide quick and easy access to pandas data structures across a wide range
   of use cases. This makes interactive work intuitive, as there's little new
   to learn if you already know how to deal with Python dictionaries and NumPy
   arrays. However, since the type of the data to be accessed isn't known in
   advance, directly using standard operators has some optimization limits. For
   production code, we recommended that you take advantage of the optimized
   pandas data access methods exposed in this chapter.

.. warning::

   Whether a copy or a reference is returned for a setting operation, may
   depend on the context.  This is sometimes called ``chained assignment`` and
   should be avoided.  See :ref:`Returning a View versus Copy
   <indexing.view_versus_copy>`

See the :ref:`cookbook<cookbook.selection>` for some advanced strategies

Different Choices for Indexing (``loc``, ``iloc``, and ``ix``)
--------------------------------------------------------------

.. versionadded:: 0.11.0

Object selection has had a number of user-requested additions in order to
support more explicit location based indexing. pandas now supports three types
of multi-axis indexing.

- ``.loc`` is strictly label based, will raise ``KeyError`` when the items are
  not found, allowed inputs are:

  - A single label, e.g. ``5`` or ``'a'``, (note that ``5`` is interpreted as a
    *label* of the index. This use is **not** an integer position along the
    index)
  - A list or array of labels ``['a', 'b', 'c']``
  - A slice object with labels ``'a':'f'``, (note that contrary to usual python
    slices, **both** the start and the stop are included!)
  - A boolean array

  See more at :ref:`Selection by Label <indexing.label>`

- ``.iloc`` is strictly integer position based (from ``0`` to ``length-1`` of
  the axis), will raise ``IndexError`` if an indexer is requested and it
  is out-of-bounds, except *slice* indexers which allow out-of-bounds indexing.
  (this conforms with python/numpy *slice* semantics). Allowed inputs are:

  - An integer e.g. ``5``
  - A list or array of integers ``[4, 3, 0]``
  - A slice object with ints ``1:7``

  See more at :ref:`Selection by Position <indexing.integer>`

- ``.ix`` supports mixed integer and label based access. It is primarily label
  based, but will fallback to integer positional access. ``.ix`` is the most
  general and will support any of the inputs to ``.loc`` and ``.iloc``, as well
  as support for floating point label schemes. ``.ix`` is especially useful
  when dealing with mixed positional and label based hierarchial indexes.
  As using integer slices with ``.ix`` have different behavior depending on
  whether the slice is interpreted as position based or label based, it's
  usually better to be explicit and use ``.iloc`` or ``.loc``.

  See more at :ref:`Advanced Indexing <indexing.advanced>`, :ref:`Advanced
  Hierarchical <indexing.advanced_hierarchical>` and :ref:`Fallback Indexing
  <indexing.fallback>`

Getting values from an object with multi-axes selection uses the following
notation (using ``.loc`` as an example, but applies to ``.iloc`` and ``.ix`` as
well). Any of the axes accessors may be the null slice ``:``. Axes left out of
the specification are assumed to be ``:``. (e.g. ``p.loc['a']`` is equiv to
``p.loc['a', :, :]``)

.. csv-table::
    :header: "Object Type", "Indexers"
    :widths: 30, 50
    :delim: ;

    Series; ``s.loc[indexer]``
    DataFrame; ``df.loc[row_indexer,column_indexer]``
    Panel; ``p.loc[item_indexer,major_indexer,minor_indexer]``

Deprecations
------------

Beginning with version 0.11.0, it's recommended that you transition away from
the following methods as they *may* be deprecated in future versions.

  - ``irow``
  - ``icol``
  - ``iget_value``

See the section :ref:`Selection by Position <indexing.integer>` for substitutes.

.. _indexing.basics:

Basics
------

As mentioned when introducing the data structures in the :ref:`last section
<basics>`, the primary function of indexing with ``[]`` (a.k.a. ``__getitem__``
for those familiar with implementing class behavior in Python) is selecting out
lower-dimensional slices. Thus,

.. csv-table::
    :header: "Object Type", "Selection", "Return Value Type"
    :widths: 30, 30, 60
    :delim: ;

    Series; ``series[label]``; scalar value
    DataFrame; ``frame[colname]``; ``Series`` corresponding to colname
    Panel; ``panel[itemname]``; ``DataFrame`` corresponing to the itemname

Here we construct a simple time series data set to use for illustrating the
indexing functionality:

.. ipython:: python

   dates = date_range('1/1/2000', periods=8)
   df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
   df
   panel = Panel({'one' : df, 'two' : df - df.mean()})
   panel

.. note::

   None of the indexing functionality is time series specific unless
   specifically stated.

Thus, as per above, we have the most basic indexing using ``[]``:

.. ipython:: python

   s = df['A']
   s[dates[5]]
   panel['two']

You can pass a list of columns to ``[]`` to select columns in that order.
If a column is not contained in the DataFrame, an exception will be
raised. Multiple columns can also be set in this manner:

.. ipython:: python

   df
   df[['B', 'A']] = df[['A', 'B']]
   df

You may find this useful for applying a transform (in-place) to a subset of the
columns.

Attribute Access
----------------

.. _indexing.columns.multiple:

.. _indexing.df_cols:

.. _indexing.attribute_access:

You may access an index on a ``Series``, column on a ``DataFrame``, and a item on a ``Panel`` directly
as an attribute:

.. ipython:: python

   sa = Series([1,2,3],index=list('abc'))
   dfa = df.copy()

.. ipython:: python

   sa.b
   dfa.A
   panel.one

You can use attribute access to modify an existing element of a Series or column of a DataFrame, but be careful;
if you try to use attribute access to create a new column, it fails silently, creating a new attribute rather than a
new column.

.. ipython:: python

   sa.a = 5
   sa
   dfa.A = list(range(len(dfa.index)))       # ok if A already exists
   dfa
   dfa['A'] = list(range(len(dfa.index)))    # use this form to create a new column
   dfa

.. warning::

   - You can use this access only if the index element is a valid python identifier, e.g. ``s.1`` is not allowed.
     see `here for an explanation of valid identifiers
     <http://docs.python.org/2.7/reference/lexical_analysis.html#identifiers>`__.

   - The attribute will not be available if it conflicts with an existing method name, e.g. ``s.min`` is not allowed.

   - The ``Series/Panel`` accesses are available starting in 0.13.0.

If you are using the IPython environment, you may also use tab-completion to
see these accessable attributes.

Slicing ranges
--------------

The most robust and consistent way of slicing ranges along arbitrary axes is
described in the :ref:`Selection by Position <indexing.integer>` section
detailing the ``.iloc`` method. For now, we explain the semantics of slicing using the ``[]`` operator.

With Series, the syntax works exactly as with an ndarray, returning a slice of
the values and the corresponding labels:

.. ipython:: python

   s[:5]
   s[::2]
   s[::-1]

Note that setting works as well:

.. ipython:: python

   s2 = s.copy()
   s2[:5] = 0
   s2

With DataFrame, slicing inside of ``[]`` **slices the rows**. This is provided
largely as a convenience since it is such a common operation.

.. ipython:: python

   df[:3]
   df[::-1]

.. _indexing.label:

Selection By Label
------------------

.. warning::

   Whether a copy or a reference is returned for a setting operation, may depend on the context.
   This is sometimes called ``chained assignment`` and should be avoided.
   See :ref:`Returning a View versus Copy <indexing.view_versus_copy>`

pandas provides a suite of methods in order to have **purely label based indexing**. This is a strict inclusion based protocol.
**ALL** of the labels for which you ask, must be in the index or a ``KeyError`` will be raised! When slicing, the start bound is *included*, **AND** the stop bound is *included*. Integers are valid labels, but they refer to the label **and not the position**.

The ``.loc`` attribute is the primary access method. The following are valid inputs:

- A single label, e.g. ``5`` or ``'a'``, (note that ``5`` is interpreted as a *label* of the index. This use is **not** an integer position along the index)
- A list or array of labels ``['a', 'b', 'c']``
- A slice object with labels ``'a':'f'`` (note that contrary to usual python slices, **both** the start and the stop are included!)
- A boolean array

.. ipython:: python

   s1 = Series(np.random.randn(6),index=list('abcdef'))
   s1
   s1.loc['c':]
   s1.loc['b']

Note that setting works as well:

.. ipython:: python

   s1.loc['c':] = 0
   s1

With a DataFrame

.. ipython:: python

   df1 = DataFrame(np.random.randn(6,4),
                   index=list('abcdef'),
                   columns=list('ABCD'))
   df1
   df1.loc[['a','b','d'],:]

Accessing via label slices

.. ipython:: python

   df1.loc['d':,'A':'C']

For getting a cross section using a label (equiv to ``df.xs('a')``)

.. ipython:: python

   df1.loc['a']

For getting values with a boolean array

.. ipython:: python

   df1.loc['a']>0
   df1.loc[:,df1.loc['a']>0]

For getting a value explicity (equiv to deprecated ``df.get_value('a','A')``)

.. ipython:: python

   # this is also equivalent to ``df1.at['a','A']``
   df1.loc['a','A']

.. _indexing.integer:

Selection By Position
---------------------

.. warning::

   Whether a copy or a reference is returned for a setting operation, may depend on the context.
   This is sometimes called ``chained assignment`` and should be avoided.
   See :ref:`Returning a View versus Copy <indexing.view_versus_copy>`

pandas provides a suite of methods in order to get **purely integer based indexing**. The semantics follow closely python and numpy slicing. These are ``0-based`` indexing. When slicing, the start bounds is *included*, while the upper bound is *excluded*. Trying to use a non-integer, even a **valid** label will raise a ``IndexError``.

The ``.iloc`` attribute is the primary access method. The following are valid inputs:

- An integer e.g. ``5``
- A list or array of integers ``[4, 3, 0]``
- A slice object with ints ``1:7``

.. ipython:: python

   s1 = Series(np.random.randn(5),index=list(range(0,10,2)))
   s1
   s1.iloc[:3]
   s1.iloc[3]

Note that setting works as well:

.. ipython:: python

   s1.iloc[:3] = 0
   s1

With a DataFrame

.. ipython:: python

   df1 = DataFrame(np.random.randn(6,4),
                   index=list(range(0,12,2)),
                   columns=list(range(0,8,2)))
   df1

Select via integer slicing

.. ipython:: python

   df1.iloc[:3]
   df1.iloc[1:5,2:4]

Select via integer list

.. ipython:: python

   df1.iloc[[1,3,5],[1,3]]

For slicing rows explicitly (equiv to deprecated ``df.irow(slice(1,3))``).

.. ipython:: python

   df1.iloc[1:3,:]

For slicing columns explicitly (equiv to deprecated ``df.icol(slice(1,3))``).

.. ipython:: python

   df1.iloc[:,1:3]

For getting a scalar via integer position (equiv to deprecated ``df.get_value(1,1)``)

.. ipython:: python

   # this is also equivalent to ``df1.iat[1,1]``
   df1.iloc[1,1]

For getting a cross section using an integer position (equiv to ``df.xs(1)``)

.. ipython:: python

   df1.iloc[1]

There is one signficant departure from standard python/numpy slicing semantics.
python/numpy allow slicing past the end of an array without an associated error.

.. ipython:: python

    # these are allowed in python/numpy.
    x = list('abcdef')
    x[4:10]
    x[8:10]

- as of v0.14.0, ``iloc`` will now accept out-of-bounds indexers for slices, e.g. a value that exceeds the length of the object being
  indexed. These will be excluded. This will make pandas conform more with pandas/numpy indexing of out-of-bounds
  values. A single indexer / list of indexers that is out-of-bounds will still raise
  ``IndexError`` (:issue:`6296`, :issue:`6299`). This could result in an empty axis (e.g. an empty DataFrame being returned)

.. ipython:: python

   dfl = DataFrame(np.random.randn(5,2),columns=list('AB'))
   dfl
   dfl.iloc[:,2:3]
   dfl.iloc[:,1:3]
   dfl.iloc[4:6]

These are out-of-bounds selections

.. code-block:: python

   dfl.iloc[[4,5,6]]
   IndexError: positional indexers are out-of-bounds

   dfl.iloc[:,4]
   IndexError: single positional indexer is out-of-bounds

.. _indexing.basics.partial_setting:

Setting With Enlargement
------------------------

.. versionadded:: 0.13

The ``.loc/.ix/[]`` operations can perform enlargement when setting a non-existant key for that axis.

In the ``Series`` case this is effectively an appending operation

.. ipython:: python

   se = Series([1,2,3])
   se
   se[5] = 5.
   se

A ``DataFrame`` can be enlarged on either axis via ``.loc``

.. ipython:: python

   dfi = DataFrame(np.arange(6).reshape(3,2),
                   columns=['A','B'])
   dfi
   dfi.loc[:,'C'] = dfi.loc[:,'A']
   dfi

This is like an ``append`` operation on the ``DataFrame``.

.. ipython:: python

   dfi.loc[3] = 5
   dfi

.. _indexing.basics.get_value:

Fast scalar value getting and setting
-------------------------------------

Since indexing with ``[]`` must handle a lot of cases (single-label access,
slicing, boolean indexing, etc.), it has a bit of overhead in order to figure
out what you're asking for. If you only want to access a scalar value, the
fastest way is to use the ``at`` and ``iat`` methods, which are implemented on
all of the data structures.

Similary to ``loc``, ``at`` provides **label** based scalar lookups, while, ``iat`` provides **integer** based lookups analagously to ``iloc``

.. ipython:: python

   s.iat[5]
   df.at[dates[5], 'A']
   df.iat[3, 0]

You can also set using these same indexers.

.. ipython:: python

   df.at[dates[5], 'E'] = 7
   df.iat[3, 0] = 7

``at`` may enlarge the object in-place as above if the indexer is missing.

.. ipython:: python

   df.at[dates[-1]+1, 0] = 7
   df

Boolean indexing
----------------

.. _indexing.boolean:

Another common operation is the use of boolean vectors to filter the data.
The operators are: ``|`` for ``or``, ``&`` for ``and``, and ``~`` for ``not``. These **must** be grouped by using parentheses.

Using a boolean vector to index a Series works exactly as in a numpy ndarray:

.. ipython:: python

   s[s > 0]
   s[(s < 0) & (s > -0.5)]
   s[(s < -1) | (s > 1 )]
   s[~(s < 0)]

You may select rows from a DataFrame using a boolean vector the same length as
the DataFrame's index (for example, something derived from one of the columns
of the DataFrame):

.. ipython:: python

   df[df['A'] > 0]

List comprehensions and ``map`` method of Series can also be used to produce
more complex criteria:

.. ipython:: python

   df2 = DataFrame({'a' : ['one', 'one', 'two', 'three', 'two', 'one', 'six'],
                    'b' : ['x', 'y', 'y', 'x', 'y', 'x', 'x'],
                    'c' : randn(7)})

   # only want 'two' or 'three'
   criterion = df2['a'].map(lambda x: x.startswith('t'))

   df2[criterion]

   # equivalent but slower
   df2[[x.startswith('t') for x in df2['a']]]

   # Multiple criteria
   df2[criterion & (df2['b'] == 'x')]

Note, with the choice methods :ref:`Selection by Label <indexing.label>`, :ref:`Selection by Position <indexing.integer>`,
and :ref:`Advanced Indexing <indexing.advanced>` you may select along more than one axis using boolean vectors combined with other indexing expressions.

.. ipython:: python

   df2.loc[criterion & (df2['b'] == 'x'),'b':'c']

.. _indexing.basics.indexing_isin:

Indexing with isin
~~~~~~~~~~~~~~~~~~

Consider the ``isin`` method of Series, which returns a boolean vector that is
true wherever the Series elements exist in the passed list. This allows you to
select rows where one or more columns have values you want:

.. ipython:: python

   s = Series(np.arange(5),index=np.arange(5)[::-1],dtype='int64')

   s

   s.isin([2, 4])

   s[s.isin([2, 4])]


DataFrame also has an ``isin`` method.  When calling ``isin``, pass a set of
values as either an array or dict.  If values is an array, ``isin`` returns
a DataFrame of booleans that is the same shape as the original DataFrame, with True
wherever the element is in the sequence of values.

.. ipython:: python

   df = DataFrame({'vals': [1, 2, 3, 4], 'ids': ['a', 'b', 'f', 'n'],
                   'ids2': ['a', 'n', 'c', 'n']})

   values = ['a', 'b', 1, 3]

   df.isin(values)

Oftentimes you'll want to match certain values with certain columns.
Just make values a ``dict`` where the key is the column, and the value is
a list of items you want to check for.

.. ipython:: python

   values = {'ids': ['a', 'b'], 'vals': [1, 3]}

   df.isin(values)

Combine DataFrame's ``isin`` with the ``any()`` and ``all()`` methods to
quickly select subsets of your data that meet a given criteria.
To select a row where each column meets its own criterion:

.. ipython:: python

  values = {'ids': ['a', 'b'], 'ids2': ['a', 'c'], 'vals': [1, 3]}

  row_mask = df.isin(values).all(1)

  df[row_mask]

The :meth:`~pandas.DataFrame.where` Method and Masking
------------------------------------------------------

Selecting values from a Series with a boolean vector generally returns a
subset of the data. To guarantee that selection output has the same shape as
the original data, you can use the ``where`` method in ``Series`` and ``DataFrame``.

To return only the selected rows

.. ipython:: python

   s[s > 0]

To return a Series of the same shape as the original

.. ipython:: python

   s.where(s > 0)

Selecting values from a DataFrame with a boolean critierion now also preserves
input data shape. ``where`` is used under the hood as the implementation.
Equivalent is ``df.where(df < 0)``

.. ipython:: python
   :suppress:

   dates = date_range('1/1/2000', periods=8)
   df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])

.. ipython:: python

   df[df < 0]

In addition, ``where`` takes an optional ``other`` argument for replacement of
values where the condition is False, in the returned copy.

.. ipython:: python

   df.where(df < 0, -df)

You may wish to set values based on some boolean criteria.
This can be done intuitively like so:

.. ipython:: python

   s2 = s.copy()
   s2[s2 < 0] = 0
   s2

   df2 = df.copy()
   df2[df2 < 0] = 0
   df2

By default, ``where`` returns a modified copy of the data. There is an
optional parameter ``inplace`` so that the original data can be modified
without creating a copy:

.. ipython:: python

   df_orig = df.copy()
   df_orig.where(df > 0, -df, inplace=True);
   df_orig

**alignment**

Furthermore, ``where`` aligns the input boolean condition (ndarray or DataFrame),
such that partial selection with setting is possible. This is analagous to
partial setting via ``.ix`` (but on the contents rather than the axis labels)

.. ipython:: python

   df2 = df.copy()
   df2[ df2[1:4] > 0 ] = 3
   df2

.. versionadded:: 0.13

Where can also accept ``axis`` and ``level`` parameters to align the input when
performing the ``where``.

.. ipython:: python

   df2 = df.copy()
   df2.where(df2>0,df2['A'],axis='index')

This is equivalent (but faster than) the following.

.. ipython:: python

   df2 = df.copy()
   df.apply(lambda x, y: x.where(x>0,y), y=df['A'])

**mask**

``mask`` is the inverse boolean operation of ``where``.

.. ipython:: python

   s.mask(s >= 0)
   df.mask(df >= 0)

.. _indexing.query:

The :meth:`~pandas.DataFrame.query` Method (Experimental)
---------------------------------------------------------

.. versionadded:: 0.13

:class:`~pandas.DataFrame` objects have a :meth:`~pandas.DataFrame.query`
method that allows selection using an expression.

You can get the value of the frame where column ``b`` has values
between the values of columns ``a`` and ``c``. For example:

.. ipython:: python
   :suppress:

   from numpy.random import randint, rand
   np.random.seed(1234)

.. ipython:: python

   n = 10
   df = DataFrame(rand(n, 3), columns=list('abc'))
   df

   # pure python
   df[(df.a < df.b) & (df.b < df.c)]

   # query
   df.query('(a < b) & (b < c)')

Do the same thing but fallback on a named index if there is no column
with the name ``a``.

.. ipython:: python

   df = DataFrame(randint(n / 2, size=(n, 2)), columns=list('bc'))
   df.index.name = 'a'
   df
   df.query('a < b and b < c')

If instead you don't want to or cannot name your index, you can use the name
``index`` in your query expression:

.. ipython:: python
   :suppress:

   old_index = index
   del index

.. ipython:: python

   df = DataFrame(randint(n, size=(n, 2)), columns=list('bc'))
   df
   df.query('index < b < c')

.. ipython:: python
   :suppress:

   index = old_index
   del old_index


.. note::

   If the name of your index overlaps with a column name, the column name is
   given precedence. For example,

   .. ipython:: python

      df = DataFrame({'a': randint(5, size=5)})
      df.index.name = 'a'
      df.query('a > 2') # uses the column 'a', not the index

   You can still use the index in a query expression by using the special
   identifier 'index':

   .. ipython:: python

      df.query('index > 2')

   If for some reason you have a column named ``index``, then you can refer to
   the index as ``ilevel_0`` as well, but at this point you should consider
   renaming your columns to something less ambiguous.


:class:`~pandas.MultiIndex` :meth:`~pandas.DataFrame.query` Syntax
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also use the levels of a ``DataFrame`` with a
:class:`~pandas.MultiIndex` as if they were columns in the frame:

.. ipython:: python

   import pandas.util.testing as tm

   n = 10
   colors = tm.choice(['red', 'green'], size=n)
   foods = tm.choice(['eggs', 'ham'], size=n)
   colors
   foods

   index = MultiIndex.from_arrays([colors, foods], names=['color', 'food'])
   df = DataFrame(randn(n, 2), index=index)
   df
   df.query('color == "red"')

If the levels of the ``MultiIndex`` are unnamed, you can refer to them using
special names:


.. ipython:: python

   df.index.names = [None, None]
   df
   df.query('ilevel_0 == "red"')


The convention is ``ilevel_0``, which means "index level 0" for the 0th level
of the ``index``.


:meth:`~pandas.DataFrame.query` Use Cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A use case for :meth:`~pandas.DataFrame.query` is when you have a collection of
:class:`~pandas.DataFrame` objects that have a subset of column names (or index
levels/names) in common. You can pass the same query to both frames *without*
having to specify which frame you're interested in querying

.. ipython:: python

   df = DataFrame(rand(n, 3), columns=list('abc'))
   df
   df2 = DataFrame(rand(n + 2, 3), columns=df.columns)
   df2
   expr = '0.0 <= a <= c <= 0.5'
   map(lambda frame: frame.query(expr), [df, df2])

:meth:`~pandas.DataFrame.query` Python versus pandas Syntax Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full numpy-like syntax

.. ipython:: python

   df = DataFrame(randint(n, size=(n, 3)), columns=list('abc'))
   df
   df.query('(a < b) & (b < c)')
   df[(df.a < df.b) & (df.b < df.c)]

Slightly nicer by removing the parentheses (by binding making comparison
operators bind tighter than ``&``/``|``)

.. ipython:: python

   df.query('a < b & b < c')

Use English instead of symbols

.. ipython:: python

   df.query('a < b and b < c')

Pretty close to how you might write it on paper

.. ipython:: python

   df.query('a < b < c')

The ``in`` and ``not in`` operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~pandas.DataFrame.query` also supports special use of Python's ``in`` and
``not in`` comparison operators, providing a succint syntax for calling the
``isin`` method of a ``Series`` or ``DataFrame``.

.. ipython:: python
   :suppress:

   try:
       old_d = d
       del d
   except NameError:
       pass

.. ipython:: python

   # get all rows where columns "a" and "b" have overlapping values
   df = DataFrame({'a': list('aabbccddeeff'), 'b': list('aaaabbbbcccc'),
                   'c': randint(5, size=12), 'd': randint(9, size=12)})
   df
   df.query('a in b')

   # How you'd do it in pure Python
   df[df.a.isin(df.b)]

   df.query('a not in b')

   # pure Python
   df[~df.a.isin(df.b)]


You can combine this with other expressions for very succinct queries:


.. ipython:: python

   # rows where cols a and b have overlapping values and col c's values are less than col d's
   df.query('a in b and c < d')

   # pure Python
   df[df.b.isin(df.a) & (df.c < df.d)]


.. note::

   Note that ``in`` and ``not in`` are evaluated in Python, since ``numexpr``
   has no equivalent of this operation. However, **only the** ``in``/``not in``
   **expression itself** is evaluated in vanilla Python. For example, in the
   expression

   .. code-block:: python

      df.query('a in b + c + d')

   ``(b + c + d)`` is evaluated by ``numexpr`` and *then* the ``in``
   operation is evaluated in plain Python. In general, any operations that can
   be evaluated using ``numexpr`` will be.

Special use of the ``==`` operator with ``list`` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Comparing a ``list`` of values to a column using ``==``/``!=`` works similarly
to ``in``/``not in``

.. ipython:: python

   df.query('b == ["a", "b", "c"]')

   # pure Python
   df[df.b.isin(["a", "b", "c"])]

   df.query('c == [1, 2]')

   df.query('c != [1, 2]')

   # using in/not in
   df.query('[1, 2] in c')

   df.query('[1, 2] not in c')

   # pure Python
   df[df.c.isin([1, 2])]


Boolean Operators
~~~~~~~~~~~~~~~~~

You can negate boolean expressions with the word ``not`` or the ``~`` operator.

.. ipython:: python

   df = DataFrame(rand(n, 3), columns=list('abc'))
   df['bools'] = rand(len(df)) > 0.5
   df.query('~bools')
   df.query('not bools')
   df.query('not bools') == df[~df.bools]

Of course, expressions can be arbitrarily complex too

.. ipython:: python

   # short query syntax
   shorter = df.query('a < b < c and (not bools) or bools > 2')

   # equivalent in pure Python
   longer = df[(df.a < df.b) & (df.b < df.c) & (~df.bools) | (df.bools > 2)]

   shorter
   longer

   shorter == longer

.. ipython:: python
   :suppress:

   try:
       d = old_d
       del old_d
   except NameError:
       pass


Performance of :meth:`~pandas.DataFrame.query`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DataFrame.query()`` using ``numexpr`` is slightly faster than Python for
large frames

.. image:: _static/query-perf.png

.. note::

   You will only see the performance benefits of using the ``numexpr`` engine
   with ``DataFrame.query()`` if your frame has more than approximately 200,000
   rows

      .. image:: _static/query-perf-small.png

This plot was created using a ``DataFrame`` with 3 columns each containing
floating point values generated using ``numpy.random.randn()``.

.. ipython:: python
   :suppress:

   df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
   df2 = df.copy()

Take Methods
------------

.. _indexing.take:

Similar to numpy ndarrays, pandas Index, Series, and DataFrame also provides
the ``take`` method that retrieves elements along a given axis at the given
indices. The given indices must be either a list or an ndarray of integer
index positions. ``take`` will also accept negative integers as relative positions to the end of the object.

.. ipython:: python

   index = Index(randint(0, 1000, 10))
   index

   positions = [0, 9, 3]

   index[positions]
   index.take(positions)

   ser = Series(randn(10))

   ser.ix[positions]
   ser.take(positions)

For DataFrames, the given indices should be a 1d list or ndarray that specifies
row or column positions.

.. ipython:: python

   frm = DataFrame(randn(5, 3))

   frm.take([1, 4, 3])

   frm.take([0, 2], axis=1)

It is important to note that the ``take`` method on pandas objects are not
intended to work on boolean indices and may return unexpected results.

.. ipython:: python

   arr = randn(10)
   arr.take([False, False, True, True])
   arr[[0, 1]]

   ser = Series(randn(10))
   ser.take([False, False, True, True])
   ser.ix[[0, 1]]

Finally, as a small note on performance, because the ``take`` method handles
a narrower range of inputs, it can offer performance that is a good deal
faster than fancy indexing.

.. ipython::

   arr = randn(10000, 5)
   indexer = np.arange(10000)
   random.shuffle(indexer)

   timeit arr[indexer]
   timeit arr.take(indexer, axis=0)

   ser = Series(arr[:, 0])
   timeit ser.ix[indexer]
   timeit ser.take(indexer)

Duplicate Data
--------------

.. _indexing.duplicate:

If you want to identify and remove duplicate rows in a DataFrame,  there are
two methods that will help: ``duplicated`` and ``drop_duplicates``. Each
takes as an argument the columns to use to identify duplicated rows.

- ``duplicated`` returns a boolean vector whose length is the number of rows, and which indicates whether a row is duplicated.
- ``drop_duplicates`` removes duplicate rows.

By default, the first observed row of a duplicate set is considered unique, but
each method has a ``take_last`` parameter that indicates the last observed row
should be taken instead.

.. ipython:: python

   df2 = DataFrame({'a' : ['one', 'one', 'two', 'three', 'two', 'one', 'six'],
                    'b' : ['x', 'y', 'y', 'x', 'y', 'x', 'x'],
                    'c' : np.random.randn(7)})
   df2.duplicated(['a','b'])
   df2.drop_duplicates(['a','b'])
   df2.drop_duplicates(['a','b'], take_last=True)

.. _indexing.dictionarylike:

Dictionary-like :meth:`~pandas.DataFrame.get` method
----------------------------------------------------

Each of Series, DataFrame, and Panel have a ``get`` method which can return a
default value.

.. ipython:: python

   s = Series([1,2,3], index=['a','b','c'])
   s.get('a')               # equivalent to s['a']
   s.get('x', default=-1)

.. _indexing.advanced:

Advanced Indexing with ``.ix``
------------------------------

.. note::

   The recent addition of ``.loc`` and ``.iloc`` have enabled users to be quite
   explicit about indexing choices. ``.ix`` allows a great flexibility to
   specify indexing locations by *label* and/or *integer position*. pandas will
   attempt to use any passed *integer* as *label* locations first (like what
   ``.loc`` would do, then to fall back on *positional* indexing, like what
   ``.iloc``  would do). See :ref:`Fallback Indexing <indexing.fallback>` for
   an example.

The syntax of using ``.ix`` is identical to ``.loc``, in :ref:`Selection by
Label <indexing.label>`, and ``.iloc`` in :ref:`Selection by Position <indexing.integer>`.

The ``.ix`` attribute takes the following inputs:

- An integer or single label, e.g. ``5`` or ``'a'``
- A list or array of labels ``['a', 'b', 'c']`` or integers ``[4, 3, 0]``
- A slice object with ints ``1:7`` or labels ``'a':'f'``
- A boolean array

We'll illustrate all of these methods. First, note that this provides a concise
way of reindexing on multiple axes at once:

.. ipython:: python

   subindex = dates[[3,4,5]]
   df.reindex(index=subindex, columns=['C', 'B'])
   df.ix[subindex, ['C', 'B']]

Assignment / setting values is possible when using ``ix``:

.. ipython:: python

   df2 = df.copy()
   df2.ix[subindex, ['C', 'B']] = 0
   df2

Indexing with an array of integers can also be done:

.. ipython:: python

   df.ix[[4,3,1]]
   df.ix[dates[[4,3,1]]]

**Slicing** has standard Python semantics for integer slices:

.. ipython:: python

   df.ix[1:7, :2]

Slicing with labels is semantically slightly different because the slice start
and stop are **inclusive** in the label-based case:

.. ipython:: python

   date1, date2 = dates[[2, 4]]
   print(date1, date2)
   df.ix[date1:date2]
   df['A'].ix[date1:date2]

Getting and setting rows in a DataFrame, especially by their location, is much
easier:

.. ipython:: python

   df2 = df[:5].copy()
   df2.ix[3]
   df2.ix[3] = np.arange(len(df2.columns))
   df2

Column or row selection can be combined as you would expect with arrays of
labels or even boolean vectors:

.. ipython:: python

   df.ix[df['A'] > 0, 'B']
   df.ix[date1:date2, 'B']
   df.ix[date1, 'B']

Slicing with labels is closely related to the ``truncate`` method which does
precisely ``.ix[start:stop]`` but returns a copy (for legacy reasons).

The :meth:`~pandas.DataFrame.select` Method
-------------------------------------------

Another way to extract slices from an object is with the ``select`` method of
Series, DataFrame, and Panel. This method should be used only when there is no
more direct way.  ``select`` takes a function which operates on labels along
``axis`` and returns a boolean.  For instance:

.. ipython:: python

   df.select(lambda x: x == 'A', axis=1)

The :meth:`~pandas.DataFrame.lookup` Method
-------------------------------------------

Sometimes you want to extract a set of values given a sequence of row labels
and column labels, and the ``lookup`` method allows for this and returns a
numpy array.  For instance,

.. ipython:: python

  dflookup = DataFrame(np.random.rand(20,4), columns = ['A','B','C','D'])
  dflookup.lookup(list(range(0,10,2)), ['B','C','A','B','D'])

.. _indexing.float64index:

Float64Index
------------

.. note::

   As of 0.14.0, ``Float64Index`` is backed by a native ``float64`` dtype
   array. Prior to 0.14.0, ``Float64Index`` was backed by an ``object`` dtype
   array. Using a ``float64`` dtype in the backend speeds up arithmetic
   operations by about 30x and boolean indexing operations on the
   ``Float64Index`` itself are about 2x as fast.


.. versionadded:: 0.13.0

By default a ``Float64Index`` will be automatically created when passing floating, or mixed-integer-floating values in index creation.
This enables a pure label-based slicing paradigm that makes ``[],ix,loc`` for scalar indexing and slicing work exactly the
same.

.. ipython:: python

   indexf = Index([1.5, 2, 3, 4.5, 5])
   indexf
   sf = Series(range(5),index=indexf)
   sf

Scalar selection for ``[],.ix,.loc`` will always be label based. An integer will match an equal float index (e.g. ``3`` is equivalent to ``3.0``)

.. ipython:: python

   sf[3]
   sf[3.0]
   sf.ix[3]
   sf.ix[3.0]
   sf.loc[3]
   sf.loc[3.0]

The only positional indexing is via ``iloc``

.. ipython:: python

   sf.iloc[3]

A scalar index that is not found will raise ``KeyError``

Slicing is ALWAYS on the values of the index, for ``[],ix,loc`` and ALWAYS positional with ``iloc``

.. ipython:: python

   sf[2:4]
   sf.ix[2:4]
   sf.loc[2:4]
   sf.iloc[2:4]

In float indexes, slicing using floats is allowed

.. ipython:: python

   sf[2.1:4.6]
   sf.loc[2.1:4.6]

In non-float indexes, slicing using floats will raise a ``TypeError``

.. code-block:: python

   In [1]: Series(range(5))[3.5]
   TypeError: the label [3.5] is not a proper indexer for this index type (Int64Index)

   In [1]: Series(range(5))[3.5:4.5]
   TypeError: the slice start [3.5] is not a proper indexer for this index type (Int64Index)

Using a scalar float indexer will be deprecated in a future version, but is allowed for now.

.. code-block:: python

   In [3]: Series(range(5))[3.0]
   Out[3]: 3

Here is a typical use-case for using this type of indexing. Imagine that you have a somewhat
irregular timedelta-like indexing scheme, but the data is recorded as floats. This could for
example be millisecond offsets.

.. ipython:: python

   dfir = concat([DataFrame(randn(5,2),
                     index=np.arange(5) * 250.0,
                     columns=list('AB')),
                  DataFrame(randn(6,2),
                     index=np.arange(4,10) * 250.1,
                     columns=list('AB'))])
   dfir

Selection operations then will always work on a value basis, for all selection operators.

.. ipython:: python

   dfir[0:1000.4]
   dfir.loc[0:1001,'A']
   dfir.loc[1000.4]

You could then easily pick out the first 1 second (1000 ms) of data then.

.. ipython:: python

   dfir[0:1000]

Of course if you need integer based selection, then use ``iloc``

.. ipython:: python

   dfir.iloc[0:5]

.. _indexing.view_versus_copy:

Returning a view versus a copy
------------------------------

When setting values in a pandas object, care must be taken to avoid what is called
``chained indexing``. Here is an example.

.. ipython:: python

   dfmi = DataFrame([list('abcd'),
                     list('efgh'),
                     list('ijkl'),
                     list('mnop')],
                    columns=MultiIndex.from_product([['one','two'],
                                                     ['first','second']]))
   dfmi

Compare these two access methods:

.. ipython:: python

   dfmi['one']['second']

.. ipython:: python

   dfmi.loc[:,('one','second')]

These both yield the same results, so which should you use? It is instructive to understand the order
of operations on these and why method 2 (``.loc``) is much preferred over method 1 (chained ``[]``)

``dfmi['one']`` selects the first level of the columns and returns a data frame that is singly-indexed.
Then another python operation ``dfmi_with_one['second']`` selects the series indexed by ``'second'`` happens.
This is indicated by the variable ``dfmi_with_one`` because pandas sees these operations as separate events.
e.g. separate calls to ``__getitem__``, so it has to treat them as linear operations, they happen one after another.

Contrast this to ``df.loc[:,('one','second')]`` which passes a nested tuple of ``(slice(None),('one','second'))`` to a single call to
``__getitem__``. This allows pandas to deal with this as a single entity. Furthermore this order of operations *can* be significantly
faster, and allows one to index *both* axes if so desired.

Why does the assignment when using chained indexing fail!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

So, why does this show the ``SettingWithCopy`` warning / and possibly not work when you do chained indexing and assignement:

.. code-block:: python

   dfmi['one']['second'] = value

Since the chained indexing is 2 calls, it is possible that either call may return a **copy** of the data because of the way it is sliced.
Thus when setting, you are actually setting a **copy**, and not the original frame data. It is impossible for pandas to figure this out because their are 2 separate python operations that are not connected.

The ``SettingWithCopy`` warning is a 'heuristic' to detect this (meaning it tends to catch most cases but is simply a lightweight check). Figuring this out for real is way complicated.

The ``.loc`` operation is a single python operation, and thus can select a slice (which still may be a copy), but allows pandas to assign that slice back into the frame after it is modified, thus setting the values as you would think.

The reason for having the ``SettingWithCopy`` warning is this. Sometimes when you slice an array you will simply get a view back, which means you can set it no problem. However, even a single dtyped array can generate a copy if it is sliced in a particular way. A multi-dtyped DataFrame (meaning it has say ``float`` and ``object`` data), will almost always yield a copy. Whether a view is created is dependent on the memory layout of the array.

Evaluation order matters
~~~~~~~~~~~~~~~~~~~~~~~~

Furthermore, in chained expressions, the order may determine whether a copy is returned or not.
If an expression will set values on a copy of a slice, then a ``SettingWithCopy``
exception will be raised (this raise/warn behavior is new starting in 0.13.0)

You can control the action of a chained assignment via the option ``mode.chained_assignment``,
which can take the values ``['raise','warn',None]``, where showing a warning is the default.

.. ipython:: python

   dfb = DataFrame({'a' : ['one', 'one', 'two',
                           'three', 'two', 'one', 'six'],
                    'c' : np.arange(7)})

   # passed via reference (will stay)
   dfb['c'][dfb.a.str.startswith('o')] = 42

This however is operating on a copy and will not work.

::

   >>> pd.set_option('mode.chained_assignment','warn')
   >>> dfb[dfb.a.str.startswith('o')]['c'] = 42
   Traceback (most recent call last)
        ...
   SettingWithCopyWarning:
        A value is trying to be set on a copy of a slice from a DataFrame.
        Try using .loc[row_index,col_indexer] = value instead

A chained assignment can also crop up in setting in a mixed dtype frame.

.. note::

   These setting rules apply to all of ``.loc/.iloc/.ix``

This is the correct access method

.. ipython:: python

   dfc = DataFrame({'A':['aaa','bbb','ccc'],'B':[1,2,3]})
   dfc.loc[0,'A'] = 11
   dfc

This *can* work at times, but is not guaranteed, and so should be avoided

.. ipython:: python

   dfc = dfc.copy()
   dfc['A'][0] = 111
   dfc

This will **not** work at all, and so should be avoided

::

   >>> pd.set_option('mode.chained_assignment','raise')
   >>> dfc.loc[0]['A'] = 1111
   Traceback (most recent call last)
        ...
   SettingWithCopyException:
        A value is trying to be set on a copy of a slice from a DataFrame.
        Try using .loc[row_index,col_indexer] = value instead

.. warning::

   The chained assignment warnings / exceptions are aiming to inform the user of a possibly invalid
   assignment. There may be false positives; situations where a chained assignment is inadvertantly
   reported.


Fallback indexing
-----------------

.. _indexing.fallback:

Float indexes should be used only with caution. If you have a float indexed
``DataFrame`` and try to select using an integer, the row that pandas returns
might not be what you expect. pandas first attempts to use the *integer*
as a *label* location, but fails to find a match (because the types
are not equal). pandas then falls back to back to positional indexing.

.. ipython:: python

    df = pd.DataFrame(np.random.randn(4,4),
        columns=list('ABCD'), index=[1.0, 2.0, 3.0, 4.0])
    df
    df.ix[1]

To select the row you do expect, instead use a float label or
use ``iloc``.

.. ipython:: python

    df.ix[1.0]
    df.iloc[0]

Instead of using a float index, it is often better to
convert to an integer index:

.. ipython:: python

    df_new = df.reset_index()
    df_new[df_new['index'] == 1.0]
    # now you can also do "float selection"
    df_new[(df_new['index'] >= 1.0) & (df_new['index'] < 2)]


.. _indexing.class:

Index objects
-------------

The pandas :class:`~pandas.Index` class and its subclasses can be viewed as
implementing an *ordered multiset*. Duplicates are allowed. However, if you try
to convert an :class:`~pandas.Index` object with duplicate entries into a
``set``, an exception will be raised.

:class:`~pandas.Index` also provides the infrastructure necessary for
lookups, data alignment, and reindexing. The easiest way to create an
:class:`~pandas.Index` directly is to pass a ``list`` or other sequence to
:class:`~pandas.Index`:

.. ipython:: python

   index = Index(['e', 'd', 'a', 'b'])
   index
   'd' in index

You can also pass a ``name`` to be stored in the index:


.. ipython:: python

   index = Index(['e', 'd', 'a', 'b'], name='something')
   index.name

Starting with pandas 0.5, the name, if set, will be shown in the console
display:

.. ipython:: python

   index = Index(list(range(5)), name='rows')
   columns = Index(['A', 'B', 'C'], name='cols')
   df = DataFrame(np.random.randn(5, 3), index=index, columns=columns)
   df
   df['A']


Set operations on Index objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.set_ops:

The three main operations are ``union (|)``, ``intersection (&)``, and ``diff
(-)``. These can be directly called as instance methods or used via overloaded
operators:

.. ipython:: python

   a = Index(['c', 'b', 'a'])
   b = Index(['c', 'e', 'd'])
   a.union(b)
   a | b
   a & b
   a - b

Also available is the ``sym_diff (^)`` operation, which returns elements
that appear in either ``idx1`` or ``idx2`` but not both. This is
equivalent to the Index created by ``(idx1 - idx2) + (idx2 - idx1)``,
with duplicates dropped.

.. ipython:: python

   idx1 = Index([1, 2, 3, 4])
   idx2 = Index([2, 3, 4, 5])
   idx1.sym_diff(idx2)
   idx1 ^ idx2

The ``isin`` method of Index objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One additional operation is the ``isin`` method that works analogously to the
``Series.isin`` method found :ref:`here <indexing.boolean>`.

.. _indexing.hierarchical:

Hierarchical indexing (MultiIndex)
----------------------------------

Hierarchical indexing (also referred to as "multi-level" indexing) is brand new
in the pandas 0.4 release. It is very exciting as it opens the door to some
quite sophisticated data analysis and manipulation, especially for working with
higher dimensional data. In essence, it enables you to store and manipulate
data with an arbitrary number of dimensions in lower dimensional data
structures like Series (1d) and DataFrame (2d).

In this section, we will show what exactly we mean by "hierarchical" indexing
and how it integrates with the all of the pandas indexing functionality
described above and in prior sections. Later, when discussing :ref:`group by
<groupby>` and :ref:`pivoting and reshaping data <reshaping>`, we'll show
non-trivial applications to illustrate how it aids in structuring data for
analysis.

See the :ref:`cookbook<cookbook.multi_index>` for some advanced strategies

.. note::

   Given that hierarchical indexing is so new to the library, it is definitely
   "bleeding-edge" functionality but is certainly suitable for production. But,
   there may inevitably be some minor API changes as more use cases are
   explored and any weaknesses in the design / implementation are identified.
   pandas aims to be "eminently usable" so any feedback about new
   functionality like this is extremely helpful.

Creating a MultiIndex (hierarchical index) object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``MultiIndex`` object is the hierarchical analogue of the standard
``Index`` object which typically stores the axis labels in pandas objects. You
can think of ``MultiIndex`` an array of tuples where each tuple is unique. A
``MultiIndex`` can be created from a list of arrays (using
``MultiIndex.from_arrays``), an array of tuples (using
``MultiIndex.from_tuples``), or a crossed set of iterables (using
``MultiIndex.from_product``).  The ``Index`` constructor will attempt to return
a ``MultiIndex`` when it is passed a list of tuples.  The following examples
demo different ways to initialize MultiIndexes.


.. ipython:: python

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   tuples = list(zip(*arrays))
   tuples

   index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
   index

   s = Series(randn(8), index=index)
   s

When you want every pairing of the elements in two iterables, it can be easier
to use the ``MultiIndex.from_product`` function:

.. ipython:: python

   iterables = [['bar', 'baz', 'foo', 'qux'], ['one', 'two']]
   MultiIndex.from_product(iterables, names=['first', 'second'])

As a convenience, you can pass a list of arrays directly into Series or
DataFrame to construct a MultiIndex automatically:

.. ipython:: python

   arrays = [np.array(['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'])
   ,
             np.array(['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two'])
             ]
   s = Series(randn(8), index=arrays)
   s
   df = DataFrame(randn(8, 4), index=arrays)
   df

All of the ``MultiIndex`` constructors accept a ``names`` argument which stores
string names for the levels themselves. If no names are provided, ``None`` will
be assigned:

.. ipython:: python

   df.index.names

This index can back any axis of a pandas object, and the number of **levels**
of the index is up to you:

.. ipython:: python

   df = DataFrame(randn(3, 8), index=['A', 'B', 'C'], columns=index)
   df
   DataFrame(randn(6, 6), index=index[:6], columns=index[:6])

We've "sparsified" the higher levels of the indexes to make the console output a
bit easier on the eyes.

It's worth keeping in mind that there's nothing preventing you from using
tuples as atomic labels on an axis:

.. ipython:: python

   Series(randn(8), index=tuples)

The reason that the ``MultiIndex`` matters is that it can allow you to do
grouping, selection, and reshaping operations as we will describe below and in
subsequent areas of the documentation. As you will see in later sections, you
can find yourself working with hierarchically-indexed data without creating a
``MultiIndex`` explicitly yourself. However, when loading data from a file, you
may wish to generate your own ``MultiIndex`` when preparing the data set.

Note that how the index is displayed by be controlled using the
``multi_sparse`` option in ``pandas.set_printoptions``:

.. ipython:: python

   pd.set_option('display.multi_sparse', False)
   df
   pd.set_option('display.multi_sparse', True)

Reconstructing the level labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.get_level_values:

The method ``get_level_values`` will return a vector of the labels for each
location at a particular level:

.. ipython:: python

   index.get_level_values(0)
   index.get_level_values('second')


Basic indexing on axis with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the important features of hierarchical indexing is that you can select
data by a "partial" label identifying a subgroup in the data. **Partial**
selection "drops" levels of the hierarchical index in the result in a
completely analogous way to selecting a column in a regular DataFrame:

.. ipython:: python

   df['bar']
   df['bar', 'one']
   df['bar']['one']
   s['qux']

See :ref:`Cross-section with hierarchical index <indexing.xs>` for how to select
on a deeper level.


Data alignment and using ``reindex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Operations between differently-indexed objects having ``MultiIndex`` on the
axes will work as you expect; data alignment will work the same as an Index of
tuples:

.. ipython:: python

   s + s[:-2]
   s + s[::2]

``reindex`` can be called with another ``MultiIndex`` or even a list or array
of tuples:

.. ipython:: python

   s.reindex(index[:3])
   s.reindex([('foo', 'two'), ('bar', 'one'), ('qux', 'one'), ('baz', 'one')])

.. _indexing.advanced_hierarchical:

Advanced indexing with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Syntactically integrating ``MultiIndex`` in advanced indexing with ``.loc/.ix`` is a
bit challenging, but we've made every effort to do so. for example the
following works as you would expect:

.. ipython:: python

   df = df.T
   df
   df.loc['bar']
   df.loc['bar', 'two']

"Partial" slicing also works quite nicely.

.. ipython:: python

   df.loc['baz':'foo']

You can slice with a 'range' of values, by providing a slice of tuples.

.. ipython:: python

   df.loc[('baz', 'two'):('qux', 'one')]
   df.loc[('baz', 'two'):'foo']

Passing a list of labels or tuples works similar to reindexing:

.. ipython:: python

   df.ix[[('bar', 'two'), ('qux', 'one')]]

.. _indexing.mi_slicers:

Multiindexing using slicers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.14.0

In 0.14.0 we added a new way to slice multi-indexed objects.
You can slice a multi-index by providing multiple indexers.

You can provide any of the selectors as if you are indexing by label, see :ref:`Selection by Label <indexing.label>`,
including slices, lists of labels, labels, and boolean indexers.

You can use ``slice(None)`` to select all the contents of *that* level. You do not need to specify all the
*deeper* levels, they will be implied as ``slice(None)``.

As usual, **both sides** of the slicers are included as this is label indexing.

.. warning::

   You should specify all axes in the ``.loc`` specifier, meaning the indexer for the **index** and
   for the **columns**. Their are some ambiguous cases where the passed indexer could be mis-interpreted
   as indexing *both* axes, rather than into say the MuliIndex for the rows.

   You should do this:

   .. code-block:: python

      df.loc[(slice('A1','A3'),.....),:]

   rather than this:

   .. code-block:: python

      df.loc[(slice('A1','A3'),.....)]

.. warning::

   You will need to make sure that the selection axes are fully lexsorted!

.. ipython:: python

   def mklbl(prefix,n):
       return ["%s%s" % (prefix,i)  for i in range(n)]

   miindex = MultiIndex.from_product([mklbl('A',4),
                                      mklbl('B',2),
                                      mklbl('C',4),
                                      mklbl('D',2)])
   micolumns = MultiIndex.from_tuples([('a','foo'),('a','bar'),
                                       ('b','foo'),('b','bah')],
                                        names=['lvl0', 'lvl1'])
   dfmi = DataFrame(np.arange(len(miindex)*len(micolumns)).reshape((len(miindex),len(micolumns))),
                    index=miindex,
                    columns=micolumns).sortlevel().sortlevel(axis=1)
   dfmi

Basic multi-index slicing using slices, lists, and labels.

.. ipython:: python

   dfmi.loc[(slice('A1','A3'),slice(None), ['C1','C3']),:]

You can use a ``pd.IndexSlice`` to shortcut the creation of these slices

.. ipython:: python

   idx = pd.IndexSlice
   dfmi.loc[idx[:,:,['C1','C3']],idx[:,'foo']]

It is possible to perform quite complicated selections using this method on multiple
axes at the same time.

.. ipython:: python

   dfmi.loc['A1',(slice(None),'foo')]
   dfmi.loc[idx[:,:,['C1','C3']],idx[:,'foo']]

Using a boolean indexer you can provide selection related to the *values*.

.. ipython:: python

   mask = dfmi[('a','foo')]>200
   dfmi.loc[idx[mask,:,['C1','C3']],idx[:,'foo']]

You can also specify the ``axis`` argument to ``.loc`` to interpret the passed
slicers on a single axis.

.. ipython:: python

   dfmi.loc(axis=0)[:,:,['C1','C3']]

Furthermore you can *set* the values using these methods

.. ipython:: python

   df2 = dfmi.copy()
   df2.loc(axis=0)[:,:,['C1','C3']] = -10
   df2

You can use a right-hand-side of an alignable object as well.

.. ipython:: python

   df2 = dfmi.copy()
   df2.loc[idx[:,:,['C1','C3']],:] = df2*1000
   df2

.. _indexing.xs:

Cross-section with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``xs`` method of ``DataFrame`` additionally takes a level argument to make
selecting data at a particular level of a MultiIndex easier.

.. ipython:: python

    df.xs('one', level='second')

.. ipython:: python

   # using the slicers (new in 0.14.0)
   df.loc[(slice(None),'one'),:]

You can also select on the columns with :meth:`~pandas.MultiIndex.xs`, by
providing the axis argument

.. ipython:: python

   df = df.T
   df.xs('one', level='second', axis=1)

.. ipython:: python

   # using the slicers (new in 0.14.0)
   df.loc[:,(slice(None),'one')]

:meth:`~pandas.MultiIndex.xs` also allows selection with multiple keys

.. ipython:: python

   df.xs(('one', 'bar'), level=('second', 'first'), axis=1)

.. ipython:: python

   # using the slicers (new in 0.14.0)
   df.loc[:,('bar','one')]

.. versionadded:: 0.13.0

You can pass ``drop_level=False`` to :meth:`~pandas.MultiIndex.xs` to retain
the level that was selected

.. ipython:: python

   df.xs('one', level='second', axis=1, drop_level=False)

versus the result with ``drop_level=True`` (the default value)

.. ipython:: python

   df.xs('one', level='second', axis=1, drop_level=True)

.. ipython:: python
   :suppress:

   df = df.T

.. _indexing.advanced_reindex:

Advanced reindexing and alignment with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameter ``level`` has been added to the ``reindex`` and ``align`` methods
of pandas objects. This is useful to broadcast values across a level. For
instance:

.. ipython:: python

   midx = MultiIndex(levels=[['zero', 'one'], ['x','y']],
                     labels=[[1,1,0,0],[1,0,1,0]])
   df = DataFrame(randn(4,2), index=midx)
   print(df)
   df2 = df.mean(level=0)
   print(df2)
   print(df2.reindex(df.index, level=0))
   df_aligned, df2_aligned = df.align(df2, level=0)
   print(df_aligned)
   print(df2_aligned)


The need for sortedness with :class:`~pandas.MultiIndex`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Caveat emptor**: the present implementation of ``MultiIndex`` requires that
the labels be sorted for some of the slicing / indexing routines to work
correctly. You can think about breaking the axis into unique groups, where at
the hierarchical level of interest, each distinct group shares a label, but no
two have the same label. However, the ``MultiIndex`` does not enforce this:
**you are responsible for ensuring that things are properly sorted**. There is
an important new method ``sortlevel`` to sort an axis within a ``MultiIndex``
so that its labels are grouped and sorted by the original ordering of the
associated factor at that level. Note that this does not necessarily mean the
labels will be sorted lexicographically!

.. ipython:: python

   import random; random.shuffle(tuples)
   s = Series(randn(8), index=MultiIndex.from_tuples(tuples))
   s
   s.sortlevel(0)
   s.sortlevel(1)

.. _indexing.sortlevel_byname:

Note, you may also pass a level name to ``sortlevel`` if the MultiIndex levels
are named.

.. ipython:: python

   s.index.set_names(['L1', 'L2'], inplace=True)
   s.sortlevel(level='L1')
   s.sortlevel(level='L2')

Some indexing will work even if the data are not sorted, but will be rather
inefficient and will also return a copy of the data rather than a view:

.. ipython:: python

   s['qux']
   s.sortlevel(1)['qux']

On higher dimensional objects, you can sort any of the other axes by level if
they have a MultiIndex:

.. ipython:: python

   df.T.sortlevel(1, axis=1)

The ``MultiIndex`` object has code to **explicity check the sort depth**. Thus,
if you try to index at a depth at which the index is not sorted, it will raise
an exception. Here is a concrete example to illustrate this:

.. ipython:: python

   tuples = [('a', 'a'), ('a', 'b'), ('b', 'a'), ('b', 'b')]
   idx = MultiIndex.from_tuples(tuples)
   idx.lexsort_depth

   reordered = idx[[1, 0, 3, 2]]
   reordered.lexsort_depth

   s = Series(randn(4), index=reordered)
   s.ix['a':'a']

However:

::

   >>> s.ix[('a', 'b'):('b', 'a')]
   Traceback (most recent call last)
        ...
   KeyError: Key length (3) was greater than MultiIndex lexsort depth (2)

Swapping levels with :meth:`~pandas.MultiIndex.swaplevel`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``swaplevel`` function can switch the order of two levels:

.. ipython:: python

   df[:5]
   df[:5].swaplevel(0, 1, axis=0)

.. _indexing.reorderlevels:

Reordering levels with :meth:`~pandas.MultiIndex.reorder_levels`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``reorder_levels`` function generalizes the ``swaplevel`` function,
allowing you to permute the hierarchical index levels in one step:

.. ipython:: python

   df[:5].reorder_levels([1,0], axis=0)


Some gory internal details
~~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, the ``MultiIndex`` consists of a few things: the **levels**, the
integer **labels**, and the level **names**:

.. ipython:: python

   index
   index.levels
   index.labels
   index.names

You can probably guess that the labels determine which unique element is
identified with that location at each layer of the index. It's important to
note that sortedness is determined **solely** from the integer labels and does
not check (or care) whether the levels themselves are sorted. Fortunately, the
constructors ``from_tuples`` and ``from_arrays`` ensure that this is true, but
if you compute the levels and labels yourself, please be careful.


Setting index metadata (``name(s)``, ``levels``, ``labels``)
------------------------------------------------------------

.. versionadded:: 0.13.0

.. _indexing.set_metadata:

Indexes are "mostly immutable", but it is possible to set and change their
metadata, like the index ``name`` (or, for ``MultiIndex``, ``levels`` and
``labels``).

You can use the ``rename``, ``set_names``, ``set_levels``, and ``set_labels``
to set these attributes directly. They default to returning a copy; however,
you can specify ``inplace=True`` to have the data change inplace.

.. ipython:: python

  ind = Index([1, 2, 3])
  ind.rename("apple")
  ind
  ind.set_names(["apple"], inplace=True)
  ind.name = "bob"
  ind

Adding an index to an existing DataFrame
----------------------------------------

Occasionally you will load or create a data set into a DataFrame and want to
add an index after you've already done so. There are a couple of different
ways.

Add an index using DataFrame columns
------------------------------------

.. _indexing.set_index:

DataFrame has a ``set_index`` method which takes a column name (for a regular
``Index``) or a list of column names (for a ``MultiIndex``), to create a new,
indexed DataFrame:

.. ipython:: python
   :suppress:

   data = DataFrame({'a' : ['bar', 'bar', 'foo', 'foo'],
                     'b' : ['one', 'two', 'one', 'two'],
                     'c' : ['z', 'y', 'x', 'w'],
                     'd' : [1., 2., 3, 4]})

.. ipython:: python

   data
   indexed1 = data.set_index('c')
   indexed1
   indexed2 = data.set_index(['a', 'b'])
   indexed2

The ``append`` keyword option allow you to keep the existing index and append
the given columns to a MultiIndex:

.. ipython:: python

   frame = data.set_index('c', drop=False)
   frame = frame.set_index(['a', 'b'], append=True)
   frame

Other options in ``set_index`` allow you not drop the index columns or to add
the index in-place (without creating a new object):

.. ipython:: python

   data.set_index('c', drop=False)
   data.set_index(['a', 'b'], inplace=True)
   data

Remove / reset the index,  ``reset_index``
------------------------------------------

As a convenience, there is a new function on DataFrame called ``reset_index``
which transfers the index values into the DataFrame's columns and sets a simple
integer index. This is the inverse operation to ``set_index``

.. ipython:: python

   data
   data.reset_index()

The output is more similar to a SQL table or a record array. The names for the
columns derived from the index are the ones stored in the ``names`` attribute.

You can use the ``level`` keyword to remove only a portion of the index:

.. ipython:: python

   frame
   frame.reset_index(level=1)


``reset_index`` takes an optional parameter ``drop`` which if true simply
discards the index, instead of putting index values in the DataFrame's columns.

.. note::

   The ``reset_index`` method used to be called ``delevel`` which is now
   deprecated.

Adding an ad hoc index
----------------------

If you create an index yourself, you can just assign it to the ``index`` field:

.. code-block:: python

   data.index = index

Indexing internal details
-------------------------

.. note::

   The following is largely relevant for those actually working on the pandas
   codebase. The source code is still the best place to look at the specifics
   of how things are implemented.

In pandas there are a few objects implemented which can serve as valid
containers for the axis labels:

  - ``Index``: the generic "ordered set" object, an ndarray of object dtype
    assuming nothing about its contents. The labels must be hashable (and
    likely immutable) and unique. Populates a dict of label to location in
    Cython to do :math:`O(1)` lookups.
  - ``Int64Index``: a version of ``Index`` highly optimized for 64-bit integer
    data, such as time stamps
  - ``MultiIndex``: the standard hierarchical index object
  - ``PeriodIndex``: An Index object with Period elements
  - ``DatetimeIndex``: An Index object with Timestamp elements
  - ``date_range``: fixed frequency date range generated from a time rule or
    DateOffset. An ndarray of Python datetime objects

The motivation for having an ``Index`` class in the first place was to enable
different implementations of indexing. This means that it's possible for you,
the user, to implement a custom ``Index`` subclass that may be better suited to
a particular application than the ones provided in pandas.

From an internal implementation point of view, the relevant methods that an
``Index`` must define are one or more of the following (depending on how
incompatible the new object internals are with the ``Index`` functions):

  - ``get_loc``: returns an "indexer" (an integer, or in some cases a
    slice object) for a label
  - ``slice_locs``: returns the "range" to slice between two labels
  - ``get_indexer``: Computes the indexing vector for reindexing / data
    alignment purposes. See the source / docstrings for more on this
  - ``get_indexer_non_unique``: Computes the indexing vector for reindexing / data
    alignment purposes when the index is non-unique. See the source / docstrings
    for more on this
  - ``reindex``: Does any pre-conversion of the input index then calls
    ``get_indexer``
  - ``union``, ``intersection``: computes the union or intersection of two
    Index objects
  - ``insert``: Inserts a new label into an Index, yielding a new object
  - ``delete``: Delete a label, yielding a new object
  - ``drop``: Deletes a set of labels
  - ``take``: Analogous to ndarray.take
