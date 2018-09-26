.. currentmodule:: pandas
.. _text:

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import lrange
   pd.options.display.max_rows=15

======================
Working with Text Data
======================

.. _text.string_methods:

Series and Index are equipped with a set of string processing methods
that make it easy to operate on each element of the array. Perhaps most
importantly, these methods exclude missing/NA values automatically. These are
accessed via the ``str`` attribute and generally have names matching
the equivalent (scalar) built-in string methods:

.. ipython:: python

   s = pd.Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
   s.str.lower()
   s.str.upper()
   s.str.len()

.. ipython:: python

   idx = pd.Index([' jack', 'jill ', ' jesse ', 'frank'])
   idx.str.strip()
   idx.str.lstrip()
   idx.str.rstrip()

The string methods on Index are especially useful for cleaning up or
transforming DataFrame columns. For instance, you may have columns with
leading or trailing whitespace:

.. ipython:: python

   df = pd.DataFrame(randn(3, 2), columns=[' Column A ', ' Column B '],
                     index=range(3))
   df

Since ``df.columns`` is an Index object, we can use the ``.str`` accessor

.. ipython:: python

   df.columns.str.strip()
   df.columns.str.lower()

These string methods can then be used to clean up the columns as needed.
Here we are removing leading and trailing white spaces, lower casing all names,
and replacing any remaining white spaces with underscores:

.. ipython:: python

   df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')
   df

.. note::

    If you have a ``Series`` where lots of elements are repeated
    (i.e. the number of unique elements in the ``Series`` is a lot smaller than the length of the
    ``Series``), it can be faster to convert the original ``Series`` to one of type
    ``category`` and then use ``.str.<method>`` or ``.dt.<property>`` on that.
    The performance difference comes from the fact that, for ``Series`` of type ``category``, the
    string operations are done on the ``.categories`` and not on each element of the
    ``Series``.

    Please note that a ``Series`` of type ``category`` with string ``.categories`` has
    some limitations in comparison of ``Series`` of type string (e.g. you can't add strings to
    each other: ``s + " " + s`` won't work if ``s`` is a ``Series`` of type ``category``). Also,
    ``.str`` methods which operate on elements of type ``list`` are not available on such a
    ``Series``.


Splitting and Replacing Strings
-------------------------------

.. _text.split:

Methods like ``split`` return a Series of lists:

.. ipython:: python

   s2 = pd.Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])
   s2.str.split('_')

Elements in the split lists can be accessed using ``get`` or ``[]`` notation:

.. ipython:: python

   s2.str.split('_').str.get(1)
   s2.str.split('_').str[1]

It is easy to expand this to return a DataFrame using ``expand``.

.. ipython:: python

   s2.str.split('_', expand=True)

It is also possible to limit the number of splits:

.. ipython:: python

   s2.str.split('_', expand=True, n=1)

``rsplit`` is similar to ``split`` except it works in the reverse direction,
i.e., from the end of the string to the beginning of the string:

.. ipython:: python

   s2.str.rsplit('_', expand=True, n=1)

``replace`` by default replaces `regular expressions
<https://docs.python.org/3/library/re.html>`__:

.. ipython:: python

   s3 = pd.Series(['A', 'B', 'C', 'Aaba', 'Baca',
                  '', np.nan, 'CABA', 'dog', 'cat'])
   s3
   s3.str.replace('^.a|dog', 'XX-XX ', case=False)

Some caution must be taken to keep regular expressions in mind! For example, the
following code will cause trouble because of the regular expression meaning of
`$`:

.. ipython:: python

   # Consider the following badly formatted financial data
   dollars = pd.Series(['12', '-$10', '$10,000'])

   # This does what you'd naively expect:
   dollars.str.replace('$', '')

   # But this doesn't:
   dollars.str.replace('-$', '-')

   # We need to escape the special character (for >1 len patterns)
   dollars.str.replace(r'-\$', '-')

.. versionadded:: 0.23.0

If you do want literal replacement of a string (equivalent to
:meth:`str.replace`), you can set the optional ``regex`` parameter to
``False``, rather than escaping each character. In this case both ``pat``
and ``repl`` must be strings:

.. ipython:: python

    # These lines are equivalent
    dollars.str.replace(r'-\$', '-')
    dollars.str.replace('-$', '-', regex=False)

.. versionadded:: 0.20.0

The ``replace`` method can also take a callable as replacement. It is called
on every ``pat`` using :func:`re.sub`. The callable should expect one
positional argument (a regex object) and return a string.

.. ipython:: python

   # Reverse every lowercase alphabetic word
   pat = r'[a-z]+'
   repl = lambda m: m.group(0)[::-1]
   pd.Series(['foo 123', 'bar baz', np.nan]).str.replace(pat, repl)

   # Using regex groups
   pat = r"(?P<one>\w+) (?P<two>\w+) (?P<three>\w+)"
   repl = lambda m: m.group('two').swapcase()
   pd.Series(['Foo Bar Baz', np.nan]).str.replace(pat, repl)

.. versionadded:: 0.20.0

The ``replace`` method also accepts a compiled regular expression object
from :func:`re.compile` as a pattern. All flags should be included in the
compiled regular expression object.

.. ipython:: python

   import re
   regex_pat = re.compile(r'^.a|dog', flags=re.IGNORECASE)
   s3.str.replace(regex_pat, 'XX-XX ')

Including a ``flags`` argument when calling ``replace`` with a compiled
regular expression object will raise a ``ValueError``.

.. ipython::

    @verbatim
    In [1]: s3.str.replace(regex_pat, 'XX-XX ', flags=re.IGNORECASE)
    ---------------------------------------------------------------------------
    ValueError: case and flags cannot be set when pat is a compiled regex

.. _text.concatenate:

Concatenation
-------------

There are several ways to concatenate a ``Series`` or ``Index``, either with itself or others, all based on :meth:`~Series.str.cat`,
resp. ``Index.str.cat``.

Concatenating a single Series into a string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The content of a ``Series`` (or ``Index``) can be concatenated:

.. ipython:: python

    s = pd.Series(['a', 'b', 'c', 'd'])
    s.str.cat(sep=',')
    
If not specified, the keyword ``sep`` for the separator defaults to the empty string, ``sep=''``:

.. ipython:: python

    s.str.cat()

By default, missing values are ignored. Using ``na_rep``, they can be given a representation:

.. ipython:: python

    t = pd.Series(['a', 'b', np.nan, 'd'])
    t.str.cat(sep=',')
    t.str.cat(sep=',', na_rep='-')

Concatenating a Series and something list-like into a Series
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first argument to :meth:`~Series.str.cat` can be a list-like object, provided that it matches the length of the calling ``Series`` (or ``Index``).

.. ipython:: python

    s.str.cat(['A', 'B', 'C', 'D'])
    
Missing values on either side will result in missing values in the result as well, *unless* ``na_rep`` is specified:

.. ipython:: python

    s.str.cat(t)
    s.str.cat(t, na_rep='-')

Concatenating a Series and something array-like into a Series
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 0.23.0

The parameter ``others`` can also be two-dimensional. In this case, the number or rows must match the lengths of the calling ``Series`` (or ``Index``).

.. ipython:: python

    d = pd.concat([t, s], axis=1)
    s
    d
    s.str.cat(d, na_rep='-')
    
Concatenating a Series and an indexed object into a Series, with alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 0.23.0

For concatenation with a ``Series`` or ``DataFrame``, it is possible to align the indexes before concatenation by setting
the ``join``-keyword.

.. ipython:: python
   :okwarning:

   u = pd.Series(['b', 'd', 'a', 'c'], index=[1, 3, 0, 2])
   s
   u
   s.str.cat(u)
   s.str.cat(u, join='left')

.. warning::

    If the ``join`` keyword is not passed, the method :meth:`~Series.str.cat` will currently fall back to the behavior before version 0.23.0 (i.e. no alignment),
    but a ``FutureWarning`` will be raised if any of the involved indexes differ, since this default will change to ``join='left'`` in a future version.

The usual options are available for ``join`` (one of ``'left', 'outer', 'inner', 'right'``).
In particular, alignment also means that the different lengths do not need to coincide anymore.

.. ipython:: python

    v = pd.Series(['z', 'a', 'b', 'd', 'e'], index=[-1, 0, 1, 3, 4])
    s
    v
    s.str.cat(v, join='left', na_rep='-')
    s.str.cat(v, join='outer', na_rep='-')

The same alignment can be used when ``others`` is a ``DataFrame``:

.. ipython:: python

    f = d.loc[[3, 2, 1, 0], :]
    s
    f
    s.str.cat(f, join='left', na_rep='-')

Concatenating a Series and many objects into a Series
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All one-dimensional list-likes can be combined in a list-like container (including iterators, ``dict``-views, etc.):

.. ipython:: python

    s
    u
    s.str.cat([u.values, ['A', 'B', 'C', 'D'], map(str, u.index)], na_rep='-')

All elements must match in length to the calling ``Series`` (or ``Index``), except those having an index if ``join`` is not None:

.. ipython:: python

    v
    s.str.cat([u, v, ['A', 'B', 'C', 'D']], join='outer', na_rep='-')

If using ``join='right'`` on a list of ``others`` that contains different indexes,
the union of these indexes will be used as the basis for the final concatenation:

.. ipython:: python

    u.loc[[3]]
    v.loc[[-1, 0]]
    s.str.cat([u.loc[[3]], v.loc[[-1, 0]]], join='right', na_rep='-')

Indexing with ``.str``
----------------------

.. _text.indexing:

You can use ``[]`` notation to directly index by position locations. If you index past the end
of the string, the result will be a ``NaN``.


.. ipython:: python

   s = pd.Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan,
                  'CABA', 'dog', 'cat'])

   s.str[0]
   s.str[1]

Extracting Substrings
---------------------

.. _text.extract:

Extract first match in each subject (extract)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   In version 0.18.0, ``extract`` gained the ``expand`` argument. When
   ``expand=False`` it returns a ``Series``, ``Index``, or
   ``DataFrame``, depending on the subject and regular expression
   pattern (same behavior as pre-0.18.0). When ``expand=True`` it
   always returns a ``DataFrame``, which is more consistent and less
   confusing from the perspective of a user. ``expand=True`` is the
   default since version 0.23.0.

The ``extract`` method accepts a `regular expression
<https://docs.python.org/3/library/re.html>`__ with at least one
capture group.

Extracting a regular expression with more than one group returns a
DataFrame with one column per group.

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('([ab])(\d)', expand=False)

Elements that do not match return a row filled with ``NaN``. Thus, a
Series of messy strings can be "converted" into a like-indexed Series
or DataFrame of cleaned-up or more useful strings, without
necessitating ``get()`` to access tuples or ``re.match`` objects. The
dtype of the result is always object, even if no match is found and
the result only contains ``NaN``.

Named groups like

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('(?P<letter>[ab])(?P<digit>\d)', expand=False)

and optional groups like

.. ipython:: python

   pd.Series(['a1', 'b2', '3']).str.extract('([ab])?(\d)', expand=False)

can also be used. Note that any capture group names in the regular
expression will be used for column names; otherwise capture group
numbers will be used.

Extracting a regular expression with one group returns a ``DataFrame``
with one column if ``expand=True``.

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('[ab](\d)', expand=True)

It returns a Series if ``expand=False``.

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('[ab](\d)', expand=False)

Calling on an ``Index`` with a regex with exactly one capture group
returns a ``DataFrame`` with one column if ``expand=True``.

.. ipython:: python

   s = pd.Series(["a1", "b2", "c3"], ["A11", "B22", "C33"])
   s
   s.index.str.extract("(?P<letter>[a-zA-Z])", expand=True)

It returns an ``Index`` if ``expand=False``.

.. ipython:: python

   s.index.str.extract("(?P<letter>[a-zA-Z])", expand=False)

Calling on an ``Index`` with a regex with more than one capture group
returns a ``DataFrame`` if ``expand=True``.

.. ipython:: python

   s.index.str.extract("(?P<letter>[a-zA-Z])([0-9]+)", expand=True)

It raises ``ValueError`` if ``expand=False``.

.. code-block:: python

    >>> s.index.str.extract("(?P<letter>[a-zA-Z])([0-9]+)", expand=False)
    ValueError: only one regex group is supported with Index

The table below summarizes the behavior of ``extract(expand=False)``
(input subject in first column, number of groups in regex in
first row)

+--------+---------+------------+
|        | 1 group | >1 group   |
+--------+---------+------------+
| Index  | Index   | ValueError |
+--------+---------+------------+
| Series | Series  | DataFrame  |
+--------+---------+------------+

Extract all matches in each subject (extractall)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _text.extractall:

.. versionadded:: 0.18.0

Unlike ``extract`` (which returns only the first match),

.. ipython:: python

   s = pd.Series(["a1a2", "b1", "c1"], index=["A", "B", "C"])
   s
   two_groups = '(?P<letter>[a-z])(?P<digit>[0-9])'
   s.str.extract(two_groups, expand=True)

the ``extractall`` method returns every match. The result of
``extractall`` is always a ``DataFrame`` with a ``MultiIndex`` on its
rows. The last level of the ``MultiIndex`` is named ``match`` and
indicates the order in the subject.

.. ipython:: python

   s.str.extractall(two_groups)

When each subject string in the Series has exactly one match,

.. ipython:: python

   s = pd.Series(['a3', 'b3', 'c2'])
   s

then ``extractall(pat).xs(0, level='match')`` gives the same result as
``extract(pat)``.

.. ipython:: python

   extract_result = s.str.extract(two_groups, expand=True)
   extract_result
   extractall_result = s.str.extractall(two_groups)
   extractall_result
   extractall_result.xs(0, level="match")

``Index`` also supports ``.str.extractall``. It returns a ``DataFrame`` which has the
same result as a ``Series.str.extractall`` with a default index (starts from 0).

.. versionadded:: 0.19.0

.. ipython:: python

   pd.Index(["a1a2", "b1", "c1"]).str.extractall(two_groups)

   pd.Series(["a1a2", "b1", "c1"]).str.extractall(two_groups)


Testing for Strings that Match or Contain a Pattern
---------------------------------------------------

You can check whether elements contain a pattern:

.. ipython:: python

   pattern = r'[0-9][a-z]'
   pd.Series(['1', '2', '3a', '3b', '03c']).str.contains(pattern)

Or whether elements match a pattern:

.. ipython:: python

   pd.Series(['1', '2', '3a', '3b', '03c']).str.match(pattern)

The distinction between ``match`` and ``contains`` is strictness: ``match``
relies on strict ``re.match``, while ``contains`` relies on ``re.search``.

Methods like ``match``, ``contains``, ``startswith``, and ``endswith`` take
an extra ``na`` argument so missing values can be considered True or False:

.. ipython:: python

   s4 = pd.Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
   s4.str.contains('A', na=False)

.. _text.indicator:

Creating Indicator Variables
----------------------------

You can extract dummy variables from string columns.
For example if they are separated by a ``'|'``:

.. ipython:: python

    s = pd.Series(['a', 'a|b', np.nan, 'a|c'])
    s.str.get_dummies(sep='|')

String ``Index`` also supports ``get_dummies`` which returns a ``MultiIndex``.

.. versionadded:: 0.18.1

.. ipython:: python

    idx = pd.Index(['a', 'a|b', np.nan, 'a|c'])
    idx.str.get_dummies(sep='|')

See also :func:`~pandas.get_dummies`.

Method Summary
--------------

.. _text.summary:

.. csv-table::
    :header: "Method", "Description"
    :widths: 20, 80
    :delim: ;

    :meth:`~Series.str.cat`;Concatenate strings
    :meth:`~Series.str.split`;Split strings on delimiter
    :meth:`~Series.str.rsplit`;Split strings on delimiter working from the end of the string
    :meth:`~Series.str.get`;Index into each element (retrieve i-th element)
    :meth:`~Series.str.join`;Join strings in each element of the Series with passed separator
    :meth:`~Series.str.get_dummies`;Split strings on the delimiter returning DataFrame of dummy variables
    :meth:`~Series.str.contains`;Return boolean array if each string contains pattern/regex
    :meth:`~Series.str.replace`;Replace occurrences of pattern/regex/string with some other string or the return value of a callable given the occurrence
    :meth:`~Series.str.repeat`;Duplicate values (``s.str.repeat(3)`` equivalent to ``x * 3``)
    :meth:`~Series.str.pad`;"Add whitespace to left, right, or both sides of strings"
    :meth:`~Series.str.center`;Equivalent to ``str.center``
    :meth:`~Series.str.ljust`;Equivalent to ``str.ljust``
    :meth:`~Series.str.rjust`;Equivalent to ``str.rjust``
    :meth:`~Series.str.zfill`;Equivalent to ``str.zfill``
    :meth:`~Series.str.wrap`;Split long strings into lines with length less than a given width
    :meth:`~Series.str.slice`;Slice each string in the Series
    :meth:`~Series.str.slice_replace`;Replace slice in each string with passed value
    :meth:`~Series.str.count`;Count occurrences of pattern
    :meth:`~Series.str.startswith`;Equivalent to ``str.startswith(pat)`` for each element
    :meth:`~Series.str.endswith`;Equivalent to ``str.endswith(pat)`` for each element
    :meth:`~Series.str.findall`;Compute list of all occurrences of pattern/regex for each string
    :meth:`~Series.str.match`;"Call ``re.match`` on each element, returning matched groups as list"
    :meth:`~Series.str.extract`;"Call ``re.search`` on each element, returning DataFrame with one row for each element and one column for each regex capture group"
    :meth:`~Series.str.extractall`;"Call ``re.findall`` on each element, returning DataFrame with one row for each match and one column for each regex capture group"
    :meth:`~Series.str.len`;Compute string lengths
    :meth:`~Series.str.strip`;Equivalent to ``str.strip``
    :meth:`~Series.str.rstrip`;Equivalent to ``str.rstrip``
    :meth:`~Series.str.lstrip`;Equivalent to ``str.lstrip``
    :meth:`~Series.str.partition`;Equivalent to ``str.partition``
    :meth:`~Series.str.rpartition`;Equivalent to ``str.rpartition``
    :meth:`~Series.str.lower`;Equivalent to ``str.lower``
    :meth:`~Series.str.upper`;Equivalent to ``str.upper``
    :meth:`~Series.str.find`;Equivalent to ``str.find``
    :meth:`~Series.str.rfind`;Equivalent to ``str.rfind``
    :meth:`~Series.str.index`;Equivalent to ``str.index``
    :meth:`~Series.str.rindex`;Equivalent to ``str.rindex``
    :meth:`~Series.str.capitalize`;Equivalent to ``str.capitalize``
    :meth:`~Series.str.swapcase`;Equivalent to ``str.swapcase``
    :meth:`~Series.str.normalize`;Return Unicode normal form. Equivalent to ``unicodedata.normalize``
    :meth:`~Series.str.translate`;Equivalent to ``str.translate``
    :meth:`~Series.str.isalnum`;Equivalent to ``str.isalnum``
    :meth:`~Series.str.isalpha`;Equivalent to ``str.isalpha``
    :meth:`~Series.str.isdigit`;Equivalent to ``str.isdigit``
    :meth:`~Series.str.isspace`;Equivalent to ``str.isspace``
    :meth:`~Series.str.islower`;Equivalent to ``str.islower``
    :meth:`~Series.str.isupper`;Equivalent to ``str.isupper``
    :meth:`~Series.str.istitle`;Equivalent to ``str.istitle``
    :meth:`~Series.str.isnumeric`;Equivalent to ``str.isnumeric``
    :meth:`~Series.str.isdecimal`;Equivalent to ``str.isdecimal``
