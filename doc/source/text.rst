.. currentmodule:: pandas
.. _text:

.. ipython:: python
   :suppress:

   import numpy as np
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import lrange
   options.display.max_rows=15

======================
Working with Text Data
======================

.. _text.string_methods:

Series is equipped with a set of string processing methods
that make it easy to operate on each element of the array. Perhaps most
importantly, these methods exclude missing/NA values automatically. These are
accessed via the Series's ``str`` attribute and generally have names matching
the equivalent (scalar) built-in string methods:

.. ipython:: python

   s = Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
   s.str.lower()
   s.str.upper()
   s.str.len()

Splitting and Replacing Strings
-------------------------------

.. _text.split:

Methods like ``split`` return a Series of lists:

.. ipython:: python

   s2 = Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])
   s2.str.split('_')

Easy to expand this to return a DataFrame

.. ipython:: python

   s2.str.split('_').apply(Series)

Elements in the split lists can be accessed using ``get`` or ``[]`` notation:

.. ipython:: python

   s2.str.split('_').str.get(1)
   s2.str.split('_').str[1]

Methods like ``replace`` and ``findall`` take `regular expressions
<https://docs.python.org/2/library/re.html>`__, too:

.. ipython:: python

   s3 = Series(['A', 'B', 'C', 'Aaba', 'Baca',
               '', np.nan, 'CABA', 'dog', 'cat'])
   s3
   s3.str.replace('^.a|dog', 'XX-XX ', case=False)

Some caution must be taken to keep regular expressions in mind! For example, the
following code will cause trouble because of the regular expression meaning of
`$`:

.. ipython:: python

   # Consider the following badly formatted financial data
   dollars = Series(['12', '-$10', '$10,000'])

   # This does what you'd naively expect:
   dollars.str.replace('$', '')

   # But this doesn't:
   dollars.str.replace('-$', '-')

   # We need to escape the special character (for >1 len patterns)
   dollars.str.replace(r'-\$', '-')

Indexing with ``.str``
----------------------

.. _text.indexing:

You can use ``[]`` notation to directly index by position locations. If you index past the end
of the string, the result will be a ``NaN``.


.. ipython:: python

   s = Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan,
               'CABA', 'dog', 'cat'])

   s.str[0]
   s.str[1]

Extracting Substrings
---------------------

.. _text.extract:

The method ``extract`` (introduced in version 0.13) accepts `regular expressions
<https://docs.python.org/2/library/re.html>`__ with match groups. Extracting a
regular expression with one group returns a Series of strings.

.. ipython:: python

   Series(['a1', 'b2', 'c3']).str.extract('[ab](\d)')

Elements that do not match return ``NaN``. Extracting a regular expression
with more than one group returns a DataFrame with one column per group.

.. ipython:: python

   Series(['a1', 'b2', 'c3']).str.extract('([ab])(\d)')

Elements that do not match return a row filled with ``NaN``.
Thus, a Series of messy strings can be "converted" into a
like-indexed Series or DataFrame of cleaned-up or more useful strings,
without necessitating ``get()`` to access tuples or ``re.match`` objects.

The results dtype always is object, even if no match is found and the result
only contains ``NaN``.

Named groups like

.. ipython:: python

   Series(['a1', 'b2', 'c3']).str.extract('(?P<letter>[ab])(?P<digit>\d)')

and optional groups like

.. ipython:: python

   Series(['a1', 'b2', '3']).str.extract('(?P<letter>[ab])?(?P<digit>\d)')

can also be used.

Testing for Strings that Match or Contain a Pattern
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can check whether elements contain a pattern:

.. ipython:: python

   pattern = r'[a-z][0-9]'
   Series(['1', '2', '3a', '3b', '03c']).str.contains(pattern)

or match a pattern:


.. ipython:: python

   Series(['1', '2', '3a', '3b', '03c']).str.match(pattern, as_indexer=True)

The distinction between ``match`` and ``contains`` is strictness: ``match``
relies on strict ``re.match``, while ``contains`` relies on ``re.search``.

.. warning::

   In previous versions, ``match`` was for *extracting* groups,
   returning a not-so-convenient Series of tuples. The new method ``extract``
   (described in the previous section) is now preferred.

   This old, deprecated behavior of ``match`` is still the default. As
   demonstrated above, use the new behavior by setting ``as_indexer=True``.
   In this mode, ``match`` is analogous to ``contains``, returning a boolean
   Series. The new behavior will become the default behavior in a future
   release.

Methods like ``match``, ``contains``, ``startswith``, and ``endswith`` take
 an extra ``na`` argument so missing values can be considered True or False:

.. ipython:: python

   s4 = Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
   s4.str.contains('A', na=False)

Creating Indicator Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can extract dummy variables from string columns.
For example if they are separated by a ``'|'``:

  .. ipython:: python

      s = Series(['a', 'a|b', np.nan, 'a|c'])
      s.str.get_dummies(sep='|')

See also :func:`~pandas.get_dummies`.

Method Summary
--------------

.. _text.summary:

.. csv-table::
    :header: "Method", "Description"
    :widths: 20, 80

    :meth:`~core.strings.StringMethods.cat`,Concatenate strings
    :meth:`~core.strings.StringMethods.split`,Split strings on delimiter
    :meth:`~core.strings.StringMethods.get`,Index into each element (retrieve i-th element)
    :meth:`~core.strings.StringMethods.join`,Join strings in each element of the Series with passed separator
    :meth:`~core.strings.StringMethods.contains`,Return boolean array if each string contains pattern/regex
    :meth:`~core.strings.StringMethods.replace`,Replace occurrences of pattern/regex with some other string
    :meth:`~core.strings.StringMethods.repeat`,Duplicate values (``s.str.repeat(3)`` equivalent to ``x * 3``)
    :meth:`~core.strings.StringMethods.pad`,"Add whitespace to left, right, or both sides of strings"
    :meth:`~core.strings.StringMethods.center`,Equivalent to ``pad(side='both')``
    :meth:`~core.strings.StringMethods.wrap`,Split long strings into lines with length less than a given width
    :meth:`~core.strings.StringMethods.slice`,Slice each string in the Series
    :meth:`~core.strings.StringMethods.slice_replace`,Replace slice in each string with passed value
    :meth:`~core.strings.StringMethods.count`,Count occurrences of pattern
    :meth:`~core.strings.StringMethods.startswith`,Equivalent to ``str.startswith(pat)`` for each element
    :meth:`~core.strings.StringMethods.endswith`,Equivalent to ``str.endswith(pat)`` for each element
    :meth:`~core.strings.StringMethods.findall`,Compute list of all occurrences of pattern/regex for each string
    :meth:`~core.strings.StringMethods.match`,"Call ``re.match`` on each element, returning matched groups as list"
    :meth:`~core.strings.StringMethods.extract`,"Call ``re.match`` on each element, as ``match`` does, but return matched groups as strings for convenience."
    :meth:`~core.strings.StringMethods.len`,Compute string lengths
    :meth:`~core.strings.StringMethods.strip`,Equivalent to ``str.strip``
    :meth:`~core.strings.StringMethods.rstrip`,Equivalent to ``str.rstrip``
    :meth:`~core.strings.StringMethods.lstrip`,Equivalent to ``str.lstrip``
    :meth:`~core.strings.StringMethods.lower`,Equivalent to ``str.lower``
    :meth:`~core.strings.StringMethods.upper`,Equivalent to ``str.upper``
