.. currentmodule:: pandas
.. _text:

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import lrange
   options.display.max_rows=15

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
Here we are removing leading and trailing whitespaces, lowercasing all names,
and replacing any remaining whitespaces with underscores:

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

Easy to expand this to return a DataFrame using ``expand``.

.. ipython:: python

   s2.str.split('_', expand=True)

It is also possible to limit the number of splits:

.. ipython:: python

   s2.str.split('_', expand=True, n=1)

``rsplit`` is similar to ``split`` except it works in the reverse direction,
i.e., from the end of the string to the beginning of the string:

.. ipython:: python

   s2.str.rsplit('_', expand=True, n=1)

Methods like ``replace`` and ``findall`` take `regular expressions
<https://docs.python.org/2/library/re.html>`__, too:

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

The method ``extract`` (introduced in version 0.13) accepts `regular expressions
<https://docs.python.org/2/library/re.html>`__ with match groups. Extracting a
regular expression with one group returns a Series of strings.

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('[ab](\d)')

Elements that do not match return ``NaN``. Extracting a regular expression
with more than one group returns a DataFrame with one column per group.

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('([ab])(\d)')

Elements that do not match return a row filled with ``NaN``.
Thus, a Series of messy strings can be "converted" into a
like-indexed Series or DataFrame of cleaned-up or more useful strings,
without necessitating ``get()`` to access tuples or ``re.match`` objects.

The results dtype always is object, even if no match is found and the result
only contains ``NaN``.

Named groups like

.. ipython:: python

   pd.Series(['a1', 'b2', 'c3']).str.extract('(?P<letter>[ab])(?P<digit>\d)')

and optional groups like

.. ipython:: python

   pd.Series(['a1', 'b2', '3']).str.extract('(?P<letter>[ab])?(?P<digit>\d)')

can also be used.

Testing for Strings that Match or Contain a Pattern
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can check whether elements contain a pattern:

.. ipython:: python

   pattern = r'[a-z][0-9]'
   pd.Series(['1', '2', '3a', '3b', '03c']).str.contains(pattern)

or match a pattern:


.. ipython:: python

   pd.Series(['1', '2', '3a', '3b', '03c']).str.match(pattern, as_indexer=True)

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

   s4 = pd.Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
   s4.str.contains('A', na=False)

Creating Indicator Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can extract dummy variables from string columns.
For example if they are separated by a ``'|'``:

  .. ipython:: python

      s = pd.Series(['a', 'a|b', np.nan, 'a|c'])
      s.str.get_dummies(sep='|')

See also :func:`~pandas.get_dummies`.

Method Summary
--------------

.. _text.summary:

.. csv-table::
    :header: "Method", "Description"
    :widths: 20, 80

    :meth:`~Series.str.cat`,Concatenate strings
    :meth:`~Series.str.split`,Split strings on delimiter
    :meth:`~Series.str.rsplit`,Split strings on delimiter working from the end of the string
    :meth:`~Series.str.get`,Index into each element (retrieve i-th element)
    :meth:`~Series.str.join`,Join strings in each element of the Series with passed separator
    :meth:`~Series.str.contains`,Return boolean array if each string contains pattern/regex
    :meth:`~Series.str.replace`,Replace occurrences of pattern/regex with some other string
    :meth:`~Series.str.repeat`,Duplicate values (``s.str.repeat(3)`` equivalent to ``x * 3``)
    :meth:`~Series.str.pad`,"Add whitespace to left, right, or both sides of strings"
    :meth:`~Series.str.center`,Equivalent to ``str.center``
    :meth:`~Series.str.ljust`,Equivalent to ``str.ljust``
    :meth:`~Series.str.rjust`,Equivalent to ``str.rjust``
    :meth:`~Series.str.zfill`,Equivalent to ``str.zfill``
    :meth:`~Series.str.wrap`,Split long strings into lines with length less than a given width
    :meth:`~Series.str.slice`,Slice each string in the Series
    :meth:`~Series.str.slice_replace`,Replace slice in each string with passed value
    :meth:`~Series.str.count`,Count occurrences of pattern
    :meth:`~Series.str.startswith`,Equivalent to ``str.startswith(pat)`` for each element
    :meth:`~Series.str.endswith`,Equivalent to ``str.endswith(pat)`` for each element
    :meth:`~Series.str.findall`,Compute list of all occurrences of pattern/regex for each string
    :meth:`~Series.str.match`,"Call ``re.match`` on each element, returning matched groups as list"
    :meth:`~Series.str.extract`,"Call ``re.match`` on each element, as ``match`` does, but return matched groups as strings for convenience."
    :meth:`~Series.str.len`,Compute string lengths
    :meth:`~Series.str.strip`,Equivalent to ``str.strip``
    :meth:`~Series.str.rstrip`,Equivalent to ``str.rstrip``
    :meth:`~Series.str.lstrip`,Equivalent to ``str.lstrip``
    :meth:`~Series.str.partition`,Equivalent to ``str.partition``
    :meth:`~Series.str.rpartition`,Equivalent to ``str.rpartition``
    :meth:`~Series.str.lower`,Equivalent to ``str.lower``
    :meth:`~Series.str.upper`,Equivalent to ``str.upper``
    :meth:`~Series.str.find`,Equivalent to ``str.find``
    :meth:`~Series.str.rfind`,Equivalent to ``str.rfind``
    :meth:`~Series.str.index`,Equivalent to ``str.index``
    :meth:`~Series.str.rindex`,Equivalent to ``str.rindex``
    :meth:`~Series.str.capitalize`,Equivalent to ``str.capitalize``
    :meth:`~Series.str.swapcase`,Equivalent to ``str.swapcase``
    :meth:`~Series.str.normalize`,Return Unicode normal form. Equivalent to ``unicodedata.normalize``
    :meth:`~Series.str.translate`,Equivalent to ``str.translate``
    :meth:`~Series.str.isalnum`,Equivalent to ``str.isalnum``
    :meth:`~Series.str.isalpha`,Equivalent to ``str.isalpha``
    :meth:`~Series.str.isdigit`,Equivalent to ``str.isdigit``
    :meth:`~Series.str.isspace`,Equivalent to ``str.isspace``
    :meth:`~Series.str.islower`,Equivalent to ``str.islower``
    :meth:`~Series.str.isupper`,Equivalent to ``str.isupper``
    :meth:`~Series.str.istitle`,Equivalent to ``str.istitle``
    :meth:`~Series.str.isnumeric`,Equivalent to ``str.isnumeric``
    :meth:`~Series.str.isdecimal`,Equivalent to ``str.isdecimal``
