.. _duplicates:

****************
Duplicate Labels
****************

:class:`Index` objects are not required to be unique; you can have duplicate row
or column labels. This may be a bit confusing at first. If you're familiar with
SQL, you know that row labels are similar to a primary key on a table, and you
would never want duplicates in a SQL table. But one of pandas' roles is to clean
messy, real-world data before it goes to some downstream system. And real-world
data has duplicates, even in fields that are supposed to be unique.

This section describes how duplicate labels change the behavior of certain
operations, and how prevent duplicates from arising during operations, or to
detect them if they do.

.. ipython:: python

   import pandas as pd
   import numpy as np

Consequences of Duplicate Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some pandas methods (:meth:`Series.reindex` for example) just don't work with
duplicates present. The output can't be determined, and so pandas raises.

.. ipython:: python
   :okexcept:

   s1 = pd.Series([0, 1, 2], index=['a', 'b', 'b'])
   s1.reindex(['a', 'b', 'c'])

Other methods, like indexing, can give very surprising results. Typically
indexing with a scalar will *reduce dimensionality*. Slicing a ``DataFrame``
with a scalar will return a ``Series``. Slicing a ``Series`` with a scalar will
return a scalar. But with duplicates, this isn't the case.

.. ipython:: python

   df1 = pd.DataFrame([[0, 1, 2], [3, 4, 5]], columns=['A', 'A', 'B'])
   df1

We have duplicates in the columns. If we slice ``'B'``, we get back a ``Series``

.. ipython:: python

   df1['B']  # a series

But slicing ``'A'`` returns a ``DataFrame``


.. ipython:: python

   df1['A']  # a DataFrame

This applies to row labels as well

.. ipython:: python

   df2 = pd.DataFrame({"A": [0, 1, 2]}, index=['a', 'a', 'b'])
   df2
   df2.loc['b', 'A']  # a scalar
   df2.loc['a', 'A']  # a Series

Duplicate Label Detection
~~~~~~~~~~~~~~~~~~~~~~~~~

You can check whether an :class:`Index` (storing the row or column labels) is
unique with :attr:`Index.is_unique`:

.. ipython:: python

   df2
   df2.index.is_unique
   df2.columns.is_unique

.. note::

   Checking whether an index is unique is somewhat expensive for large datasets.
   Pandas does cache this result, so re-checking on the same index is very fast.

:meth:`Index.duplicated` will return a boolean ndarray indicating whether a
label is a repeat.

.. ipython:: python

   df2.index.duplicated()

Which can be used as a boolean filter to drop duplicate rows.

.. ipython:: python

   df2.loc[~df2.index.duplicated(), :]

If you need additional logic to handle duplicate labels, rather than just
dropping the repeats, using :meth:`~DataFrame.groupby` on the index is a common
trick. For example, we'll resolve duplicates by taking the average of all rows
with the same label.

.. ipython:: python

   df2.groupby(level=0).mean()

.. _duplicates.disallow:

Disallowing Duplicate Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 1.1.0

As noted above, handling duplicates is an important feature when reading in raw
data. That said, you may want to avoid introducing duplicates as part of a data
processing pipeline (from methods like :meth:`pandas.concat`,
:meth:`~DataFrame.rename`, etc.). Both :class:`Series` and :class:`DataFrame`
can be created with the argument ``allows_duplicate_labels=False`` to *disallow*
duplicate labels (the default is to allow them). If there are duplicate labels,
an exception will be raised.

.. ipython:: python
   :okexcept:

   pd.Series([0, 1, 2], index=['a', 'b', 'b'], allows_duplicate_labels=False)

This applies to both row and column labels for a :class:`DataFrame`

.. ipython:: python
   :okexcept:

   pd.DataFrame([[0, 1, 2], [3, 4, 5]], columns=["A", "B", "C"],
                allows_duplicate_labels=False)

This attribute can be checked or set with :attr:`~DataFrame.allows_duplicate_labels`,
which indicates whether that object can have duplicate labels.

.. ipython:: python

   df = pd.DataFrame({"A": [0, 1, 2, 3]},
                     index=['x', 'y', 'X', 'Y'],
                     allows_duplicate_labels=False)
   df
   df.allows_duplicate_labels

Performing an operation that introduces duplicate labels on a ``Series`` or
``DataFrame`` that disallows duplicates will raise an
:class:`errors.DuplicateLabelError`.

.. ipython:: python
   :okexcept:

   df.rename(str.upper)

Duplicate Label Propagation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In general, disallowing duplicates is "sticky". It's preserved through
operations.

.. ipython:: python
   :okexcept:

   s1 = pd.Series(0, index=['a', 'b'], allows_duplicate_labels=False)
   s1
   s1.head().rename({"a": "b"})

.. warning::

   This is an experimental feature. Currently, many methods fail to
   propagate the ``allows_duplicate_labels`` value. In future versions
   it is expected that every method taking or returning one or more
   DataFrame or Series objects will propagate ``allows_duplicate_labels``.
