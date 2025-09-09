.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import pandas as pd
   import numpy as np

.. _boolean:

**************************
Nullable Boolean data type
**************************

.. note::

   BooleanArray is currently experimental. Its API or implementation may
   change without warning.

.. _boolean.indexing:

Indexing with NA values
-----------------------

pandas allows indexing with ``NA`` values in a boolean array, which are treated as ``False``.

.. ipython:: python
   :okexcept:

   s = pd.Series([1, 2, 3])
   mask = pd.array([True, False, pd.NA], dtype="boolean")
   s[mask]

If you would prefer to keep the ``NA`` values you can manually fill them with ``fillna(True)``.

.. ipython:: python

   s[mask.fillna(True)]

If you create a column of ``NA`` values (for example to fill them later)
with ``df['new_col'] = pd.NA``, the ``dtype`` would be set to ``object`` in the
new column. The performance on this column will be worse than with
the appropriate type. It's better to use
``df['new_col'] = pd.Series(pd.NA, dtype="boolean")``
(or another ``dtype`` that supports ``NA``).

.. ipython:: python

   df = pd.DataFrame()
   df['objects'] = pd.NA
   df.dtypes

.. _boolean.kleene:

Kleene logical operations
-------------------------

:class:`arrays.BooleanArray` implements `Kleene Logic`_ (sometimes called three-value logic) for
logical operations like ``&`` (and), ``|`` (or) and ``^`` (exclusive-or).

This table demonstrates the results for every combination. These operations are symmetrical,
so flipping the left- and right-hand side makes no difference in the result.

================= =========
Expression        Result
================= =========
``True & True``   ``True``
``True & False``  ``False``
``True & NA``     ``NA``
``False & False`` ``False``
``False & NA``    ``False``
``NA & NA``       ``NA``
``True | True``   ``True``
``True | False``  ``True``
``True | NA``     ``True``
``False | False`` ``False``
``False | NA``    ``NA``
``NA | NA``       ``NA``
``True ^ True``   ``False``
``True ^ False``  ``True``
``True ^ NA``     ``NA``
``False ^ False`` ``False``
``False ^ NA``    ``NA``
``NA ^ NA``       ``NA``
================= =========

When an ``NA`` is present in an operation, the output value is ``NA`` only if
the result cannot be determined solely based on the other input. For example,
``True | NA`` is ``True``, because both ``True | True`` and ``True | False``
are ``True``. In that case, we don't actually need to consider the value
of the ``NA``.

On the other hand, ``True & NA`` is ``NA``. The result depends on whether
the ``NA`` really is ``True`` or ``False``, since ``True & True`` is ``True``,
but ``True & False`` is ``False``, so we can't determine the output.


This differs from how ``np.nan`` behaves in logical operations. pandas treated
``np.nan`` is *always false in the output*.

In ``or``

.. ipython:: python

   pd.Series([True, False, np.nan], dtype="object") | True
   pd.Series([True, False, np.nan], dtype="boolean") | True

In ``and``

.. ipython:: python

   pd.Series([True, False, np.nan], dtype="object") & True
   pd.Series([True, False, np.nan], dtype="boolean") & True

.. _Kleene Logic: https://en.wikipedia.org/wiki/Three-valued_logic#Kleene_and_Priest_logics
