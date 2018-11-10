.. currentmodule:: pandas

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd

 .. _integer_na:

********************************
Integer Data with Missing Values
********************************

.. versionadded:: 0.24.0

In :ref:`missing_data`, we say that pandas primarily uses ``NaN`` to represent
missing data. Because ``NaN`` is a float, this forces an array of integers with
any missing values to become floating point. In some cases, this may not matter
much. But if your integer column is, say, and identifier, casting to float can
lead to bad outcomes.

Pandas can represent integer data with missing values with the
:class:`arrays.IntegerArray` array. This is an :ref:`extension types <extending.extension-types>`
implemented within pandas. It is not the default dtype and will not be inferred,
you must explicitly create an :class:`api.extensions.IntegerArray` using :func:`integer_array`.

.. ipython:: python

   arr = integer_array([1, 2, np.nan])
   arr

This array can be stored in a :class:`DataFrame` or :class:`Series` like any
NumPy array.

.. ipython:: python

   pd.Series(arr)

Alternatively, you can instruct pandas to treat an array-like as an
:class:`api.extensions.IntegerArray` by specifying a dtype with a capital "I".

.. ipython:: python

   s = pd.Series([1, 2, np.nan], dtype="Int64")
   s

Note that by default (if you don't specify `dtype`), NumPy is used, and you'll end
up with a ``float64`` dtype Series:

.. ipython:: python

   pd.Series([1, 2, np.nan])


Operations involving an integer array will behave similar to NumPy arrays.
Missing values will be propagated, and and the data will be coerced to another
dtype if needed.

.. ipython:: python

   # arithmetic
   s + 1

   # comparison
   s == 1

   # indexing
   s.iloc[1:3]

   # operate with other dtypes
   s + s.iloc[1:3].astype('Int8')

   # coerce when needed
   s + 0.01

Reduction and groupby operations such as 'sum' work as well.

.. ipython:: python

   df.sum()
   df.groupby('B').A.sum()
