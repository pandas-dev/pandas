.. currentmodule:: pandas
.. _missing_data:

*************************
Working with missing data
*************************

.. ipython:: python
   :suppress:

   import numpy as np; randn = np.random.randn
   from pandas import *

.. note::

    The choice of using ``NaN`` for missing data was largely for simplicity and
    performance reasons. It differs from the MaskedArray approach of, for
    example, :mod:`scikits.timeseries`. For a discussion of the issues with the
    various approaches, :ref:`see here <missing_data>`. We are hopeful that
    NumPy will be able to natively provide a NA dtype solution performant enough
    to be used in pandas.

Missing data basics
-------------------

When / why does data become missing?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some might quibble over our usage of *missing*. By "missing" we simply mean
**null** or "not present for whatever reason". Many data sets simply arrive with
missing data, either because it exists and was not collected or it never
existed. For example, in a collection of financial time series, some of the time
series might start on different dates. Thus, values prior to the start date
would generally be marked as missing.

In pandas, one of the most common ways that missing data is **introduced** into
a data set is by reindexing. For example

.. ipython:: python

   df = DataFrame(randn(5, 3), index=['a', 'c', 'e', 'f', 'h'],
                  columns=['one', 'two', 'three'])
   df['four'] = 'bar'
   df['five'] = df['one'] > 0
   df
   df2 = df.reindex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
   df2

Values considered "missing"
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As data comes in many shapes and forms, pandas aims to be flexible with regard
to handling missing data. While ``NaN`` is the default missing value marker for
reasons of computational speed and convenience, we need to be able to easily
detect this value with data of different types: floating point, integer,
boolean, and general object. In many cases, however, the Python ``None`` will
arise and we wish to also consider that "missing" or "null". Lastly, for legacy
reasons ``inf`` and ``-inf`` are also considered to be "null" in
computations. Since in NumPy divide-by-zero generates ``inf`` or ``-inf`` and
not ``NaN``, I think you will find this is a worthwhile trade-off (Zen of
Python: "practicality beats purity").

To make detecting missing values easier (and across different array dtypes),
pandas provides the ``isnull`` and ``notnull`` functions:

.. ipython:: python

   df2['one']
   isnull(df2['one'])
   notnull(df2['four'])

**Summary:** ``NaN``, ``inf``, ``-inf``, and ``None`` (in object arrays) are all
considered missing by the ``isnull`` and ``notnull`` functions.

Calculations with missing data
------------------------------

Cleaning / replacing missing data
---------------------------------

Dropping missing rows / columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Filling missing values
~~~~~~~~~~~~~~~~~~~~~~

Missing data casting and indexing rules
---------------------------------------
