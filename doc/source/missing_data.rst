.. currentmodule:: pandas
.. _missing_data:

*************************
Working with missing data
*************************

In this section, we will discuss missing (also referred to as NA) values in
pandas.

.. ipython:: python
   :suppress:

   import numpy as np; randn = np.random.randn
   from pandas import *
   import matplotlib.pyplot as plt

.. note::

    The choice of using ``NaN`` internally to denote missing data was largely
    for simplicity and performance reasons. It differs from the MaskedArray
    approach of, for example, :mod:`scikits.timeseries`. We are hopeful that
    NumPy will soon be able to provide a native NA type solution (similar to R)
    performant enough to be used in pandas.

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

.. _missing.isnull:

To make detecting missing values easier (and across different array dtypes),
pandas provides the :func:`~pandas.core.common.isnull` and
:func:`~pandas.core.common.notnull` functions, which are also methods on
``Series`` objects:

.. ipython:: python

   df2['one']
   isnull(df2['one'])
   df2['four'].notnull()

**Summary:** ``NaN``, ``inf``, ``-inf``, and ``None`` (in object arrays) are
all considered missing by the ``isnull`` and ``notnull`` functions.

Calculations with missing data
------------------------------

Missing values propagate naturally through arithmetic operations between pandas
objects.

.. ipython:: python
   :suppress:

   df = df2.ix[:, ['one', 'two', 'three']]
   a = df2.ix[:5, ['one', 'two']].fillna(method='pad')
   b = df2.ix[:5, ['one', 'two', 'three']]

.. ipython:: python

   a
   b
   a + b

The descriptive statistics and computational methods discussed in the
:ref:`data structure overview <basics.stats>` (and listed :ref:`here
<api.series.stats>` and :ref:`here <api.dataframe.stats>`) are all written to
account for missing data. For example:

  * When summing data, NA (missing) values will be treated as zero
  * If the data are all NA, the result will be NA
  * Methods like **cumsum** and **cumprod** ignore NA values, but preserve them
    in the resulting arrays

.. ipython:: python

   df
   df['one'].sum()
   df.mean(1)
   df.cumsum()

NA values in GroupBy
~~~~~~~~~~~~~~~~~~~~

NA groups in GroupBy are automatically excluded. This behavior is consistent
with R, for example.



Cleaning / filling missing data
--------------------------------

pandas objects are equipped with various data manipulation methods for dealing
with missing data.

.. _missing_data.fillna:

Filling missing values: fillna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The **fillna** function can "fill in" NA values with non-null data in a couple
of ways, which we illustrate:

**Replace NA with a scalar value**

.. ipython:: python

   df2
   df2.fillna(0)
   df2['four'].fillna('missing')

**Fill gaps forward or backward**

Using the same filling arguments as :ref:`reindexing <basics.reindexing>`, we
can propagate non-null values forward or backward:

.. ipython:: python

   df
   df.fillna(method='pad')

To remind you, these are the available filling methods:

.. csv-table::
    :header: "Method", "Action"
    :widths: 30, 50

    pad / ffill, Fill values forward
    bfill / backfill, Fill values backward

With time series data, using pad/ffill is extremely common so that the "last
known value" is available at every time point.

.. _missing_data.dropna:

Dropping axis labels with missing data: dropna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may wish to simply exclude labels from a data set which refer to missing
data. To do this, use the **dropna** method:

.. ipython:: python
   :suppress:

   df['two'] = df['two'].fillna(0)
   df['three'] = df['three'].fillna(0)

.. ipython:: python

   df
   df.dropna(axis=0)
   df.dropna(axis=1)
   df['one'].dropna()

**dropna** is presently only implemented for Series and DataFrame, but will be
eventually added to Panel. Series.dropna is a simpler method as it only has one
axis to consider. DataFrame.dropna has considerably more options, which can be
examined :ref:`in the API <api.dataframe.missing>`.

Interpolation
~~~~~~~~~~~~~

A basic linear **interpolate** method has been implemented on Series with
intended use for time series data. There has not been a great deal of demand
for interpolation methods outside of the filling methods described above.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   ts = Series(randn(100), index=DateRange('1/1/2000', periods=100,
                                           timeRule='EOM'))
   ts[20:40] = np.nan
   ts[60:80] = np.nan
   ts = ts.cumsum()

.. ipython:: python

   fig, axes = plt.subplots(ncols=2, figsize=(8, 4))
   ts.plot(ax=axes[0])
   ts.interpolate().plot(ax=axes[1])
   axes[0].set_title('Not interpolated')
   @savefig series_interpolate.png width=6in
   axes[1].set_title('Interpolated')

   plt.close('all')

Missing data casting rules and indexing
---------------------------------------

While pandas supports storing arrays of integer and boolean type, these types
are not capable of storing missing data. Until we can switch to using a native
NA type in NumPy, we've established some "casting rules" when reindexing will
cause missing data to be introduced into, say, a Series or DataFrame. Here they
are:

.. csv-table::
    :header: "data type", "Cast to"
    :widths: 40, 40

	integer, float
    boolean, object
    float, no cast
    object, no cast

For example:

.. ipython:: python

   s = Series(randn(5), index=[0, 2, 4, 6, 7])
   s > 0
   (s > 0).dtype
   crit = (s > 0).reindex(range(8))
   crit
   crit.dtype

Ordinarily NumPy will complain if you try to use an object array (even if it
contains boolean values) instead of a boolean array to get or set values from
an ndarray (e.g. selecting values based on some criteria). If a boolean vector
contains NAs, an exception will be generated:

.. ipython:: python
   :okexcept:

   reindexed = s.reindex(range(8)).fillna(0)
   reindexed[crit]

However, these can be filled in using **fillna** and it will work fine:

.. ipython:: python

   reindexed[crit.fillna(False)]
   reindexed[crit.fillna(True)]

