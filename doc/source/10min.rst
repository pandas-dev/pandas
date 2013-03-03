.. _10min:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   import os
   np.random.seed(123456)
   from pandas import *
   import pandas as pd
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)

   #### portions of this were borrowed from the
   #### Pandas cheatsheet  
   #### created during the PyData Workshop-Sprint 2012 
   #### Hannah Chen, Henry Chow, Eric Cox, Robert Mauriello 


********************
10 Minutes to Pandas
********************

This is a short introduction to pandas, geared mainly for new users.

Customarily, we import as follows

.. ipython:: python

   import pandas as pd
   import numpy as np

Object Creation
---------------

See the :ref:`Data Structure Intro section <dsintro>`

Creating a ``Series`` by passing a list of values, letting pandas create a default 
integer index

.. ipython:: python

   s = pd.Series([1,3,5,np.nan,6,8])
   s

Creating a ``DataFrame`` by passing a numpy array, with a datetime index and labeled columns.

.. ipython:: python

   dates = pd.date_range('20130101',periods=6)
   dates
   df = pd.DataFrame(np.random.randn(6,4),index=dates,columns=list('ABCD'))
   df

Creating a ``DataFrame`` by passing a dict of objects that can be converted to series-like.

.. ipython:: python

   df2 = pd.DataFrame({ 'A' : 1., 
                        'B' : pd.Timestamp('20130102'), 
                        'C' : pd.Series(1,index=range(4),dtype='float32'),
                        'D' : np.array([3] * 4,dtype='int32'), 
                        'E' : 'foo' })
   df2

Having specific dtypes

.. ipython:: python

   df2.dtypes

Viewing Data
------------

See the :ref:`Basics section <basics>`

See the top & bottom rows of the frame

.. ipython:: python

   df.head()
   df.tail()

Display the index,columns, and the underlying numpy data

.. ipython:: python

   df.index
   df.columns
   df.values

Describe shows a quick statistic summary of your data

.. ipython:: python

   df.describe()

Selection
---------

See the :ref:`Indexing section <indexing>`


Getting
~~~~~~~

Selecting a single column, which yields a ``Series``

.. ipython:: python

   df['A']

Selecting via ``[]``, which slices the rows.

.. ipython:: python

   df[0:3]
   df['20130102':'20130104']

Selection by Label
~~~~~~~~~~~~~~~~~~

For getting a cross section using a label

.. ipython:: python

   df.loc[dates[0]]

Selecting on a multi-axis by label

.. ipython:: python

   df.loc[:,['A','B']]

Showing label slicing, both endpoints are *included*

.. ipython:: python

   df.loc['20130102':'20130104',['A','B']]

Reduction in the dimensions of the returned object

.. ipython:: python

   df.loc['20130102',['A','B']]

For getting a scalar value

.. ipython:: python

   df.loc[dates[0],'A']

For getting fast access to a scalar (equiv to the prior method)

.. ipython:: python

   df.at[dates[0],'A']

Selection by Position
~~~~~~~~~~~~~~~~~~~~~

Select via the position of the passed integers

.. ipython:: python

   # this is a cross-section of the object
   df.iloc[3]

By integer slices, acting similar to numpy/python

.. ipython:: python

   df.iloc[3:5,0:2]

By lists of integer position locations, similar to the numpy/python style

.. ipython:: python

   df.iloc[[1,2,4],[0,2]]

For slicing rows explicitly

.. ipython:: python

   df.iloc[1:3,:]

For slicing columns explicitly

.. ipython:: python

   df.iloc[:,1:3]

For getting a value explicity

.. ipython:: python

   df.iloc[1,1]

For getting fast access to a scalar (equiv to the prior method)

.. ipython:: python

   df.iat[1,1]

There is one signficant departure from standard python/numpy slicing semantics.
python/numpy allow slicing past the end of an array without an associated error.

.. ipython:: python

    # these are allowed in python/numpy.
    x = list('abcdef')
    x[4:10]
    x[8:10]

Pandas will detect this and raise ``IndexError``, rather than return an empty structure.

::

    >>> df.iloc[:,3:6]
    IndexError: out-of-bounds on slice (end)

Boolean Indexing
~~~~~~~~~~~~~~~~

Using a single column's values to select data.

.. ipython:: python

   df[df.A > 0]

A ``where`` operation.

.. ipython:: python

   df[df > 0]


Setting
~~~~~~~

Setting a new column automatically aligns the data
by the indexes

.. ipython:: python

   s1 = pd.Series([1,2,3,4,5,6],index=date_range('20130102',periods=6))
   s1
   df['F'] = s1

Setting values by label

.. ipython:: python

   df.at[dates[0],'A'] = 0

Setting values by position

.. ipython:: python

   df.iat[0,1] = 0

Setting by assigning with a numpy array

.. ipython:: python

   df.loc[:,'D'] = np.array([5] * len(df))
   df

Missing Data
------------

Pandas primarily uses the value ``np.nan`` to represent missing data. It
is by default not included in computations. See the :ref:`Missing Data section <missing_data>`

Reindexing allows you to change/add/delete the index on a specified axis. This
returns a copy of the data.

.. ipython:: python

   df1 = df.reindex(index=dates[0:4],columns=list(df.columns) + ['E'])
   df1.loc[dates[0]:dates[1],'E'] = 1
   df1

To drop any rows that have missing data.

.. ipython:: python

   df1.dropna(how='any')

Filling missing data

.. ipython:: python

   df1.fillna(value=5)


Operations
----------

See the :ref:`Basic section on Binary Ops <basics.binop>`

Stats
~~~~~

Performing a descriptive statistic

.. ipython:: python

   df.mean()

Same operation on the other axis

.. ipython:: python

   df.mean(1)

Operations on missing data, exclude the data

.. ipython:: python

  df1.mean()

Apply
~~~~~

Applying functions to the data

.. ipython:: python

   df.apply(np.cumsum)
   df.apply(lambda x: x.max() - x.min())

Merge
-----

Concat
~~~~~~

Pandas provides various facilities for easily combining together Series,
DataFrame, and Panel objects with various kinds of set logic for the indexes
and relational algebra functionality in the case of join / merge-type
operations.

See the :ref:`Merging section <merging>`

Concatenating pandas objects together 

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 4))
   df

   # break it into pieces
   pieces = [df[:3], df[3:7], df[7:]]

   concat(pieces)

Join
~~~~

SQL style merges. See the :ref:`Database style joining <merging.join>`

.. ipython:: python

   left = pd.DataFrame({'key': ['foo', 'foo'], 'lval': [1, 2]})
   right = pd.DataFrame({'key': ['foo', 'foo'], 'rval': [4, 5]})
   left
   right
   merge(left, right, on='key')

Append
~~~~~~

Append rows to a dataframe. See the :ref:`Appending <merging.concatenation>`

.. ipython:: python

   df = pd.DataFrame(np.random.randn(8, 4), columns=['A','B','C','D'])
   df
   s = df.iloc[3]
   df.append(s, ignore_index=True)
   df


Grouping
--------

By "group by" we are referring to a process involving one or more of the following
steps

 - **Splitting** the data into groups based on some criteria
 - **Applying** a function to each group independently
 - **Combining** the results into a data structure

See the :ref:`Grouping section <groupby>`

.. ipython:: python

   df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                            'foo', 'bar', 'foo', 'foo'],
                      'B' : ['one', 'one', 'two', 'three',
                            'two', 'two', 'one', 'three'],
                      'C' : randn(8), 'D' : randn(8)})
   df

Grouping and then applying a function ``sum`` to the resulting groups.

.. ipython:: python

   df.groupby('A').sum()

Grouping by multiple columns forms a hierarchical index, which we then apply the function.

.. ipython:: python

   df.groupby(['A','B']).sum()

Reshaping
---------

See the section on :ref:`Hierarchical Indexing <indexing.hierarchical>` and
see the section on :ref:`Reshaping <reshaping.stacking>`).

.. ipython:: python

   tuples = zip(*[['bar', 'bar', 'baz', 'baz',
                   'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two',
                   'one', 'two', 'one', 'two']])
   index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])
   df = pd.DataFrame(randn(8, 2), index=index, columns=['A', 'B'])
   df2 = df[:4]
   df2

The ``stack`` function "compresses" a level in the DataFrame's columns. to

.. ipython:: python

   stacked = df2.stack()
   stacked

With a "stacked" DataFrame or Series (having a ``MultiIndex`` as the
``index``), the inverse operation of ``stack`` is ``unstack``, which by default
unstacks the **last level**:

.. ipython:: python

   stacked.unstack()
   stacked.unstack(1)
   stacked.unstack(0)

Time Series
-----------

Pandas has simple, powerful, and efficient functionality for
performing resampling operations during frequency conversion (e.g., converting
secondly data into 5-minutely data). This is extremely common in, but not
limited to, financial applications. See the :ref:`Time Series section <timeseries>`

.. ipython:: python

   rng = pd.date_range('1/1/2012', periods=100, freq='S')
   ts = pd.Series(randint(0, 500, len(rng)), index=rng)
   ts.resample('5Min', how='sum')

Time zone representation

.. ipython:: python

   rng = pd.date_range('3/6/2012 00:00', periods=5, freq='D')
   ts = pd.Series(randn(len(rng)), rng)
   ts_utc = ts.tz_localize('UTC')
   ts_utc

Convert to another time zone

.. ipython:: python

   ts_utc.tz_convert('US/Eastern')

Converting between time span representations

.. ipython:: python

   rng = pd.date_range('1/1/2012', periods=5, freq='M')
   ts = pd.Series(randn(len(rng)), index=rng)
   ts
   ps = ts.to_period()
   ps
   ps.to_timestamp()

Converting between period and timestamp enables some convenient arithmetic
functions to be used. In the following example, we convert a quarterly
frequency with year ending in November to 9am of the end of the month following
the quarter end:

.. ipython:: python

   prng = period_range('1990Q1', '2000Q4', freq='Q-NOV')
   ts = Series(randn(len(prng)), prng)
   ts.index = (prng.asfreq('M', 'e') + 1).asfreq('H', 's') + 9
   ts.head()


Plotting
--------

.. ipython:: python
   :suppress:

   import matplotlib.pyplot as plt
   plt.close('all')

.. ipython:: python

   ts = pd.Series(randn(1000), index=pd.date_range('1/1/2000', periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png width=4.5in
   ts.plot()

On DataFrame, ``plot`` is a convenience to plot all of the columns with labels:

.. ipython:: python

   df = pd.DataFrame(randn(1000, 4), index=ts.index,
                     columns=['A', 'B', 'C', 'D'])
   df = df.cumsum()

   @savefig frame_plot_basic.png width=4.5in
   plt.figure(); df.plot(); plt.legend(loc='best')

Getting Data In/Out
-------------------

CSV
~~~

:ref:`Writing to a csv file <io.store_in_csv>`

.. ipython:: python

   df.to_csv('foo.csv')

:ref:`Reading from a csv file <io.read_csv_table>`

.. ipython:: python

   pd.read_csv('foo.csv')

.. ipython:: python
   :suppress:

   os.remove('foo.csv')

HDF5
~~~~

Reading and writing to :ref:`HDFStores <io.hdf5>`

Writing to a HDF5 Store

.. ipython:: python

   store = pd.HDFStore('foo.h5')
   store['df'] = df

Reading from a HDF5 Store

.. ipython:: python

   store['df']

.. ipython:: python
   :suppress:

   store.close()
   os.remove('foo.h5')

