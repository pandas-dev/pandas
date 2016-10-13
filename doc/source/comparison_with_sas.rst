.. currentmodule:: pandas
.. _compare_with_sas:

Comparison with SAS
********************
For potential users coming from `SAS <https://en.wikipedia.org/wiki/SAS_(software)>`__
this page is meant to demonstrate how different SAS operations would be
performed in pandas.

If you're new to pandas, you might want to first read through :ref:`10 Minutes to pandas<10min>`
to familiarize yourself with the library.

As is customary, we import pandas and numpy as follows:

.. ipython:: python

    import pandas as pd
    import numpy as np


.. note::

   Throughout this tutorial, the pandas ``DataFrame`` will be displayed by calling
   ``df.head()``, which displays the first N (default 5) rows of the ``DataFrame``.
   This is often used in interactive work (e.g. `Jupyter notebook
   <https://jupyter.org/>`_ or terminal) - the equivalent in SAS would be:

   .. code-block:: none

      proc print data=df(obs=5);
      run;

Data Structures
---------------

General Terminology Translation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "pandas", "SAS"
    :widths: 20, 20

    ``DataFrame``, data set
    column, variable
    row, observation
    groupby, BY-group
    ``NaN``, ``.``


``DataFrame`` / ``Series``
~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``DataFrame`` in pandas is analogous to a SAS data set - a two-dimensional
data source with labeled columns that can be of different types. As will be
shown in this document, almost any operation that can be applied to a data set
using SAS's ``DATA`` step, can also be accomplished in pandas.

A ``Series`` is the data structure that represents one column of a
``DataFrame``. SAS doesn't have a separate data structure for a single column,
but in general, working with a ``Series`` is analogous to referencing a column
in the ``DATA`` step.

``Index``
~~~~~~~~~

Every ``DataFrame`` and ``Series`` has an ``Index`` - which are labels on the
*rows* of the data. SAS does not have an exactly analogous concept. A data set's
row are essentially unlabeled, other than an implicit integer index that can be
accessed during the ``DATA`` step (``_N_``).

In pandas, if no index is specified, an integer index is also used by default
(first row = 0, second row = 1, and so on). While using a labeled ``Index`` or
``MultiIndex`` can enable sophisticated analyses and is ultimately an important
part of pandas to understand, for this comparison we will essentially ignore the
``Index`` and just treat the ``DataFrame`` as a collection of columns. Please
see the :ref:`indexing documentation<indexing>` for much more on how to use an
``Index`` effectively.


Data Input / Output
-------------------

Constructing a DataFrame from Values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A SAS data set can be built from specified values by
placing the data after a ``datalines`` statement and
specifying the column names.

.. code-block:: none

   data df;
       input x y;
       datalines;
       1 2
       3 4
       5 6
       ;
   run;

A pandas ``DataFrame`` can be constructed in many different ways,
but for a small number of values, it is often convenient to specify it as
a python dictionary, where the keys are the column names
and the values are the data.

.. ipython:: python

   df = pd.DataFrame({
             'x': [1, 3, 5],
             'y': [2, 4, 6]})
   df


Reading External Data
~~~~~~~~~~~~~~~~~~~~~

Like SAS, pandas provides utilities for reading in data from
many formats.  The ``tips`` dataset, found within the pandas
tests (`csv <https://raw.github.com/pandas-dev/pandas/master/pandas/tests/data/tips.csv>`_)
will be used in many of the following examples.

SAS provides ``PROC IMPORT`` to read csv data into a data set.

.. code-block:: none

   proc import datafile='tips.csv' dbms=csv out=tips replace;
       getnames=yes;
   run;

The pandas method is :func:`read_csv`, which works similarly.

.. ipython:: python

   url = 'https://raw.github.com/pandas-dev/pandas/master/pandas/tests/data/tips.csv'
   tips = pd.read_csv(url)
   tips.head()


Like ``PROC IMPORT``, ``read_csv`` can take a number of parameters to specify
how the data should be parsed.  For example, if the data was instead tab delimited,
and did not have column names, the pandas command would be:

.. code-block:: python

   tips = pd.read_csv('tips.csv', sep='\t', header=None)

   # alternatively, read_table is an alias to read_csv with tab delimiter
   tips = pd.read_table('tips.csv', header=None)

In addition to text/csv, pandas supports a variety of other data formats
such as Excel, HDF5, and SQL databases.  These are all read via a ``pd.read_*``
function.  See the :ref:`IO documentation<io>` for more details.

Exporting Data
~~~~~~~~~~~~~~

The inverse of ``PROC IMPORT`` in SAS is ``PROC EXPORT``

.. code-block:: none

   proc export data=tips outfile='tips2.csv' dbms=csv;
   run;

Similarly in pandas, the opposite of ``read_csv`` is :meth:`~DataFrame.to_csv`,
and other data formats follow a similar api.

.. code-block:: python

   tips.to_csv('tips2.csv')


Data Operations
---------------

Operations on Columns
~~~~~~~~~~~~~~~~~~~~~

In the ``DATA`` step, arbitrary math expressions can
be used on new or existing columns.

.. code-block:: none

   data tips;
       set tips;
       total_bill = total_bill - 2;
       new_bill = total_bill / 2;
   run;

pandas provides similar vectorized operations by
specifying the individual ``Series`` in the ``DataFrame``.
New columns can be assigned in the same way.

.. ipython:: python

   tips['total_bill'] = tips['total_bill'] - 2
   tips['new_bill'] = tips['total_bill'] / 2.0
   tips.head()

.. ipython:: python
   :suppress:

   tips = tips.drop('new_bill', axis=1)

Filtering
~~~~~~~~~

Filtering in SAS is done with an ``if`` or ``where`` statement, on one
or more columns.

.. code-block:: none

   data tips;
       set tips;
       if total_bill > 10;
   run;

   data tips;
       set tips;
       where total_bill > 10;
       /* equivalent in this case - where happens before the
          DATA step begins and can also be used in PROC statements */
   run;

DataFrames can be filtered in multiple ways; the most intuitive of which is using
:ref:`boolean indexing <indexing.boolean>`

.. ipython:: python

   tips[tips['total_bill'] > 10].head()

If/Then Logic
~~~~~~~~~~~~~

In SAS, if/then logic can be used to create new columns.

.. code-block:: none

   data tips;
       set tips;
       format bucket $4.;

       if total_bill < 10 then bucket = 'low';
       else bucket = 'high';
   run;

The same operation in pandas can be accomplished using
the ``where`` method from ``numpy``.

.. ipython:: python

   tips['bucket'] = np.where(tips['total_bill'] < 10, 'low', 'high')
   tips.head()

.. ipython:: python
   :suppress:

   tips = tips.drop('bucket', axis=1)

Date Functionality
~~~~~~~~~~~~~~~~~~

SAS provides a variety of functions to do operations on
date/datetime columns.

.. code-block:: none

   data tips;
       set tips;
       format date1 date2 date1_plusmonth mmddyy10.;
       date1 = mdy(1, 15, 2013);
       date2 = mdy(2, 15, 2015);
       date1_year = year(date1);
       date2_month = month(date2);
       * shift date to beginning of next interval;
       date1_next = intnx('MONTH', date1, 1);
       * count intervals between dates;
       months_between = intck('MONTH', date1, date2);
   run;

The equivalent pandas operations are shown below.  In addition to these
functions pandas supports other Time Series features
not available in Base SAS (such as resampling and and custom offsets) -
see the :ref:`timeseries documentation<timeseries>` for more details.

.. ipython:: python

   tips['date1'] = pd.Timestamp('2013-01-15')
   tips['date2'] = pd.Timestamp('2015-02-15')
   tips['date1_year'] = tips['date1'].dt.year
   tips['date2_month'] = tips['date2'].dt.month
   tips['date1_next'] = tips['date1'] + pd.offsets.MonthBegin()
   tips['months_between'] = (tips['date2'].dt.to_period('M') -
                             tips['date1'].dt.to_period('M'))

   tips[['date1','date2','date1_year','date2_month',
         'date1_next','months_between']].head()

.. ipython:: python
   :suppress:

   tips = tips.drop(['date1','date2','date1_year',
      'date2_month','date1_next','months_between'], axis=1)

Selection of Columns
~~~~~~~~~~~~~~~~~~~~

SAS provides keywords in the ``DATA`` step to select,
drop, and rename columns.

.. code-block:: none

   data tips;
       set tips;
       keep sex total_bill tip;
   run;

   data tips;
       set tips;
       drop sex;
   run;

   data tips;
       set tips;
       rename total_bill=total_bill_2;
   run;

The same operations are expressed in pandas below.

.. ipython:: python

   # keep
   tips[['sex', 'total_bill', 'tip']].head()

   # drop
   tips.drop('sex', axis=1).head()

   # rename
   tips.rename(columns={'total_bill':'total_bill_2'}).head()


Sorting by Values
~~~~~~~~~~~~~~~~~

Sorting in SAS is accomplished via ``PROC SORT``

.. code-block:: none

   proc sort data=tips;
       by sex total_bill;
   run;

pandas objects have a :meth:`~DataFrame.sort_values` method, which
takes a list of columns to sort by.

.. ipython:: python

   tips = tips.sort_values(['sex', 'total_bill'])
   tips.head()

Merging
-------

The following tables will be used in the merge examples

.. ipython:: python

   df1 = pd.DataFrame({'key': ['A', 'B', 'C', 'D'],
                       'value': np.random.randn(4)})
   df1
   df2 = pd.DataFrame({'key': ['B', 'D', 'D', 'E'],
                        'value': np.random.randn(4)})
   df2

In SAS, data must be explicitly sorted before merging.  Different
types of joins are accomplished using the ``in=`` dummy
variables to track whether a match was found in one or both
input frames.

.. code-block:: none

   proc sort data=df1;
       by key;
   run;

   proc sort data=df2;
       by key;
   run;

   data left_join inner_join right_join outer_join;
       merge df1(in=a) df2(in=b);

       if a and b then output inner_join;
       if a then output left_join;
       if b then output right_join;
       if a or b then output outer_join;
   run;

pandas DataFrames have a :meth:`~DataFrame.merge` method, which provides
similar functionality.  Note that the data does not have
to be sorted ahead of time, and different join
types are accomplished via the ``how`` keyword.

.. ipython:: python

   inner_join = df1.merge(df2, on=['key'], how='inner')
   inner_join

   left_join = df1.merge(df2, on=['key'], how='left')
   left_join

   right_join = df1.merge(df2, on=['key'], how='right')
   right_join

   outer_join = df1.merge(df2, on=['key'], how='outer')
   outer_join


Missing Data
------------

Like SAS, pandas has a representation for missing data - which is the
special float value ``NaN`` (not a number).  Many of the semantics
are the same, for example missing data propagates through numeric
operations, and is ignored by default for aggregations.

.. ipython:: python

   outer_join
   outer_join['value_x'] + outer_join['value_y']
   outer_join['value_x'].sum()

One difference is that missing data cannot be compared to its sentinel value.
For example, in SAS you could do this to filter missing values.

.. code-block:: none

   data outer_join_nulls;
       set outer_join;
       if value_x = .;
   run;

   data outer_join_no_nulls;
       set outer_join;
       if value_x ^= .;
   run;

Which doesn't work in in pandas.  Instead, the ``pd.isnull`` or ``pd.notnull`` functions
should be used for comparisons.

.. ipython:: python

   outer_join[pd.isnull(outer_join['value_x'])]
   outer_join[pd.notnull(outer_join['value_x'])]

pandas also provides a variety of methods to work with missing data - some of
which would be challenging to express in SAS. For example, there are methods to
drop all rows with any missing values, replacing missing values with a specified
value, like the mean, or forward filling from previous rows. See the
:ref:`missing data documentation<missing_data>` for more.

.. ipython:: python

   outer_join.dropna()
   outer_join.fillna(method='ffill')
   outer_join['value_x'].fillna(outer_join['value_x'].mean())


GroupBy
-------

Aggregation
~~~~~~~~~~~

SAS's PROC SUMMARY can be used to group by one or
more key variables and compute aggregations on
numeric columns.

.. code-block:: none

   proc summary data=tips nway;
       class sex smoker;
       var total_bill tip;
       output out=tips_summed sum=;
   run;

pandas provides a flexible ``groupby`` mechanism that
allows similar aggregations.  See the :ref:`groupby documentation<groupby>`
for more details and examples.

.. ipython:: python

   tips_summed = tips.groupby(['sex', 'smoker'])['total_bill', 'tip'].sum()
   tips_summed.head()


Transformation
~~~~~~~~~~~~~~

In SAS, if the group aggregations need to be used with
the original frame, it must be merged back together.  For
example, to subtract the mean for each observation by smoker group.

.. code-block:: none

   proc summary data=tips missing nway;
       class smoker;
       var total_bill;
       output out=smoker_means mean(total_bill)=group_bill;
   run;

   proc sort data=tips;
       by smoker;
   run;

   data tips;
       merge tips(in=a) smoker_means(in=b);
       by smoker;
       adj_total_bill = total_bill - group_bill;
       if a and b;
   run;


pandas ``groubpy`` provides a ``transform`` mechanism that allows
these type of operations to be succinctly expressed in one
operation.

.. ipython:: python

   gb = tips.groupby('smoker')['total_bill']
   tips['adj_total_bill'] = tips['total_bill'] - gb.transform('mean')
   tips.head()


By Group Processing
~~~~~~~~~~~~~~~~~~~

In addition to aggregation, pandas ``groupby`` can be used to
replicate most other by group processing from SAS. For example,
this ``DATA`` step reads the data by sex/smoker group and filters to
the first entry for each.

.. code-block:: none

   proc sort data=tips;
      by sex smoker;
   run;

   data tips_first;
       set tips;
       by sex smoker;
       if FIRST.sex or FIRST.smoker then output;
   run;

In pandas this would be written as:

.. ipython:: python

   tips.groupby(['sex','smoker']).first()


Other Considerations
--------------------

Disk vs Memory
~~~~~~~~~~~~~~

pandas operates exclusively in memory, where a SAS data set exists on disk.
This means that the size of data able to be loaded in pandas is limited by your
machine's memory, but also that the operations on that data may be faster.

If out of core processing is needed, one possibility is the
`dask.dataframe <http://dask.pydata.org/en/latest/dataframe.html>`_
library (currently in development) which 
provides a subset of pandas functionality for an on-disk ``DataFrame``

Data Interop
~~~~~~~~~~~~

pandas provides a :func:`read_sas` method that can read SAS data saved in
the XPORT format.  The ability to read SAS's binary format is planned for a
future release.

.. code-block:: none

   libname xportout xport 'transport-file.xpt';
   data xportout.tips;
       set tips(rename=(total_bill=tbill));
       * xport variable names limited to 6 characters;
   run;

.. code-block:: python

   df = pd.read_sas('transport-file.xpt')

XPORT is a relatively limited format and the parsing of it is not as
optimized as some of the other pandas readers. An alternative way
to interop data between SAS and pandas is to serialize to csv.

.. code-block:: python

   # version 0.17, 10M rows

   In [8]: %time df = pd.read_sas('big.xpt')
   Wall time: 14.6 s

   In [9]: %time df = pd.read_csv('big.csv')
   Wall time: 4.86 s
