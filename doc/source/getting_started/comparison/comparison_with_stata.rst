.. _compare_with_stata:

{{ header }}

Comparison with Stata
*********************
For potential users coming from `Stata <https://en.wikipedia.org/wiki/Stata>`__
this page is meant to demonstrate how different Stata operations would be
performed in pandas.

If you're new to pandas, you might want to first read through :ref:`10 Minutes to pandas<10min>`
to familiarize yourself with the library.

As is customary, we import pandas and NumPy as follows. This means that we can refer to the
libraries as ``pd`` and ``np``, respectively, for the rest of the document.

.. ipython:: python

    import pandas as pd
    import numpy as np


.. note::

   Throughout this tutorial, the pandas ``DataFrame`` will be displayed by calling
   ``df.head()``, which displays the first N (default 5) rows of the ``DataFrame``.
   This is often used in interactive work (e.g. `Jupyter notebook
   <https://jupyter.org/>`_ or terminal) -- the equivalent in Stata would be:

   .. code-block:: stata

      list in 1/5

Data structures
---------------

General terminology translation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "pandas", "Stata"
    :widths: 20, 20

    ``DataFrame``, data set
    column, variable
    row, observation
    groupby, bysort
    ``NaN``, ``.``


``DataFrame`` / ``Series``
~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``DataFrame`` in pandas is analogous to a Stata data set -- a two-dimensional
data source with labeled columns that can be of different types. As will be
shown in this document, almost any operation that can be applied to a data set
in Stata can also be accomplished in pandas.

A ``Series`` is the data structure that represents one column of a
``DataFrame``. Stata doesn't have a separate data structure for a single column,
but in general, working with a ``Series`` is analogous to referencing a column
of a data set in Stata.

``Index``
~~~~~~~~~

Every ``DataFrame`` and ``Series`` has an ``Index`` -- labels on the
*rows* of the data. Stata does not have an exactly analogous concept. In Stata, a data set's
rows are essentially unlabeled, other than an implicit integer index that can be
accessed with ``_n``.

In pandas, if no index is specified, an integer index is also used by default
(first row = 0, second row = 1, and so on). While using a labeled ``Index`` or
``MultiIndex`` can enable sophisticated analyses and is ultimately an important
part of pandas to understand, for this comparison we will essentially ignore the
``Index`` and just treat the ``DataFrame`` as a collection of columns. Please
see the :ref:`indexing documentation<indexing>` for much more on how to use an
``Index`` effectively.


Data input / output
-------------------

Constructing a DataFrame from values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Stata data set can be built from specified values by
placing the data after an ``input`` statement and
specifying the column names.

.. code-block:: stata

   input x y
   1 2
   3 4
   5 6
   end

A pandas ``DataFrame`` can be constructed in many different ways,
but for a small number of values, it is often convenient to specify it as
a Python dictionary, where the keys are the column names
and the values are the data.

.. ipython:: python

   df = pd.DataFrame({'x': [1, 3, 5], 'y': [2, 4, 6]})
   df


Reading external data
~~~~~~~~~~~~~~~~~~~~~

Like Stata, pandas provides utilities for reading in data from
many formats.  The ``tips`` data set, found within the pandas
tests (`csv <https://raw.github.com/pandas-dev/pandas/master/pandas/tests/data/tips.csv>`_)
will be used in many of the following examples.

Stata provides ``import delimited`` to read csv data into a data set in memory.
If the ``tips.csv`` file is in the current working directory, we can import it as follows.

.. code-block:: stata

   import delimited tips.csv

The pandas method is :func:`read_csv`, which works similarly. Additionally, it will automatically download
the data set if presented with a url.

.. ipython:: python

   url = ('https://raw.github.com/pandas-dev'
          '/pandas/master/pandas/tests/data/tips.csv')
   tips = pd.read_csv(url)
   tips.head()

Like ``import delimited``, :func:`read_csv` can take a number of parameters to specify
how the data should be parsed.  For example, if the data were instead tab delimited,
did not have column names, and existed in the current working directory,
the pandas command would be:

.. code-block:: python

   tips = pd.read_csv('tips.csv', sep='\t', header=None)

   # alternatively, read_table is an alias to read_csv with tab delimiter
   tips = pd.read_table('tips.csv', header=None)

Pandas can also read Stata data sets in ``.dta`` format with the :func:`read_stata` function.

.. code-block:: python

   df = pd.read_stata('data.dta')

In addition to text/csv and Stata files, pandas supports a variety of other data formats
such as Excel, SAS, HDF5, Parquet, and SQL databases.  These are all read via a ``pd.read_*``
function.  See the :ref:`IO documentation<io>` for more details.


Exporting data
~~~~~~~~~~~~~~

The inverse of ``import delimited`` in Stata is ``export delimited``

.. code-block:: stata

   export delimited tips2.csv

Similarly in pandas, the opposite of ``read_csv`` is :meth:`DataFrame.to_csv`.

.. code-block:: python

   tips.to_csv('tips2.csv')

Pandas can also export to Stata file format with the :meth:`DataFrame.to_stata` method.

.. code-block:: python

   tips.to_stata('tips2.dta')


Data operations
---------------

Operations on columns
~~~~~~~~~~~~~~~~~~~~~

In Stata, arbitrary math expressions can be used with the ``generate`` and
``replace`` commands on new or existing columns. The ``drop`` command drops
the column from the data set.

.. code-block:: stata

   replace total_bill = total_bill - 2
   generate new_bill = total_bill / 2
   drop new_bill

pandas provides similar vectorized operations by
specifying the individual ``Series`` in the ``DataFrame``.
New columns can be assigned in the same way. The :meth:`DataFrame.drop` method
drops a column from the ``DataFrame``.

.. ipython:: python

   tips['total_bill'] = tips['total_bill'] - 2
   tips['new_bill'] = tips['total_bill'] / 2
   tips.head()

   tips = tips.drop('new_bill', axis=1)

Filtering
~~~~~~~~~

Filtering in Stata is done with an ``if`` clause on one or more columns.

.. code-block:: stata

   list if total_bill > 10

DataFrames can be filtered in multiple ways; the most intuitive of which is using
:ref:`boolean indexing <indexing.boolean>`.

.. ipython:: python

   tips[tips['total_bill'] > 10].head()

If/then logic
~~~~~~~~~~~~~

In Stata, an ``if`` clause can also be used to create new columns.

.. code-block:: stata

   generate bucket = "low" if total_bill < 10
   replace bucket = "high" if total_bill >= 10

The same operation in pandas can be accomplished using
the ``where`` method from ``numpy``.

.. ipython:: python

   tips['bucket'] = np.where(tips['total_bill'] < 10, 'low', 'high')
   tips.head()

.. ipython:: python
   :suppress:

   tips = tips.drop('bucket', axis=1)

Date functionality
~~~~~~~~~~~~~~~~~~

Stata provides a variety of functions to do operations on
date/datetime columns.

.. code-block:: stata

   generate date1 = mdy(1, 15, 2013)
   generate date2 = date("Feb152015", "MDY")

   generate date1_year = year(date1)
   generate date2_month = month(date2)

   * shift date to beginning of next month
   generate date1_next = mdy(month(date1) + 1, 1, year(date1)) if month(date1) != 12
   replace date1_next = mdy(1, 1, year(date1) + 1) if month(date1) == 12
   generate months_between = mofd(date2) - mofd(date1)

   list date1 date2 date1_year date2_month date1_next months_between

The equivalent pandas operations are shown below.  In addition to these
functions, pandas supports other Time Series features
not available in Stata (such as time zone handling and custom offsets) --
see the :ref:`timeseries documentation<timeseries>` for more details.

.. ipython:: python

   tips['date1'] = pd.Timestamp('2013-01-15')
   tips['date2'] = pd.Timestamp('2015-02-15')
   tips['date1_year'] = tips['date1'].dt.year
   tips['date2_month'] = tips['date2'].dt.month
   tips['date1_next'] = tips['date1'] + pd.offsets.MonthBegin()
   tips['months_between'] = (tips['date2'].dt.to_period('M')
                             - tips['date1'].dt.to_period('M'))

   tips[['date1', 'date2', 'date1_year', 'date2_month', 'date1_next',
         'months_between']].head()

.. ipython:: python
   :suppress:

   tips = tips.drop(['date1', 'date2', 'date1_year', 'date2_month',
                     'date1_next', 'months_between'], axis=1)

Selection of columns
~~~~~~~~~~~~~~~~~~~~

Stata provides keywords to select, drop, and rename columns.

.. code-block:: stata

   keep sex total_bill tip

   drop sex

   rename total_bill total_bill_2

The same operations are expressed in pandas below. Note that in contrast to Stata, these
operations do not happen in place. To make these changes persist, assign the operation back
to a variable.

.. ipython:: python

   # keep
   tips[['sex', 'total_bill', 'tip']].head()

   # drop
   tips.drop('sex', axis=1).head()

   # rename
   tips.rename(columns={'total_bill': 'total_bill_2'}).head()


Sorting by values
~~~~~~~~~~~~~~~~~

Sorting in Stata is accomplished via ``sort``

.. code-block:: stata

   sort sex total_bill

pandas objects have a :meth:`DataFrame.sort_values` method, which
takes a list of columns to sort by.

.. ipython:: python

   tips = tips.sort_values(['sex', 'total_bill'])
   tips.head()


String processing
-----------------

Finding length of string
~~~~~~~~~~~~~~~~~~~~~~~~

Stata determines the length of a character string with the :func:`strlen` and
:func:`ustrlen` functions for ASCII and Unicode strings, respectively.

.. code-block:: stata

   generate strlen_time = strlen(time)
   generate ustrlen_time = ustrlen(time)

Python determines the length of a character string with the ``len`` function.
In Python 3, all strings are Unicode strings. ``len`` includes trailing blanks.
Use ``len`` and ``rstrip`` to exclude trailing blanks.

.. ipython:: python

   tips['time'].str.len().head()
   tips['time'].str.rstrip().str.len().head()


Finding position of substring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stata determines the position of a character in a string with the :func:`strpos` function.
This takes the string defined by the first argument and searches for the
first position of the substring you supply as the second argument.

.. code-block:: stata

   generate str_position = strpos(sex, "ale")

Python determines the position of a character in a string with the
:func:`find` function.  ``find`` searches for the first position of the
substring.  If the substring is found, the function returns its
position.  Keep in mind that Python indexes are zero-based and
the function will return -1 if it fails to find the substring.

.. ipython:: python

   tips['sex'].str.find("ale").head()


Extracting substring by position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stata extracts a substring from a string based on its position with the :func:`substr` function.

.. code-block:: stata

   generate short_sex = substr(sex, 1, 1)

With pandas you can use ``[]`` notation to extract a substring
from a string by position locations.  Keep in mind that Python
indexes are zero-based.

.. ipython:: python

   tips['sex'].str[0:1].head()


Extracting nth word
~~~~~~~~~~~~~~~~~~~

The Stata :func:`word` function returns the nth word from a string.
The first argument is the string you want to parse and the
second argument specifies which word you want to extract.

.. code-block:: stata

   clear
   input str20 string
   "John Smith"
   "Jane Cook"
   end

   generate first_name = word(name, 1)
   generate last_name = word(name, -1)

Python extracts a substring from a string based on its text
by using regular expressions. There are much more powerful
approaches, but this just shows a simple approach.

.. ipython:: python

   firstlast = pd.DataFrame({'string': ['John Smith', 'Jane Cook']})
   firstlast['First_Name'] = firstlast['string'].str.split(" ", expand=True)[0]
   firstlast['Last_Name'] = firstlast['string'].str.rsplit(" ", expand=True)[0]
   firstlast


Changing case
~~~~~~~~~~~~~

The Stata :func:`strupper`, :func:`strlower`, :func:`strproper`,
:func:`ustrupper`, :func:`ustrlower`, and :func:`ustrtitle` functions
change the case of ASCII and Unicode strings, respectively.

.. code-block:: stata

   clear
   input str20 string
   "John Smith"
   "Jane Cook"
   end

   generate upper = strupper(string)
   generate lower = strlower(string)
   generate title = strproper(string)
   list

The equivalent Python functions are ``upper``, ``lower``, and ``title``.

.. ipython:: python

   firstlast = pd.DataFrame({'string': ['John Smith', 'Jane Cook']})
   firstlast['upper'] = firstlast['string'].str.upper()
   firstlast['lower'] = firstlast['string'].str.lower()
   firstlast['title'] = firstlast['string'].str.title()
   firstlast

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

In Stata, to perform a merge, one data set must be in memory
and the other must be referenced as a file name on disk. In
contrast, Python must have both ``DataFrames`` already in memory.

By default, Stata performs an outer join, where all observations
from both data sets are left in memory after the merge. One can
keep only observations from the initial data set, the merged data set,
or the intersection of the two by using the values created in the
``_merge`` variable.

.. code-block:: stata

   * First create df2 and save to disk
   clear
   input str1 key
   B
   D
   D
   E
   end
   generate value = rnormal()
   save df2.dta

   * Now create df1 in memory
   clear
   input str1 key
   A
   B
   C
   D
   end
   generate value = rnormal()

   preserve

   * Left join
   merge 1:n key using df2.dta
   keep if _merge == 1

   * Right join
   restore, preserve
   merge 1:n key using df2.dta
   keep if _merge == 2

   * Inner join
   restore, preserve
   merge 1:n key using df2.dta
   keep if _merge == 3

   * Outer join
   restore
   merge 1:n key using df2.dta

pandas DataFrames have a :meth:`DataFrame.merge` method, which provides
similar functionality. Note that different join
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


Missing data
------------

Like Stata, pandas has a representation for missing data -- the
special float value ``NaN`` (not a number).  Many of the semantics
are the same; for example missing data propagates through numeric
operations, and is ignored by default for aggregations.

.. ipython:: python

   outer_join
   outer_join['value_x'] + outer_join['value_y']
   outer_join['value_x'].sum()

One difference is that missing data cannot be compared to its sentinel value.
For example, in Stata you could do this to filter missing values.

.. code-block:: stata

   * Keep missing values
   list if value_x == .
   * Keep non-missing values
   list if value_x != .

This doesn't work in pandas.  Instead, the :func:`pd.isna` or :func:`pd.notna` functions
should be used for comparisons.

.. ipython:: python

   outer_join[pd.isna(outer_join['value_x'])]
   outer_join[pd.notna(outer_join['value_x'])]

Pandas also provides a variety of methods to work with missing data -- some of
which would be challenging to express in Stata. For example, there are methods to
drop all rows with any missing values, replacing missing values with a specified
value, like the mean, or forward filling from previous rows. See the
:ref:`missing data documentation<missing_data>` for more.

.. ipython:: python

   # Drop rows with any missing value
   outer_join.dropna()

   # Fill forwards
   outer_join.fillna(method='ffill')

   # Impute missing values with the mean
   outer_join['value_x'].fillna(outer_join['value_x'].mean())


GroupBy
-------

Aggregation
~~~~~~~~~~~

Stata's ``collapse`` can be used to group by one or
more key variables and compute aggregations on
numeric columns.

.. code-block:: stata

   collapse (sum) total_bill tip, by(sex smoker)

pandas provides a flexible ``groupby`` mechanism that
allows similar aggregations.  See the :ref:`groupby documentation<groupby>`
for more details and examples.

.. ipython:: python

   tips_summed = tips.groupby(['sex', 'smoker'])[['total_bill', 'tip']].sum()
   tips_summed.head()


Transformation
~~~~~~~~~~~~~~

In Stata, if the group aggregations need to be used with the
original data set, one would usually use ``bysort`` with :func:`egen`.
For example, to subtract the mean for each observation by smoker group.

.. code-block:: stata

   bysort sex smoker: egen group_bill = mean(total_bill)
   generate adj_total_bill = total_bill - group_bill


pandas ``groupby`` provides a ``transform`` mechanism that allows
these type of operations to be succinctly expressed in one
operation.

.. ipython:: python

   gb = tips.groupby('smoker')['total_bill']
   tips['adj_total_bill'] = tips['total_bill'] - gb.transform('mean')
   tips.head()


By group processing
~~~~~~~~~~~~~~~~~~~

In addition to aggregation, pandas ``groupby`` can be used to
replicate most other ``bysort`` processing from Stata. For example,
the following example lists the first observation in the current
sort order by sex/smoker group.

.. code-block:: stata

   bysort sex smoker: list if _n == 1

In pandas this would be written as:

.. ipython:: python

   tips.groupby(['sex', 'smoker']).first()


Other considerations
--------------------

Disk vs memory
~~~~~~~~~~~~~~

Pandas and Stata both operate exclusively in memory. This means that the size of
data able to be loaded in pandas is limited by your machine's memory.
If out of core processing is needed, one possibility is the
`dask.dataframe <https://dask.pydata.org/en/latest/dataframe.html>`_
library, which provides a subset of pandas functionality for an
on-disk ``DataFrame``.
