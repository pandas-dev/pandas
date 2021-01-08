.. _compare_with_stata:

{{ header }}

Comparison with Stata
*********************
For potential users coming from `Stata <https://en.wikipedia.org/wiki/Stata>`__
this page is meant to demonstrate how different Stata operations would be
performed in pandas.

.. include:: includes/introduction.rst


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


``DataFrame``
~~~~~~~~~~~~~

A ``DataFrame`` in pandas is analogous to a Stata data set -- a two-dimensional
data source with labeled columns that can be of different types. As will be
shown in this document, almost any operation that can be applied to a data set
in Stata can also be accomplished in pandas.

``Series``
~~~~~~~~~~

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


Copies vs. in place operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: includes/copies.rst


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

.. include:: includes/construct_dataframe.rst

Reading external data
~~~~~~~~~~~~~~~~~~~~~

Like Stata, pandas provides utilities for reading in data from
many formats.  The ``tips`` data set, found within the pandas
tests (`csv <https://raw.github.com/pandas-dev/pandas/master/pandas/tests/io/data/csv/tips.csv>`_)
will be used in many of the following examples.

Stata provides ``import delimited`` to read csv data into a data set in memory.
If the ``tips.csv`` file is in the current working directory, we can import it as follows.

.. code-block:: stata

   import delimited tips.csv

The pandas method is :func:`read_csv`, which works similarly. Additionally, it will automatically download
the data set if presented with a url.

.. ipython:: python

   url = (
       "https://raw.github.com/pandas-dev"
       "/pandas/master/pandas/tests/io/data/csv/tips.csv"
   )
   tips = pd.read_csv(url)
   tips

Like ``import delimited``, :func:`read_csv` can take a number of parameters to specify
how the data should be parsed.  For example, if the data were instead tab delimited,
did not have column names, and existed in the current working directory,
the pandas command would be:

.. code-block:: python

   tips = pd.read_csv("tips.csv", sep="\t", header=None)

   # alternatively, read_table is an alias to read_csv with tab delimiter
   tips = pd.read_table("tips.csv", header=None)

pandas can also read Stata data sets in ``.dta`` format with the :func:`read_stata` function.

.. code-block:: python

   df = pd.read_stata("data.dta")

In addition to text/csv and Stata files, pandas supports a variety of other data formats
such as Excel, SAS, HDF5, Parquet, and SQL databases.  These are all read via a ``pd.read_*``
function.  See the :ref:`IO documentation<io>` for more details.


Limiting output
~~~~~~~~~~~~~~~

.. include:: includes/limit.rst

The equivalent in Stata would be:

.. code-block:: stata

   list in 1/5


Exporting data
~~~~~~~~~~~~~~

The inverse of ``import delimited`` in Stata is ``export delimited``

.. code-block:: stata

   export delimited tips2.csv

Similarly in pandas, the opposite of ``read_csv`` is :meth:`DataFrame.to_csv`.

.. code-block:: python

   tips.to_csv("tips2.csv")

pandas can also export to Stata file format with the :meth:`DataFrame.to_stata` method.

.. code-block:: python

   tips.to_stata("tips2.dta")


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

.. include:: includes/column_operations.rst


Filtering
~~~~~~~~~

Filtering in Stata is done with an ``if`` clause on one or more columns.

.. code-block:: stata

   list if total_bill > 10

.. include:: includes/filtering.rst

If/then logic
~~~~~~~~~~~~~

In Stata, an ``if`` clause can also be used to create new columns.

.. code-block:: stata

   generate bucket = "low" if total_bill < 10
   replace bucket = "high" if total_bill >= 10

.. include:: includes/if_then.rst

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

.. include:: includes/time_date.rst

Selection of columns
~~~~~~~~~~~~~~~~~~~~

Stata provides keywords to select, drop, and rename columns.

.. code-block:: stata

   keep sex total_bill tip

   drop sex

   rename total_bill total_bill_2

.. include:: includes/column_selection.rst


Sorting by values
~~~~~~~~~~~~~~~~~

Sorting in Stata is accomplished via ``sort``

.. code-block:: stata

   sort sex total_bill

.. include:: includes/sorting.rst

String processing
-----------------

Finding length of string
~~~~~~~~~~~~~~~~~~~~~~~~

Stata determines the length of a character string with the :func:`strlen` and
:func:`ustrlen` functions for ASCII and Unicode strings, respectively.

.. code-block:: stata

   generate strlen_time = strlen(time)
   generate ustrlen_time = ustrlen(time)

.. include:: includes/length.rst


Finding position of substring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stata determines the position of a character in a string with the :func:`strpos` function.
This takes the string defined by the first argument and searches for the
first position of the substring you supply as the second argument.

.. code-block:: stata

   generate str_position = strpos(sex, "ale")

.. include:: includes/find_substring.rst


Extracting substring by position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stata extracts a substring from a string based on its position with the :func:`substr` function.

.. code-block:: stata

   generate short_sex = substr(sex, 1, 1)

.. include:: includes/extract_substring.rst


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

.. include:: includes/nth_word.rst


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

.. include:: includes/case.rst


Merging
-------

.. include:: includes/merge_setup.rst

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

.. include:: includes/merge.rst


Missing data
------------

Both pandas and Stata have a representation for missing data.

.. include:: includes/missing_intro.rst

One difference is that missing data cannot be compared to its sentinel value.
For example, in Stata you could do this to filter missing values.

.. code-block:: stata

   * Keep missing values
   list if value_x == .
   * Keep non-missing values
   list if value_x != .

.. include:: includes/missing.rst


GroupBy
-------

Aggregation
~~~~~~~~~~~

Stata's ``collapse`` can be used to group by one or
more key variables and compute aggregations on
numeric columns.

.. code-block:: stata

   collapse (sum) total_bill tip, by(sex smoker)

.. include:: includes/groupby.rst


Transformation
~~~~~~~~~~~~~~

In Stata, if the group aggregations need to be used with the
original data set, one would usually use ``bysort`` with :func:`egen`.
For example, to subtract the mean for each observation by smoker group.

.. code-block:: stata

   bysort sex smoker: egen group_bill = mean(total_bill)
   generate adj_total_bill = total_bill - group_bill

.. include:: includes/transform.rst


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

   tips.groupby(["sex", "smoker"]).first()


Other considerations
--------------------

Disk vs memory
~~~~~~~~~~~~~~

pandas and Stata both operate exclusively in memory. This means that the size of
data able to be loaded in pandas is limited by your machine's memory.
If out of core processing is needed, one possibility is the
`dask.dataframe <https://dask.pydata.org/en/latest/dataframe.html>`_
library, which provides a subset of pandas functionality for an
on-disk ``DataFrame``.
