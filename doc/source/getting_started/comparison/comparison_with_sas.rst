.. _compare_with_sas:

{{ header }}

Comparison with SAS
********************

For potential users coming from `SAS <https://en.wikipedia.org/wiki/SAS_(software)>`__
this page is meant to demonstrate how different SAS operations would be
performed in pandas.

.. include:: includes/introduction.rst


Data structures
---------------

General terminology translation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "pandas", "SAS"
    :widths: 20, 20

    ``DataFrame``, data set
    column, variable
    row, observation
    groupby, BY-group
    ``NaN``, ``.``


``DataFrame``
~~~~~~~~~~~~~

A ``DataFrame`` in pandas is analogous to a SAS data set - a two-dimensional
data source with labeled columns that can be of different types. As will be
shown in this document, almost any operation that can be applied to a data set
using SAS's ``DATA`` step, can also be accomplished in pandas.

``Series``
~~~~~~~~~~

A ``Series`` is the data structure that represents one column of a
``DataFrame``. SAS doesn't have a separate data structure for a single column,
but in general, working with a ``Series`` is analogous to referencing a column
in the ``DATA`` step.

``Index``
~~~~~~~~~

Every ``DataFrame`` and ``Series`` has an ``Index`` - which are labels on the
*rows* of the data. SAS does not have an exactly analogous concept. A data set's
rows are essentially unlabeled, other than an implicit integer index that can be
accessed during the ``DATA`` step (``_N_``).

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

A SAS data set can be built from specified values by
placing the data after a ``datalines`` statement and
specifying the column names.

.. code-block:: sas

   data df;
       input x y;
       datalines;
       1 2
       3 4
       5 6
       ;
   run;

.. include:: includes/construct_dataframe.rst

Reading external data
~~~~~~~~~~~~~~~~~~~~~

Like SAS, pandas provides utilities for reading in data from
many formats.  The ``tips`` dataset, found within the pandas
tests (`csv <https://raw.github.com/pandas-dev/pandas/master/pandas/tests/io/data/csv/tips.csv>`_)
will be used in many of the following examples.

SAS provides ``PROC IMPORT`` to read csv data into a data set.

.. code-block:: sas

   proc import datafile='tips.csv' dbms=csv out=tips replace;
       getnames=yes;
   run;

The pandas method is :func:`read_csv`, which works similarly.

.. ipython:: python

   url = (
       "https://raw.github.com/pandas-dev/"
       "pandas/master/pandas/tests/io/data/csv/tips.csv"
   )
   tips = pd.read_csv(url)
   tips


Like ``PROC IMPORT``, ``read_csv`` can take a number of parameters to specify
how the data should be parsed.  For example, if the data was instead tab delimited,
and did not have column names, the pandas command would be:

.. code-block:: python

   tips = pd.read_csv("tips.csv", sep="\t", header=None)

   # alternatively, read_table is an alias to read_csv with tab delimiter
   tips = pd.read_table("tips.csv", header=None)

In addition to text/csv, pandas supports a variety of other data formats
such as Excel, HDF5, and SQL databases.  These are all read via a ``pd.read_*``
function.  See the :ref:`IO documentation<io>` for more details.

Limiting output
~~~~~~~~~~~~~~~

.. include:: includes/limit.rst

The equivalent in SAS would be:

.. code-block:: sas

   proc print data=df(obs=5);
   run;


Exporting data
~~~~~~~~~~~~~~

The inverse of ``PROC IMPORT`` in SAS is ``PROC EXPORT``

.. code-block:: sas

   proc export data=tips outfile='tips2.csv' dbms=csv;
   run;

Similarly in pandas, the opposite of ``read_csv`` is :meth:`~DataFrame.to_csv`,
and other data formats follow a similar api.

.. code-block:: python

   tips.to_csv("tips2.csv")


Data operations
---------------

Operations on columns
~~~~~~~~~~~~~~~~~~~~~

In the ``DATA`` step, arbitrary math expressions can
be used on new or existing columns.

.. code-block:: sas

   data tips;
       set tips;
       total_bill = total_bill - 2;
       new_bill = total_bill / 2;
   run;

.. include:: includes/column_operations.rst


Filtering
~~~~~~~~~

Filtering in SAS is done with an ``if`` or ``where`` statement, on one
or more columns.

.. code-block:: sas

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

.. include:: includes/filtering.rst

If/then logic
~~~~~~~~~~~~~

In SAS, if/then logic can be used to create new columns.

.. code-block:: sas

   data tips;
       set tips;
       format bucket $4.;

       if total_bill < 10 then bucket = 'low';
       else bucket = 'high';
   run;

.. include:: includes/if_then.rst

Date functionality
~~~~~~~~~~~~~~~~~~

SAS provides a variety of functions to do operations on
date/datetime columns.

.. code-block:: sas

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
not available in Base SAS (such as resampling and custom offsets) -
see the :ref:`timeseries documentation<timeseries>` for more details.

.. include:: includes/time_date.rst

Selection of columns
~~~~~~~~~~~~~~~~~~~~

SAS provides keywords in the ``DATA`` step to select,
drop, and rename columns.

.. code-block:: sas

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

.. include:: includes/column_selection.rst


Sorting by values
~~~~~~~~~~~~~~~~~

Sorting in SAS is accomplished via ``PROC SORT``

.. code-block:: sas

   proc sort data=tips;
       by sex total_bill;
   run;

.. include:: includes/sorting.rst

String processing
-----------------

Finding length of string
~~~~~~~~~~~~~~~~~~~~~~~~

SAS determines the length of a character string with the
`LENGTHN <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/viewer.htm#a002284668.htm>`__
and `LENGTHC <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/viewer.htm#a002283942.htm>`__
functions. ``LENGTHN`` excludes trailing blanks and ``LENGTHC`` includes trailing blanks.

.. code-block:: sas

   data _null_;
   set tips;
   put(LENGTHN(time));
   put(LENGTHC(time));
   run;

.. include:: includes/length.rst


Finding position of substring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SAS determines the position of a character in a string with the
`FINDW <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/viewer.htm#a002978282.htm>`__ function.
``FINDW`` takes the string defined by the first argument and searches for the first position of the substring
you supply as the second argument.

.. code-block:: sas

   data _null_;
   set tips;
   put(FINDW(sex,'ale'));
   run;

.. include:: includes/find_substring.rst


Extracting substring by position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SAS extracts a substring from a string based on its position with the
`SUBSTR <https://www2.sas.com/proceedings/sugi25/25/cc/25p088.pdf>`__ function.

.. code-block:: sas

   data _null_;
   set tips;
   put(substr(sex,1,1));
   run;

.. include:: includes/extract_substring.rst


Extracting nth word
~~~~~~~~~~~~~~~~~~~

The SAS `SCAN <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/viewer.htm#a000214639.htm>`__
function returns the nth word from a string. The first argument is the string you want to parse and the
second argument specifies which word you want to extract.

.. code-block:: sas

   data firstlast;
   input String $60.;
   First_Name = scan(string, 1);
   Last_Name = scan(string, -1);
   datalines2;
   John Smith;
   Jane Cook;
   ;;;
   run;

.. include:: includes/nth_word.rst


Changing case
~~~~~~~~~~~~~

The SAS `UPCASE <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/viewer.htm#a000245965.htm>`__
`LOWCASE <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/viewer.htm#a000245912.htm>`__ and
`PROPCASE <https://support.sas.com/documentation/cdl/en/lrdict/64316/HTML/default/a002598106.htm>`__
functions change the case of the argument.

.. code-block:: sas

   data firstlast;
   input String $60.;
   string_up = UPCASE(string);
   string_low = LOWCASE(string);
   string_prop = PROPCASE(string);
   datalines2;
   John Smith;
   Jane Cook;
   ;;;
   run;

.. include:: includes/case.rst


Merging
-------

.. include:: includes/merge_setup.rst

In SAS, data must be explicitly sorted before merging.  Different
types of joins are accomplished using the ``in=`` dummy
variables to track whether a match was found in one or both
input frames.

.. code-block:: sas

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

.. include:: includes/merge.rst


Missing data
------------

Both pandas and SAS have a representation for missing data.

.. include:: includes/missing_intro.rst

One difference is that missing data cannot be compared to its sentinel value.
For example, in SAS you could do this to filter missing values.

.. code-block:: sas

   data outer_join_nulls;
       set outer_join;
       if value_x = .;
   run;

   data outer_join_no_nulls;
       set outer_join;
       if value_x ^= .;
   run;

.. include:: includes/missing.rst


GroupBy
-------

Aggregation
~~~~~~~~~~~

SAS's ``PROC SUMMARY`` can be used to group by one or
more key variables and compute aggregations on
numeric columns.

.. code-block:: sas

   proc summary data=tips nway;
       class sex smoker;
       var total_bill tip;
       output out=tips_summed sum=;
   run;

.. include:: includes/groupby.rst


Transformation
~~~~~~~~~~~~~~

In SAS, if the group aggregations need to be used with
the original frame, it must be merged back together.  For
example, to subtract the mean for each observation by smoker group.

.. code-block:: sas

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

.. include:: includes/transform.rst


By group processing
~~~~~~~~~~~~~~~~~~~

In addition to aggregation, pandas ``groupby`` can be used to
replicate most other by group processing from SAS. For example,
this ``DATA`` step reads the data by sex/smoker group and filters to
the first entry for each.

.. code-block:: sas

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

   tips.groupby(["sex", "smoker"]).first()


Other considerations
--------------------

Disk vs memory
~~~~~~~~~~~~~~

pandas operates exclusively in memory, where a SAS data set exists on disk.
This means that the size of data able to be loaded in pandas is limited by your
machine's memory, but also that the operations on that data may be faster.

If out of core processing is needed, one possibility is the
`dask.dataframe <https://dask.pydata.org/en/latest/dataframe.html>`_
library (currently in development) which
provides a subset of pandas functionality for an on-disk ``DataFrame``

Data interop
~~~~~~~~~~~~

pandas provides a :func:`read_sas` method that can read SAS data saved in
the XPORT or SAS7BDAT binary format.

.. code-block:: sas

   libname xportout xport 'transport-file.xpt';
   data xportout.tips;
       set tips(rename=(total_bill=tbill));
       * xport variable names limited to 6 characters;
   run;

.. code-block:: python

   df = pd.read_sas("transport-file.xpt")
   df = pd.read_sas("binary-file.sas7bdat")

You can also specify the file format directly. By default, pandas will try
to infer the file format based on its extension.

.. code-block:: python

   df = pd.read_sas("transport-file.xpt", format="xport")
   df = pd.read_sas("binary-file.sas7bdat", format="sas7bdat")

XPORT is a relatively limited format and the parsing of it is not as
optimized as some of the other pandas readers. An alternative way
to interop data between SAS and pandas is to serialize to csv.

.. code-block:: ipython

   # version 0.17, 10M rows

   In [8]: %time df = pd.read_sas('big.xpt')
   Wall time: 14.6 s

   In [9]: %time df = pd.read_csv('big.csv')
   Wall time: 4.86 s
