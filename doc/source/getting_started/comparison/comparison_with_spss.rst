.. _compare_with_spss:

{{ header }}

Comparison with SPSS
********************
For potential users coming from `SPSS <https://www.ibm.com/spss>`__, this page is meant to demonstrate
how various SPSS operations would be performed using pandas.

.. include:: includes/introduction.rst

Data structures
---------------

General terminology translation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "pandas", "SPSS"
    :widths: 20, 20

    :class:`DataFrame`, data file
    column, variable
    row, case
    groupby, split file
    :class:`NaN`, system-missing

:class:`DataFrame`
~~~~~~~~~~~~~~~~~~

A :class:`DataFrame` in pandas is analogous to an SPSS data file - a two-dimensional
data source with labeled columns that can be of different types. As will be shown in this
document, almost any operation that can be performed in SPSS can also be accomplished in pandas.

:class:`Series`
~~~~~~~~~~~~~~~

A :class:`Series` is the data structure that represents one column of a :class:`DataFrame`. SPSS doesn't have a
separate data structure for a single variable, but in general, working with a :class:`Series` is analogous
to working with a variable in SPSS.

:class:`Index`
~~~~~~~~~~~~~~

Every :class:`DataFrame` and :class:`Series` has an :class:`Index` -- labels on the *rows* of the data. SPSS does not
have an exact analogue, as cases are simply numbered sequentially from 1. In pandas, if no index is
specified, a :class:`RangeIndex` is used by default (first row = 0, second row = 1, and so on).

While using a labeled :class:`Index` or :class:`MultiIndex` can enable sophisticated analyses and is ultimately an
important part of pandas to understand, for this comparison we will essentially ignore the :class:`Index` and
just treat the :class:`DataFrame` as a collection of columns. Please see the :ref:`indexing documentation<indexing>`
for much more on how to use an :class:`Index` effectively.


Copies vs. in place operations
------------------------------

.. include:: includes/copies.rst


Data input / output
-------------------

Reading external data
~~~~~~~~~~~~~~~~~~~~~

Like SPSS, pandas provides utilities for reading in data from many formats. The ``tips`` dataset, found within
the pandas tests (`csv <https://raw.githubusercontent.com/pandas-dev/pandas/main/pandas/tests/io/data/csv/tips.csv>`_)
will be used in many of the following examples.

In SPSS, you would use File > Open > Data to import a CSV file:

.. code-block:: text

    FILE > OPEN > DATA
    /TYPE=CSV
    /FILE='tips.csv'
    /DELIMITERS=","
    /FIRSTCASE=2
    /VARIABLES=col1 col2 col3.

The pandas equivalent would use :func:`read_csv`:

.. code-block:: python

   url = (
       "https://raw.githubusercontent.com/pandas-dev"
       "/pandas/main/pandas/tests/io/data/csv/tips.csv"
   )
   tips = pd.read_csv(url)
   tips

Like SPSS's data import wizard, ``read_csv`` can take a number of parameters to specify how the data should be parsed.
For example, if the data was instead tab delimited, and did not have column names, the pandas command would be:

.. code-block:: python

   tips = pd.read_csv("tips.csv", sep="\t", header=None)

   # alternatively, read_table is an alias to read_csv with tab delimiter
   tips = pd.read_table("tips.csv", header=None)


Data operations
---------------

Filtering
~~~~~~~~~

In SPSS, filtering is done through Data > Select Cases:

.. code-block:: text

    SELECT IF (total_bill > 10).
    EXECUTE.

In pandas, boolean indexing can be used:

.. code-block:: python

    tips[tips["total_bill"] > 10]


Sorting
~~~~~~~

In SPSS, sorting is done through Data > Sort Cases:

.. code-block:: text

    SORT CASES BY sex total_bill.
    EXECUTE.

In pandas, this would be written as:

.. code-block:: python

    tips.sort_values(["sex", "total_bill"])


String processing
-----------------

Finding length of string
~~~~~~~~~~~~~~~~~~~~~~~~

In SPSS:

.. code-block:: text

    COMPUTE length = LENGTH(time).
    EXECUTE.

.. include:: includes/length.rst


Changing case
~~~~~~~~~~~~~

In SPSS:

.. code-block:: text

    COMPUTE upper = UPCASE(time).
    COMPUTE lower = LOWER(time).
    EXECUTE.

.. include:: includes/case.rst


Merging
-------

In SPSS, merging data files is done through Data > Merge Files.

.. include:: includes/merge_setup.rst
.. include:: includes/merge.rst


GroupBy operations
------------------

Split-file processing
~~~~~~~~~~~~~~~~~~~~~

In SPSS, split-file analysis is done through Data > Split File:

.. code-block:: text

    SORT CASES BY sex.
    SPLIT FILE BY sex.
    DESCRIPTIVES VARIABLES=total_bill tip
      /STATISTICS=MEAN STDDEV MIN MAX.

The pandas equivalent would be:

.. code-block:: python

    tips.groupby("sex")[["total_bill", "tip"]].agg(["mean", "std", "min", "max"])


Missing data
------------

SPSS uses the period (``.``) for numeric missing values and blank spaces for string missing values.
pandas uses ``NaN`` (Not a Number) for numeric missing values and ``None`` or ``NaN`` for string
missing values.

.. include:: includes/missing.rst


Other considerations
--------------------

Output management
-----------------

While pandas does not have a direct equivalent to SPSS's Output Management System (OMS), you can
capture and export results in various ways:

.. code-block:: python

   # Save summary statistics to CSV
   tips.groupby('sex')[['total_bill', 'tip']].mean().to_csv('summary.csv')

   # Save multiple results to Excel sheets
   with pd.ExcelWriter('results.xlsx') as writer:
       tips.describe().to_excel(writer, sheet_name='Descriptives')
       tips.groupby('sex').mean().to_excel(writer, sheet_name='Means by Gender')
