.. _overview:

.. currentmodule:: pandas

****************
Package overview
****************

:mod:`pandas` is an open source, BSD-licensed library providing high-performance, 
easy-to-use data structures and data analysis tools for the `Python <https://www.python.org/>`__
programming language.

:mod:`pandas` consists of the following elements:

* A set of labeled array data structures, the primary of which are
  Series and DataFrame.
* Index objects enabling both simple axis indexing and multi-level /
  hierarchical axis indexing.
* An integrated group by engine for aggregating and transforming data sets.
* Date range generation (date_range) and custom date offsets enabling the
  implementation of customized frequencies.
* Input/Output tools: loading tabular data from flat files (CSV, delimited,
  Excel 2003), and saving and loading pandas objects from the fast and
  efficient PyTables/HDF5 format.
* Memory-efficient "sparse" versions of the standard data structures for storing
  data that is mostly missing or mostly constant (some fixed value).
* Moving window statistics (rolling mean, rolling standard deviation, etc.).

Data Structures
---------------

.. csv-table::
    :header: "Dimensions", "Name", "Description"
    :widths: 15, 20, 50

    1, "Series", "1D labeled homogeneously-typed array"
    2, "DataFrame", "General 2D labeled, size-mutable tabular structure with potentially heterogeneously-typed column"

Why more than one data structure?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best way to think about the pandas data structures is as flexible
containers for lower dimensional data. For example, DataFrame is a container
for Series, and Series is a container for scalars. We would like to be
able to insert and remove objects from these containers in a dictionary-like
fashion.

Also, we would like sensible default behaviors for the common API functions
which take into account the typical orientation of time series and
cross-sectional data sets. When using ndarrays to store 2- and 3-dimensional
data, a burden is placed on the user to consider the orientation of the data
set when writing functions; axes are considered more or less equivalent (except
when C- or Fortran-contiguousness matters for performance). In pandas, the axes
are intended to lend more semantic meaning to the data; i.e., for a particular
data set there is likely to be a "right" way to orient the data. The goal,
then, is to reduce the amount of mental effort required to code up data
transformations in downstream functions.

For example, with tabular data (DataFrame) it is more semantically helpful to
think of the **index** (the rows) and the **columns** rather than axis 0 and
axis 1. Iterating through the columns of the DataFrame thus results in more
readable code:

::

    for col in df.columns:
        series = df[col]
        # do something with series

Mutability and copying of data
------------------------------

All pandas data structures are value-mutable (the values they contain can be
altered) but not always size-mutable. The length of a Series cannot be
changed, but, for example, columns can be inserted into a DataFrame. However,
the vast majority of methods produce new objects and leave the input data
untouched. In general we like to **favor immutability** where sensible.

Getting Support
---------------

The first stop for pandas issues and ideas is the `Github Issue Tracker
<https://github.com/pandas-dev/pandas/issues>`__. If you have a general question,
pandas community experts can answer through `Stack Overflow
<http://stackoverflow.com/questions/tagged/pandas>`__.

Community
---------

pandas is actively supported today by a community of like-minded individuals around 
the world who contribute their valuable time and energy to help make open source 
pandas possible. Thanks to `all of our contributors <https://github.com/pandas-dev/pandas/graphs/contributors>`__.

If you're interested in contributing, please
visit `Contributing to pandas webpage <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`__.

pandas is a `NumFOCUS <https://www.numfocus.org/open-source-projects/>`__ sponsored project.
This will help ensure the success of development of pandas as a world-class open-source
project, and makes it possible to `donate <https://pandas.pydata.org/donate.html>`__ to the project.

Project Governance
------------------

The governance process that pandas project has used informally since its inception in 2008 is formalized in `Project Governance documents <https://github.com/pandas-dev/pandas-governance>`__.
The documents clarify how decisions are made and how the various elements of our community interact, including the relationship between open source collaborative development and work that may be funded by for-profit or non-profit entities.

Wes McKinney is the Benevolent Dictator for Life (BDFL).

Development Team
-----------------

The list of the Core Team members and more detailed information can be found on the `people’s page <https://github.com/pandas-dev/pandas-governance/blob/master/people.md>`__ of the governance repo.
 

Institutional Partners
----------------------

The information about current institutional partners can be found on `pandas website page <https://pandas.pydata.org/about.html>`__.

Modules Privacy
----------------

The following modules in pandas are considered PRIVATE:

::

  * pandas.compat
  
  * pandas.util
      -decorators
      -print_versions
      -doctools
      -validators
      -depr_module
      
  * pandas.core
      -computation
      -indexes
      -sparse
      -reshape
      -dtypes
  	      

License
-------

.. literalinclude:: ../../LICENSE

