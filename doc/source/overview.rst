.. _overview:

.. currentmodule:: pandas

****************
Package overview
****************

:mod:`pandas` consists of the following things

 * A set of labeled array data structures, the primary of which are
   Series/TimeSeries and DataFrame
 * Index objects enabling both simple axis indexing and multi-level /
   hierarchical axis indexing
 * An integrated group by engine for aggregating and transforming data sets
 * Date range generation (DateRange) and custom date offsets enabling the
   implementation of customized frequencies
 * Input/Output tools: loading tabular data from flat files (CSV, delimited,
   Excel 2003), and saving and loading pandas objects from the fast and
   efficient PyTables/HDF5 format.
 * Memory-efficent "sparse" versions of the standard data structures for storing
   data that is mostly missing or mostly constant (some fixed value)
 * Moving window statistics (rolling mean, rolling standard deviation, etc.)
 * Static and moving window linear and `panel regression
   <http://en.wikipedia.org/wiki/Panel_data>`__

License
-------

pandas is released under a standard 3-clause BSD license

Data structures at a glance
---------------------------

.. csv-table::
    :header: "Dimensions", "Name", "Description"
    :widths: 15, 20, 50

    1, Series, "1D labeled homogeneously-typed array"
    1, TimeSeries, "Series with index containing datetimes"
    2, DataFrame, "General 2D labeled, size-mutable tabular structure with
    potentially heterogeneously-typed columns"
    3, Panel, "General 3D labeled, also size-mutable array"

Why more than 1 data structure?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best way to think about the pandas data structures is as flexible
containers for lower dimensional data. For example, DataFrame is a container
for Series, and Panel is a container for DataFrame objects. We would like to be
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
axis 1. And iterating through the columns of the DataFrame thus results in more
readable code:

::

    for col in df.columns:
        series = df[col]
        # do something with series

Mutability and copying of data
------------------------------

All pandas data structures are value-mutable (the values they contain can be
altered) but not always size-mutable. The length of a Series cannot be changed,
but, for example, columns can be inserted into a DataFrame. However, the vast
majority of methods produce new objects and leave the input data untouched. In
general, though, we like to **favor immutability** where sensible.


History
-------

pandas development began at `AQR Capital Management <http://www.aqr.com>`__ in
April 2008. It was open-sourced at the end of 2009 and continues to be actively
used and maintained.

Contact
-------

Please feel free to send comments or questions directly to Wes McKinney at
wesmckinn (-at-) gmail (-dot-) com or the pystatsmodels mailing list.
