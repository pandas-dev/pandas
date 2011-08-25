.. _overview:

.. currentmodule:: pandas

****************
Package overview
****************

:mod:`pandas` is a library providing, among other things, a set of
convenient and powerful data structures for working with labeled
statistical (financial, economic, econometric) data sets. We will
refer to this data as *time series* and *cross-sectional* (or
*longitudinal*) which are common terms in statistics and
econometrics. pandas has multiple target audiences:

 * Users of R or MATLAB who wish to switch to Python for interactive
   data analysis and implementation of statistical models

 * NumPy users who are looking for richer data structures for working
   with time series and cross-sectional data.

 * System developers who wish to have a robust and well-tested library
   for building production applications involving such data sets.

Data structures at a glance
---------------------------

.. csv-table::
    :header: "Dimensions", "Name", "Description"
    :widths: 10, 15, 50

    1, Series, "Most generic 1D structure"
    1, TimeSeries, "Series indexed by datetimes"
    2, DataFrame, "General 2D indexed tabular structure"
    3, WidePanel, "General 3D panel data"
    3, LongPanel, "Stacked (2D) format panel data"

Why more than 1 data structure?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best way to think about the pandas data structures is as flexible containers
for lower dimensional data. For example, DataFrame is a container for Series,
and WidePanel is a container for DataFrame objects. We would like to be able to
insert and remove objects from these containers in a dictionary-like fashion.

Also, we would like sensible default behaviors for the common API functions
which take into account the typical orientation of time series and
cross-sectional data sets. When using ndarrays to store 2- and 3-dimensional
data, a burden is placed on the user to consider the orientation of the data set
when writing functions; axes are considered more or less equivalent (except when
C- or Fortran-contiguousness matters for performance). In pandas, the axes are
intended to lend more semantic meaning to the data; i.e., for a particular data
set there is likely to be a "right" way to orient the data. The goal, then, is
to reduce the amount of thought required to code up data transformations in
downstream functions.

Lest we be too hand-wavy, here are some common use cases to
illustrate:

 (A) :ref:`DataFrame <dataframe>` containing multiple related time series

  * **columns**: "data type" associated with each time series
  * **index**:  dates shared by time series

 (B) :ref:`DataFrame <dataframe>` containing multiple cross-sections

  * **columns**: "data type" associated with each cross-section
  * **index**:  individual / entity labels common to cross-sections

 (C) :ref:`WidePanel <panel>` containing panel data

  * **items**: "data type" associated with each collection of time series
  * **major_axis**: dates shared by time series
  * **minor_axis**: individual / entity labels common to time series

Lastly, particularly if you don't buy the above explanation, having a
specialized vocabulary to refer to types of data sets often serves as
a benefit when discussing a dataset with other users (or reading their
code).

A quick note on mutation
~~~~~~~~~~~~~~~~~~~~~~~~

Most instance methods on the pandas data structures return a new
object, rather than updating the original object in-place. However,
when working with the contents (e.g. a column in a DataFrame),
mutations **will** be reflected in the original structure. In general,
though, we like to "favor immutability" where sensible.

What else is in the package?
----------------------------
