.. _cookbook:

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

********
Cookbook
********

This is a respository for *short and sweet* examples and links for useful pandas recipes.
We encourage users to add to this documentation.

This is a great *First Pull Request* (to add interesting links and/or put short code inline
for existing links)

Idioms
------

.. _cookbook.idioms:

These are some neat pandas ``idioms``

`How to do if-then-else?
<http://stackoverflow.com/questions/17128302/python-pandas-idiom-for-if-then-else>`__

`How to split a frame with a boolean criterion?
<http://stackoverflow.com/questions/14957116/how-to-split-a-dataframe-according-to-a-boolean-criterion>`__

`How to select from a frame with complex criteria?
<http://stackoverflow.com/questions/15315452/selecting-with-complex-criteria-from-pandas-dataframe>`__

`Select rows closest to a user defined numer
<http://stackoverflow.com/questions/17758023/return-rows-in-a-dataframe-closest-to-a-user-defined-number>`__

.. _cookbook.selection:

Selection
---------

The :ref:`indexing <indexing>` docs.

Indexing using both row labels and conditionals, see
`here
<http://stackoverflow.com/questions/14725068/pandas-using-row-labels-in-boolean-indexing>`__

Use loc for label-oriented slicing and iloc positional slicing, see
`here
<https://github.com/pydata/pandas/issues/2904>`__

Extend a panel frame by transposing, adding a new dimension, and transposing back to the original dimensions, see
`here
<http://stackoverflow.com/questions/15364050/extending-a-pandas-panel-frame-along-the-minor-axis>`__

Mask a panel by using ``np.where`` and then reconstructing the panel with the new masked values
`here
<http://stackoverflow.com/questions/14650341/boolean-mask-in-pandas-panel>`__

Using ``~`` to take the complement of a boolean array, see
`here
<http://stackoverflow.com/questions/14986510/picking-out-elements-based-on-complement-of-indices-in-python-pandas>`__

`Efficiently creating columns using applymap
<http://stackoverflow.com/questions/16575868/efficiently-creating-additional-columns-in-a-pandas-dataframe-using-map>`__

.. _cookbook.multi_index:

MultiIndexing
-------------

The :ref:`multindexing <indexing.hierarchical>` docs.

`Creating a multi-index from a labeled frame
<http://stackoverflow.com/questions/14916358/reshaping-dataframes-in-pandas-based-on-column-labels>`__

Slicing
~~~~~~~

`Slicing a multi-index with xs
<http://stackoverflow.com/questions/12590131/how-to-slice-multindex-columns-in-pandas-dataframes>`__

`Slicing a multi-index with xs #2
<http://stackoverflow.com/questions/14964493/multiindex-based-indexing-in-pandas>`__

Sorting
~~~~~~~

`Multi-index sorting
<http://stackoverflow.com/questions/14733871/mutli-index-sorting-in-pandas>`__

`Partial Selection, the need for sortedness
<https://github.com/pydata/pandas/issues/2995>`__

Levels
~~~~~~

`Prepending a level to a multiindex
<http://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex>`__

`Flatten Hierarchical columns
<http://stackoverflow.com/questions/14507794/python-pandas-how-to-flatten-a-hierarchical-index-in-columns>`__

.. _cookbook.missing_data:

Missing Data
------------

The :ref:`missing data <missing_data>` docs.

Replace
~~~~~~~

`Using replace with backrefs
<http://stackoverflow.com/questions/16818871/extracting-value-and-creating-new-column-out-of-it>`__

.. _cookbook.grouping:

Grouping
--------

The :ref:`grouping <groupby>` docs.

`Basic grouping with apply
<http://stackoverflow.com/questions/15322632/python-pandas-df-groupy-agg-column-reference-in-agg>`__

`Using get_group
<http://stackoverflow.com/questions/14734533/how-to-access-pandas-groupby-dataframe-by-key>`__

`Apply to different items in a group
<http://stackoverflow.com/questions/15262134/apply-different-functions-to-different-items-in-group-object-python-pandas>`__

`Expanding Apply
<http://stackoverflow.com/questions/14542145/reductions-down-a-column-in-pandas>`__

`Replacing values with groupby means
<http://stackoverflow.com/questions/14760757/replacing-values-with-groupby-means>`__

`Sort by group with aggregation
<http://stackoverflow.com/questions/14941366/pandas-sort-by-group-aggregate-and-column>`__

`Create multiple aggregated columns
<http://stackoverflow.com/questions/14897100/create-multiple-columns-in-pandas-aggregation-function>`__

`Create a value counts column and reassign back to the DataFrame
<http://stackoverflow.com/questions/17709270/i-want-to-create-a-column-of-value-counts-in-my-pandas-dataframe>`__

Expanding Data
~~~~~~~~~~~~~~

`Alignment and to-date
<http://stackoverflow.com/questions/15489011/python-time-series-alignment-and-to-date-functions>`__

`Rolling Computation window based on values instead of counts
<http://stackoverflow.com/questions/14300768/pandas-rolling-computation-with-window-based-on-values-instead-of-counts>`__

`Rolling Mean by Time Interval
<http://stackoverflow.com/questions/15771472/pandas-rolling-mean-by-time-interval>`__

Splitting
~~~~~~~~~

`Splitting a frame
<http://stackoverflow.com/questions/13353233/best-way-to-split-a-dataframe-given-an-edge/15449992#15449992>`__

.. _cookbook.pivot:

Pivot
~~~~~
The :ref:`Pivot <reshaping.pivot>` docs.

`Partial sums and subtotals
<http://stackoverflow.com/questions/15570099/pandas-pivot-tables-row-subtotals/15574875#15574875>`__

`Frequency table like plyr in R
<http://stackoverflow.com/questions/15589354/frequency-tables-in-pandas-like-plyr-in-r>`__

Apply
~~~~~

`Turning embeded lists into a multi-index frame
<http://stackoverflow.com/questions/17349981/converting-pandas-dataframe-with-categorical-values-into-binary-values>`__

Timeseries
----------

`Between times
<http://stackoverflow.com/questions/14539992/pandas-drop-rows-outside-of-time-range>`__

`Using indexer between time
<http://stackoverflow.com/questions/17559885/pandas-dataframe-mask-based-on-index>`__

`Vectorized Lookup
<http://stackoverflow.com/questions/13893227/vectorized-look-up-of-values-in-pandas-dataframe>`__

Turn a matrix with hours in columns and days in rows into a continous row sequence in the form of a time series.
`How to rearrange a python pandas dataframe?
<http://stackoverflow.com/questions/15432659/how-to-rearrange-a-python-pandas-dataframe>`__

.. _cookbook.resample:

Resampling
~~~~~~~~~~

The :ref:`Resample <timeseries.resampling>` docs.

`TimeGrouping of values grouped across time
<http://stackoverflow.com/questions/15297053/how-can-i-divide-single-values-of-a-dataframe-by-monthly-averages>`__

`TimeGrouping #2
<http://stackoverflow.com/questions/14569223/timegrouper-pandas>`__

`Using TimeGrouper and another grouping to create subgroups, then apply a custom function
<https://github.com/pydata/pandas/issues/3791>`__

`Resampling with custom periods
<http://stackoverflow.com/questions/15408156/resampling-with-custom-periods>`__

`Resample intraday frame without adding new days
<http://stackoverflow.com/questions/14898574/resample-intrday-pandas-dataframe-without-add-new-days>`__

`Resample minute data
<http://stackoverflow.com/questions/14861023/resampling-minute-data>`__

.. _cookbook.merge:

Merge
-----

The :ref:`Concat <merging.concatenation>` docs. The :ref:`Join <merging.join>` docs.

`emulate R rbind
<http://stackoverflow.com/questions/14988480/pandas-version-of-rbind>`__

`Self Join
<https://github.com/pydata/pandas/issues/2996>`__

`How to set the index and join
<http://stackoverflow.com/questions/14341805/pandas-merge-pd-merge-how-to-set-the-index-and-join>`__

`KDB like asof join
<http://stackoverflow.com/questions/12322289/kdb-like-asof-join-for-timeseries-data-in-pandas/12336039#12336039>`__

`Join with a criteria based on the values
<http://stackoverflow.com/questions/15581829/how-to-perform-an-inner-or-outer-join-of-dataframes-with-pandas-on-non-simplisti>`__

.. _cookbook.plotting:

Plotting
--------

The :ref:`Plotting <visualization>` docs.

`Make Matplotlib look like R
<http://stackoverflow.com/questions/14349055/making-matplotlib-graphs-look-like-r-by-default>`__

`Setting x-axis major and minor labels
<http://stackoverflow.com/questions/12945971/pandas-timeseries-plot-setting-x-axis-major-and-minor-ticks-and-labels>`__

`Plotting multiple charts in an ipython notebook
<http://stackoverflow.com/questions/16392921/make-more-than-one-chart-in-same-ipython-notebook-cell>`__

`Creating a multi-line plot
<http://stackoverflow.com/questions/16568964/make-a-multiline-plot-from-csv-file-in-matplotlib>`__

`Plotting a heatmap
<http://stackoverflow.com/questions/17050202/plot-timeseries-of-histograms-in-python>`__

`Annotate a time-series plot
<http://stackoverflow.com/questions/11067368/annotate-time-series-plot-in-matplotlib>`__

Data In/Out
-----------

`Performance comparison of SQL vs HDF5
<http://stackoverflow.com/questions/16628329/hdf5-and-sqlite-concurrency-compression-i-o-performance>`__

.. _cookbook.csv:

CSV
~~~

The :ref:`CSV <io.read_csv_table>` docs

`read_csv in action
<http://wesmckinney.com/blog/?p=635>`__

`appending to a csv
<http://stackoverflow.com/questions/17134942/pandas-dataframe-output-end-of-csv>`__

`Reading a csv chunk-by-chunk
<http://stackoverflow.com/questions/11622652/large-persistent-dataframe-in-pandas/12193309#12193309>`__

`Reading the first few lines of a frame
<http://stackoverflow.com/questions/15008970/way-to-read-first-few-lines-for-pandas-dataframe>`__

Reading a file that is compressed but not by ``gzip/bz2`` (the native compresed formats which ``read_csv`` understands).
This example shows a ``WinZipped`` file, but is a general application of opening the file within a context manager and
using that handle to read.
`See here
<http://stackoverflow.com/questions/17789907/pandas-convert-winzipped-csv-file-to-data-frame>`__

`Inferring dtypes from a file
<http://stackoverflow.com/questions/15555005/get-inferred-dataframe-types-iteratively-using-chunksize>`__

`Dealing with bad lines
<https://github.com/pydata/pandas/issues/2886>`__

`Dealing with bad lines II
<http://nipunbatra.wordpress.com/2013/06/06/reading-unclean-data-csv-using-pandas/>`__

`Reading CSV with Unix timestamps and converting to local timezone
<http://nipunbatra.wordpress.com/2013/06/07/pandas-reading-csv-with-unix-timestamps-and-converting-to-local-timezone/>`__

`Write a multi-row index CSV without writing duplicates
<http://stackoverflow.com/questions/17349574/pandas-write-multiindex-rows-with-to-csv>`__

.. _cookbook.sql:

SQL
~~~

The :ref:`SQL <io.sql>` docs

`Reading from databases with SQL
<http://stackoverflow.com/questions/10065051/python-pandas-and-databases-like-mysql>`__

.. _cookbook.excel:

Excel
~~~~~

The :ref:`Excel <io.excel>` docs

`Reading from a filelike handle
<http://stackoverflow.com/questions/15588713/sheets-of-excel-workbook-from-a-url-into-a-pandas-dataframe>`__

.. _cookbook.hdf:

HDFStore
~~~~~~~~

The :ref:`HDFStores <io.hdf5>` docs

`Simple Queries with a Timestamp Index
<http://stackoverflow.com/questions/13926089/selecting-columns-from-pandas-hdfstore-table>`__

`Managing heteregenous data using a linked multiple table hierarchy
<https://github.com/pydata/pandas/issues/3032>`__

`Merging on-disk tables with millions of rows
<http://stackoverflow.com/questions/14614512/merging-two-tables-with-millions-of-rows-in-python/14617925#14617925>`__

Deduplicating a large store by chunks, essentially a recusive reduction operation. Shows a function for taking in data from
csv file and creating a store by chunks, with date parsing as well.
`See here
<http://stackoverflow.com/questions/16110252/need-to-compare-very-large-files-around-1-5gb-in-python/16110391#16110391>`__

`Large Data work flows
<http://stackoverflow.com/questions/14262433/large-data-work-flows-using-pandas>`__

`Reading in a sequence of files, then providing a global unique index to a store while appending
<http://stackoverflow.com/questions/16997048/how-does-one-append-large-amounts-of-data-to-a-pandas-hdfstore-and-get-a-natural>`__

`Groupby on a HDFStore
<http://stackoverflow.com/questions/15798209/pandas-group-by-query-on-large-data-in-hdfstore>`__

`Troubleshoot HDFStore exceptions
<http://stackoverflow.com/questions/15488809/how-to-trouble-shoot-hdfstore-exception-cannot-find-the-correct-atom-type>`__

`Setting min_itemsize with strings
<http://stackoverflow.com/questions/15988871/hdfstore-appendstring-dataframe-fails-when-string-column-contents-are-longer>`__

Storing Attributes to a group node

.. ipython:: python

    df = DataFrame(np.random.randn(8,3))
    store = HDFStore('test.h5')
    store.put('df',df)

    # you can store an arbitrary python object via pickle
    store.get_storer('df').attrs.my_attribute = dict(A = 10)
    store.get_storer('df').attrs.my_attribute

.. ipython:: python
   :suppress:

    store.close()
    os.remove('test.h5')

Computation
-----------

`Numerical integration (sample-based) of a time series
<http://nbviewer.ipython.org/5720498>`__

Miscellaneous
-------------

The :ref:`Timedeltas <timeseries.timedeltas>` docs.

`Operating with timedeltas
<https://github.com/pydata/pandas/pull/2899>`__

`Create timedeltas with date differences
<http://stackoverflow.com/questions/15683588/iterating-through-a-pandas-dataframe>`__

`Adding days to dates in a dataframe
<http://stackoverflow.com/questions/16385785/add-days-to-dates-in-dataframe>`__

Aliasing Axis Names
-------------------

To globally provide aliases for axis names, one can define these 2 functions:

.. ipython:: python

   def set_axis_alias(cls, axis, alias):
        if axis not in cls._AXIS_NUMBERS:
            raise Exception("invalid axis [%s] for alias [%s]" % (axis, alias))
        cls._AXIS_ALIASES[alias] = axis

   def clear_axis_alias(cls, axis, alias):
        if axis not in cls._AXIS_NUMBERS:
            raise Exception("invalid axis [%s] for alias [%s]" % (axis, alias))
        cls._AXIS_ALIASES.pop(alias,None)


   set_axis_alias(DataFrame,'columns', 'myaxis2')
   df2 = DataFrame(randn(3,2),columns=['c1','c2'],index=['i1','i2','i3'])
   df2.sum(axis='myaxis2')
   clear_axis_alias(DataFrame,'columns', 'myaxis2')
