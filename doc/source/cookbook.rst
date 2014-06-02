.. _cookbook:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   import os
   np.random.seed(123456)
   from pandas import *
   options.display.max_rows=15
   options.display.mpl_style='default'
   import pandas as pd
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)

********
Cookbook
********

This is a repository for *short and sweet* examples and links for useful pandas recipes.
We encourage users to add to this documentation.

This is a great *First Pull Request* (to add interesting links and/or put short code inline
for existing links)

Idioms
------

.. _cookbook.idioms:

These are some neat pandas ``idioms``

`How to do if-then-else?
<http://stackoverflow.com/questions/17128302/python-pandas-idiom-for-if-then-else>`__

`How to do if-then-else #2
<http://stackoverflow.com/questions/19913659/pandas-conditional-creation-of-a-series-dataframe-column>`__

`How to split a frame with a boolean criterion?
<http://stackoverflow.com/questions/14957116/how-to-split-a-dataframe-according-to-a-boolean-criterion>`__

`How to select from a frame with complex criteria?
<http://stackoverflow.com/questions/15315452/selecting-with-complex-criteria-from-pandas-dataframe>`__

`Select rows closest to a user-defined number
<http://stackoverflow.com/questions/17758023/return-rows-in-a-dataframe-closest-to-a-user-defined-number>`__

`How to reduce a sequence (e.g. of Series) using a binary operator
<http://stackoverflow.com/questions/21058254/pandas-boolean-operation-in-a-python-list/21058331>`__


.. _cookbook.selection:

Selection
---------

The :ref:`indexing <indexing>` docs.

`Indexing using both row labels and conditionals
<http://stackoverflow.com/questions/14725068/pandas-using-row-labels-in-boolean-indexing>`__

`Use loc for label-oriented slicing and iloc positional slicing
<https://github.com/pydata/pandas/issues/2904>`__

`Extend a panel frame by transposing, adding a new dimension, and transposing back to the original dimensions
<http://stackoverflow.com/questions/15364050/extending-a-pandas-panel-frame-along-the-minor-axis>`__

`Mask a panel by using np.where and then reconstructing the panel with the new masked values
<http://stackoverflow.com/questions/14650341/boolean-mask-in-pandas-panel>`__

`Using ~ to take the complement of a boolean array, see
<http://stackoverflow.com/questions/14986510/picking-out-elements-based-on-complement-of-indices-in-python-pandas>`__

`Efficiently creating columns using applymap
<http://stackoverflow.com/questions/16575868/efficiently-creating-additional-columns-in-a-pandas-dataframe-using-map>`__

`Keep other columns when using min() with groupby
<http://stackoverflow.com/questions/23394476/keep-other-columns-when-using-min-with-groupby>`__

.. _cookbook.multi_index:

MultiIndexing
-------------

The :ref:`multindexing <indexing.hierarchical>` docs.

`Creating a multi-index from a labeled frame
<http://stackoverflow.com/questions/14916358/reshaping-dataframes-in-pandas-based-on-column-labels>`__

Arithmetic
~~~~~~~~~~

`Performing arithmetic with a multi-index that needs broadcasting
<http://stackoverflow.com/questions/19501510/divide-entire-pandas-multiindex-dataframe-by-dataframe-variable/19502176#19502176>`__

Slicing
~~~~~~~

`Slicing a multi-index with xs
<http://stackoverflow.com/questions/12590131/how-to-slice-multindex-columns-in-pandas-dataframes>`__

`Slicing a multi-index with xs #2
<http://stackoverflow.com/questions/14964493/multiindex-based-indexing-in-pandas>`__

`Setting portions of a multi-index with xs
<http://stackoverflow.com/questions/19319432/pandas-selecting-a-lower-level-in-a-dataframe-to-do-a-ffill>`__

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

panelnd
~~~~~~~

The :ref:`panelnd<dsintro.panelnd>` docs.

`Construct a 5D panelnd
<http://stackoverflow.com/questions/18748598/why-my-panelnd-factory-throwing-a-keyerror>`__

.. _cookbook.missing_data:

Missing Data
------------

The :ref:`missing data<missing_data>` docs.

Fill forward a reversed timeseries

.. ipython:: python

   df = pd.DataFrame(np.random.randn(6,1), index=pd.date_range('2013-08-01', periods=6, freq='B'), columns=list('A'))
   df.ix[3,'A'] = np.nan
   df
   df.reindex(df.index[::-1]).ffill()

`cumsum reset at NaN values
<http://stackoverflow.com/questions/18196811/cumsum-reset-at-nan>`__

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

`Shift groups of the values in a column based on the index
<http://stackoverflow.com/q/23198053/190597>`__

.. ipython:: python

   df = pd.DataFrame(
        {u'line_race': [10L, 10L, 8L, 10L, 10L, 8L],
         u'beyer': [99L, 102L, 103L, 103L, 88L, 100L]},
        index=[u'Last Gunfighter', u'Last Gunfighter', u'Last Gunfighter',
               u'Paynter', u'Paynter', u'Paynter']); df

   df['beyer_shifted'] = df.groupby(level=0)['beyer'].shift(1)
   df

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

`Turning embedded lists into a multi-index frame
<http://stackoverflow.com/questions/17349981/converting-pandas-dataframe-with-categorical-values-into-binary-values>`__

`Rolling apply with a DataFrame returning a Series
<http://stackoverflow.com/questions/19121854/using-rolling-apply-on-a-dataframe-object>`__

`Rolling apply with a DataFrame returning a Scalar
<http://stackoverflow.com/questions/21040766/python-pandas-rolling-apply-two-column-input-into-function/21045831#21045831>`__

Timeseries
----------

`Between times
<http://stackoverflow.com/questions/14539992/pandas-drop-rows-outside-of-time-range>`__

`Using indexer between time
<http://stackoverflow.com/questions/17559885/pandas-dataframe-mask-based-on-index>`__

`Constructing a datetime range that excludes weekends and includes only certain times
<http://stackoverflow.com/questions/24010830/pandas-generate-sequential-timestamp-with-jump/24014440#24014440?>`__

`Vectorized Lookup
<http://stackoverflow.com/questions/13893227/vectorized-look-up-of-values-in-pandas-dataframe>`__

Turn a matrix with hours in columns and days in rows into a continuous row sequence in the form of a time series.
`How to rearrange a python pandas DataFrame?
<http://stackoverflow.com/questions/15432659/how-to-rearrange-a-python-pandas-dataframe>`__

`Dealing with duplicates when reindexing a timeseries to a specified frequency
<http://stackoverflow.com/questions/22244383/pandas-df-refill-adding-two-columns-of-different-shape>`__

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

`Resample with groupby <http://stackoverflow.com/q/18677271/564538>`__

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

`Annotate a time-series plot #2
<http://stackoverflow.com/questions/17891493/annotating-points-from-a-pandas-dataframe-in-matplotlib-plot>`__

`Generate Embedded plots in excel files using Pandas, Vincent and xlsxwriter
<http://pandas-xlsxwriter-charts.readthedocs.org/en/latest/introduction.html>`__

`Boxplot for each quartile of a stratifying variable
<http://stackoverflow.com/questions/23232989/boxplot-stratified-by-column-in-python-pandas>`__

.. ipython:: python

    df = pd.DataFrame(
        {u'stratifying_var': np.random.uniform(0, 100, 20),
         u'price': np.random.normal(100, 5, 20)}
    )
    df[u'quartiles'] = pd.qcut(
        df[u'stratifying_var'],
        4,
        labels=[u'0-25%', u'25-50%', u'50-75%', u'75-100%']
    )

    @savefig quartile_boxplot.png
    df.boxplot(column=u'price', by=u'quartiles')


Data In/Out
-----------

`Performance comparison of SQL vs HDF5
<http://stackoverflow.com/questions/16628329/hdf5-and-sqlite-concurrency-compression-i-o-performance>`__

.. _cookbook.csv:

CSV
~~~

The :ref:`CSV <io.read_csv_table>` docs

`read_csv in action <http://wesmckinney.com/blog/?p=635>`__

`appending to a csv
<http://stackoverflow.com/questions/17134942/pandas-dataframe-output-end-of-csv>`__

`Reading a csv chunk-by-chunk
<http://stackoverflow.com/questions/11622652/large-persistent-dataframe-in-pandas/12193309#12193309>`__

`Reading only certain rows of a csv chunk-by-chunk
<http://stackoverflow.com/questions/19674212/pandas-data-frame-select-rows-and-clear-memory>`__

`Reading the first few lines of a frame
<http://stackoverflow.com/questions/15008970/way-to-read-first-few-lines-for-pandas-dataframe>`__

Reading a file that is compressed but not by ``gzip/bz2`` (the native compressed formats which ``read_csv`` understands).
This example shows a ``WinZipped`` file, but is a general application of opening the file within a context manager and
using that handle to read.
`See here
<http://stackoverflow.com/questions/17789907/pandas-convert-winzipped-csv-file-to-data-frame>`__

`Inferring dtypes from a file
<http://stackoverflow.com/questions/15555005/get-inferred-dataframe-types-iteratively-using-chunksize>`__

`Dealing with bad lines
<http://github.com/pydata/pandas/issues/2886>`__

`Dealing with bad lines II
<http://nipunbatra.github.io/2013/06/reading-unclean-data-csv-using-pandas/>`__

`Reading CSV with Unix timestamps and converting to local timezone
<http://nipunbatra.github.io/2013/06/pandas-reading-csv-with-unix-timestamps-and-converting-to-local-timezone/>`__

`Write a multi-row index CSV without writing duplicates
<http://stackoverflow.com/questions/17349574/pandas-write-multiindex-rows-with-to-csv>`__

Parsing date components in multi-columns is faster with a format

.. code-block:: python

    In [30]: i = pd.date_range('20000101',periods=10000)

    In [31]: df = pd.DataFrame(dict(year = i.year, month = i.month, day = i.day))

    In [32]: df.head()
    Out[32]:
       day  month  year
    0    1      1  2000
    1    2      1  2000
    2    3      1  2000
    3    4      1  2000
    4    5      1  2000

    In [33]: %timeit pd.to_datetime(df.year*10000+df.month*100+df.day,format='%Y%m%d')
    100 loops, best of 3: 7.08 ms per loop

    # simulate combinging into a string, then parsing
    In [34]: ds = df.apply(lambda x: "%04d%02d%02d" % (x['year'],x['month'],x['day']),axis=1)

    In [35]: ds.head()
    Out[35]:
    0    20000101
    1    20000102
    2    20000103
    3    20000104
    4    20000105
    dtype: object

    In [36]: %timeit pd.to_datetime(ds)
    1 loops, best of 3: 488 ms per loop

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

.. _cookbook.html:

`Reading HTML tables from a server that cannot handle the default request
header <http://stackoverflow.com/a/18939272/564538>`__

.. _cookbook.hdf:

HDFStore
~~~~~~~~

The :ref:`HDFStores <io.hdf5>` docs

`Simple Queries with a Timestamp Index
<http://stackoverflow.com/questions/13926089/selecting-columns-from-pandas-hdfstore-table>`__

`Managing heterogeneous data using a linked multiple table hierarchy
<http://github.com/pydata/pandas/issues/3032>`__

`Merging on-disk tables with millions of rows
<http://stackoverflow.com/questions/14614512/merging-two-tables-with-millions-of-rows-in-python/14617925#14617925>`__

Deduplicating a large store by chunks, essentially a recursive reduction operation. Shows a function for taking in data from
csv file and creating a store by chunks, with date parsing as well.
`See here
<http://stackoverflow.com/questions/16110252/need-to-compare-very-large-files-around-1-5gb-in-python/16110391#16110391>`__

`Creating a store chunk-by-chunk from a csv file
<http://stackoverflow.com/questions/20428355/appending-column-to-frame-of-hdf-file-in-pandas/20428786#20428786>`__

`Appending to a store, while creating a unique index
<http://stackoverflow.com/questions/16997048/how-does-one-append-large-amounts-of-data-to-a-pandas-hdfstore-and-get-a-natural/16999397#16999397>`__

`Large Data work flows
<http://stackoverflow.com/questions/14262433/large-data-work-flows-using-pandas>`__

`Reading in a sequence of files, then providing a global unique index to a store while appending
<http://stackoverflow.com/questions/16997048/how-does-one-append-large-amounts-of-data-to-a-pandas-hdfstore-and-get-a-natural>`__

`Groupby on a HDFStore
<http://stackoverflow.com/questions/15798209/pandas-group-by-query-on-large-data-in-hdfstore>`__

`Hierarchical queries on a HDFStore
<http://stackoverflow.com/questions/22777284/improve-query-performance-from-a-large-hdfstore-table-with-pandas/22820780#22820780>`__

`Counting with a HDFStore
<http://stackoverflow.com/questions/20497897/converting-dict-of-dicts-into-pandas-dataframe-memory-issues>`__

`Troubleshoot HDFStore exceptions
<http://stackoverflow.com/questions/15488809/how-to-trouble-shoot-hdfstore-exception-cannot-find-the-correct-atom-type>`__

`Setting min_itemsize with strings
<http://stackoverflow.com/questions/15988871/hdfstore-appendstring-dataframe-fails-when-string-column-contents-are-longer>`__

`Using ptrepack to create a completely-sorted-index on a store
<http://stackoverflow.com/questions/17893370/ptrepack-sortby-needs-full-index>`__

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


.. _cookbook.binary:

Binary Files
~~~~~~~~~~~~

pandas readily accepts numpy record arrays, if you need to read in a binary
file consisting of an array of C structs. For example, given this C program
in a file called ``main.c`` compiled with ``gcc main.c -std=gnu99`` on a
64-bit machine,

.. code-block:: c

   #include <stdio.h>
   #include <stdint.h>

   typedef struct _Data
   {
       int32_t count;
       double avg;
       float scale;
   } Data;

   int main(int argc, const char *argv[])
   {
       size_t n = 10;
       Data d[n];

       for (int i = 0; i < n; ++i)
       {
           d[i].count = i;
           d[i].avg = i + 1.0;
           d[i].scale = (float) i + 2.0f;
       }

       FILE *file = fopen("binary.dat", "wb");
       fwrite(&d, sizeof(Data), n, file);
       fclose(file);

       return 0;
   }

the following Python code will read the binary file ``'binary.dat'`` into a
pandas ``DataFrame``, where each element of the struct corresponds to a column
in the frame:

.. code-block:: python

   import numpy as np
   from pandas import DataFrame

   names = 'count', 'avg', 'scale'

   # note that the offsets are larger than the size of the type because of
   # struct padding
   offsets = 0, 8, 16
   formats = 'i4', 'f8', 'f4'
   dt = np.dtype({'names': names, 'offsets': offsets, 'formats': formats},
                 align=True)
   df = DataFrame(np.fromfile('binary.dat', dt))

.. note::

   The offsets of the structure elements may be different depending on the
   architecture of the machine on which the file was created. Using a raw
   binary file format like this for general data storage is not recommended, as
   it is not cross platform. We recommended either HDF5 or msgpack, both of
   which are supported by pandas' IO facilities.

Computation
-----------

`Numerical integration (sample-based) of a time series
<http://nbviewer.ipython.org/5720498>`__

Miscellaneous
-------------

The :ref:`Timedeltas <timeseries.timedeltas>` docs.

`Operating with timedeltas
<http://github.com/pydata/pandas/pull/2899>`__

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
