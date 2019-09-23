.. _scale:

*************************
Scaling to large datasets
*************************

Pandas provide data structures for in-memory analytics. This makes using pandas
to analyze larger than memory datasets somewhat tricky.

This document provides a few recommendations for scaling to larger datasets.
It's a complement to :ref:`enhancingperf`, which focuses on speeding up analysis
for datasets that fit in memory.

But first, it's worth considering *not using pandas*. Pandas isn't the right
tool for all situations. If you're working with very large datasets and a tool
like PostgreSQL fits your needs, then you should probably be using that.
Assuming you want or need the expressivity and power of pandas, let's carry on.

.. ipython:: python

   import pandas as pd
   import numpy as np
   from pandas.util.testing import make_timeseries


Use more efficient file formats
-------------------------------

Depending on your workload, data loading may be a bottleneck. In these case you
might consider switching from a slow format like CSV to a faster format like
Parquet. Loading from a file format like Parquet will also require less memory
usage, letting you load larger datasets into pandas before running out of
memory.

.. ipython:: python

   # Make a random in-memory dataset
   ts = make_timeseries(freq="30S", seed=0)
   ts


We'll now write and read the file using CSV and parquet.


.. ipython:: python

   %time ts.to_csv("timeseries.csv")

.. ipython:: python

   %time ts2 = pd.read_csv("timeseries.csv", index_col="timestamp", parse_dates=["timestamp"])

.. ipython:: python

   %time ts.to_parquet("timeseries.parquet")

.. ipython:: python

   %time _ = pd.read_parquet("timeseries.parquet")

Notice that parquet gives much higher performance for reading and writing, both
in terms of speed and lower peak memory usage. See :ref:`io` for more.

Load less data
--------------

Suppose our raw dataset on disk has many columns, but we need just a subset
for our analysis. To get those columns, we can either

1. Load the entire dataset then select those columns.
2. Just load the columns we need.

Loading just the columns you need can be much faster and requires less memory.

.. ipython:: python

   # make a similar dataset with many columns
   timeseries = [
       make_timeseries(freq="1T", seed=i).rename(columns=lambda x: f"{x}_{i}")
       for i in range(10)
   ]
   ts_wide = pd.concat(timeseries, axis=1)
   ts_wide.head()
   ts_wide.to_parquet("timeseries_wide.parquet")


Option 1 loads in all the data and then filters to what we need.

.. ipython:: python

   columns = ['id_0', 'name_0', 'x_0', 'y_0']
   
   %time _ = pd.read_parquet("timeseries_wide.parquet")[columns]

Option 2 only loads the columns we request. This is faster and has a lower peak
memory usage, since the entire dataset isn't in memory at once.

.. ipython:: python

   %time _ = pd.read_parquet("timeseries_wide.parquet", columns=columns)


With :func:`pandas.read_csv`, you can specify ``usecols`` to limit the columns
read into memory.


Use efficient datatypes
-----------------------

The default pandas data types are not the most memory efficient. This is
especially true for high-cardinality text data (columns with relatively few
unique values). By using more efficient data types you can store larger datasets
in memory.

.. ipython:: python

   ts.dtypes

.. ipython:: python

   ts.memory_usage(deep=True)  # memory usage in bytes


The ``name`` column is taking up much more memory than any other. It has just a
few unique values, so it's a good candidate for converting to a
:class:`Categorical`. With a Categorical, we store each unique name once and use
space-efficient integers to know which specific name is used in each row.


.. ipython:: python

   ts2 = ts.copy()
   ts2['name'] = ts2['name'].astype('category')
   ts2.memory_usage(deep=True)

We can go a bit further and downcast the numeric columns to their smallest types
using :func:`pandas.to_numeric`.

.. ipython:: python

   ts2['id'] = pd.to_numeric(ts2['id'], downcast='unsigned')
   ts2[['x', 'y']] = ts2[['x', 'y']].apply(pd.to_numeric, downcast='float')
   ts2.dtypes

.. ipython:: python

   ts2.memory_usage(deep=True)

.. ipython:: python

   reduction = ts2.memory_usage(deep=True).sum() / ts.memory_usage(deep=True).sum()
   print(f"{reduction:0.2f}")

In all, we've reduced the in-memory footprint of this dataset to 1/5 of its
original size.

See :ref:`categorical` for more on ``Categorical`` and :ref:`basics.dtypes`
for an overview of all of pandas' dtypes.

Use Other libraries
-------------------

Pandas is just one library offering a DataFrame API. Because of its popularity,
pandas' API has become something of a standard that other libraries implement.

For example, `Dask`_, a parallel computing library, has `dask.dataframe`_, a
pandas-like API for working with larger than memory datasets in parallel. Dask
can use multiple threads or processes on a single machine, or a cluster of
machines to process data in parallel.

Let's make a larger dataset on disk (as parquet files) that's split into chunks,
one per year.

.. ipython:: python
             
   import pathlib
   
   N = 12
   starts = [f'20{i:>02d}-01-01' for i in range(N)]
   ends = [f'20{i:>02d}-12-13' for i in range(N)]
   
   pathlib.Path("data/timeseries").mkdir(exist_ok=True)
   
   for i, (start, end) in enumerate(zip(starts, ends)):
       ts = make_timeseries(start=start, end=end, freq='1T', seed=i)
       ts.to_parquet(f"data/timeseries/ts-{i}.parquet")
   
We'll import ``dask.dataframe`` and notice that the API feels similar to pandas.
We can use Dask's ``read_parquet`` function, but provide a globstring of files to read in.

.. ipython:: python

   import dask.dataframe as dd

   ddf = dd.read_parquet("data/timeseries/ts*.parquet", engine="pyarrow")
   ddf

Inspecting the ``ddf`` object, we see a few things

* There are familiar attributes like ``.columns`` and ``.dtypes``
* There are familiar methods like ``.groupby``, ``.sum``, etc.
* There are new attributes like ``.npartitions`` and ``.divisions``

The partitions and divisions are how Dask parallizes computation. A **Dask**
DataFrame is made up of many **Pandas** DataFrames. A single method call on a
Dask DataFrame ends up making many pandas method calls, and Dask knows how to
coordinate everything to get the result.

.. ipython:: python

   ddf.columns
   ddf.dtypes
   ddf.npartitions

One major difference: the ``dask.dataframe`` API is *lazy*. If you look at the
repr above, you'll notice that the values aren't actually printed out; just the
column names and dtypes. That's because Dask hasn't actually read the data yet.
Rather than executing immediately, doing operations build up a **task graph**.

.. ipython:: python

   ddf
   ddf['name']
   ddf['name'].value_counts()

Each of these calls is instant because the result isn't being computed yet.
We're just building up a list of computation to do when someone needs the
result. Dask knows that the return type of a ``pandas.Series.value_counts``
is a pandas Series with a certain dtype and a certain name. So the Dask version
returns a Dask Series with the same dtype and the same name.

To get the actual result you can call ``.compute()``.
             
.. ipython:: python
             
   %time ddf['name'].value_counts().compute()

At that point, you get back the same thing you'd get with pandas, in this case
a concrete pandas Series with the count of each ``name``.

Calling ``.compute`` causes the full task graph to be executed. This includes
reading the data, selecting the columns, and doing the ``value_counts``. The
execution is done *in parallel* where possible, and Dask tries to keep the
overall memory footprint small. You can work with datasets that are much larger
than memory, as long as each partition (a regular pandas DataFrame) fits in memory.

By default, ``dask.dataframe`` operations use a threadpool to do operations in
parallel. We can also connect to a cluster to distribute the work on many
machines. In this case we'll connect to a local "cluster" made up of several
processes on this single machine.

.. ipython:: python

   from dask.distributed import Client, LocalCluster
   
   cluster = LocalCluster()
   client = Client(cluster)
   client

Once this ``client`` is created, all of Dask's computation will take place on
the cluster (which is just processes in this case).

Dask implements the most used parts of the pandas API. For example, we can do
a familiar groupby aggregation.

.. ipython:: python

   %time ddf.groupby('name')[['x', 'y']].mean().compute().head()

The grouping and aggregation is done out-of-core and in parallel.

When Dask knows the ``divisions`` of a dataset, certain optimizations are
possible. When reading parquet datasets written by dask, the divisions will be
known automatically. In this case, since we created the parquet files manually,
we need to supply the divisions manually.

.. ipython:: python

   divisions = tuple(pd.to_datetime(starts)) + (pd.Timestamp(ends[-1]),)
   ddf.divisions = divisions
   ddf

Now we can do things like fast random access with ``.loc``.

.. ipython:: python

   ddf.loc['2002-01-01 12:01':'2002-01-01 12:05'].compute()

Dask knows to just look in the 3rd partition for selecting values in `2002`. It
doesn't need to look at any other data.

Many workflows involve a large amount of data and processing it in a way that
reduces the size to something that fits in memory. In this case, we'll resample
to daily frequency and take the mean. Once we've taken the mean, we know the
results will fit in memory, so we can safely call ``compute`` without running
out of memory. At that point it's just a regular pandas object.

.. ipython:: python

   @savefig dask_resample.png
   ddf[['x', 'y']].resample("1D").mean().cumsum().compute().plot()

These Dask examples have all be done using multiple processes on a single
machine. Dask can be `deployed on a cluster
<https://docs.dask.org/en/latest/setup.html>`_ to scale up to even larger
datasets.

You see more dask examples at https://examples.dask.org.

Use chunking
------------

If using another library like Dask is not an option, you can achieve similar
results with a bit of work.

For example, we can recreate the out-of-core ``value_counts`` we did earlier
with Dask. The peak memory usage of this will be the size of the single largest
DataFrame.

.. ipython:: python

   files = list(pathlib.Path("data/timeseries/").glob("ts*.parquet"))
   files

.. ipython:: python

   %%time
   counts = pd.Series(dtype=int)
   for path in files:
       df = pd.read_parquet(path)
       counts = counts.add(df['name'].value_counts(), fill_value=0)
   counts.astype(int)

This matches the counts we saw above with Dask.

Some readers, like :meth:`pandas.read_csv` offer parameters to control the
``chunksize``. Manually chunking is an OK option for workflows that don't
require too sophisticated of operations. Some operations, like ``groupby``, are
much harder to do chunkwise. In these cases, you may be better switching to a
library like Dask, which implements these chunked algorithms for you.

.. ipython:: python

   del client, cluster

.. _Dask: https://dask.org
.. _dask.dataframe: https://docs.dask.org/en/latest/dataframe.html
