.. _scale:

*************************
Scaling to large datasets
*************************

pandas provides data structures for in-memory analytics, which makes using pandas
to analyze datasets that are larger than memory datasets somewhat tricky. Even datasets
that are a sizable fraction of memory become unwieldy, as some pandas operations need
to make intermediate copies.

This document provides a few recommendations for scaling your analysis to larger datasets.
It's a complement to :ref:`enhancingperf`, which focuses on speeding up analysis
for datasets that fit in memory.

Load less data
--------------

Suppose our raw dataset on disk has many columns.

.. ipython:: python
   :okwarning:

   import pandas as pd
   import numpy as np

   def make_timeseries(start="2000-01-01", end="2000-12-31", freq="1D", seed=None):
       index = pd.date_range(start=start, end=end, freq=freq, name="timestamp")
       n = len(index)
       state = np.random.RandomState(seed)
       columns = {
           "name": state.choice(["Alice", "Bob", "Charlie"], size=n),
           "id": state.poisson(1000, size=n),
           "x": state.rand(n) * 2 - 1,
           "y": state.rand(n) * 2 - 1,
       }
       df = pd.DataFrame(columns, index=index, columns=sorted(columns))
       if df.index[-1] == end:
           df = df.iloc[:-1]
       return df

   timeseries = [
       make_timeseries(freq="1min", seed=i).rename(columns=lambda x: f"{x}_{i}")
       for i in range(10)
   ]
   ts_wide = pd.concat(timeseries, axis=1)
   ts_wide.head()
   ts_wide.to_parquet("timeseries_wide.parquet")

To load the columns we want, we have two options.
Option 1 loads in all the data and then filters to what we need.

.. ipython:: python
   :okwarning:

   columns = ["id_0", "name_0", "x_0", "y_0"]

   pd.read_parquet("timeseries_wide.parquet")[columns]

Option 2 only loads the columns we request.

.. ipython:: python
   :okwarning:

   pd.read_parquet("timeseries_wide.parquet", columns=columns)

.. ipython:: python
   :suppress:

   import os

   os.remove("timeseries_wide.parquet")

If we were to measure the memory usage of the two calls, we'd see that specifying
``columns`` uses about 1/10th the memory in this case.

With :func:`pandas.read_csv`, you can specify ``usecols`` to limit the columns
read into memory. Not all file formats that can be read by pandas provide an option
to read a subset of columns.

Use efficient datatypes
-----------------------

The default pandas data types are not the most memory efficient. This is
especially true for text data columns with relatively few unique values (commonly
referred to as "low-cardinality" data). By using more efficient data types, you
can store larger datasets in memory.

.. ipython:: python
   :okwarning:

   ts = make_timeseries(freq="30s", seed=0)
   ts.to_parquet("timeseries.parquet")
   ts = pd.read_parquet("timeseries.parquet")
   ts

.. ipython:: python
   :suppress:

   os.remove("timeseries.parquet")

Now, let's inspect the data types and memory usage to see where we should focus our
attention.

.. ipython:: python

   ts.dtypes

.. ipython:: python

   ts.memory_usage(deep=True)  # memory usage in bytes


The ``name`` column is taking up much more memory than any other. It has just a
few unique values, so it's a good candidate for converting to a
:class:`pandas.Categorical`. With a :class:`pandas.Categorical`, we store each unique name once and use
space-efficient integers to know which specific name is used in each row.


.. ipython:: python

   ts2 = ts.copy()
   ts2["name"] = ts2["name"].astype("category")
   ts2.memory_usage(deep=True)

We can go a bit further and downcast the numeric columns to their smallest types
using :func:`pandas.to_numeric`.

.. ipython:: python

   ts2["id"] = pd.to_numeric(ts2["id"], downcast="unsigned")
   ts2[["x", "y"]] = ts2[["x", "y"]].apply(pd.to_numeric, downcast="float")
   ts2.dtypes

.. ipython:: python

   ts2.memory_usage(deep=True)

.. ipython:: python

   reduction = ts2.memory_usage(deep=True).sum() / ts.memory_usage(deep=True).sum()
   print(f"{reduction:0.2f}")

In all, we've reduced the in-memory footprint of this dataset to 1/5 of its
original size.

See :ref:`categorical` for more on :class:`pandas.Categorical` and :ref:`basics.dtypes`
for an overview of all of pandas' dtypes.

Use chunking
------------

Some workloads can be achieved with chunking by splitting a large problem into a bunch of small problems. For example,
converting an individual CSV file into a Parquet file and repeating that for each file in a directory. As long as each chunk
fits in memory, you can work with datasets that are much larger than memory.

.. note::

   Chunking works well when the operation you're performing requires zero or minimal
   coordination between chunks. For more complicated workflows, you're better off
   :ref:`using another library <scale.other_libraries>`.

Suppose we have an even larger "logical dataset" on disk that's a directory of parquet
files. Each file in the directory represents a different year of the entire dataset.

.. ipython:: python
   :okwarning:

   import pathlib

   N = 12
   starts = [f"20{i:>02d}-01-01" for i in range(N)]
   ends = [f"20{i:>02d}-12-13" for i in range(N)]

   pathlib.Path("data/timeseries").mkdir(exist_ok=True)

   for i, (start, end) in enumerate(zip(starts, ends)):
       ts = make_timeseries(start=start, end=end, freq="1min", seed=i)
       ts.to_parquet(f"data/timeseries/ts-{i:0>2d}.parquet")


::

   data
   └── timeseries
       ├── ts-00.parquet
       ├── ts-01.parquet
       ├── ts-02.parquet
       ├── ts-03.parquet
       ├── ts-04.parquet
       ├── ts-05.parquet
       ├── ts-06.parquet
       ├── ts-07.parquet
       ├── ts-08.parquet
       ├── ts-09.parquet
       ├── ts-10.parquet
       └── ts-11.parquet

Now we'll implement an out-of-core :meth:`pandas.Series.value_counts`. The peak memory usage of this
workflow is the single largest chunk, plus a small series storing the unique value
counts up to this point. As long as each individual file fits in memory, this will
work for arbitrary-sized datasets.

.. ipython:: python
   :okwarning:

   %%time
   files = pathlib.Path("data/timeseries/").glob("ts*.parquet")
   counts = pd.Series(dtype=int)
   for path in files:
       df = pd.read_parquet(path)
       counts = counts.add(df["name"].value_counts(), fill_value=0)
   counts.astype(int)

Some readers, like :meth:`pandas.read_csv`, offer parameters to control the
``chunksize`` when reading a single file.

Manually chunking is an OK option for workflows that don't
require too sophisticated of operations. Some operations, like :meth:`pandas.DataFrame.groupby`, are
much harder to do chunkwise. In these cases, you may be better switching to a
different library that implements these out-of-core algorithms for you.

.. _scale.other_libraries:

Use Dask
--------

pandas is just one library offering a DataFrame API. Because of its popularity,
pandas' API has become something of a standard that other libraries implement.
The pandas documentation maintains a list of libraries implementing a DataFrame API
in `the ecosystem page <https://pandas.pydata.org/community/ecosystem.html>`_.

For example, `Dask`_, a parallel computing library, has `dask.dataframe`_, a
pandas-like API for working with larger than memory datasets in parallel. Dask
can use multiple threads or processes on a single machine, or a cluster of
machines to process data in parallel.


We'll import ``dask.dataframe`` and notice that the API feels similar to pandas.
We can use Dask's ``read_parquet`` function, but provide a globstring of files to read in.

.. ipython:: python
   :okwarning:

   import dask.dataframe as dd

   ddf = dd.read_parquet("data/timeseries/ts*.parquet", engine="pyarrow")
   ddf

Inspecting the ``ddf`` object, we see a few things

* There are familiar attributes like ``.columns`` and ``.dtypes``
* There are familiar methods like ``.groupby``, ``.sum``, etc.
* There are new attributes like ``.npartitions`` and ``.divisions``

The partitions and divisions are how Dask parallelizes computation. A **Dask**
DataFrame is made up of many pandas :class:`pandas.DataFrame`. A single method call on a
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
   :okwarning:

   ddf
   ddf["name"]
   ddf["name"].value_counts()

Each of these calls is instant because the result isn't being computed yet.
We're just building up a list of computation to do when someone needs the
result. Dask knows that the return type of a :class:`pandas.Series.value_counts`
is a pandas :class:`pandas.Series` with a certain dtype and a certain name. So the Dask version
returns a Dask Series with the same dtype and the same name.

To get the actual result you can call ``.compute()``.

.. ipython:: python
   :okwarning:

   %time ddf["name"].value_counts().compute()

At that point, you get back the same thing you'd get with pandas, in this case
a concrete pandas :class:`pandas.Series` with the count of each ``name``.

Calling ``.compute`` causes the full task graph to be executed. This includes
reading the data, selecting the columns, and doing the ``value_counts``. The
execution is done *in parallel* where possible, and Dask tries to keep the
overall memory footprint small. You can work with datasets that are much larger
than memory, as long as each partition (a regular pandas :class:`pandas.DataFrame`) fits in memory.

By default, ``dask.dataframe`` operations use a threadpool to do operations in
parallel. We can also connect to a cluster to distribute the work on many
machines. In this case we'll connect to a local "cluster" made up of several
processes on this single machine.

.. code-block:: python

   >>> from dask.distributed import Client, LocalCluster

   >>> cluster = LocalCluster()
   >>> client = Client(cluster)
   >>> client
   <Client: 'tcp://127.0.0.1:53349' processes=4 threads=8, memory=17.18 GB>

Once this ``client`` is created, all of Dask's computation will take place on
the cluster (which is just processes in this case).

Dask implements the most used parts of the pandas API. For example, we can do
a familiar groupby aggregation.

.. ipython:: python
   :okwarning:

   %time ddf.groupby("name")[["x", "y"]].mean().compute().head()

The grouping and aggregation is done out-of-core and in parallel.

When Dask knows the ``divisions`` of a dataset, certain optimizations are
possible. When reading parquet datasets written by dask, the divisions will be
known automatically. In this case, since we created the parquet files manually,
we need to supply the divisions manually.

.. ipython:: python
   :okwarning:

   N = 12
   starts = [f"20{i:>02d}-01-01" for i in range(N)]
   ends = [f"20{i:>02d}-12-13" for i in range(N)]

   divisions = tuple(pd.to_datetime(starts)) + (pd.Timestamp(ends[-1]),)
   ddf.divisions = divisions
   ddf

Now we can do things like fast random access with ``.loc``.

.. ipython:: python
   :okwarning:

   ddf.loc["2002-01-01 12:01":"2002-01-01 12:05"].compute()

Dask knows to just look in the 3rd partition for selecting values in 2002. It
doesn't need to look at any other data.

Many workflows involve a large amount of data and processing it in a way that
reduces the size to something that fits in memory. In this case, we'll resample
to daily frequency and take the mean. Once we've taken the mean, we know the
results will fit in memory, so we can safely call ``compute`` without running
out of memory. At that point it's just a regular pandas object.

.. ipython:: python
   :okwarning:

   @savefig dask_resample.png
   ddf[["x", "y"]].resample("1D").mean().cumsum().compute().plot()

.. ipython:: python
   :suppress:

   import shutil

   shutil.rmtree("data/timeseries")

These Dask examples have all be done using multiple processes on a single
machine. Dask can be `deployed on a cluster
<https://docs.dask.org/en/latest/setup.html>`_ to scale up to even larger
datasets.

You see more dask examples at https://examples.dask.org.

.. _Dask: https://dask.org
.. _dask.dataframe: https://docs.dask.org/en/latest/dataframe.html
