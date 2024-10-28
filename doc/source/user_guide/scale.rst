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

   columns = ["id_0", "name_0", "x_0", "y_0"]

   pd.read_parquet("timeseries_wide.parquet")[columns]

Option 2 only loads the columns we request.

.. ipython:: python

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
   :ref:`using other libraries <scale.other_libraries>`.

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

Use Other Libraries
-------------------

There are other libraries which provide similar APIs to pandas and work nicely with pandas DataFrame,
and can give you the ability to scale your large dataset processing and analytics
by parallel runtime, distributed memory, clustering, etc. You can find more information
in `the ecosystem page <https://pandas.pydata.org/community/ecosystem.html#out-of-core>`_.
