.. _io.parquet:

=======
Parquet
=======

`Apache Parquet <https://parquet.apache.org/>`__ provides a partitioned binary columnar serialization for data frames. It is designed to
make reading and writing data frames efficient, and to make sharing data across data analysis
languages easy. Parquet can use a variety of compression techniques to shrink the file size as much as possible
while still maintaining good read performance.

Parquet is designed to faithfully serialize and de-serialize ``DataFrame`` s, supporting all of the pandas
dtypes, including extension dtypes such as datetime with tz.

Several caveats.

* Duplicate column names and non-string columns names are not supported.
* The ``pyarrow`` engine always writes the index to the output, but ``fastparquet`` only writes non-default
  indexes. This extra column can cause problems for non-pandas consumers that are not expecting it. You can
  force including or omitting indexes with the ``index`` argument, regardless of the underlying engine.
* Index level names, if specified, must be strings.
* In the ``pyarrow`` engine, categorical dtypes for non-string types can be serialized to parquet, but will de-serialize as their primitive dtype.
* The ``pyarrow`` engine preserves the ``ordered`` flag of categorical dtypes with string types. ``fastparquet`` does not preserve the ``ordered`` flag.
* Non supported types include ``Interval`` and actual Python object types. These will raise a helpful error message
  on an attempt at serialization. ``Period`` type is supported with pyarrow >= 0.16.0.
* The ``pyarrow`` engine preserves extension data types such as the nullable integer and string data
  type (requiring pyarrow >= 0.16.0, and requiring the extension type to implement the needed protocols,
  see the :ref:`extension types documentation <extending.extension.arrow>`).

You can specify an ``engine`` to direct the serialization. This can be one of ``pyarrow``, or ``fastparquet``, or ``auto``.
If the engine is NOT specified, then the ``pd.options.io.parquet.engine`` option is checked; if this is also ``auto``,
then ``pyarrow`` is tried, and falling back to ``fastparquet``.

See the documentation for `pyarrow <https://arrow.apache.org/docs/python/>`__ and `fastparquet <https://fastparquet.readthedocs.io/en/latest/>`__.

.. note::

   These engines are very similar and should read/write nearly identical parquet format files.
   ``pyarrow>=8.0.0`` supports timedelta data, ``fastparquet>=0.1.4`` supports timezone aware datetimes.
   These libraries differ by having different underlying dependencies (``fastparquet`` by using ``numba``, while ``pyarrow`` uses a c-library).

.. ipython:: python

   df = pd.DataFrame(
       {
           "a": list("abc"),
           "b": list(range(1, 4)),
           "c": np.arange(3, 6).astype("u1"),
           "d": np.arange(4.0, 7.0, dtype="float64"),
           "e": [True, False, True],
           "f": pd.date_range("20130101", periods=3),
           "g": pd.date_range("20130101", periods=3, tz="US/Eastern"),
           "h": pd.Categorical(list("abc")),
           "i": pd.Categorical(list("abc"), ordered=True),
       }
   )

   df
   df.dtypes

Write to a parquet file.

.. ipython:: python

   df.to_parquet("example_pa.parquet", engine="pyarrow")
   df.to_parquet("example_fp.parquet", engine="fastparquet")

Read from a parquet file.

.. ipython:: python

   result = pd.read_parquet("example_fp.parquet", engine="fastparquet")
   result = pd.read_parquet("example_pa.parquet", engine="pyarrow")

   result.dtypes

By setting the ``dtype_backend`` argument you can control the default dtypes used for the resulting DataFrame.

.. ipython:: python

   result = pd.read_parquet("example_pa.parquet", engine="pyarrow", dtype_backend="pyarrow")

   result.dtypes

.. note::

   Note that this is not supported for ``fastparquet``.


Read only certain columns of a parquet file.

.. ipython:: python

   result = pd.read_parquet(
       "example_fp.parquet",
       engine="fastparquet",
       columns=["a", "b"],
   )
   result = pd.read_parquet(
       "example_pa.parquet",
       engine="pyarrow",
       columns=["a", "b"],
   )
   result.dtypes


.. ipython:: python
   :suppress:

   os.remove("example_pa.parquet")
   os.remove("example_fp.parquet")


Handling indexes
''''''''''''''''

Serializing a ``DataFrame`` to parquet may include the implicit index as one or
more columns in the output file. Thus, this code:

.. ipython:: python

    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    df.to_parquet("test.parquet", engine="pyarrow")

creates a parquet file with *three* columns if you use ``pyarrow`` for serialization:
``a``, ``b``, and ``__index_level_0__``. If you're using ``fastparquet``, the
index `may or may not <https://fastparquet.readthedocs.io/en/latest/api.html#fastparquet.write>`_
be written to the file.

This unexpected extra column causes some databases like Amazon Redshift to reject
the file, because that column doesn't exist in the target table.

If you want to omit a dataframe's indexes when writing, pass ``index=False`` to
:func:`~pandas.DataFrame.to_parquet`:

.. ipython:: python

    df.to_parquet("test.parquet", index=False)

This creates a parquet file with just the two expected columns, ``a`` and ``b``.
If your ``DataFrame`` has a custom index, you won't get it back when you load
this file into a ``DataFrame``.

Passing ``index=True`` will *always* write the index, even if that's not the
underlying engine's default behavior.

.. ipython:: python
   :suppress:

   os.remove("test.parquet")


Partitioning Parquet files
''''''''''''''''''''''''''

Parquet supports partitioning of data based on the values of one or more columns.

.. ipython:: python

    df = pd.DataFrame({"a": [0, 0, 1, 1], "b": [0, 1, 0, 1]})
    df.to_parquet(path="test", engine="pyarrow", partition_cols=["a"], compression=None)

The ``path`` specifies the parent directory to which data will be saved.
The ``partition_cols`` are the column names by which the dataset will be partitioned.
Columns are partitioned in the order they are given. The partition splits are
determined by the unique values in the partition columns.
The above example creates a partitioned dataset that may look like:

.. code-block:: text

    test
    ├── a=0
    │   ├── 0bac803e32dc42ae83fddfd029cbdebc.parquet
    │   └──  ...
    └── a=1
        ├── e6ab24a4f45147b49b54a662f0c412a3.parquet
        └── ...

.. ipython:: python
   :suppress:

   from shutil import rmtree

   try:
       rmtree("test")
   except OSError:
       pass
