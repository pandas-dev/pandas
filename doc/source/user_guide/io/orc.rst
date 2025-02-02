.. _io.orc:

===
ORC
===

Similar to the :ref:`parquet <io.parquet>` format, the `ORC Format <https://orc.apache.org/>`__ is a binary columnar serialization
for data frames. It is designed to make reading data frames efficient. pandas provides both the reader and the writer for the
ORC format, :func:`~pandas.read_orc` and :func:`~pandas.DataFrame.to_orc`. This requires the `pyarrow <https://arrow.apache.org/docs/python/>`__ library.

.. warning::

   * It is *highly recommended* to install pyarrow using conda due to some issues occurred by pyarrow.
   * :func:`~pandas.DataFrame.to_orc` requires pyarrow>=7.0.0.
   * :func:`~pandas.read_orc` and :func:`~pandas.DataFrame.to_orc` are not supported on Windows yet, you can find valid environments on :ref:`install optional dependencies <install.warn_orc>`.
   * For supported dtypes please refer to `supported ORC features in Arrow <https://arrow.apache.org/docs/cpp/orc.html#data-types>`__.
   * Currently timezones in datetime columns are not preserved when a dataframe is converted into ORC files.

.. ipython:: python

   df = pd.DataFrame(
       {
           "a": list("abc"),
           "b": list(range(1, 4)),
           "c": np.arange(4.0, 7.0, dtype="float64"),
           "d": [True, False, True],
           "e": pd.date_range("20130101", periods=3),
       }
   )

   df
   df.dtypes

Write to an orc file.

.. ipython:: python

   df.to_orc("example_pa.orc", engine="pyarrow")

Read from an orc file.

.. ipython:: python

   result = pd.read_orc("example_pa.orc")

   result.dtypes

Read only certain columns of an orc file.

.. ipython:: python

   result = pd.read_orc(
       "example_pa.orc",
       columns=["a", "b"],
   )
   result.dtypes


.. ipython:: python
   :suppress:

   os.remove("example_pa.orc")
