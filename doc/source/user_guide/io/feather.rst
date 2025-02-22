.. _io.feather:

=======
Feather
=======

Feather provides binary columnar serialization for data frames. It is designed to make reading and writing data
frames efficient, and to make sharing data across data analysis languages easy.

Feather is designed to faithfully serialize and de-serialize DataFrames, supporting all of the pandas
dtypes, including extension dtypes such as categorical and datetime with tz.

Several caveats:

* The format will NOT write an ``Index``, or ``MultiIndex`` for the
  ``DataFrame`` and will raise an error if a non-default one is provided. You
  can ``.reset_index()`` to store the index or ``.reset_index(drop=True)`` to
  ignore it.
* Duplicate column names and non-string columns names are not supported
* Actual Python objects in object dtype columns are not supported. These will
  raise a helpful error message on an attempt at serialization.

See the `Full Documentation <https://github.com/wesm/feather>`__.

.. ipython:: python

   import pytz

   df = pd.DataFrame(
       {
           "a": list("abc"),
           "b": list(range(1, 4)),
           "c": np.arange(3, 6).astype("u1"),
           "d": np.arange(4.0, 7.0, dtype="float64"),
           "e": [True, False, True],
           "f": pd.Categorical(list("abc")),
           "g": pd.date_range("20130101", periods=3),
           "h": pd.date_range("20130101", periods=3, tz=pytz.timezone("US/Eastern")),
           "i": pd.date_range("20130101", periods=3, freq="ns"),
       }
   )

   df
   df.dtypes

Write to a feather file.

.. ipython:: python
   :okwarning:

   df.to_feather("example.feather")

Read from a feather file.

.. ipython:: python
   :okwarning:

   result = pd.read_feather("example.feather")
   result

   # we preserve dtypes
   result.dtypes

.. ipython:: python
   :suppress:

   import os
   os.remove("example.feather")
