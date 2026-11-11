.. _pyarrow:

{{ header }}

*********************
PyArrow Functionality
*********************

pandas can utilize `PyArrow <https://arrow.apache.org/docs/python/index.html>`__ to extend functionality and improve the performance
of various APIs. This includes:

* More extensive `data types <https://arrow.apache.org/docs/python/api/datatypes.html>`__ compared to NumPy
* Missing data support (NA) for all data types
* Performant IO reader integration
* Facilitate interoperability with other dataframe libraries based on the Apache Arrow specification (e.g. polars, cuDF)

To use this functionality, please ensure you have :ref:`installed the minimum supported PyArrow version. <install.optional_dependencies>`


Data Structure Integration
--------------------------

A :class:`Series`, :class:`Index`, or the columns of a :class:`DataFrame` can be directly backed by a :external+pyarrow:py:class:`pyarrow.ChunkedArray`
which is similar to a NumPy array. To construct these from the main pandas data structures, you can pass in a string of the type followed by
``[pyarrow]``, e.g. ``"int64[pyarrow]"`` into the ``dtype`` parameter

.. ipython:: python

   ser = pd.Series([-1.5, 0.2, None], dtype="float32[pyarrow]")
   ser

   idx = pd.Index([True, None], dtype="bool[pyarrow]")
   idx

   df = pd.DataFrame([[1, 2], [3, 4]], dtype="uint64[pyarrow]")
   df

.. note::

   The string alias ``"string[pyarrow]"`` maps to ``pd.StringDtype("pyarrow")`` which is not equivalent to
   specifying ``dtype=pd.ArrowDtype(pa.string())``. Generally, operations on the data will behave similarly
   except ``pd.StringDtype("pyarrow")`` can return NumPy-backed nullable types while ``pd.ArrowDtype(pa.string())``
   will return :class:`ArrowDtype`.

   .. ipython:: python

      import pyarrow as pa
      data = list("abc")
      ser_sd = pd.Series(data, dtype="string[pyarrow]")
      ser_ad = pd.Series(data, dtype=pd.ArrowDtype(pa.string()))
      ser_ad.dtype == ser_sd.dtype
      ser_sd.str.contains("a")
      ser_ad.str.contains("a")

For PyArrow types that accept parameters, you can pass in a PyArrow type with those parameters
into :class:`ArrowDtype` to use in the ``dtype`` parameter.

.. ipython:: python

   import pyarrow as pa
   list_str_type = pa.list_(pa.string())
   ser = pd.Series([["hello"], ["there"]], dtype=pd.ArrowDtype(list_str_type))
   ser

.. ipython:: python

   from datetime import time
   idx = pd.Index([time(12, 30), None], dtype=pd.ArrowDtype(pa.time64("us")))
   idx

.. ipython:: python

   from decimal import Decimal
   decimal_type = pd.ArrowDtype(pa.decimal128(3, scale=2))
   data = [[Decimal("3.19"), None], [None, Decimal("-1.23")]]
   df = pd.DataFrame(data, dtype=decimal_type)
   df

If you already have an :external+pyarrow:py:class:`pyarrow.Array` or :external+pyarrow:py:class:`pyarrow.ChunkedArray`,
you can pass it into :class:`.arrays.ArrowExtensionArray` to construct the associated :class:`Series`, :class:`Index`
or :class:`DataFrame` object.

.. ipython:: python

   pa_array = pa.array(
       [{"1": "2"}, {"10": "20"}, None],
       type=pa.map_(pa.string(), pa.string()),
   )
   ser = pd.Series(pd.arrays.ArrowExtensionArray(pa_array))
   ser

To retrieve a pyarrow :external+pyarrow:py:class:`pyarrow.ChunkedArray` from a :class:`Series` or :class:`Index`, you can call
the pyarrow array constructor on the :class:`Series` or :class:`Index`.

.. ipython:: python

   ser = pd.Series([1, 2, None], dtype="uint8[pyarrow]")
   pa.array(ser)

   idx = pd.Index(ser)
   pa.array(idx)

To convert a :external+pyarrow:py:class:`pyarrow.Table` to a :class:`DataFrame`, you can call the
:external+pyarrow:py:meth:`pyarrow.Table.to_pandas` method with ``types_mapper=pd.ArrowDtype``.

.. ipython:: python

   table = pa.table([pa.array([1, 2, 3], type=pa.int64())], names=["a"])

   df = table.to_pandas(types_mapper=pd.ArrowDtype)
   df
   df.dtypes


Operations
----------

PyArrow data structure integration is implemented through pandas' :class:`~pandas.api.extensions.ExtensionArray` :ref:`interface <extending.extension-types>`;
therefore, supported functionality exists where this interface is integrated within the pandas API. Additionally, this functionality
is accelerated with PyArrow `compute functions <https://arrow.apache.org/docs/python/api/compute.html>`__ where available. This includes:

* Numeric aggregations
* Numeric arithmetic
* Numeric rounding
* Logical and comparison functions
* String functionality
* Datetime functionality

The following are just some examples of operations that are accelerated by native PyArrow compute functions.

.. ipython:: python

   import pyarrow as pa
   ser = pd.Series([-1.545, 0.211, None], dtype="float32[pyarrow]")
   ser.mean()
   ser + ser
   ser > (ser + 1)

   ser.dropna()
   ser.isna()
   ser.fillna(0)

.. ipython:: python

   ser_str = pd.Series(["a", "b", None], dtype=pd.ArrowDtype(pa.string()))
   ser_str.str.startswith("a")

.. ipython:: python

   from datetime import datetime
   pa_type = pd.ArrowDtype(pa.timestamp("ns"))
   ser_dt = pd.Series([datetime(2022, 1, 1), None], dtype=pa_type)
   ser_dt.dt.strftime("%Y-%m")

.. _pyarrow.deferred_take:

Deferred row selection of string and binary data
------------------------------------------------

.. versionadded:: 3.2.0

Selecting rows of string or binary data copies every byte of every selected
value. In merges and filters of wide frames, that copy often dominates the
run time, and many of the copied columns are dropped or never used afterwards.
For PyArrow-backed string and binary data, pandas therefore defers it: the
result records the original array and which of its rows are selected, and
copies the values (gathers them) the first time they are needed.

The result behaves exactly like one gathered right away. Deferral changes only
when the copy happens, and so when its time and memory are spent.

**Which data.** Arrays with a PyArrow ``string``, ``large_string``, ``binary``
or ``large_binary`` type: the default ``str`` dtype, ``"string[pyarrow]"``, and
:class:`ArrowDtype` with one of those types. All other data, including
NumPy-backed, ``object`` and ``StringDtype("python")`` columns and other PyArrow
types, is copied right away.

**What is deferred.** Any selection of existing rows by position or by boolean
mask. This includes:

* :func:`merge` and :meth:`DataFrame.join`, for the rows that have a match
* boolean indexing, such as ``df[mask]``, ``ser[mask]``, ``df.loc[mask]`` and
  :meth:`DataFrame.query`
* :meth:`DataFrame.take`, and ``iloc`` and ``loc`` with a list of positions or
  existing labels
* :meth:`DataFrame.reindex` to existing labels
* methods that reorder or drop rows, such as :meth:`DataFrame.sort_values`,
  :meth:`DataFrame.sort_index`, :meth:`DataFrame.sample`,
  :meth:`DataFrame.nlargest`, :meth:`DataFrame.dropna` and
  :meth:`DataFrame.drop_duplicates` (for the columns not used to find the rows)

Selecting rows of a deferred result, by position, mask or slice, is also
deferred and is combined with the earlier selection, so the values are copied
once, straight from the original array.

A selection that introduces missing values, such as the unmatched rows of a
``"left"``, ``"right"`` or ``"outer"`` merge or :meth:`DataFrame.reindex` to new
labels, copies right away. So does a selection of no rows.

**When the values are gathered.** The first time an operation reads them. For
example:

* computations: comparisons, :meth:`Series.isna`, the ``.str`` methods,
  :meth:`Series.unique`, and grouping by or aggregating the column
* setting values, :func:`concat`, and :meth:`DataFrame.memory_usage` or
  :meth:`DataFrame.info`
* converting to other objects: :meth:`Series.to_numpy`, iteration,
  ``pyarrow.array(ser)``, writing to a file, and pickling

These operations do not gather:

* the length, shape and dtype
* selecting columns, and adding, renaming or dropping other columns
* reading a single value, such as ``ser.iloc[0]`` or ``df.at[label, "col"]``
* selecting rows again
* displaying a large frame, which gathers only the rows it shows

Views, such as ``df["col"]``, share the deferred selection, so they gather once
between them and share the result. :meth:`DataFrame.copy` gathers the values
into new memory and leaves the original object deferred.

**Memory.** Until its values are gathered, a deferred result keeps the original
array in memory, as a slice does. If you select a few rows of a large frame
and delete the frame, its string data stays in memory until the selected
columns are used. Call :meth:`DataFrame.copy` on the selection to release it.

**Errors.** An error raised by the gather, such as a ``MemoryError``, is raised
by the operation that first reads the values, not by the selection.

I/O Reading
-----------

PyArrow also provides IO reading functionality that has been integrated into several pandas IO readers. The following
functions provide an ``engine`` keyword that can dispatch to PyArrow to accelerate reading from an IO source.

* :func:`read_csv`
* :func:`read_feather`
* :func:`read_json`
* :func:`read_orc`
* :func:`read_parquet`
* :func:`read_table` (experimental)

.. ipython:: python

   import io
   data = io.StringIO("""a,b,c
      1,2.5,True
      3,4.5,False
   """)
   df = pd.read_csv(data, engine="pyarrow")
   df

By default, these functions and all other IO reader functions return NumPy-backed data. These readers can return
PyArrow-backed data by specifying the parameter ``dtype_backend="pyarrow"``. A reader does not need to set
``engine="pyarrow"`` to necessarily return PyArrow-backed data.

.. ipython:: python

    import io
    data = io.StringIO("""a,b,c,d,e,f,g,h,i
        1,2.5,True,a,,,,,
        3,4.5,False,b,6,7.5,True,a,
    """)
    df_pyarrow = pd.read_csv(data, dtype_backend="pyarrow")
    df_pyarrow.dtypes

Several non-IO reader functions can also use the ``dtype_backend`` argument to return PyArrow-backed data including:

* :func:`to_numeric`
* :meth:`DataFrame.convert_dtypes`
* :meth:`Series.convert_dtypes`
