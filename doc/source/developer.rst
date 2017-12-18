.. _developer:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.options.display.max_rows = 15

*********
Developer
*********

This section will focus on downstream applications of pandas.

.. _apache.parquet:

Storing pandas DataFrame objects in Apache Parquet format
---------------------------------------------------------

The `Apache Parquet <https://github.com/apache/parquet-format>`__ format
provides key-value metadata at the file and column level, stored in the footer
of the Parquet file:

.. code-block:: shell

  5: optional list<KeyValue> key_value_metadata

where ``KeyValue`` is

.. code-block:: shell

   struct KeyValue {
     1: required string key
     2: optional string value
   }

So that a ``pandas.DataFrame`` can be faithfully reconstructed, we store a
``pandas`` metadata key in the ``FileMetaData`` with the value stored as :

.. code-block:: text

   {'index_columns': ['__index_level_0__', '__index_level_1__', ...],
    'column_indexes': [<ci0>, <ci1>, ..., <ciN>],
    'columns': [<c0>, <c1>, ...],
    'pandas_version': $VERSION}

Here, ``<c0>``/``<ci0>`` and so forth are dictionaries containing the metadata
for each column. This has JSON form:

.. code-block:: text

   {'name': column_name,
    'pandas_type': pandas_type,
    'numpy_type': numpy_type,
    'metadata': metadata}

``pandas_type`` is the logical type of the column, and is one of:

* Boolean: ``'bool'``
* Integers: ``'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'``
* Floats: ``'float16', 'float32', 'float64'``
* Date and Time Types: ``'datetime', 'datetimetz'``, ``'timedelta'``
* String: ``'unicode', 'bytes'``
* Categorical: ``'categorical'``
* Other Python objects: ``'object'``

The ``numpy_type`` is the physical storage type of the column, which is the
result of ``str(dtype)`` for the underlying NumPy array that holds the data. So
for ``datetimetz`` this is ``datetime64[ns]`` and for categorical, it may be
any of the supported integer categorical types.

The ``metadata`` field is ``None`` except for:

* ``datetimetz``: ``{'timezone': zone, 'unit': 'ns'}``, e.g. ``{'timezone',
  'America/New_York', 'unit': 'ns'}``. The ``'unit'`` is optional, and if
  omitted it is assumed to be nanoseconds.
* ``categorical``: ``{'num_categories': K, 'ordered': is_ordered, 'type': $TYPE}``

  * Here ``'type'`` is optional, and can be a nested pandas type specification
    here (but not categorical)

* ``unicode``: ``{'encoding': encoding}``

  * The encoding is optional, and if not present is UTF-8

* ``object``: ``{'encoding': encoding}``. Objects can be serialized and stored
  in ``BYTE_ARRAY`` Parquet columns. The encoding can be one of:

  * ``'pickle'``
  * ``'msgpack'``
  * ``'bson'``
  * ``'json'``

* ``timedelta``: ``{'unit': 'ns'}``. The ``'unit'`` is optional, and if omitted
  it is assumed to be nanoseconds. This metadata is optional altogether

For types other than these, the ``'metadata'`` key can be
omitted. Implementations can assume ``None`` if the key is not present.

As an example of fully-formed metadata:

.. code-block:: text

   {'index_columns': ['__index_level_0__'],
    'column_indexes': [
        {'name': None,
         'pandas_type': 'string',
         'numpy_type': 'object',
         'metadata': None}
    ],
    'columns': [
        {'name': 'c0',
         'pandas_type': 'int8',
         'numpy_type': 'int8',
         'metadata': None},
        {'name': 'c1',
         'pandas_type': 'bytes',
         'numpy_type': 'object',
         'metadata': None},
        {'name': 'c2',
         'pandas_type': 'categorical',
         'numpy_type': 'int16',
         'metadata': {'num_categories': 1000, 'ordered': False}},
        {'name': 'c3',
         'pandas_type': 'datetimetz',
         'numpy_type': 'datetime64[ns]',
         'metadata': {'timezone': 'America/Los_Angeles'}},
        {'name': 'c4',
         'pandas_type': 'object',
         'numpy_type': 'object',
         'metadata': {'encoding': 'pickle'}},
        {'name': '__index_level_0__',
         'pandas_type': 'int64',
         'numpy_type': 'int64',
         'metadata': None}
    ],
    'pandas_version': '0.20.0'}

.. _developer.custom-array-types:

Custom Array Types
------------------

.. versionadded:: 0.23.0

.. warning::
   Support for custom array types is experimental.

Sometimes the NumPy type system isn't rich enough for your needs. Pandas has
made a few extensions internally (e.g. ``Categorical``). While this has worked
well for pandas, not all custom data types belong in pandas itself.

Pandas defines an interface for custom arrays. Arrays implementing this
interface will be stored correctly in ``Series`` or ``DataFrame``. The ABCs
that must be implemented are

1. :class:`ExtensionDtype` A class describing your data type itself. This is
   similar to a ``numpy.dtype``.
2. :class:`ExtensionArray`: A container for your data.

Throughout this document, we'll use the example of storing IPv6 addresses. An
IPv6 address is 128 bits, so NumPy doesn't have a native data type for it. We'll
model it as a structured array with two ``uint64`` fields, which together
represent the 128-bit integer that is the IP Address.

Extension Dtype
'''''''''''''''

This should describe your data type. The most important fields are ``name`` and
``base``:

.. code-block:: python

   class IPv6Type(ExtensionDtype):
       name = 'IPv6'
       base = np.dtype([('hi', '>u8'), ('lo', '>u8')])
       type = IPTypeType
       kind = 'O'
       fill_value = np.array([(0, 0)], dtype=base)

``base`` describe the underlying storage of individual items in your array.
TODO: is this true? Or does ``.base`` refer to the original memory this
is a view on? Different meanings for ``np.dtype.base`` vs. ``np.ndarray.base``?

In our IPAddress case, we're using a NumPy structured array with two fields.

Extension Array
'''''''''''''''

This is the actual array container for your data, similar to a ``Categorical``,
and requires the most work to implement correctly. *pandas makes no assumptions
about how you store the data*. You're free to use NumPy arrays or PyArrow
arrays, or even just Python lists. That said, several of the methods required by
the interface expect NumPy arrays as the return value.

* ``dtype``: Should be an *instance* of your custom ``ExtensionType``
* ``formtting_values(self)``: Used for printing Series and DataFrame
* ``concat_same_type(concat)``: Used in :func:`pd.concat`
