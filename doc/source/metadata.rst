.. _metadata:

.. currentmodule:: pandas

**********************************************
Storing pandas Objects in Various File Formats
**********************************************

This document provides specifications for metadata to assist with reading and
writing pandas objects to different third party file formats.

Apache Parquet
--------------

The Apache Parquet format provides key-value metadata at the file and column
level, stored in the footer of the Parquet file:

.. code-block:: shell

  5: optional list<KeyValue> key_value_metadata

where ``KeyValue`` is

.. code-block:: shell

   struct KeyValue {
     1: required string key
     2: optional string value
   }

So that a ``pandas.DataFrame`` can be faithfully reconstructed, we store a
``pandas`` metadata key in the ``FileMetaData`` with the the value stored as :

.. code-block:: text

   {'index_columns': ['__index_level_0__', '__index_level_1__', ...],
    'columns': [<c0>, <c1>, ...],
    'pandas_version': '0.20.0'}

Here, ``<c0>`` and so forth are dictionaries containing the metadata for each
column. This has JSON form:

.. code-block:: text

   {'name': column_name,
    'type': pandas_type,
    'numpy_dtype': numpy_type,
    'metadata': type_metadata}

``pandas_type`` is the logical type of the column, and is one of:

* Boolean: ``'bool'``
* Integers: ``'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'``
* Floats: ``'float16', 'float32', 'float64'``
* Datetime: ``'datetime', 'datetimetz'``
* String: ``'unicode', 'bytes'``
* Categorical: ``'categorical'``

The ``numpy_type`` is the physical storage type of the column, which is the
result of ``str(dtype)`` for the underlying NumPy array that holds the data. So
for ``datetimetz`` this is ``datetime64[ns]`` and for categorical, it may be
any of the supported integer categorical types.

The ``type_metadata`` is ``None`` except for:

* ``datetimetz``: ``{'timezone': zone}``, e.g. ``{'timezone', 'America/New_York'}``
* ``categorical``: ``{'num_categories': K}``

For types other than these, the ``'metadata'`` key can be
omitted. Implementations can assume ``None`` if the key is not present.

As an example of fully-formed metadata:

.. code-block:: text

   {'index_columns': ['__index_level_0__'],
    'columns': [
        {'name': 'c0',
         'type': 'int8',
         'numpy_type': 'int8',
         'metadata': None},
        {'name': 'c1',
         'type': 'bytes',
         'numpy_type': 'object',
         'metadata': None},
        {'name': 'c2',
         'type': 'categorical',
         'numpy_type': 'int16',
         'metadata': {'num_categories': 1000}},
        {'name': 'c3',
         'type': 'datetimetz',
         'numpy_type': 'datetime64[ns]',
         'metadata': {'timezone': 'America/Los_Angeles'}},
        {'name': '__index_level_0__',
         'type': 'int64',
         'numpy_type': 'int64',
         'metadata': None}
    ],
    'pandas_version': '0.20.0'}
