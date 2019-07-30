.. _developer:

{{ header }}

.. currentmodule:: pandas

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
for each column, *including the index columns*. This has JSON form:

.. code-block:: text

   {'name': column_name,
    'field_name': parquet_column_name,
    'pandas_type': pandas_type,
    'numpy_type': numpy_type,
    'metadata': metadata}

.. note::

   Every index column is stored with a name matching the pattern
   ``__index_level_\d+__`` and its corresponding column information is can be
   found with the following code snippet.

   Following this naming convention isn't strictly necessary, but strongly
   suggested for compatibility with Arrow.

   Here's an example of how the index metadata is structured in pyarrow:

    .. code-block:: python

       # assuming there's at least 3 levels in the index
       index_columns = metadata['index_columns']  # noqa: F821
       columns = metadata['columns']  # noqa: F821
       ith_index = 2
       assert index_columns[ith_index] == '__index_level_2__'
       ith_index_info = columns[-len(index_columns):][ith_index]
       ith_index_level_name = ith_index_info['name']

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
         'field_name': 'None',
         'pandas_type': 'unicode',
         'numpy_type': 'object',
         'metadata': {'encoding': 'UTF-8'}}
    ],
    'columns': [
        {'name': 'c0',
         'field_name': 'c0',
         'pandas_type': 'int8',
         'numpy_type': 'int8',
         'metadata': None},
        {'name': 'c1',
         'field_name': 'c1',
         'pandas_type': 'bytes',
         'numpy_type': 'object',
         'metadata': None},
        {'name': 'c2',
         'field_name': 'c2',
         'pandas_type': 'categorical',
         'numpy_type': 'int16',
         'metadata': {'num_categories': 1000, 'ordered': False}},
        {'name': 'c3',
         'field_name': 'c3',
         'pandas_type': 'datetimetz',
         'numpy_type': 'datetime64[ns]',
         'metadata': {'timezone': 'America/Los_Angeles'}},
        {'name': 'c4',
         'field_name': 'c4',
         'pandas_type': 'object',
         'numpy_type': 'object',
         'metadata': {'encoding': 'pickle'}},
        {'name': None,
         'field_name': '__index_level_0__',
         'pandas_type': 'int64',
         'numpy_type': 'int64',
         'metadata': None}
    ],
    'pandas_version': '0.20.0'}
