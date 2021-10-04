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

   {'index_columns': [<descr0>, <descr1>, ...],
    'column_indexes': [<ci0>, <ci1>, ..., <ciN>],
    'columns': [<c0>, <c1>, ...],
    'pandas_version': $VERSION,
    'creator': {
      'library': $LIBRARY,
      'version': $LIBRARY_VERSION
    }}

The "descriptor" values ``<descr0>`` in the ``'index_columns'`` field are
strings (referring to a column) or dictionaries with values as described below.

The ``<c0>``/``<ci0>`` and so forth are dictionaries containing the metadata
for each column, *including the index columns*. This has JSON form:

.. code-block:: text

   {'name': column_name,
    'field_name': parquet_column_name,
    'pandas_type': pandas_type,
    'numpy_type': numpy_type,
    'metadata': metadata}

See below for the detailed specification for these.

Index metadata descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~

``RangeIndex`` can be stored as metadata only, not requiring serialization. The
descriptor format for these as is follows:

.. code-block:: python

   index = pd.RangeIndex(0, 10, 2)
   {
       "kind": "range",
       "name": index.name,
       "start": index.start,
       "stop": index.stop,
       "step": index.step,
   }

Other index types must be serialized as data columns along with the other
DataFrame columns. The metadata for these is a string indicating the name of
the field in the data columns, for example ``'__index_level_0__'``.

If an index has a non-None ``name`` attribute, and there is no other column
with a name matching that value, then the ``index.name`` value can be used as
the descriptor. Otherwise (for unnamed indexes and ones with names colliding
with other column names) a disambiguating name with pattern matching
``__index_level_\d+__`` should be used. In cases of named indexes as data
columns, ``name`` attribute is always stored in the column descriptors as
above.

Column metadata
~~~~~~~~~~~~~~~

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
    'pandas_version': '0.20.0',
    'creator': {
      'library': 'pyarrow',
      'version': '0.13.0'
    }}
