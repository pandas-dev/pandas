{{ header }}

.. _api.general_utility_functions:

=========================
General utility functions
=========================
.. currentmodule:: pandas

Working with options
--------------------
.. autosummary::
   :toctree: api/

   describe_option
   reset_option
   get_option
   set_option
   option_context

.. _api.general.testing:

Testing functions
-----------------
.. autosummary::
   :toctree: api/

   testing.assert_frame_equal
   testing.assert_series_equal
   testing.assert_index_equal
   testing.assert_equal
   testing.assert_almost_equal
   testing.assert_categorical_equal
   testing.assert_datetime_array_equal
   testing.assert_extension_array_equal
   testing.assert_interval_array_equal
   testing.assert_numpy_array_equal
   testing.assert_period_array_equal
   testing.assert_sp_array_equal
   testing.assert_timedelta_array_equal

Exceptions and warnings
-----------------------
.. autosummary::
   :toctree: api/

   errors.DtypeWarning
   errors.EmptyDataError
   errors.OutOfBoundsDatetime
   errors.ParserError
   errors.ParserWarning
   errors.PerformanceWarning
   errors.UnsortedIndexError
   errors.UnsupportedFunctionCall

Data types related functionality
--------------------------------
.. autosummary::
   :toctree: api/

   api.types.union_categoricals
   api.types.infer_dtype
   api.types.pandas_dtype

Dtype introspection
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

    api.types.is_bool_dtype
    api.types.is_categorical_dtype
    api.types.is_complex_dtype
    api.types.is_datetime64_any_dtype
    api.types.is_datetime64_dtype
    api.types.is_datetime64_ns_dtype
    api.types.is_datetime64tz_dtype
    api.types.is_extension_type
    api.types.is_extension_array_dtype
    api.types.is_float_dtype
    api.types.is_int64_dtype
    api.types.is_integer_dtype
    api.types.is_interval_dtype
    api.types.is_numeric_dtype
    api.types.is_object_dtype
    api.types.is_period_dtype
    api.types.is_signed_integer_dtype
    api.types.is_string_dtype
    api.types.is_timedelta64_dtype
    api.types.is_timedelta64_ns_dtype
    api.types.is_unsigned_integer_dtype
    api.types.is_sparse

Iterable introspection
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

    api.types.is_dict_like
    api.types.is_file_like
    api.types.is_list_like
    api.types.is_named_tuple
    api.types.is_iterator

Scalar introspection
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

    api.types.is_bool
    api.types.is_categorical
    api.types.is_complex
    api.types.is_float
    api.types.is_hashable
    api.types.is_integer
    api.types.is_interval
    api.types.is_number
    api.types.is_re
    api.types.is_re_compilable
    api.types.is_scalar
