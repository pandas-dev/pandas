{{ header }}

.. _api.extensions:

==========
Extensions
==========
.. currentmodule:: pandas

These are primarily intended for library authors looking to extend pandas
objects.

.. autosummary::
   :toctree: api/

   api.extensions.register_extension_dtype
   api.extensions.register_dataframe_accessor
   api.extensions.register_series_accessor
   api.extensions.register_index_accessor
   api.extensions.ExtensionDtype

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   api.extensions.ExtensionArray
   api.extensions.NDArrayBackedExtensionArray
   arrays.NumpyExtensionArray

.. We need this autosummary so that methods and attributes are generated.
.. Separate block, since they aren't classes.

   .. autosummary::
      :toctree: api/

      api.extensions.ExtensionArray._accumulate
      api.extensions.ExtensionArray._concat_same_type
      api.extensions.ExtensionArray._explode
      api.extensions.ExtensionArray._formatter
      api.extensions.ExtensionArray._from_factorized
      api.extensions.ExtensionArray._from_sequence
      api.extensions.ExtensionArray._from_sequence_of_strings
      api.extensions.ExtensionArray._hash_pandas_object
      api.extensions.ExtensionArray._pad_or_backfill
      api.extensions.ExtensionArray._reduce
      api.extensions.ExtensionArray._values_for_argsort
      api.extensions.ExtensionArray._values_for_factorize
      api.extensions.ExtensionArray.argsort
      api.extensions.ExtensionArray.astype
      api.extensions.ExtensionArray.copy
      api.extensions.ExtensionArray.view
      api.extensions.ExtensionArray.dropna
      api.extensions.ExtensionArray.duplicated
      api.extensions.ExtensionArray.equals
      api.extensions.ExtensionArray.factorize
      api.extensions.ExtensionArray.fillna
      api.extensions.ExtensionArray.insert
      api.extensions.ExtensionArray.interpolate
      api.extensions.ExtensionArray.isin
      api.extensions.ExtensionArray.isna
      api.extensions.ExtensionArray.ravel
      api.extensions.ExtensionArray.repeat
      api.extensions.ExtensionArray.searchsorted
      api.extensions.ExtensionArray.shift
      api.extensions.ExtensionArray.take
      api.extensions.ExtensionArray.unique
      api.extensions.ExtensionArray.dtype
      api.extensions.ExtensionArray.nbytes
      api.extensions.ExtensionArray.ndim
      api.extensions.ExtensionArray.shape
      api.extensions.ExtensionArray.tolist
      api.extensions.NDArrayBackedExtensionArray.dtype
      api.extensions.NDArrayBackedExtensionArray.argmax
      api.extensions.NDArrayBackedExtensionArray.argmin
      api.extensions.NDArrayBackedExtensionArray.argsort
      api.extensions.NDArrayBackedExtensionArray.astype
      api.extensions.NDArrayBackedExtensionArray.dropna
      api.extensions.NDArrayBackedExtensionArray.equals
      api.extensions.NDArrayBackedExtensionArray.factorize
      api.extensions.NDArrayBackedExtensionArray.fillna
      api.extensions.NDArrayBackedExtensionArray.insert
      api.extensions.NDArrayBackedExtensionArray.isin
      api.extensions.NDArrayBackedExtensionArray.isna
      api.extensions.NDArrayBackedExtensionArray.searchsorted
      api.extensions.NDArrayBackedExtensionArray.shift
      api.extensions.NDArrayBackedExtensionArray.take
      api.extensions.NDArrayBackedExtensionArray.to_numpy
      api.extensions.NDArrayBackedExtensionArray.tolist
      api.extensions.NDArrayBackedExtensionArray.unique
      api.extensions.NDArrayBackedExtensionArray.value_counts
      api.extensions.NDArrayBackedExtensionArray.view

Additionally, we have some utility methods for ensuring your object
behaves correctly.

.. autosummary::
  :toctree: api/

  api.indexers.check_array_indexer


The sentinel ``pandas.api.extensions.no_default`` is used as the default
value in some methods. Use an ``is`` comparison to check if the user
provides a non-default value.
