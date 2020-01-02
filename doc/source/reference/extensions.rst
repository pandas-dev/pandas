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
   arrays.PandasArray

.. We need this autosummary so that methods and attributes are generated.
.. Separate block, since they aren't classes.

   .. autosummary::
      :toctree: api/

      api.extensions.ExtensionArray._concat_same_type
      api.extensions.ExtensionArray._formatter
      api.extensions.ExtensionArray._from_factorized
      api.extensions.ExtensionArray._from_sequence
      api.extensions.ExtensionArray._from_sequence_of_strings
      api.extensions.ExtensionArray._ndarray_values
      api.extensions.ExtensionArray._reduce
      api.extensions.ExtensionArray._values_for_argsort
      api.extensions.ExtensionArray._values_for_factorize
      api.extensions.ExtensionArray.argsort
      api.extensions.ExtensionArray.astype
      api.extensions.ExtensionArray.copy
      api.extensions.ExtensionArray.view
      api.extensions.ExtensionArray.dropna
      api.extensions.ExtensionArray.factorize
      api.extensions.ExtensionArray.fillna
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

Additionally, we have some utility methods for ensuring your object
behaves correctly.

.. autosummary::
  :toctree: api/

  api.indexers.check_bool_array_indexer
