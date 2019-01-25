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
   api.extensions.ExtensionArray
   arrays.PandasArray
