{{ header }}

.. _api:

=============
API Reference
=============

This page gives an overview of all public pandas objects, functions and
methods. All classes and functions exposed in ``pandas.*`` namespace are public.

Some subpackages are public which include ``pandas.errors``,
``pandas.plotting``, and ``pandas.testing``. Public functions in
``pandas.io`` and ``pandas.tseries`` submodules are mentioned in
the documentation. ``pandas.api.types`` subpackage holds some
public functions related to data types in pandas.

.. warning::

    The ``pandas.core``, ``pandas.compat``, and ``pandas.util`` top-level modules are PRIVATE. Stable functionality in such modules is not guaranteed.

.. toctree::
   :maxdepth: 2

   io
   general_functions
   series
   frame
   panel
   indexing
   scalars
   offset_frequency
   window
   groupby
   resampling
   style
   plotting
   general_utility_functions
   extensions

.. This is to prevent warnings in the doc build. We don't want to encourage
.. these methods.

.. toctree::
   :hidden:

   generated/pandas.DataFrame.blocks
   generated/pandas.DataFrame.as_matrix
   generated/pandas.DataFrame.ix
   generated/pandas.Index.asi8
   generated/pandas.Index.data
   generated/pandas.Index.flags
   generated/pandas.Index.holds_integer
   generated/pandas.Index.is_type_compatible
   generated/pandas.Index.nlevels
   generated/pandas.Index.sort
   generated/pandas.Panel.agg
   generated/pandas.Panel.aggregate
   generated/pandas.Panel.blocks
   generated/pandas.Panel.empty
   generated/pandas.Panel.is_copy
   generated/pandas.Panel.items
   generated/pandas.Panel.ix
   generated/pandas.Panel.major_axis
   generated/pandas.Panel.minor_axis
   generated/pandas.Series.asobject
   generated/pandas.Series.blocks
   generated/pandas.Series.from_array
   generated/pandas.Series.ix
   generated/pandas.Series.imag
   generated/pandas.Series.real


.. Can't convince sphinx to generate toctree for this class attribute.
.. So we do it manually to avoid a warning

.. toctree::
   :hidden:

   generated/pandas.api.extensions.ExtensionDtype.na_value
