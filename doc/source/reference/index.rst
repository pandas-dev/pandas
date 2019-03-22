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

.. If you update this toctree, also update the manual toctree in the
   main index.rst.template

.. toctree::
   :maxdepth: 2

   io
   general_functions
   series
   frame
   arrays
   panel
   indexing
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

..
    .. toctree::

        api/pandas.DataFrame.blocks
        api/pandas.DataFrame.as_matrix
        api/pandas.DataFrame.ix
        api/pandas.Index.asi8
        api/pandas.Index.data
        api/pandas.Index.flags
        api/pandas.Index.holds_integer
        api/pandas.Index.is_type_compatible
        api/pandas.Index.nlevels
        api/pandas.Index.sort
        api/pandas.Series.asobject
        api/pandas.Series.blocks
        api/pandas.Series.from_array
        api/pandas.Series.ix
        api/pandas.Series.imag
        api/pandas.Series.real


.. Can't convince sphinx to generate toctree for this class attribute.
.. So we do it manually to avoid a warning

..
    .. toctree::

        api/pandas.api.extensions.ExtensionDtype.na_value
