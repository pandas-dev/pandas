{{ header }}

.. _api:

=============
API reference
=============

This page gives an overview of all public pandas objects, functions and
methods. All classes and functions exposed in ``pandas.*`` namespace are public.

The following subpackages are public.

- ``pandas.errors``: Custom exception and warnings classes that are raised by pandas.
- ``pandas.plotting``: Plotting public API.
- ``pandas.testing``: Functions that are useful for writing tests involving pandas objects.
- ``pandas.api.extensions``: Functions and classes for extending pandas objects.
- ``pandas.api.indexers``: Functions and classes for rolling window indexers.
- ``pandas.api.interchange``: DataFrame interchange protocol.
- ``pandas.api.types``: Datatype classes and functions.
- ``pandas.api.typing``: Classes that may be necessary for type-hinting.
  These are classes that are encountered as intermediate results but should not be instantiated
  directly by users. These classes are not to be confused with classes from the
  `pandas-stubs <https://github.com/pandas-dev/pandas-stubs>`_ package
  which has classes in addition to those that occur in pandas for type-hinting.

In addition, public functions in ``pandas.io``, ``pandas.tseries``, ``pandas.util`` submodules
are explicitly mentioned in the documentation. Further APIs in these modules are not guaranteed
to be stable.


.. warning::

    The ``pandas.core``, ``pandas.compat`` top-level modules are PRIVATE. Stable functionality in such modules is not guaranteed.

.. If you update this toctree, also update the manual toctree in the
.. main index.rst.template

.. toctree::
   :maxdepth: 2

   io
   general_functions
   series
   frame
   arrays
   indexing
   offset_frequency
   window
   groupby
   resampling
   style
   plotting
   options
   extensions
   testing
   missing_value

.. This is to prevent warnings in the doc build. We don't want to encourage
.. these methods.

..
    .. toctree::

        api/pandas.Index.nlevels
        api/pandas.Index.sort


.. Can't convince sphinx to generate toctree for this class attribute.
.. So we do it manually to avoid a warning

..
    .. toctree::

        api/pandas.api.extensions.ExtensionDtype.na_value
