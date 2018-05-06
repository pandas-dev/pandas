.. _extending:

****************
Extending Pandas
****************

While pandas provides a rich set of methods, containers, and data types, your
needs may not be fully satisfied. Pandas offers a few options for extending
pandas.

.. _extending.register-accessors:

Registering Custom Accessors
----------------------------

Libraries can use the decorators
:func:`pandas.api.extensions.register_dataframe_accessor`,
:func:`pandas.api.extensions.register_series_accessor`, and
:func:`pandas.api.extensions.register_index_accessor`, to add additional
"namespaces" to pandas objects. All of these follow a similar convention: you
decorate a class, providing the name of attribute to add. The class's
``__init__`` method gets the object being decorated. For example:

.. code-block:: python

   @pd.api.extensions.register_dataframe_accessor("geo")
   class GeoAccessor(object):
       def __init__(self, pandas_obj):
           self._obj = pandas_obj

       @property
       def center(self):
           # return the geographic center point of this DataFrame
           lat = self._obj.latitude
           lon = self._obj.longitude
           return (float(lon.mean()), float(lat.mean()))

       def plot(self):
           # plot this array's data on a map, e.g., using Cartopy
           pass

Now users can access your methods using the ``geo`` namespace:

      >>> ds = pd.DataFrame({'longitude': np.linspace(0, 10),
      ...                    'latitude': np.linspace(0, 20)})
      >>> ds.geo.center
      (5.0, 10.0)
      >>> ds.geo.plot()
      # plots data on a map

This can be a convenient way to extend pandas objects without subclassing them.
If you write a custom accessor, make a pull request adding it to our
:ref:`ecosystem` page.

.. _extending.extension-types:

Extension Types
---------------

.. versionadded:: 0.23.0

.. warning::

   The :class:`pandas.api.extension.ExtensionDtype` and :class:`pandas.api.extension.ExtensionArray` APIs are new and
   experimental. They may change between versions without warning.

Pandas defines an interface for implementing data types and arrays that *extend*
NumPy's type system. Pandas itself uses the extension system for some types
that aren't built into NumPy (categorical, period, interval, datetime with
timezone).

Libraries can define a custom array and data type. When pandas encounters these
objects, they will be handled properly (i.e. not converted to an ndarray of
objects). Many methods like :func:`pandas.isna` will dispatch to the extension
type's implementation.

If you're building a library that implements the interface, please publicize it
on :ref:`ecosystem.extensions`.

The interface consists of two classes.

:class:`~pandas.api.extension.ExtensionDtype`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A :class:`pandas.api.extension.ExtensionDtype` is similar to a ``numpy.dtype`` object. It describes the
data type. Implementors are responsible for a few unique items like the name.

One particularly important item is the ``type`` property. This should be the
class that is the scalar type for your data. For example, if you were writing an
extension array for IP Address data, this might be ``ipaddress.IPv4Address``.

See the `extension dtype source`_ for interface definition.

:class:`~pandas.api.extension.ExtensionArray`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class provides all the array-like functionality. ExtensionArrays are
limited to 1 dimension. An ExtensionArray is linked to an ExtensionDtype via the
``dtype`` attribute.

Pandas makes no restrictions on how an extension array is created via its
``__new__`` or ``__init__``, and puts no restrictions on how you store your
data. We do require that your array be convertible to a NumPy array, even if
this is relatively expensive (as it is for ``Categorical``).

They may be backed by none, one, or many NumPy arrays. For example,
``pandas.Categorical`` is an extension array backed by two arrays,
one for codes and one for categories. An array of IPv6 addresses may
be backed by a NumPy structured array with two fields, one for the
lower 64 bits and one for the upper 64 bits. Or they may be backed
by some other storage type, like Python lists.

See the `extension array source`_ for the interface definition. The docstrings
and comments contain guidance for properly implementing the interface.

We provide a test suite for ensuring that your extension arrays satisfy the expected
behavior. To use the test suite, you must provide several pytest fixtures and inherit
from the base test class. The required fixtures are found in
https://github.com/pandas-dev/pandas/blob/master/pandas/tests/extension/conftest.py.

To use a test, subclass it:

.. code-block:: python

   from pandas.tests.extension import base

   class TestConstructors(base.BaseConstructorsTests):
       pass


See https://github.com/pandas-dev/pandas/blob/master/pandas/tests/extension/base/__init__.py
for a list of all the tests available.

.. _extension dtype source: https://github.com/pandas-dev/pandas/blob/master/pandas/core/dtypes/base.py
.. _extension array source: https://github.com/pandas-dev/pandas/blob/master/pandas/core/arrays/base.py

.. _extending.subclassing-pandas:

Subclassing pandas Data Structures
----------------------------------

.. warning:: There are some easier alternatives before considering subclassing ``pandas`` data structures.

  1. Extensible method chains with :ref:`pipe <basics.pipe>`

  2. Use *composition*. See `here <http://en.wikipedia.org/wiki/Composition_over_inheritance>`_.

  3. Extending by :ref:`registering an accessor <extending.register-accessors>`

  4. Extending by :ref:`extension type <extending.extension-types>`

This section describes how to subclass ``pandas`` data structures to meet more specific needs. There are two points that need attention:

1. Override constructor properties.
2. Define original properties

.. note::

   You can find a nice example in `geopandas <https://github.com/geopandas/geopandas>`_ project.

Override Constructor Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each data structure has several *constructor properties* for returning a new
data structure as the result of an operation. By overriding these properties,
you can retain subclasses through ``pandas`` data manipulations.

There are 3 constructor properties to be defined:

- ``_constructor``: Used when a manipulation result has the same dimesions as the original.
- ``_constructor_sliced``: Used when a manipulation result has one lower dimension(s) as the original, such as ``DataFrame`` single columns slicing.
- ``_constructor_expanddim``: Used when a manipulation result has one higher dimension as the original, such as ``Series.to_frame()`` and ``DataFrame.to_panel()``.

Following table shows how ``pandas`` data structures define constructor properties by default.

===========================  ======================= =============
Property Attributes          ``Series``              ``DataFrame``      
===========================  ======================= =============
``_constructor``             ``Series``              ``DataFrame``      
``_constructor_sliced``      ``NotImplementedError`` ``Series``         
``_constructor_expanddim``   ``DataFrame``           ``Panel``          
===========================  ======================= =============

Below example shows how to define ``SubclassedSeries`` and ``SubclassedDataFrame`` overriding constructor properties.

.. code-block:: python

   class SubclassedSeries(Series):

       @property
       def _constructor(self):
           return SubclassedSeries

       @property
       def _constructor_expanddim(self):
           return SubclassedDataFrame

   class SubclassedDataFrame(DataFrame):

       @property
       def _constructor(self):
           return SubclassedDataFrame

       @property
       def _constructor_sliced(self):
           return SubclassedSeries

.. code-block:: python

   >>> s = SubclassedSeries([1, 2, 3])
   >>> type(s)
   <class '__main__.SubclassedSeries'>

   >>> to_framed = s.to_frame()
   >>> type(to_framed)
   <class '__main__.SubclassedDataFrame'>

   >>> df = SubclassedDataFrame({'A', [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]})
   >>> df
      A  B  C
   0  1  4  7
   1  2  5  8
   2  3  6  9

   >>> type(df)
   <class '__main__.SubclassedDataFrame'>

   >>> sliced1 = df[['A', 'B']]
   >>> sliced1
      A  B
   0  1  4
   1  2  5
   2  3  6
   >>> type(sliced1)
   <class '__main__.SubclassedDataFrame'>

   >>> sliced2 = df['A']
   >>> sliced2
   0    1
   1    2
   2    3
   Name: A, dtype: int64
   >>> type(sliced2)
   <class '__main__.SubclassedSeries'>

Define Original Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

To let original data structures have additional properties, you should let ``pandas`` know what properties are added. ``pandas`` maps unknown properties to data names overriding ``__getattribute__``. Defining original properties can be done in one of 2 ways:

1. Define ``_internal_names`` and ``_internal_names_set`` for temporary properties which WILL NOT be passed to manipulation results.
2. Define ``_metadata`` for normal properties which will be passed to manipulation results.

Below is an example to define two original properties, "internal_cache" as a temporary property and "added_property" as a normal property

.. code-block:: python

   class SubclassedDataFrame2(DataFrame):

       # temporary properties
       _internal_names = pd.DataFrame._internal_names + ['internal_cache']
       _internal_names_set = set(_internal_names)

       # normal properties
       _metadata = ['added_property']

       @property
       def _constructor(self):
           return SubclassedDataFrame2

.. code-block:: python

   >>> df = SubclassedDataFrame2({'A': [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]})
   >>> df
      A  B  C
   0  1  4  7
   1  2  5  8
   2  3  6  9

   >>> df.internal_cache = 'cached'
   >>> df.added_property = 'property'

   >>> df.internal_cache
   cached
   >>> df.added_property
   property

   # properties defined in _internal_names is reset after manipulation
   >>> df[['A', 'B']].internal_cache
   AttributeError: 'SubclassedDataFrame2' object has no attribute 'internal_cache'

   # properties defined in _metadata are retained
   >>> df[['A', 'B']].added_property
   property
