.. _extending:

{{ header }}

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
   class GeoAccessor:
       def __init__(self, pandas_obj):
           self._validate(pandas_obj)
           self._obj = pandas_obj

       @staticmethod
       def _validate(obj):
           # verify there is a column latitude and a column longitude
           if 'latitude' not in obj.columns or 'longitude' not in obj.columns:
               raise AttributeError("Must have 'latitude' and 'longitude'.")

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

We highly recommend validating the data in your accessor's `__init__`.
In our ``GeoAccessor``, we validate that the data contains the expected columns,
raising an ``AttributeError`` when the validation fails.
For a ``Series`` accessor, you should validate the ``dtype`` if the accessor
applies only to certain dtypes.


.. _extending.extension-types:

Extension Types
---------------

.. versionadded:: 0.23.0

.. warning::

   The :class:`pandas.api.extensions.ExtensionDtype` and :class:`pandas.api.extensions.ExtensionArray` APIs are new and
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

:class:`~pandas.api.extensions.ExtensionDtype`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A :class:`pandas.api.extensions.ExtensionDtype` is similar to a ``numpy.dtype`` object. It describes the
data type. Implementors are responsible for a few unique items like the name.

One particularly important item is the ``type`` property. This should be the
class that is the scalar type for your data. For example, if you were writing an
extension array for IP Address data, this might be ``ipaddress.IPv4Address``.

See the `extension dtype source`_ for interface definition.

.. versionadded:: 0.24.0

:class:`pandas.api.extension.ExtensionDtype` can be registered to pandas to allow creation via a string dtype name.
This allows one to instantiate ``Series`` and ``.astype()`` with a registered string name, for
example ``'category'`` is a registered string accessor for the ``CategoricalDtype``.

See the `extension dtype dtypes`_ for more on how to register dtypes.

:class:`~pandas.api.extensions.ExtensionArray`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. _extending.extension.operator:

:class:`~pandas.api.extensions.ExtensionArray` Operator Support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 0.24.0

By default, there are no operators defined for the class :class:`~pandas.api.extensions.ExtensionArray`.
There are two approaches for providing operator support for your ExtensionArray:

1. Define each of the operators on your ``ExtensionArray`` subclass.
2. Use an operator implementation from pandas that depends on operators that are already defined
   on the underlying elements (scalars) of the ExtensionArray.

.. note::

   Regardless of the approach, you may want to set ``__array_priority__``
   if you want your implementation to be called when involved in binary operations
   with NumPy arrays.

For the first approach, you define selected operators, e.g., ``__add__``, ``__le__``, etc. that
you want your ``ExtensionArray`` subclass to support.

The second approach assumes that the underlying elements (i.e., scalar type) of the ``ExtensionArray``
have the individual operators already defined.  In other words, if your ``ExtensionArray``
named ``MyExtensionArray`` is implemented so that each element is an instance
of the class ``MyExtensionElement``, then if the operators are defined
for ``MyExtensionElement``, the second approach will automatically
define the operators for ``MyExtensionArray``.

A mixin class, :class:`~pandas.api.extensions.ExtensionScalarOpsMixin` supports this second
approach.  If developing an ``ExtensionArray`` subclass, for example ``MyExtensionArray``,
can simply include ``ExtensionScalarOpsMixin`` as a parent class of ``MyExtensionArray``,
and then call the methods :meth:`~MyExtensionArray._add_arithmetic_ops` and/or
:meth:`~MyExtensionArray._add_comparison_ops` to hook the operators into
your ``MyExtensionArray`` class, as follows:

.. code-block:: python

    from pandas.api.extensions import ExtensionArray, ExtensionScalarOpsMixin

    class MyExtensionArray(ExtensionArray, ExtensionScalarOpsMixin):
        pass


    MyExtensionArray._add_arithmetic_ops()
    MyExtensionArray._add_comparison_ops()


.. note::

   Since ``pandas`` automatically calls the underlying operator on each
   element one-by-one, this might not be as performant as implementing your own
   version of the associated operators directly on the ``ExtensionArray``.

For arithmetic operations, this implementation will try to reconstruct a new
``ExtensionArray`` with the result of the element-wise operation. Whether
or not that succeeds depends on whether the operation returns a result
that's valid for the ``ExtensionArray``. If an ``ExtensionArray`` cannot
be reconstructed, an ndarray containing the scalars returned instead.

For ease of implementation and consistency with operations between pandas
and NumPy ndarrays, we recommend *not* handling Series and Indexes in your binary ops.
Instead, you should detect these cases and return ``NotImplemented``.
When pandas encounters an operation like ``op(Series, ExtensionArray)``, pandas
will

1. unbox the array from the ``Series`` (``Series.array``)
2. call ``result = op(values, ExtensionArray)``
3. re-box the result in a ``Series``

:class:`~pandas.api.extensions.ExtensionArray` Series Operations Support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 0.25.0

In addition to operators like `__mul__` and `__add__`, the pandas Series
namespace provides a long list of useful operations such as :meth:`Series.round`,
:meth:`Series.sum`, :meth:`Series.abs`, etc'. Some of these are handled by
pandas own algorithm implementations (via a dispatch function), while others
simply call an equivalent numpy function with data from the underlying array.
In order to support this operations in a new ExtensionArray, you must provide
an implementation for them.

.. note::

    There is a third category of operation which live on the `pandas`
    namespace, for example `:meth:pd.concat`. There is an equivalent
    numpy function `:meth:np.concatenate`, but this is not used. n
    general, these function should just work with your EA, you do not
    need to impelment more than the general EA interface.

As of 0.25.0, the list of series operations which pandas' provides its own
implementations are: :meth:`Series.any`, :meth:`Series.all`,
:meth:`Series.min`, :meth:`Series.max`, :meth:`Series.sum`,
:meth:`Series.mean`, :meth:`Series.median`, :meth:`Series.prod`
(and its alias :meth:`Series.product`), :meth:`Series.std`,
:meth:`Series.var`, :meth:`Series.sem`,
:meth:`Series.kurt`, and :meth:`Series.skew`.

In order to implement any of this functions, your ExtensionArray must include
an Implementation of :meth:`ExtensionArray._reduce`. Once you provide an
implementation of :meth:`ExtensionArray._reduce` which handles a particular
method, calling the method on the Series will invoke the implementation
on your ExtensionArray. All these methods are reduction functions, and
so are expected to return a scalr value of some type. However it is perfectly
acceptable to return some instance of an :class:`pandas.api.extensions.ExtensionArray`.

Series operations which are not handled by :meth:`ExtensionArray._reduce`,
such as :meth:`Series.round`, will generally invoke an equivalent numpy
function with your extension array as the argument. Pandas only guarantees
that your array will be passed to a numpy function, it does not dictate
how your ExtensionArray should interact with numpy's dispatch logic
in order to achieve its goal, since there are several alternative ways
of achieving similar results.

However, the details of numpy's dispatch logic are not entirely simple, and
there nuances which you should be aware of. For that reason, and in order to
make it easier to create new pandas extensions, we will now cover some of
possible approaches of dealing with numpy.

The first alternative, and the simplest, is to simply provide an `__array__`
method for your ExtensionArray. This is a standard numpy function documented
here (TBD), which must return a numpy array equivalent of your ExtensionArray.
This will usually be an array whose dtype is `object` and whose values are
instances of some class which your ExtensionArray wraps into an array. For
example, the pandas tests include an ExtensionArray example called
`DecimalArray`, if it used this method, its `__array__` method would return an
ndarray of `decimal.Decimal` objects.

Implementing `__array__`  is easy, but it usually isn't satisfactory because
it means most Series operations will return a Series of object dtype, instead
of maintaining your ExtensionArray's dtype.

The second approach is more involved, but it does a proper job of maintaining
the ExtensionArray's dtype through operations. It requires a detailed
understanding of how numpy functions operate on non ndarray objects.

Just as pandas handles some operation via :meth:`ExtensionArray._reduce`
and others by delegating to numpy, numpy makes a distinction between
between two types of opersions: ufuncs (such as `np.floor`, `np.ceil`,
and `np.abs`), and non-ufuncs (for example `np.round`, and `np.repeat`).

We will deal with ufuncs first. You can find a list of numpy's ufuncs here
(TBD). In order to support numpy ufuncs, a convenient approach is to implement
numpy's `__array_ufunc__` interface, specified in
[NEP13](https://www.numpy.org/neps/nep-0013-ufunc-overrides.html). In brief,
if your ExtensionArray implements a compliant `__array_ufunc__` interface,
when a numpy ufunc such as `np.floor` is invoked on your array, its
implementation of `__array_ufunc__`  will bec called first and given the
opportunity to compute the function. The return value needn't be a numpy
ndarray (though it can be). In general, you want the return value to be an
instance of your ExtensionArray. In some cases, your implementation can
calculate the result itself (see for example TBD), or, if your ExtensionArray
already has a numeric ndarray backing it, your implementation will itself
invoke the numpy ufunc itself on it (see for example TBD). In either case,
after computing the values, you will usually want to wrap the result as a new
ExtensionArray instance and return it to the caller. Pandas will automatically
use that Array as the backing ExtensionArray for a new Series object.

.. note::
    Before [NEP13](https://www.numpy.org/neps/nep-0013-ufunc-overrides.html),
    numpy already provides a way of wrapping ufunc functions via `__array_prepare__`
    and `__array_wrap__`, as documented in the Numpy docuemntation section
    ["Subclassing Numpy"](http://docs.python.org/doc/numpy/user/basics.subclassing.html).
    However, NEP13 seems to have largely superceded that mechanism.


With ufuncs out of the way, we turn to the remaining numpy operations, such
as `np.round`. The simplest way to support these operations is to simply
implement a compatible method on your ExtensionArray. For example, if your
ExtensionArray has a compatible `round` method on your ExtensionArray,
When python involes `ser.round()`, :meth:``Series.round` will invoke
`np.round(self.array)`, which will pass your ExtensionArray to the `np.round`
method. Numpy will detect that your EA implements a compatible `round`
and will invoke it to perform the operation. As in the ufunc case,
your implemntation will generally perform the calculaion itself,
or call numpy on its own acking numeric array, and in either case
will wrap the result as a new instance of ExtensionArray and return that
as a result. It is usually possible to write generic code to handle
most ufuncs without having to provide a special case for each. For an example, see TBD.

.. important::

  When providing implementations of numpy functions such as `np.round`,
  It essential that function signature is compatible with the numpy original.
  Otherwise,, numpy will ignore it.

  For example, the signature for `np.round` is `np.round(a, decimals=0, out=None)`.
  if you implement a round function which omits the `out` keyword:

.. code-block:: python

    def round(self, decimals=0):
        pass


  numpy will ignore it. The following will work however:

.. code-block:: python

    def round(self, decimals=0, **kwds):
        pass

An alternative approach to implementing individual functions, is to override
`__getattr__` in your ExtensionArray, and to intercept requests for method
names which you wish to support (such as `round`). For most functions,
you can return a dynamiclly generated function, which simply calls
the numpy function on your existing backing numeric array, wraps
the result in your ExtensionArray, and returns it. This approach can
reduce boilerplate significantly, but you do have to maintain a whitelist,
and may require more than one case, based on signature.


A third possible approach, is to use the `__array_function__`
mechanism introduced by [NEP18](https://www.numpy.org/neps/nep-0018-array-function-protocol.html).
This is an opt-in mechanism in numpy 1.16 (by setting an environment variable), and
is enabled by default starting with numpy 1.17. As of 1.17 it is still considered
experimental, and its design is actively being revised. We will not discuss it further
here, but it is certainly possible to make use of it to achieve the same goal. Your
mileage may vary.

.. important::
    Implementing `__array_function__` is not a substitute for implementing `__array_ufunc__`.
    The `__array_function__` mechanism complements (and to a degree copies) the`__array_ufunc__`
    mechanism, by providing the same flexibility for non-ufuncs.

.. important::
    `__array_function__` is an "all-in" solution. That means that if you cannot mix it with
    explicit implementations for some methods and using `__array_function__` for some.
    If you both `__array_function__` and also provide an implementation of `round`, numpy
    will invoke `__array_function__` for all the operations in the specification, **including**
     `round`.

With this overview in hand, you hopefully have the necessary information in order
to develop rich, full-featured ExtensionArrays that seamlessly plug in to pandas.

.. important::
    You are not required to provide implementations for the full complemnt of Series
    operations in your ExtensionArray. In fact, some of them may not even make sense
    within your context. You amay also choose to ass implementations incrementally,
    as the need arised.

    TBD: should we have a standard way of signalling not supported instead of a
    random AttributeError exception being thrown.

.. important::

The above description currently leads the state of the code considerably. Many Series
methods need to be updated to conform to this model of EA support. If you find a
bug, or something else which does not behave as described, please report it to
the pandas team by opening an issue.


Formatting Extension Arrays
^^^^^^^^^^^^^^^^^^^^^^^^

TBD

.. _extending.extension.testing:

Testing Extension Arrays
^^^^^^^^^^^^^^^^^^^^^^^^

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

.. _extension dtype dtypes: https://github.com/pandas-dev/pandas/blob/master/pandas/core/dtypes/dtypes.py
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

* ``_constructor``: Used when a manipulation result has the same dimensions as the original.
* ``_constructor_sliced``: Used when a manipulation result has one lower dimension(s) as the original, such as ``DataFrame`` single columns slicing.
* ``_constructor_expanddim``: Used when a manipulation result has one higher dimension as the original, such as ``Series.to_frame()``.

Following table shows how ``pandas`` data structures define constructor properties by default.

===========================  ======================= =============
Property Attributes          ``Series``              ``DataFrame``
===========================  ======================= =============
``_constructor``             ``Series``              ``DataFrame``
``_constructor_sliced``      ``NotImplementedError`` ``Series``
``_constructor_expanddim``   ``DataFrame``           ``NotImplementedError``
===========================  ======================= =============

Below example shows how to define ``SubclassedSeries`` and ``SubclassedDataFrame`` overriding constructor properties.

.. code-block:: python

   class SubclassedSeries(pd.Series):

       @property
       def _constructor(self):
           return SubclassedSeries

       @property
       def _constructor_expanddim(self):
           return SubclassedDataFrame


   class SubclassedDataFrame(pd.DataFrame):

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

   >>> df = SubclassedDataFrame({'A': [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]})
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

   class SubclassedDataFrame2(pd.DataFrame):

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
