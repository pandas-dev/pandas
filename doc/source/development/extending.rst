.. _extending:

{{ header }}

****************
Extending pandas
****************

While pandas provides a rich set of methods, containers, and data types, your
needs may not be fully satisfied. pandas offers a few options for extending
pandas.

.. _extending.register-accessors:

Registering custom accessors
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
           if "latitude" not in obj.columns or "longitude" not in obj.columns:
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

      >>> ds = pd.DataFrame(
      ...     {"longitude": np.linspace(0, 10), "latitude": np.linspace(0, 20)}
      ... )
      >>> ds.geo.center
      (5.0, 10.0)
      >>> ds.geo.plot()
      # plots data on a map

This can be a convenient way to extend pandas objects without subclassing them.
If you write a custom accessor, make a pull request adding it to our
`ecosystem <https://pandas.pydata.org/community/ecosystem.html>`_ page.

We highly recommend validating the data in your accessor's ``__init__``.
In our ``GeoAccessor``, we validate that the data contains the expected columns,
raising an ``AttributeError`` when the validation fails.
For a ``Series`` accessor, you should validate the ``dtype`` if the accessor
applies only to certain dtypes.


.. _extending.extension-types:

Extension types
---------------

.. note::

   The :class:`pandas.api.extensions.ExtensionDtype` and :class:`pandas.api.extensions.ExtensionArray` APIs were
   experimental prior to pandas 1.5. Starting with version 1.5, future changes will follow
   the :ref:`pandas deprecation policy <policies.version>`.

pandas defines an interface for implementing data types and arrays that *extend*
NumPy's type system. pandas itself uses the extension system for some types
that aren't built into NumPy (categorical, period, interval, datetime with
timezone).

Libraries can define a custom array and data type. When pandas encounters these
objects, they will be handled properly (i.e. not converted to an ndarray of
objects). Many methods like :func:`pandas.isna` will dispatch to the extension
type's implementation.

If you're building a library that implements the interface, please publicize it
on `the ecosystem page <https://pandas.pydata.org/community/ecosystem.html>`_.

The interface consists of two classes.

:class:`~pandas.api.extensions.ExtensionDtype`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A :class:`pandas.api.extensions.ExtensionDtype` is similar to a ``numpy.dtype`` object. It describes the
data type. Implementers are responsible for a few unique items like the name.

One particularly important item is the ``type`` property. This should be the
class that is the scalar type for your data. For example, if you were writing an
extension array for IP Address data, this might be ``ipaddress.IPv4Address``.

See the `extension dtype source`_ for interface definition.

:class:`pandas.api.extensions.ExtensionDtype` can be registered to pandas to allow creation via a string dtype name.
This allows one to instantiate ``Series`` and ``.astype()`` with a registered string name, for
example ``'category'`` is a registered string accessor for the ``CategoricalDtype``.

See the `extension dtype dtypes`_ for more on how to register dtypes.

:class:`~pandas.api.extensions.ExtensionArray`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class provides all the array-like functionality. ExtensionArrays are
limited to 1 dimension. An ExtensionArray is linked to an ExtensionDtype via the
``dtype`` attribute.

pandas makes no restrictions on how an extension array is created via its
``__new__`` or ``__init__``, and puts no restrictions on how you store your
data. We do require that your array be convertible to a NumPy array, even if
this is relatively expensive (as it is for ``Categorical``).

They may be backed by none, one, or many NumPy arrays. For example,
:class:`pandas.Categorical` is an extension array backed by two arrays,
one for codes and one for categories. An array of IPv6 addresses may
be backed by a NumPy structured array with two fields, one for the
lower 64 bits and one for the upper 64 bits. Or they may be backed
by some other storage type, like Python lists.

See the `extension array source`_ for the interface definition. The docstrings
and comments contain guidance for properly implementing the interface.

:class:`~pandas.api.extensions.NDArrayBackedExtensionArray`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For ExtensionArrays backed by a single NumPy array, the
:class:`~pandas.api.extensions.NDArrayBackedExtensionArray` class can save you
some effort. It contains a private property ``_ndarray`` with the backing NumPy
array and implements the extension array interface.

Implement the following:

``_box_func``
  Convert from array values to the type you wish to expose to users.

``_internal_fill_value``
   Scalar used to denote ``NA`` value inside our ``self._ndarray``, e.g. ``-1``
   for ``Categorical``, ``iNaT`` for ``Period``.

``_validate_scalar``
  Convert from an object to a value which can be stored in the NumPy array.

``_validate_setitem_value``
  Convert a value or values for use in setting a value or values in the backing
  NumPy array.

.. code-block:: python

    class DateArray(NDArrayBackedExtensionArray):
        _internal_fill_value = numpy.datetime64("NaT")

        def __init__(self, values):
            backing_array_dtype = "<M8[ns]"
            super().__init__(values=values, dtype=backing_array_dtype)

        def _box_func(self, value):
            if pandas.isna(x):
                return pandas.NaT
            return x.astype("datetime64[us]").item().date()

        def _validate_scalar(self, scalar):
            if pandas.isna(scalar):
                return numpy.datetime64("NaT")
            elif isinstance(scalar, datetime.date):
                return pandas.Timestamp(
                    year=scalar.year, month=scalar.month, day=scalar.day
                ).to_datetime64()
            else:
                raise TypeError("Invalid value type", scalar)

        def _validate_setitem_value(self, value):
            if pandas.api.types.is_list_like(value):
                return np.array([self._validate_scalar(v) for v in value], dtype = self.dtype)
            return self._validate_scalar(value)


To support 2D arrays, use the ``_from_backing_data`` helper function when a
method is called on multi-dimensional data of the same dtype as ``_ndarray``.

.. code-block:: python

    class CustomArray(NDArrayBackedExtensionArray):

        ...

        def min(self, *, axis: Optional[int] = None, skipna: bool = True, **kwargs):
            pandas.compat.numpy.function.validate_minnumpy_validate_min((), kwargs)
            result = pandas.core.nanops.nanmin(
                values=self._ndarray, axis=axis, mask=self.isna(), skipna=skipna
            )
            if axis is None or self.ndim == 1:
                return self._box_func(result)
            return self._from_backing_data(result)


Subclass the tests in :mod:`pandas.tests.extension.base` in your test suite to
validate your implementation.

.. code-block:: python

    @pytest.fixture
    def data():
        return CustomArray(numpy.arange(-10, 10, 1)


    class TestCustomArray(base.ExtensionTests):
        pass


.. _extending.extension.operator:

:class:`~pandas.api.extensions.ExtensionArray` operator support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. _extending.extension.ufunc:

NumPy universal functions
^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`Series` implements ``__array_ufunc__``. As part of the implementation,
pandas unboxes the ``ExtensionArray`` from the :class:`Series`, applies the ufunc,
and re-boxes it if necessary.

If applicable, we highly recommend that you implement ``__array_ufunc__`` in your
extension array to avoid coercion to an ndarray. See
`the NumPy documentation <https://numpy.org/doc/stable/reference/generated/numpy.lib.mixins.NDArrayOperatorsMixin.html>`__
for an example.

As part of your implementation, we require that you defer to pandas when a pandas
container (:class:`Series`, :class:`DataFrame`, :class:`Index`) is detected in ``inputs``.
If any of those is present, you should return ``NotImplemented``. pandas will take care of
unboxing the array from the container and re-calling the ufunc with the unwrapped input.

.. _extending.extension.testing:

Testing extension arrays
^^^^^^^^^^^^^^^^^^^^^^^^

We provide a test suite for ensuring that your extension arrays satisfy the expected
behavior. To use the test suite, you must provide several pytest fixtures and inherit
from the base test class. The required fixtures are found in
https://github.com/pandas-dev/pandas/blob/main/pandas/tests/extension/conftest.py.

To use a test, subclass it:

.. code-block:: python

   from pandas.tests.extension import base


   class TestConstructors(base.BaseConstructorsTests):
       pass


See https://github.com/pandas-dev/pandas/blob/main/pandas/tests/extension/base/__init__.py
for a list of all the tests available.

.. _extending.extension.arrow:

Compatibility with Apache Arrow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An ``ExtensionArray`` can support conversion to / from ``pyarrow`` arrays
(and thus support for example serialization to the Parquet file format)
by implementing two methods: ``ExtensionArray.__arrow_array__`` and
``ExtensionDtype.__from_arrow__``.

The ``ExtensionArray.__arrow_array__`` ensures that ``pyarrow`` knowns how
to convert the specific extension array into a ``pyarrow.Array`` (also when
included as a column in a pandas DataFrame):

.. code-block:: python

    class MyExtensionArray(ExtensionArray):
        ...

        def __arrow_array__(self, type=None):
            # convert the underlying array values to a pyarrow Array
            import pyarrow

            return pyarrow.array(..., type=type)

The ``ExtensionDtype.__from_arrow__`` method then controls the conversion
back from pyarrow to a pandas ExtensionArray. This method receives a pyarrow
``Array`` or ``ChunkedArray`` as only argument and is expected to return the
appropriate pandas ``ExtensionArray`` for this dtype and the passed values:

.. code-block:: none

    class ExtensionDtype:
        ...

        def __from_arrow__(self, array: pyarrow.Array/ChunkedArray) -> ExtensionArray:
            ...

See more in the `Arrow documentation <https://arrow.apache.org/docs/python/extending_types.html>`__.

Those methods have been implemented for the nullable integer and string extension
dtypes included in pandas, and ensure roundtrip to pyarrow and the Parquet file format.

.. _extension dtype dtypes: https://github.com/pandas-dev/pandas/blob/main/pandas/core/dtypes/dtypes.py
.. _extension dtype source: https://github.com/pandas-dev/pandas/blob/main/pandas/core/dtypes/base.py
.. _extension array source: https://github.com/pandas-dev/pandas/blob/main/pandas/core/arrays/base.py

.. _extending.subclassing-pandas:

Subclassing pandas data structures
----------------------------------

.. warning:: There are some easier alternatives before considering subclassing ``pandas`` data structures.

  1. Extensible method chains with :ref:`pipe <basics.pipe>`

  2. Use *composition*. See `here <https://en.wikipedia.org/wiki/Composition_over_inheritance>`_.

  3. Extending by :ref:`registering an accessor <extending.register-accessors>`

  4. Extending by :ref:`extension type <extending.extension-types>`

This section describes how to subclass ``pandas`` data structures to meet more specific needs. There are two points that need attention:

1. Override constructor properties.
2. Define original properties

.. note::

   You can find a nice example in `geopandas <https://github.com/geopandas/geopandas>`_ project.

Override constructor properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each data structure has several *constructor properties* for returning a new
data structure as the result of an operation. By overriding these properties,
you can retain subclasses through ``pandas`` data manipulations.

There are 3 possible constructor properties to be defined on a subclass:

* ``DataFrame/Series._constructor``: Used when a manipulation result has the same dimension as the original.
* ``DataFrame._constructor_sliced``: Used when a ``DataFrame`` (sub-)class manipulation result should be a ``Series`` (sub-)class.
* ``Series._constructor_expanddim``: Used when a ``Series`` (sub-)class manipulation result should be a ``DataFrame`` (sub-)class, e.g. ``Series.to_frame()``.

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

   >>> df = SubclassedDataFrame({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})
   >>> df
      A  B  C
   0  1  4  7
   1  2  5  8
   2  3  6  9

   >>> type(df)
   <class '__main__.SubclassedDataFrame'>

   >>> sliced1 = df[["A", "B"]]
   >>> sliced1
      A  B
   0  1  4
   1  2  5
   2  3  6

   >>> type(sliced1)
   <class '__main__.SubclassedDataFrame'>

   >>> sliced2 = df["A"]
   >>> sliced2
   0    1
   1    2
   2    3
   Name: A, dtype: int64

   >>> type(sliced2)
   <class '__main__.SubclassedSeries'>

Define original properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

To let original data structures have additional properties, you should let ``pandas`` know what properties are added. ``pandas`` maps unknown properties to data names overriding ``__getattribute__``. Defining original properties can be done in one of 2 ways:

1. Define ``_internal_names`` and ``_internal_names_set`` for temporary properties which WILL NOT be passed to manipulation results.
2. Define ``_metadata`` for normal properties which will be passed to manipulation results.

Below is an example to define two original properties, "internal_cache" as a temporary property and "added_property" as a normal property

.. code-block:: python

   class SubclassedDataFrame2(pd.DataFrame):

       # temporary properties
       _internal_names = pd.DataFrame._internal_names + ["internal_cache"]
       _internal_names_set = set(_internal_names)

       # normal properties
       _metadata = ["added_property"]

       @property
       def _constructor(self):
           return SubclassedDataFrame2

.. code-block:: python

   >>> df = SubclassedDataFrame2({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})
   >>> df
      A  B  C
   0  1  4  7
   1  2  5  8
   2  3  6  9

   >>> df.internal_cache = "cached"
   >>> df.added_property = "property"

   >>> df.internal_cache
   cached
   >>> df.added_property
   property

   # properties defined in _internal_names is reset after manipulation
   >>> df[["A", "B"]].internal_cache
   AttributeError: 'SubclassedDataFrame2' object has no attribute 'internal_cache'

   # properties defined in _metadata are retained
   >>> df[["A", "B"]].added_property
   property

.. _extending.plotting-backends:

Plotting backends
-----------------

pandas can be extended with third-party plotting backends. The
main idea is letting users select a plotting backend different than the provided
one based on Matplotlib. For example:

.. code-block:: python

    >>> pd.set_option("plotting.backend", "backend.module")
    >>> pd.Series([1, 2, 3]).plot()

This would be more or less equivalent to:

.. code-block:: python

    >>> import backend.module
    >>> backend.module.plot(pd.Series([1, 2, 3]))

The backend module can then use other visualization tools (Bokeh, Altair,...)
to generate the plots.

Libraries implementing the plotting backend should use `entry points <https://setuptools.pypa.io/en/latest/userguide/entry_point.html>`__
to make their backend discoverable to pandas. The key is ``"pandas_plotting_backends"``. For example, pandas
registers the default "matplotlib" backend as follows.

.. code-block:: python

   # in setup.py
   setup(  # noqa: F821
       ...,
       entry_points={
           "pandas_plotting_backends": [
               "matplotlib = pandas:plotting._matplotlib",
           ],
       },
   )


More information on how to implement a third-party plotting backend can be found at
https://github.com/pandas-dev/pandas/blob/main/pandas/plotting/__init__.py#L1.

.. _extending.pandas_priority:

Arithmetic with 3rd party types
-------------------------------

In order to control how arithmetic works between a custom type and a pandas type,
implement ``__pandas_priority__``.  Similar to numpy's ``__array_priority__``
semantics, arithmetic methods on :class:`DataFrame`, :class:`Series`, and :class:`Index`
objects will delegate to ``other``, if it has an attribute ``__pandas_priority__`` with a higher value.

By default, pandas objects try to operate with other objects, even if they are not types known to pandas:

.. code-block:: python

    >>> pd.Series([1, 2]) + [10, 20]
    0    11
    1    22
    dtype: int64

In the example above, if ``[10, 20]`` was a custom type that can be understood as a list, pandas objects will still operate with it in the same way.

In some cases, it is useful to delegate to the other type the operation. For example, consider I implement a
custom list object, and I want the result of adding my custom list with a pandas :class:`Series` to be an instance of my list
and not a :class:`Series` as seen in the previous example. This is now possible by defining the ``__pandas_priority__`` attribute
of my custom list, and setting it to a higher value, than the priority of the pandas objects I want to operate with.

The ``__pandas_priority__`` of :class:`DataFrame`, :class:`Series`, and :class:`Index` are ``4000``, ``3000``, and ``2000`` respectively.  The base ``ExtensionArray.__pandas_priority__`` is ``1000``.

.. code-block:: python

    class CustomList(list):
        __pandas_priority__ = 5000

        def __radd__(self, other):
            # return `self` and not the addition for simplicity
            return self

    custom = CustomList()
    series = pd.Series([1, 2, 3])

    # Series refuses to add custom, since it's an unknown type with higher priority
    assert series.__add__(custom) is NotImplemented

    # This will cause the custom class `__radd__` being used instead
    assert series + custom is custom
