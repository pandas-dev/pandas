.. _migration_guides:

{{ header }}

================
Migration Guides
================

For large changes that are difficult or impossible to deprecate in a user-friendly manner,
pandas will implement the changes under the ``future`` configuration. This section
goes into detail for each of these changes.

.. _copy_on_write:

*******************
Copy-on-Write (CoW)
*******************

.. note::

    Copy-on-Write is now the default with pandas 3.0.

Copy-on-Write was first introduced in version 1.5.0. Starting from version 2.0 most of the
optimizations that become possible through CoW are implemented and supported. All possible
optimizations are supported starting from pandas 2.1.

CoW will lead to more predictable behavior since it is not possible to update more than
one object with one statement, e.g. indexing operations or methods won't have side-effects. Additionally, through
delaying copies as long as possible, the average performance and memory usage will improve.

Previous behavior
-----------------

pandas indexing behavior is tricky to understand. Some operations return views while
other return copies. Depending on the result of the operation, mutating one object
might accidentally mutate another:

.. code-block:: ipython

    In [1]: df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    In [2]: subset = df["foo"]
    In [3]: subset.iloc[0] = 100
    In [4]: df
    Out[4]:
       foo  bar
    0  100    4
    1    2    5
    2    3    6


Mutating ``subset``, e.g. updating its values, also updated ``df``. The exact behavior was
hard to predict. Copy-on-Write solves accidentally modifying more than one object,
it explicitly disallows this. ``df`` is unchanged:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    subset = df["foo"]
    subset.iloc[0] = 100
    df

The following sections will explain what this means and how it impacts existing
applications.

.. _copy_on_write.migration_guide:

Migrating to Copy-on-Write
--------------------------

Copy-on-Write is the default and only mode in pandas 3.0. This means that users
need to migrate their code to be compliant with CoW rules.

The default mode in pandas < 3.0 raises warnings for certain cases that will actively
change behavior and thus change user intended behavior.

pandas 2.2 has a warning mode

.. code-block:: python

    pd.options.mode.copy_on_write = "warn"

that will warn for every operation that will change behavior with CoW. We expect this mode
to be very noisy, since many cases that we don't expect that they will influence users will
also emit a warning. We recommend checking this mode and analyzing the warnings, but it is
not necessary to address all of these warning. The first two items of the following lists
are the only cases that need to be addressed to make existing code work with CoW.

The following few items describe the user visible changes:

**Chained assignment will never work**

``loc`` should be used as an alternative. Check the
:ref:`chained assignment section <copy_on_write_chained_assignment>` for more details.

**Accessing the underlying array of a pandas object will return a read-only view**

.. ipython:: python

    ser = pd.Series([1, 2, 3])
    ser.to_numpy()

This example returns a NumPy array that is a view of the Series object. This view can
be modified and thus also modify the pandas object. This is not compliant with CoW
rules. The returned array is set to non-writeable to protect against this behavior.
Creating a copy of this array allows modification. You can also make the array
writeable again if you don't care about the pandas object anymore.

See the section about :ref:`read-only NumPy arrays <copy_on_write_read_only_na>`
for more details.

**Only one pandas object is updated at once**

The following code snippet updated both ``df`` and ``subset`` without CoW:

.. code-block:: ipython

    In [1]: df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    In [2]: subset = df["foo"]
    In [3]: subset.iloc[0] = 100
    In [4]: df
    Out[4]:
       foo  bar
    0  100    4
    1    2    5
    2    3    6

This is not possible anymore with CoW, since the CoW rules explicitly forbid this.
This includes updating a single column as a :class:`Series` and relying on the change
propagating back to the parent :class:`DataFrame`.
This statement can be rewritten into a single statement with ``loc`` or ``iloc`` if
this behavior is necessary. :meth:`DataFrame.where` is another suitable alternative
for this case.

Updating a column selected from a :class:`DataFrame` with an inplace method will
also not work anymore.

.. ipython:: python
    :okwarning:

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df["foo"].replace(1, 5, inplace=True)
    df

This is another form of chained assignment. This can generally be rewritten in 2
different forms:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df.replace({"foo": {1: 5}}, inplace=True)
    df

A different alternative would be to not use ``inplace``:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df["foo"] = df["foo"].replace(1, 5)
    df

**Constructors now copy NumPy arrays by default**

The Series and DataFrame constructors now copies a NumPy array by default when not
otherwise specified. This was changed to avoid mutating a pandas object when the
NumPy array is changed inplace outside of pandas. You can set ``copy=False`` to
avoid this copy.

Description
-----------

CoW means that any DataFrame or Series derived from another in any way always
behaves as a copy. As a consequence, we can only change the values of an object
through modifying the object itself. CoW disallows updating a DataFrame or a Series
that shares data with another DataFrame or Series object inplace.

This avoids side-effects when modifying values and hence, most methods can avoid
actually copying the data and only trigger a copy when necessary.

The following example will operate inplace:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df.iloc[0, 0] = 100
    df

The object ``df`` does not share any data with any other object and hence no
copy is triggered when updating the values. In contrast, the following operation
triggers a copy of the data under CoW:


.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df2 = df.reset_index(drop=True)
    df2.iloc[0, 0] = 100

    df
    df2

``reset_index`` returns a lazy copy with CoW while it copies the data without CoW.
Since both objects, ``df`` and ``df2`` share the same data, a copy is triggered
when modifying ``df2``. The object ``df`` still has the same values as initially
while ``df2`` was modified.

If the object ``df`` isn't needed anymore after performing the ``reset_index`` operation,
you can emulate an inplace-like operation through assigning the output of ``reset_index``
to the same variable:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df = df.reset_index(drop=True)
    df.iloc[0, 0] = 100
    df

The initial object gets out of scope as soon as the result of ``reset_index`` is
reassigned and hence ``df`` does not share data with any other object. No copy
is necessary when modifying the object. This is generally true for all methods
listed in :ref:`Copy-on-Write optimizations <copy_on_write.optimizations>`.

Previously, when operating on views, the view and the parent object was modified:

.. code-block:: ipython

    In [1]: df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    In [2]: subset = df["foo"]
    In [3]: subset.iloc[0] = 100
    In [4]: df
    Out[4]:
       foo  bar
    0  100    4
    1    2    5
    2    3    6

CoW triggers a copy when ``df`` is changed to avoid mutating ``view`` as well:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    view = df[:]
    df.iloc[0, 0] = 100

    df
    view

.. _copy_on_write_chained_assignment:

Chained Assignment
------------------

Chained assignment references a technique where an object is updated through
two subsequent indexing operations, e.g.

.. code-block:: ipython

    In [1]: df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    In [2]: df["foo"][df["bar"] > 5] = 100
    In [3]: df
    Out[3]:
       foo  bar
    0    1    4
    1    2    5
    2  100    6

The column ``foo`` was updated where the column ``bar`` is greater than 5.
This violated the CoW principles though, because it would have to modify the
view ``df["foo"]`` and ``df`` in one step. Hence, chained assignment will
consistently never work and raise a ``ChainedAssignmentError`` warning
with CoW enabled:

.. ipython:: python
    :okwarning:

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    df["foo"][df["bar"] > 5] = 100

With copy on write this can be done by using ``loc``.

.. ipython:: python

    df.loc[df["bar"] > 5, "foo"] = 100

.. _copy_on_write_read_only_na:

Read-only NumPy arrays
----------------------

Accessing the underlying NumPy array of a DataFrame will return a read-only array if the array
shares data with the initial DataFrame:

The array is a copy if the initial DataFrame consists of more than one array:

.. ipython:: python

    df = pd.DataFrame({"a": [1, 2], "b": [1.5, 2.5]})
    df.to_numpy()

The array shares data with the DataFrame if the DataFrame consists of only one NumPy array:

.. ipython:: python

    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    df.to_numpy()

This array is read-only, which means that it can't be modified inplace:

.. ipython:: python
    :okexcept:

    arr = df.to_numpy()
    arr[0, 0] = 100

The same holds true for a Series, since a Series always consists of a single array.

There are two potential solutions to this:

- Trigger a copy manually if you want to avoid updating DataFrames that share memory with your array.
- Make the array writeable. This is a more performant solution but circumvents Copy-on-Write rules, so
  it should be used with caution.

.. ipython:: python

    arr = df.to_numpy()
    arr.flags.writeable = True
    arr[0, 0] = 100
    arr

Patterns to avoid
-----------------

No defensive copy will be performed if two objects share the same data while
you are modifying one object inplace.

.. ipython:: python

    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df2 = df.reset_index(drop=True)
    df2.iloc[0, 0] = 100

This creates two objects that share data and thus the setitem operation will trigger a
copy. This is not necessary if the initial object ``df`` isn't needed anymore.
Simply reassigning to the same variable will invalidate the reference that is
held by the object.

.. ipython:: python

    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df = df.reset_index(drop=True)
    df.iloc[0, 0] = 100

No copy is necessary in this example.
Creating multiple references keeps unnecessary references alive
and thus will hurt performance with Copy-on-Write.

.. _copy_on_write.optimizations:

Copy-on-Write optimizations
---------------------------

A new lazy copy mechanism that defers the copy until the object in question is modified
and only if this object shares data with another object. This mechanism was added to
methods that don't require a copy of the underlying data. Popular examples are :meth:`DataFrame.drop` for ``axis=1``
and :meth:`DataFrame.rename`.

These methods return views when Copy-on-Write is enabled, which provides a significant
performance improvement compared to the regular execution.

.. _string_migration_guide:

************************
The new string data type
************************

The upcoming pandas 3.0 release introduces a new, default string data type. This
will most likely cause some work when upgrading to pandas 3.0, and this page
provides an overview of the issues you might run into and gives guidance on how
to address them.

This new dtype is already available in the pandas 2.3 release, and you can
enable it with:

.. code-block:: python

    pd.options.future.infer_string = True

This allows you to test your code before the final 3.0 release.

.. note::

   This migration guide focuses on the changes and migration steps needed when
   you are currently using ``object`` dtype for string data, which is used by
   default in pandas < 3.0. If you are already using one of the opt-in string
   dtypes, you can continue to do so without change.
   See :ref:`string_migration_guide-for_existing_users` for more details.

Background
----------

Historically, pandas has always used the NumPy ``object`` dtype as the default
to store text data. This has two primary drawbacks. First, ``object`` dtype is
not specific to strings: any Python object can be stored in an ``object``-dtype
array, not just strings, and seeing ``object`` as the dtype for a column with
strings is confusing for users. Second, this is not always very efficient (both
performance wise and for memory usage).

Since pandas 1.0, an opt-in string data type has been available, but this has
not yet been made the default, and uses the ``pd.NA`` scalar to represent
missing values.

Pandas 3.0 changes the default dtype for strings to a new string data type,
a variant of the existing optional string data type but using ``NaN`` as the
missing value indicator, to be consistent with the other default data types.

To improve performance, the new string data type will use the ``pyarrow``
package by default, if installed (and otherwise it uses object dtype under the
hood as a fallback).

See `PDEP-14: Dedicated string data type for pandas 3.0 <https://pandas.pydata.org/pdeps/0014-string-dtype.html>`__
for more background and details.

.. - brief primer on the new dtype

.. - Main characteristics:
..    - inferred by default (Default inference of a string dtype)
..    - only strings (setitem with non string fails)
..    - missing values sentinel is always NaN and uses NaN semantics

.. - Breaking changes:
..    - dtype is no longer object dtype
..    - None gets coerced to NaN
..    - setitem raises an error for non-string data

Brief introduction to the new default string dtype
--------------------------------------------------

By default, pandas will infer this new string dtype instead of object dtype for
string data (when creating pandas objects, such as in constructors or IO
functions).

Being a default dtype means that the string dtype will be used in IO methods or
constructors when the dtype is being inferred and the input is inferred to be
string data:

.. code-block:: python

   >>> pd.Series(["a", "b", None])
   0      a
   1      b
   2    NaN
   dtype: str

It can also be specified explicitly using the ``"str"`` alias:

.. code-block:: python

   >>> pd.Series(["a", "b", None], dtype="str")
   0      a
   1      b
   2    NaN
   dtype: str

Similarly, functions like :func:`read_csv`, :func:`read_parquet`, and others
will now use the new string dtype when reading string data.

In contrast to the current object dtype, the new string dtype will only store
strings. This also means that it will raise an error if you try to store a
non-string value in it (see below for more details).

Missing values with the new string dtype are always represented as ``NaN`` (``np.nan``),
and the missing value behavior is similar to other default dtypes.

This new string dtype should otherwise behave the same as the existing
``object`` dtype users are used to. For example, all string-specific methods
through the ``str`` accessor will work the same:

.. code-block:: python

   >>> ser = pd.Series(["a", "b", None], dtype="str")
   >>> ser.str.upper()
   0    A
   1    B
   2  NaN
   dtype: str

.. note::

   The new default string dtype is an instance of the :class:`pandas.StringDtype`
   class. The dtype can be constructed as ``pd.StringDtype(na_value=np.nan)``,
   but for general usage we recommend to use the shorter ``"str"`` alias.

.. _string_migration_guide-differences:

Overview of behavior differences and how to address them
---------------------------------------------------------

The dtype is no longer a numpy "object" dtype
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When inferring or reading string data, the data type of the resulting DataFrame
column or Series will silently start being the new ``"str"`` dtype instead of
the numpy ``"object"`` dtype, and this can have some impact on your code.

The new string dtype is a pandas data type ("extension dtype"), and no longer a
numpy ``np.dtype`` instance. Therefore, passing the dtype of a string column to
numpy functions will no longer work (e.g. passing it to a ``dtype=`` argument
of a numpy function, or using ``np.issubdtype`` to check the dtype).

Checking the dtype
^^^^^^^^^^^^^^^^^^

When checking the dtype, code might currently do something like:

.. code-block:: python

   >>> ser = pd.Series(["a", "b", "c"])
   >>> ser.dtype == "object"

to check for columns with string data (by checking for the dtype being
``"object"``). This will no longer work in pandas 3+, since ``ser.dtype`` will
now be ``"str"`` with the new default string dtype, and the above check will
return ``False``.

To check for columns with string data, you should instead use:

.. code-block:: python

   >>> ser.dtype == "str"

**How to write compatible code**

For code that should work on both pandas 2.x and 3.x, you can use the
:func:`pandas.api.types.is_string_dtype` function:

.. code-block:: python

   >>> pd.api.types.is_string_dtype(ser.dtype)
   True

This will return ``True`` for both the object dtype and the string dtypes.

Hardcoded use of object dtype
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have code where the dtype is hardcoded in constructors, like

.. code-block:: python

   >>> pd.Series(["a", "b", "c"], dtype="object")

this will keep using the object dtype. You will want to update this code to
ensure you get the benefits of the new string dtype.

**How to write compatible code?**

First, in many cases it can be sufficient to remove the specific data type, and
let pandas do the inference. But if you want to be specific, you can specify the
``"str"`` dtype:

.. code-block:: python

   >>> pd.Series(["a", "b", "c"], dtype="str")

This is actually compatible with pandas 2.x as well, since in pandas < 3,
``dtype="str"`` was essentially treated as an alias for object dtype.

.. attention::

   While using ``dtype="str"`` in constructors is compatible with pandas 2.x,
   specifying it as the dtype in :meth:`~Series.astype` runs into the issue
   of also stringifying missing values in pandas 2.x. See the section
   :ref:`string_migration_guide-astype_str` for more details.

.. _string_migration.select_dtypes:

For selecting string columns with :meth:`~DataFrame.select_dtypes` in a pandas
2.x and 3.x compatible way, it is not possible to use ``"str"``. While this
works for pandas 3.x, it raises an error in pandas 2.x.
As an alternative, you can select both ``object`` (for pandas 2.x) and
``"string"`` (for pandas 3.x; which will also select the default ``str`` dtype
and does not error on pandas 2.x):

.. code-block:: python

   # can use ``include=["str"]`` for pandas >= 3
   >>> df.select_dtypes(include=["object", "string"])


The missing value sentinel is now always NaN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using object dtype, multiple possible missing value sentinels are
supported, including ``None`` and ``np.nan``. With the new default string dtype,
the missing value sentinel is always NaN (``np.nan``):

.. code-block:: python

   # with object dtype, None is preserved as None and seen as missing
   >>> ser = pd.Series(["a", "b", None], dtype="object")
   >>> ser
   0       a
   1       b
   2    None
   dtype: object
   >>> print(ser[2])
   None

   # with the new string dtype, any missing value like None is coerced to NaN
   >>> ser = pd.Series(["a", "b", None], dtype="str")
   >>> ser
   0      a
   1      b
   2    NaN
   dtype: str
   >>> print(ser[2])
   nan

Generally this should be no problem when relying on missing value behavior in
pandas methods (for example, ``ser.isna()`` will give the same result as before).
But when you relied on the exact value of ``None`` being present, that can
impact your code.

**How to write compatible code?**

When checking for a missing value, instead of checking for the exact value of
``None`` or ``np.nan``, you should use the :func:`pandas.isna` function. This is
the most robust way to check for missing values, as it will work regardless of
the dtype and the exact missing value sentinel:

.. code-block:: python

   >>> pd.isna(ser[2])
   True

One caveat: this function works both on scalars and on array-likes, and in the
latter case it will return an array of bools. When using it in a Boolean context
(for example, ``if pd.isna(..): ..``) be sure to only pass a scalar to it.

"setitem" operations will now raise an error for non-string data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the new string dtype, any attempt to set a non-string value in a Series or
DataFrame will raise an error:

.. code-block:: python

   >>> ser = pd.Series(["a", "b", None], dtype="str")
   >>> ser[1] = 2.5
   ---------------------------------------------------------------------------
   TypeError                                 Traceback (most recent call last)
   ...
   TypeError: Invalid value '2.5' for dtype 'str'. Value should be a string or missing value, got 'float' instead.

If you relied on the flexible nature of object dtype being able to hold any
Python object, but your initial data was inferred as strings, your code might be
impacted by this change.

**How to write compatible code?**

You can update your code to ensure you only set string values in such columns,
or otherwise you can explicitly ensure the column has object dtype first. This
can be done by specifying the dtype explicitly in the constructor, or by using
the :meth:`~pandas.Series.astype` method:

.. code-block:: python

   >>> ser = pd.Series(["a", "b", None], dtype="str")
   >>> ser = ser.astype("object")
   >>> ser[1] = 2.5

This ``astype("object")`` call will be redundant when using pandas 2.x, but
this code will work for all versions.

Invalid unicode input
~~~~~~~~~~~~~~~~~~~~~

Python allows to have a built-in ``str`` object that represents invalid unicode
data. And since the ``object`` dtype can hold any Python object, you can have a
pandas Series with such invalid unicode data:

.. code-block:: python

   >>> ser = pd.Series(["\u2600", "\ud83d"], dtype=object)
   >>> ser
   0    ☀
   1    \ud83d
   dtype: object

However, when using the string dtype using ``pyarrow`` under the hood, this can
only store valid unicode data, and otherwise it will raise an error:

.. code-block:: python

   >>> ser = pd.Series(["\u2600", "\ud83d"])
   ---------------------------------------------------------------------------
   UnicodeEncodeError                        Traceback (most recent call last)
   ...
   UnicodeEncodeError: 'utf-8' codec can't encode character '\ud83d' in position 0: surrogates not allowed

If you want to keep the previous behaviour, you can explicitly specify
``dtype=object`` to keep working with object dtype.

When you have byte data that you want to convert to strings using ``decode()``,
the :meth:`~pandas.Series.str.decode` method now has a ``dtype`` parameter to be
able to specify object dtype instead of the default of string dtype for this use
case.

:meth:`Series.values` now returns an :class:`~pandas.api.extensions.ExtensionArray`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With object dtype, using ``.values`` on a Series will return the underlying NumPy array.

.. code-block:: python

   >>> ser = pd.Series(["a", "b", np.nan], dtype="object")
   >>> type(ser.values)
   <class 'numpy.ndarray'>

However with the new string dtype, the underlying ExtensionArray is returned instead.

.. code-block:: python

   >>> ser = pd.Series(["a", "b", pd.NA], dtype="str")
   >>> ser.values
   <ArrowStringArray>
   ['a', 'b', nan]
   Length: 3, dtype: str

If your code requires a NumPy array, you should use :meth:`Series.to_numpy`.

.. code-block:: python

   >>> ser = pd.Series(["a", "b", pd.NA], dtype="str")
   >>> ser.to_numpy()
   ['a' 'b' nan]

In general, you should always prefer :meth:`Series.to_numpy` to get a NumPy array or :meth:`Series.array` to get an ExtensionArray over using :meth:`Series.values`.

Notable bug fixes
~~~~~~~~~~~~~~~~~

.. _string_migration_guide-astype_str:

``astype(str)`` preserving missing values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The stringifying of missing values is a long standing "bug" or misfeature, as
discussed in https://github.com/pandas-dev/pandas/issues/25353, but fixing it
introduces a significant behaviour change.

With pandas < 3, when using ``astype(str)`` or ``astype("str")``, the operation
would convert every element to a string, including the missing values:

.. code-block:: python

   # OLD behavior in pandas < 3
   >>> ser = pd.Series([1.5, np.nan])
   >>> ser
   0    1.5
   1    NaN
   dtype: float64
   >>> ser.astype("str")
   0    1.5
   1    nan
   dtype: object
   >>> ser.astype("str").to_numpy()
   array(['1.5', 'nan'], dtype=object)

Note how ``NaN`` (``np.nan``) was converted to the string ``"nan"``. This was
not the intended behavior, and it was inconsistent with how other dtypes handled
missing values.

With pandas 3, this behavior has been fixed, and now ``astype("str")`` will cast
to the new string dtype, which preserves the missing values:

.. code-block:: python

   # NEW behavior in pandas 3
   >>> pd.options.future.infer_string = True
   >>> ser = pd.Series([1.5, np.nan])
   >>> ser.astype("str")
   0    1.5
   1    NaN
   dtype: str
   >>> ser.astype("str").to_numpy()
   array(['1.5', nan], dtype=object)

If you want to preserve the old behaviour of converting every object to a
string, you can use ``ser.map(str)`` instead. If you want do such conversion
while preserving the missing values in a way that works with both pandas 2.x and
3.x, you can use ``ser.map(str, na_action="ignore")`` (for pandas 3.x only, you
can do ``ser.astype("str")``).

If you want to convert to object or string dtype for pandas 2.x and 3.x,
respectively, without needing to stringify each individual element, you will
have to use a conditional check on the pandas version.
For example, to convert a categorical Series with string categories to its
dense non-categorical version with object or string dtype:

.. code-block:: python

   >>> import pandas as pd
   >>> ser = pd.Series(["a", np.nan], dtype="category")
   >>> ser.astype(object if pd.__version__ < "3" else "str")


``prod()`` raising for string data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In pandas < 3, calling the :meth:`~pandas.Series.prod` method on a Series with
string data would generally raise an error, except when the Series was empty or
contained only a single string (potentially with missing values):

.. code-block:: python

   >>> ser = pd.Series(["a", None], dtype=object)
   >>> ser.prod()
   'a'

When the Series contains multiple strings, it will raise a ``TypeError``. This
behaviour stays the same in pandas 3 when using the flexible ``object`` dtype.
But by virtue of using the new string dtype, this will generally consistently
raise an error regardless of the number of strings:

.. code-block:: python

   >>> ser = pd.Series(["a", None], dtype="str")
   >>> ser.prod()
   ---------------------------------------------------------------------------
   TypeError                                 Traceback (most recent call last)
   ...
   TypeError: Cannot perform reduction 'prod' with string dtype


.. _string_migration_guide-for_existing_users:

For existing users of the nullable ``StringDtype``
--------------------------------------------------

While pandas 3.0 introduces a new _default_ string data type, pandas had an
opt-in nullable string data type since pandas 1.0, which can be specified using
``dtype="string"``. This nullable string dtype uses ``pd.NA`` as the missing
value indicator. In addition, also through :class:`ArrowDtype` (by using
``dtypes_backend="pyarrow"``) since pandas 1.5, one could already make use of
a dedicated string dtype.

If you are already using one of the nullable string dtypes, for example by
specifying ``dtype="string"``, by using :meth:`~DataFrame.convert_dtypes`, or
by specifying the ``dtype_backend`` argument in IO functions, you can continue
to do so without change.

The migration guide above applies to code that is currently (< 3.0) using object
dtype for string data.
