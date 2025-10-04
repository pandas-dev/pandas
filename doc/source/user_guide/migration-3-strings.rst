{{ header }}

.. _string_migration_guide:

=========================================================
Migration guide for the new string data type (pandas 3.0)
=========================================================

The upcoming pandas 3.0 release introduces a new, default string data type. This
will most likely cause some work when upgrading to pandas 3.0, and this page
provides an overview of the issues you might run into and gives guidance on how
to address them.

This new dtype is already available in the pandas 2.3 release, and you can
enable it with:

.. code-block:: python

    pd.options.future.infer_string = True

This allows you to test your code before the final 3.0 release.

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

.. For existing users of the nullable ``StringDtype``
.. --------------------------------------------------

.. TODO
