.. _copy_on_write:

{{ header }}

*******************
Copy-on-Write (CoW)
*******************

Copy-on-Write was first introduced in version 1.5.0. Starting from version 2.0 most of the
optimizations that become possible through CoW are implemented and supported. All possible
optimizations are supported starting from pandas 2.1.

We expect that CoW will be enabled by default in version 3.0.

CoW will lead to more predictable behavior since it is not possible to update more than
one object with one statement, e.g. indexing operations or methods won't have side-effects. Additionally, through
delaying copies as long as possible, the average performance and memory usage will improve.

Previous behavior
-----------------

pandas indexing behavior is tricky to understand. Some operations return views while
other return copies. Depending on the result of the operation, mutation one object
might accidentally mutate another:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    subset = df["foo"]
    subset.iloc[0] = 100
    df

Mutating ``subset``, e.g. updating its values, also updates ``df``. The exact behavior is
hard to predict. Copy-on-Write solves accidentally modifying more than one object,
it explicitly disallows this. With CoW enabled, ``df`` is unchanged:

.. ipython:: python

    pd.options.mode.copy_on_write = True

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    subset = df["foo"]
    subset.iloc[0] = 100
    df

The following sections will explain what this means and how it impacts existing
applications.

Description
-----------

CoW means that any DataFrame or Series derived from another in any way always
behaves as a copy. As a consequence, we can only change the values of an object
through modifying the object itself. CoW disallows updating a DataFrame or a Series
that shares data with another DataFrame or Series object inplace.

This avoids side-effects when modifying values and hence, most methods can avoid
actually copying the data and only trigger a copy when necessary.

The following example will operate inplace with CoW:

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

.. ipython:: python

    with pd.option_context("mode.copy_on_write", False):
        df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        view = df[:]
        df.iloc[0, 0] = 100

        df
        view

CoW triggers a copy when ``df`` is changed to avoid mutating ``view`` as well:

.. ipython:: python

    df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
    view = df[:]
    df.iloc[0, 0] = 100

    df
    view

Chained Assignment
------------------

Chained assignment references a technique where an object is updated through
two subsequent indexing operations, e.g.

.. ipython:: python

    with pd.option_context("mode.copy_on_write", False):
        df = pd.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})
        df["foo"][df["bar"] > 5] = 100
        df

The column ``foo`` is updated where the column ``bar`` is greater than 5.
This violates the CoW principles though, because it would have to modify the
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

There are two potential solution to this:

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
    df2 = df.reset_index()
    df2.iloc[0, 0] = 100

This creates two objects that share data and thus the setitem operation will trigger a
copy. This is not necessary if the initial object ``df`` isn't needed anymore.
Simply reassigning to the same variable will invalidate the reference that is
held by the object.

.. ipython:: python

    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df = df.reset_index()
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

How to enable CoW
-----------------

Copy-on-Write can be enabled through the configuration option ``copy_on_write``. The option can
be turned on __globally__ through either of the following:

.. ipython:: python

    pd.set_option("mode.copy_on_write", True)

    pd.options.mode.copy_on_write = True

.. ipython:: python
    :suppress:

    pd.options.mode.copy_on_write = False
