.. _copy_on_write:

{{ header }}

*******************
Copy-on-Write (CoW)
*******************

Copy-on-Write was first introduced in version 1.5.0. Starting from version 2.0 most of the
optimizations that become possible through CoW are implemented and supported. A complete list
can be found at :ref:`Copy-on-Write optimizations <copy_on_write.optimizations>`.

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

.. _copy_on_write.optimizations:

Copy-on-Write optimizations
---------------------------

A new lazy copy mechanism that defers the copy until the object in question is modified
and only if this object shares data with another object. This mechanism was added to
following methods:

  - :meth:`DataFrame.reset_index` / :meth:`Series.reset_index`
  - :meth:`DataFrame.set_index`
  - :meth:`DataFrame.set_axis` / :meth:`Series.set_axis`
  - :meth:`DataFrame.set_flags` / :meth:`Series.set_flags`
  - :meth:`DataFrame.rename_axis` / :meth:`Series.rename_axis`
  - :meth:`DataFrame.reindex` / :meth:`Series.reindex`
  - :meth:`DataFrame.reindex_like` / :meth:`Series.reindex_like`
  - :meth:`DataFrame.assign`
  - :meth:`DataFrame.drop`
  - :meth:`DataFrame.dropna` / :meth:`Series.dropna`
  - :meth:`DataFrame.select_dtypes`
  - :meth:`DataFrame.align` / :meth:`Series.align`
  - :meth:`Series.to_frame`
  - :meth:`DataFrame.rename` / :meth:`Series.rename`
  - :meth:`DataFrame.add_prefix` / :meth:`Series.add_prefix`
  - :meth:`DataFrame.add_suffix` / :meth:`Series.add_suffix`
  - :meth:`DataFrame.drop_duplicates` / :meth:`Series.drop_duplicates`
  - :meth:`DataFrame.droplevel` / :meth:`Series.droplevel`
  - :meth:`DataFrame.reorder_levels` / :meth:`Series.reorder_levels`
  - :meth:`DataFrame.between_time` / :meth:`Series.between_time`
  - :meth:`DataFrame.filter` / :meth:`Series.filter`
  - :meth:`DataFrame.head` / :meth:`Series.head`
  - :meth:`DataFrame.tail` / :meth:`Series.tail`
  - :meth:`DataFrame.isetitem`
  - :meth:`DataFrame.pipe` / :meth:`Series.pipe`
  - :meth:`DataFrame.pop` / :meth:`Series.pop`
  - :meth:`DataFrame.replace` / :meth:`Series.replace`
  - :meth:`DataFrame.shift` / :meth:`Series.shift`
  - :meth:`DataFrame.sort_index` / :meth:`Series.sort_index`
  - :meth:`DataFrame.sort_values` / :meth:`Series.sort_values`
  - :meth:`DataFrame.squeeze` / :meth:`Series.squeeze`
  - :meth:`DataFrame.swapaxes`
  - :meth:`DataFrame.swaplevel` / :meth:`Series.swaplevel`
  - :meth:`DataFrame.take` / :meth:`Series.take`
  - :meth:`DataFrame.to_timestamp` / :meth:`Series.to_timestamp`
  - :meth:`DataFrame.to_period` / :meth:`Series.to_period`
  - :meth:`DataFrame.truncate`
  - :meth:`DataFrame.iterrows`
  - :meth:`DataFrame.tz_convert` / :meth:`Series.tz_localize`
  - :meth:`DataFrame.fillna` / :meth:`Series.fillna`
  - :meth:`DataFrame.interpolate` / :meth:`Series.interpolate`
  - :meth:`DataFrame.ffill` / :meth:`Series.ffill`
  - :meth:`DataFrame.bfill` / :meth:`Series.bfill`
  - :meth:`DataFrame.where` / :meth:`Series.where`
  - :meth:`DataFrame.infer_objects` / :meth:`Series.infer_objects`
  - :meth:`DataFrame.astype` / :meth:`Series.astype`
  - :meth:`DataFrame.convert_dtypes` / :meth:`Series.convert_dtypes`
  - :meth:`DataFrame.join`
  - :meth:`DataFrame.eval`
  - :func:`concat`
  - :func:`merge`

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
