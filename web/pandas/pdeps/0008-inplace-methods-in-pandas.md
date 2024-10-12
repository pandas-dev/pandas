# PDEP-8: In-place methods in pandas

- Created: 16 February 2023
- Status: Under discussion
- Discussion: [PR 51466](https://github.com/pandas-dev/pandas/pull/51466)
- Authors: [Thomas Li](https://github.com/lithomas1),
           [Patrick Hoefler](https://github.com/phofl),
           [Joris Van den Bossche](https://github.com/jorisvandenbossche)
- Revision: 1

## Abstract

This PDEP proposes that:

- The ``inplace`` parameter will be removed from any method which can never update the
  underlying values of a pandas object inplace or which alters the shape of the object,
  and where the `inplace=True` option is only syntactic sugar for reassigning the result
  to the calling DataFrame/Series.
- As a consequence, the `inplace` parameter is only kept for those methods that can
  modify the underlying values of a pandas object inplace, such as `fillna` or `replace`.
- With the introduction of Copy-on-Write ([PDEP-7](^1)), users don't need the `inplace`
  keyword to avoid a copy of the data.
- For those methods that will keep the `inplace=True` option:
    - the method will do an attempt to do the operation inplace but still silently copy
      when needed (for Copy-on-Write), i.e. there is no guarantee it is actually done inplace.
    - the method will return the calling object (`self`), instead of the current `None`.

## Motivation and Scope

The `inplace=True` keyword has been a controversial topic for many years. It is generally seen (at least by several
pandas maintainers and educators) as bad practice and often unnecessary, but at the same time it is also widely used,
partly because of confusion around the impact of the keyword.

Generally, we assume that people use the keyword for the following reasons:

1. Because they think it is more efficient (it is faster and/or can save memory)
2. To save the result to the same variable / update the original variable (avoid the pattern of reassigning to the same
   variable name)

For the first reason: efficiency is an important aspect. However, in practice it is not always the case
that `inplace=True` improves anything. Some of the methods with an `inplace` keyword can actually work inplace, but
others still make a copy under the hood anyway. In addition, with the introduction of Copy-on-Write ([PDEP-7](^1)), there are now other
ways to avoid making unnecessary copies by default (without needing to specify a keyword). The next section gives a
detailed overview of those different cases.

For the second reason: we are convinced that this is not worth it. While it might save some keystrokes (if you have a
long variable name), this code style also has sufficient disadvantages that we think it is not worth providing "two
ways" to achieve the same result:

- You can't use method chaining with `inplace=True`
- The ``inplace`` keyword complicates type annotations (because the return value depends on the value of `inplace`)
- Using `inplace=True` gives code that mutates the state of an object and thus has side-effects. That can introduce
  subtle bugs and is harder to debug.

Finally, there are also methods that have a `copy` keyword instead of an `inplace` keyword (which also avoids copying
the data when `copy=False`, but returns a new object referencing the same data instead of updating the calling object),
adding to the inconsistencies. This keyword is also redundant now with the introduction of Copy-on-Write.

Given the above reasons, we are convinced that there is no need for neither the `inplace` nor the `copy` keyword, except
for a small subset of methods that can actually update data inplace. Removing those keywords will give a more
consistent and less confusing API. Removing the `copy` keyword is covered by PDEP-7 about Copy-on-Write,
and this PDEP will focus on the `inplace` keyword.

Thus, in this PDEP, we aim to standardize behavior across methods to make control of inplace-ness of methods
consistent, and compatible with Copy-on-Write.

Note: there are also operations (not methods) that work inplace in pandas, such as indexing (
e.g. `df.loc[0, "col"] = val`) or inplace operators (e.g. `df += 1`). This is out of scope for this PDEP, as we focus on
the inplace behaviour of DataFrame and Series _methods_.

## Detailed description

### Status Quo

Many methods in pandas currently have the ability to perform an operation inplace. For example, some methods such
as ``DataFrame.insert`` only support inplace operations, while other methods use the `inplace` keyword to control
whether an operation is done inplace or not.

While we generally speak about "inplace" operations, this term is used in various context. Broadly speaking,
for this PDEP, we can distinguish two kinds of "inplace" operations:

* **"values-inplace"**: an operation that updates the underlying values of a Series or DataFrame columns inplace
  (without making a copy of the array).

  As illustration, an example of such a values-inplace operation without using a method:

    :::python
    # if the dtype is compatible, this setitem operation updates the underlying array inplace
    df.loc[0, "col"] = val

* **"object-inplace"**: an operation that updates a pandas DataFrame or Series _object_ inplace, but without
  updating existing column values inplace.

  As illustration, an example of such an object-inplace operation without using a method:

    :::python
    # we replace the Index on `df` inplace, but without actually
    # updating any existing array
    df.index = pd.Index(...)
    # we update the DataFrame inplace, but by completely replacing a column,
    # not by mutating the existing column's underlying array
    df["col"] = new_values

  Object-inplace operations, while not actually modifying existing column values, keep
  (a subset of) those columns and thus can avoid copying the data of those existing columns.

In addition, several methods supporting the ``inplace`` keyword cannot actually be done inplace (in either meaning)
because they make a copy as a
consequence of the operations they perform, regardless of whether ``inplace`` is ``True`` or not. This, coupled with the
fact that the ``inplace=True`` changes the return type of a method from a pandas object to ``None``, makes usage of
the ``inplace`` keyword confusing and non-intuitive.

To summarize the status quo of inplace behavior of methods, we have divided methods that can operate inplace or have
an ``inplace`` keyword into 4 groups:

**Group 1: Methods that always operate inplace (no user-control with ``inplace`` keyword)**

| Method Name   |
|:--------------|
| ``insert``    |
| ``pop``       |
| ``update``    |
| ``isetitem``  |

This group encompasses both kinds of inplace: `update` can be values-inplace, while the others are object-inplace
(for example, although ``isetitem`` operates on the original pandas object inplace,
it will not change any existing values inplace; rather it will remove the values of the column being set, and insert new values).

**Group 2: Methods that can modify the underlying data of the DataFrame/Series object ("values-inplace")**

| Method Name     |
|:----------------|
| ``where``       |
| ``fillna``      |
| ``replace``     |
| ``mask``        |
| ``interpolate`` |
| ``ffill``       |
| ``bfill``       |
| ``clip``        |

These methods don't operate inplace by default, but can be done inplace with `inplace=True` _if_ the dtypes are compatible
(e.g. the values replacing the old values can be stored in the original array without an astype). All those methods leave
the structure of the DataFrame or Series intact (shape, row/column labels), but can mutate some elements of the data of
the DataFrame or Series.

**Group 3: Methods that can modify the DataFrame/Series object, but not the pre-existing values ("object-inplace")**

| Method Name                 |
|:----------------------------|
| ``drop`` (dropping columns) |
| ``rename``                  |
| ``rename_axis``             |
| ``reset_index``             |
| ``set_index``               |

These methods can change the structure of the DataFrame or Series, such as changing the shape by adding or removing
columns, or changing the row/column labels (changing the index/columns attributes), but don't modify the existing
underlying column data of the object.

All those methods make a copy of the full data by default, but can be performed object-inplace with
avoiding copying all data (currently enabled with specifying `inplace=True`).

Note: there are also methods that have a `copy` keyword instead of an `inplace` keyword (e.g. `set_axis`). This serves
a similar purpose (avoid copying all data), but those methods don't update the original object inplace and instead
return a new object referencing the same data.

**Group 4: Methods that can never operate inplace**

| Method Name            |
|:-----------------------|
| `drop` (dropping rows) |
| `dropna`               |
| `drop_duplicates`      |
| `sort_values`          |
| `sort_index`           |
| `eval`                 |
| `query`                |

Although these methods have the `inplace` keyword, they can never operate inplace, in either meaning, because the nature of the
operation requires copying (such as reordering or dropping rows). For those methods, `inplace=True` is essentially just
syntactic sugar for reassigning the new result to the calling DataFrame/Series.

Note: in the case of a "no-op" (for example when sorting an already sorted DataFrame), some of those methods might not
need to perform a copy and could be considered as "object-inplace" in that case.
This currently happens with Copy-on-Write (regardless of `inplace`), but this is considered an
implementation detail for the purpose of this PDEP.

### Proposed changes and reasoning

The methods from **group 1 (always inplace, no keyword)** won't change behavior, and will remain always inplace.

For methods from **group 4 (never inplace)**, the `inplace` keyword has no actual effect
(except for re-assigning to the calling variable) and is effectively syntactic sugar for
manually re-assigning. For this group, we propose to remove the `inplace` keyword.

For methods from **group 3 (object-inplace)**, the `inplace=True` keyword can currently be
used to avoid a copy. However, with the introduction of Copy-on-Write, every operation
will potentially return a shallow copy of the input object by default (if the performed
operation does not require a copy of the data). This future default is therefore
equivalent to the behavior with `inplace=True` for those methods (minus the return
value).

For the above reasoning, we think there is no benefit of keeping the keyword around for
these methods. To emulate behavior of the `inplace` keyword, we can reassign the result
of an operation to the same variable:

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    df = df.reset_index()
    df.iloc[0, 1] = ...

All references to the original object will go out of scope when the result of the `reset_index` operation is re-assigned
to `df`. As a consequence, `iloc` will continue to operate inplace, and the underlying data will not be copied (with Copy-on-Write).

**Group 2 (values-inplace)** methods differ, though, since they modify the underlying
data, and therefore can actually happen inplace:

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    df.replace(to_replace=1, value=100, inplace=True)

Currently, the above updates `df` values-inplace, without requiring a copy of the data.
For this type of method, however, we can _not_ emulate the above usage of `inplace` by
re-assigning:

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    df = df.replace(to_replace=1, value=100)

If we follow the rules of Copy-on-Write[^1] where "any subset or returned
series/dataframe always behaves as a copy of the original, and thus never modifies the
original", then there is no way of doing this operation inplace by default, because the
original object `df` would be modified before the reference goes out of scope (pandas
does not know whether you will re-assign it to `df` or assign it to another variable).
That would violate the Copy-on-Write rules, and therefore the `replace()` method in the
example always needs to make a copy of the underlying data by default

For this case, an `inplace=True` option can have an actual benefit, i.e. allowing to
avoid a data copy. Therefore, we propose to keep the `inplace` argument for this
group of methods.

Summarizing for the `inplace` keyword, we propose to:

- Keep the `inplace` keyword for this subset of methods (group 2) that can update the
  underlying values inplace ("values-inplace")
- Remove the `inplace` keyword from all other methods that either can never work inplace
  (group 4) or only update the object (group 3, "object-inplace", which can be emulated
  with reassigning).

### Other design questions

#### With `inplace=True`, should we silently copy or raise an error if the data has references?

For those methods where we would keep the `inplace=True` option (group 2), there is a complication that actually operating inplace
is not always possible.

For example,

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    df.replace(to_replace=1, value=100, inplace=True)

can be performed inplace.

This is only true if `df` does not share the values it stores with another pandas object. For example, the following
operations

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    view = df[:]
    # We can't operate inplace, because view would also be modified!
    df.replace(to_replace=1, value=100, inplace=True)

would be incompatible with the Copy-on-Write rules when actually done inplace. In this case we can either

- copy the shared values before performing the operation to avoid modifying another object (i.e. follow the standard
  Copy-on-Write procedure),
- raise an error to indicate that more than one object would be changed and the inplace operation is not possible.

Raising an error here is problematic since oftentimes users do not have control over whether a method would cause a "
lazy copy" to be triggered under Copy-on-Write. It is also hard to fix, adding a `copy()` before calling a method
with `inplace=True` might actually be worse than triggering the copy under the hood. We would only copy columns that
share data with another object, not the whole object like `.copy()` would.

**Therefore, we propose to silently copy when needed.** The `inplace=True` option would thus mean "try inplace whenever possible", and not guarantee it is actually done inplace.

In the future, if there is demand for it, it could still be possible to add to option to raise a warning whenever this happens.
This would be useful in an IPython shell/Jupyter Notebook setting, where the user would have the opportunity to delete
unused references that are causing the copying to be triggered.

For example:

    :::ipython
    In [1]: import pandas as pd

    In [2]: pd.set_option("mode.copy_on_write", True)

    In [3]: ser = pd.Series([1,2,3])

    In [4]: ser_vals = ser.values # Save values to check inplace-ness

    In [5]: ser
    Out[5]:
    0    1
    1    2
    2    3
    dtype: int64

    In [6]: ser = ser[:] # Original series should go out of scope

    In [7]: ser.iloc[0] = -1 # This should be inplace

    In [8]: ser
    Out[8]:
    0   -1
    1    2
    2    3
    dtype: int64

    In [9]: ser_vals
    Out[9]: array([1, 2, 3]) # It's not modified!

    In [10]: Out[5] # IPython kept our series alive since we displayed it!

While there are ways to mitigate this[^5], it may be helpful to let the user know that an operation that they performed
was not inplace, since it is possible to go out of memory because of this.

#### Return the calling object (`self`) also when using `inplace=True`?

One of the downsides of the `inplace=True` option is that the return type of those methods
depends on the value of `inplace`, and that method chaining does not work.
Those downsides are still relevant for the cases where we keep `inplace=True`.
To address this, we can have those methods return the object that was operated on
inplace when `inplace=True`.

Advantages:

- It enables to use inplace operations in a method chain
- It simplifies type annotations

Disadvantages:

- In general, when a pandas method returns an object, this is a _new_ object, and thus following the Copy-on-Write rules
  of behaving as a copy. This would introduce a special case where an _identical_ object would be
  returned (`df2 = df.method(inplace=True); assert df2 is df`)
- It would change the behaviour of the current `inplace=True`

We generally assume that changing to return `self` should not give much problems for
existing usage (typically, the current return value of `None` is not actively used).
Further, we think the advantages of simplifing return types and enabling methods chains
outweighs the special case of returning an identical object.
**Therefore, we propose that for those methods with an `inplace=True` option, the calling object (`self`) gets returned.**

## Backward compatibility

Removing the `inplace` keyword is a breaking change, but since the affected behaviour is `inplace=True`, the default
behaviour when not specifying the keyword (i.e. `inplace=False`) will not change and the keyword itself can first be
deprecated before it is removed.

## Rejected alternatives

### Remove the `inplace` keyword altogether

In the past, it was considered to remove the `inplace` keyword entirely. This was because many methods with
the `inplace` keyword did not actually operate inplace, but made a copy and re-assigned the underlying values under
the hood, causing confusion and providing no real benefit to users.

Because a majority of the methods supporting `inplace` did not operate inplace, it was considered at the time to
deprecate and remove inplace from all methods, and add back the keyword as necessary.[^3]

For methods where the operation actually _can_ be done inplace (group 2), however, removing the `inplace`
keyword could give a significant performance regression when currently using this keyword with large
DataFrames. Therefore, we decided to keep the `inplace` keyword for this small subset of methods.

### Standardize on the `copy` keyword instead of `inplace`

It may seem more natural to standardize on the `copy` keyword instead of the `inplace` keyword, since the `copy`
keyword already returns a new object instead of None (enabling method chaining) and avoids a copy when it is set to `False`.

However, the `copy` keyword is not supported in any of the values-mutating methods listed in Group 2 above
unlike `inplace`, so semantics of future inplace mutation of values align better with the current behavior of
the `inplace` keyword, than with the current behavior of the `copy` keyword.

Furthermore, with the approved Copy-on-Write proposal, the `copy` keyword also has become superfluous. With Copy-on-Write
enabled, methods that return a new pandas object will always try to avoid a copy whenever possible, regardless of
a `copy=False` keyword. Thus, the Copy-on-Write PDEP proposes to actually remove the `copy` keyword from the methods
where it is currently used (so it would be strange to add this as a new keyword to the Group 2 methods).

Currently, when using `copy=False` in methods where it is supported, a new pandas object is returned as the result
of a method call (same as with `copy=True`), but with the values backing this object being shared with the calling
object when possible (but the calling object is never modified). With the proposed inplace behavior for Group 2 methods,
a potential `copy=False` option would return a new pandas object with identical values as the original object (that
was modified inplace, in contrast to current usage of `copy=False`), which may be confusing for users, and lead to
ambiguity with Copy on Write rules.

## History

The future of the `inplace` keyword is something that has been debated a lot over the years.

It may be helpful to review those discussions (see links) [^2] [^3] [^4] to better understand this PDEP.

## Timeline

The `inplace` keyword is widely used, and thus we need to take considerable time to
deprecate and remove this feature.

- For those methods where the `inplace` keyword will be removed, we add a
  DeprecationWarning in the first release after acceptance (2.2 if possible, otherwise
  3.0)
- Together with enabling Copy-on-Write in the pandas 3.0 major release, we already
  update those methods that will keep the `inplace` keyword with the new behaviour
  (returning `self`, working inplace when possible)
- Somewhere during the 3.x release cycle (e.g. in 3.1, depending on when the deprecation
  was started), we change the DeprecationWarning to a more visible FutureWarning.
- The deprecated keyword is removed in pandas 4.0.

When introducing the warning in 2.2 (or 3.0), users will already have the ability to
enable Copy-on-Write so they can rewrite their code in a way that avoids the deprecation
warning (remove the usage of `inplace`) while keeping the no-copy behaviour (which will
be the default with Copy-on-Write).

## PDEP History

- 16 February 2023: Initial draft

## References

[^1]: [Copy on Write Specification](https://pandas.pydata.org/pdeps/0007-copy-on-write.html)
[^2]: <https://github.com/pandas-dev/pandas/issues/48141>
[^3]: <https://github.com/pandas-dev/pandas/issues/16529>
[^4]: <https://github.com/pandas-dev/pandas/issues/50535>
[^5]: <https://stackoverflow.com/questions/37808904/prevent-ipython-from-storing-outputs-in-out-variable>
