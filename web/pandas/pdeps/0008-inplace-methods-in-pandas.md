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

- The ``inplace`` parameter will be removed from any methods that never can be done inplace
- The ``inplace`` parameter will also be removed from any methods that modify the shape of a pandas object's values or
  don't modify the internal values of a pandas object at all.
- In contrast, the ``inplace`` parameter will be kept for any methods that only modify the underlying data of a pandas
  object.
    - For example, the ``fillna`` method would retain its ``inplace`` keyword, while ``dropna`` (potentially shrinks the
      length of a ``DataFrame``/``Series``) and ``rename`` (alters labels not values) would lose their ``inplace``
      keywords
    - For those methods, since Copy-on-Write behavior will lazily copy if the result is unchanged, users should reassign
      to the same variable to imitate behavior of the ``inplace`` keyword.
      e.g. ``df = df.dropna()`` for a DataFrame with no null values.
- The ``copy`` parameter will also be removed, except in constructors and in functions/methods that convert array-likes
  to pandas objects (e.g. the ``pandas.array`` function) and functions/methods that export pandas objects to other data
  types (e.g. ``DataFrame/Series.to_numpy`` method).
- Open Questions
  (These questions are deferred to a later revision, and will not affect the acceptance process of this PDEP.)
    - Should ``inplace=True`` return the original pandas object that was operated inplace on?
    - What should happen when ``inplace=True`` but the original pandas object cannot be operated inplace on because it
      shares its values with another pandas object?

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
others still make a copy under the hood anyway. In addition, with the introduction of Copy-on-Write, there are now other
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

Given the above reasons, we are convinced that there is no need for neither the `inplace` nor the `copy` keyword (except
for a small subset of methods that can actually update data inplace). Removing those keywords will give a more
consistent and less confusing API. Removing the `copy` keyword is covered by PDEP-7 about Copy-on-Write,
and this PDEP will focus on the `inplace` keyword.

Thus, in this PDEP, we aim to standardize behavior across methods to make control of inplace-ness of operations
consistent, and compatible with Copy-on-Write.

Note: there are also operations (not methods) that work inplace in pandas, such as indexing (
e.g. `df.loc[0, "col"] = val`) or inplace operators (e.g. `df += 1`). This is out of scope for this PDEP, as we focus on
the inplace behaviour of DataFrame and Series _methods_.

## Detailed description

### Status Quo

Many methods in pandas currently have the ability to perform an operation inplace. For example, some methods such
as ``DataFrame.insert`` only support inplace operations, while other methods use the `inplace` keyword to control
whether an operation is done inplace or not.

Unfortunately, many methods supporting the ``inplace`` keyword either cannot be done inplace, or make a copy as a
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
| ``isetitem``* |

\* Although ``isetitem`` operates on the original pandas object inplace, it will not change any existing values
inplace (it will remove the values of the column being set, and insert new values).

**Group 2: Methods that modify the underlying data of the DataFrame/Series object and can be done inplace**

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

These methods don't operate inplace by default, but can be done inplace with `inplace=True` if the dtypes are compatible
(e.g. the values replacing the old values can be stored in the original array without an astype). All those methods leave
the structure of the DataFrame or Series intact (shape, row/column labels), but can mutate some elements of the data of
the DataFrame or Series.

**Group 3: Methods that modify the DataFrame/Series object, but not the pre-existing values**

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

All those methods make a copy of the full data by default, but can be performed inplace with
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

Although these methods have the `inplace` keyword, they can never operate inplace because the nature of the
operation requires copying (such as reordering or dropping rows). For those methods, `inplace=True` is essentially just
syntactic sugar for reassigning the new result to the calling DataFrame/Series.

Note: in the case of a "no-op" (for example when sorting an already sorted DataFrame), some of those methods might not
need to perform a copy. This currently happens with Copy-on-Write (regardless of `inplace`), but this is considered an
implementation detail for the purpose of this PDEP.

### Proposed changes and reasoning

The methods from group 1 won't change behavior, and will remain always inplace.

Methods in groups 3 and 4 will lose their `inplace` keyword. Under Copy-on-Write, every operation will
potentially return a shallow copy of the input object, if the performed operation does not require a copy of the data. This is
equivalent to the behavior with `inplace=True` for those methods. If users want to make a hard
copy, they can call the `copy()` method on the result of the operation.

Therefore, there is no benefit of keeping the keywords around for these methods.

To emulate behavior of the `inplace` keyword, we can reassign the result of an operation to the same variable:

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    df = df.reset_index()
    df.iloc[0, 1] = ...

All references to the original object will go out of scope when the result of the `reset_index` operation is assigned
to `df`. As a consequence, `iloc` will continue to operate inplace, and the underlying data will not be copied.

Group 2 methods differ, though, since they only modify the underlying data, and therefore can be inplace.

    :::python
    df = pd.DataFrame({"foo": [1, 2, 3]})
    df = df.replace(to_replace=1, value=100)

If we follow the rules of Copy-on-Write[^1] where "any subset or returned series/dataframe always behaves as a copy of
the original, and thus never modifies the original", then there is no way of doing this operation inplace by default.
The original object would be modified before the reference goes out of scope.

To avoid triggering a copy when a value would actually get replaced, we will keep the `inplace` argument for those
methods.

### Open Questions

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

There is another possible variant, which would be to trigger the copy (like the first option), but have an option to
raise a warning whenever this happens.
This would be useful in an IPython shell/Jupyter Notebook setting, where the user would have the opportunity to delete
unused references that are causing the copying to be triggered.

For example,

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

The downsides of keeping the `inplace=True` option for certain methods, are that the return type of those methods will
now depend on the value of `inplace`, and that method chaining will no longer work.

One way around this is to have the method return the original object that was operated on inplace when `inplace=True`.

Advantages:

- It enables to use inplace operations in a method chain
- It simplifies type annotations
- It enables to change the default for `inplace` to True under Copy-on-Write

Disadvantages:

- In general, when a pandas method returns an object, this is a _new_ object, and thus following the Copy-on-Write rules
  of behaving as a copy. This would introduce a special case where an _identical_ object would be
  returned (`df2 = df.method(inplace=True); assert df2 is df`)
- It would change the behaviour of the current `inplace=True`

Given that `inplace` is already widely used by the pandas community, we would like to collect feedback about what the
expected return type should be. Therefore, we will defer a decision on this until a later revision of this PDEP.

## Backward compatibility

Removing the `inplace` keyword is a breaking change, but since the affected behaviour is `inplace=True`, the default
behaviour when not specifying the keyword (i.e. `inplace=False`) will not change and the keyword itself can first be
deprecated before it is removed.

## Rejected ideas

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
keyword already returns a new object instead of None (enabling method chaining) and avoids a coopy when it is set to `False`.

However, the `copy` keyword is not supported in any of the values-mutating methods listed in Group 2 above
unlike `inplace`, so semantics of future inplace mutation of values align better with the current behavior of
the `inplace` keyword, than with the current behavior of the `copy` keyword.

Furthermore, with the Copy-on-Write proposal, the `copy` keyword also has become superfluous. With Copy-on-Write
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

Copy-on-Write is a relatively new feature (added in version 1.5) and some methods are missing the "lazy copy"
optimization (equivalent to `copy=False`).

Therefore, we propose deprecating the `copy` and `inplace` parameters in pandas 2.1, to
allow for bugs with Copy-on-Write to be addressed and for more optimizations to be added.

Hopefully, users will be able to switch to Copy-on-Write to keep the no-copy behavior and to silence the warnings.

The full removal of the `copy` parameter and `inplace` (where necessary) is set for pandas 3.0, which will coincide
with the enablement of Copy-on-Write for pandas by default.

## PDEP History

- 16 February 2023: Initial draft

## References

[^1]: [Copy on Write Specification](https://docs.google.com/document/d/1ZCQ9mx3LBMy-nhwRl33_jgcvWo9IWdEfxDNQ2thyTb0/edit#heading=h.iexejdstiz8u)
[^2]: <https://github.com/pandas-dev/pandas/issues/48141>
[^3]: <https://github.com/pandas-dev/pandas/issues/16529>
[^4]: <https://github.com/pandas-dev/pandas/issues/50535>
[^5]: <https://stackoverflow.com/questions/37808904/prevent-ipython-from-storing-outputs-in-out-variable>
