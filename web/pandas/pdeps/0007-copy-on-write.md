# PDEP-7: Consistent copy/view semantics in pandas with Copy-on-Write

- Created: July 2021
- Status: Implemented
- Discussion: [#36195](https://github.com/pandas-dev/pandas/issues/36195)
- Author: [Joris Van den Bossche](https://github.com/jorisvandenbossche)
- Revision: 1

[TOC]

## Abstract

Short summary of the proposal:

1. The result of _any_ indexing operation (subsetting a DataFrame or Series in any way,
   i.e. including accessing a DataFrame column as a Series) or any method returning a
   new DataFrame or Series, always _behaves as if_ it were a copy in terms of user
   API.
2. We implement Copy-on-Write (as implementation detail). This way, we can actually use
   views as much as possible under the hood, while ensuring the user API behaves as a
   copy.
3. As a consequence, if you want to modify an object (DataFrame or Series), the only way
   to do this is to directly modify that object itself .

This addresses multiple aspects: 1) a clear and consistent user API (a clear rule: _any_
subset or returned series/dataframe **always** behaves as a copy of the original, and
thus never modifies the original) and 2) improving performance by avoiding excessive
copies (e.g. a chained method workflow would no longer return an actual data copy at each
step).

Because every single indexing step behaves as a copy, this also means that with this
proposal, "chained assignment" (with multiple setitem steps) will _never_ work and
the `SettingWithCopyWarning` can be removed.

## Background

pandas' current behavior on whether indexing returns a view or copy is confusing. Even
for experienced users, it's hard to tell whether a view or copy will be returned (see
below for a summary). We'd like to provide an API that is consistent and sensible about
returning views vs. copies.

We also care about performance. Returning views from indexing operations is faster and
reduces memory usage. The same is true for several methods that don't modify the data
such as setting/resetting the index, renaming columns, etc. that can be used in a method
chaining workflow and currently return a new copy at each step.

Finally, there are API / usability issues around views. It can be challenging to know
the user's intent in operations that modify a subset of a DataFrame (column and/or row
selection), like:

```python
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
>>> df2 = df[["A", "B"]]
>>> df2.loc[df2["A"] > 1, "A"] = 1
```

Did the user intend to modify `df` when they modified `df2` (setting aside issues with
the current implementation)? In other words, if we had a perfectly consistent world
where indexing the columns always returned views or always returned a copy, does the
code above imply that the user wants to mutate `df`?

There are two possible behaviours the user might intend:

1. Case 1: I know my subset might be a view of the original and I want to modify the
   original as well.
2. Case 2: I just want to modify the subset without modifying the original.

Today, pandas' inconsistency means _neither_ of these workflows is really possible. The
first is difficult, because indexing operations often (though not always) return copies,
and even when a view is returned you sometimes get a `SettingWithCopyWarning` when
mutating. The second is somewhat possible, but requires many defensive copies (to avoid
`SettingWithCopyWarning`, or to ensure that you have a copy when a view _was_ returned).

## Proposal

For these reasons (consistency, performance, code clarity), this PDEP proposes the
following changes:

1. The result of _any_ indexing operation (subsetting a DataFrame or Series in any way,
   i.e. including accessing a DataFrame column as a Series) or any method returning a
   new DataFrame or Series, always _behaves as if_ it were a copy in terms of user
   API.
2. We implement Copy-on-Write. This way, we can actually use views as much as possible
   under the hood, while ensuring the user API behaves as a copy.

The intent is to capture the performance benefits of views as much as possible, while
providing consistent and clear behaviour to the user. This essentially makes returning
views an internal optimization, without the user needing to know if the specific
indexing operation would return a view or a copy. The new rule would be simple: any
series/dataframe derived from another series/dataframe, through an indexing operation or
a method, always behaves as a copy of the original series/dataframe.

The mechanism to ensure this consistent behaviour, Copy-on-Write, would entail the
following: the setitem operation (i.e. `df[..] = ..` or `df.loc[..] = ..` or
`df.iloc[..] = ..`, or equivalent for Series) would check if the data that is being
modified is a view on another dataframe (or is being viewed by another dataframe). If it
is, then we would copy the data before mutating.

Taking the example from above, if the user wishes to not mutate the parent, we no longer
require a defensive copy just to avoid a `SettingWithCopyWarning`.

```python
# Case 2: The user does not want mutating df2 to mutate the parent df, via CoW
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
>>> df2 = df[["A", "B"]]
>>> df2.loc[df2["A"] > 1, "A"] = 1
>>> df.iloc[1, 0]  # df was not mutated
2
```

On the other hand, if the user actually wants to modify the original df, they can no
longer rely on the fact that `df2` could be a view, as mutating a subset would now never
mutate the parent. The only way to modify the original df is by combining all indexing
steps in a single indexing operation on the original (no "chained" setitem):

```python
# Case 1: user wants mutations of df2 to be reflected in df -> no longer possible
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
>>> df2 = df[["A", "B"]]
>>> df2.loc[df2["A"] > 1, "A"] = 1  # mutating df2 will not mutate df
>>> df.loc[df["A"] > 1, "A"] = 1  # need to directly mutate df instead
```

### This proposal also extends to methods

In principle, there's nothing special about indexing when it comes to defensive copying.
_Any_ method that returns a new series/dataframe without altering existing data (rename,
set_index, assign, dropping columns, etc.) currently returns a copy by default and is a
candidate for returning a view:

```python
>>> df2 = df.rename(columns=str.lower)
>>> df3 = df2.set_index("a")
```

Now, generally, pandas users won't expect `df2` or `df3` to be a view such that mutating
`df2` or `df3` would mutate `df`. Copy-on-Write allows us to also avoid
unnecessary copies in methods such as the above (or in the variant using method chaining
like `df.rename(columns=str.lower).set_index("a")`).

### Propagating mutation forwards

Thus far we have considered the (more common) case of taking a subset, mutating the
subset, and how that should affect the parent. What about the other direction, where the
parent is mutated?

```python
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
>>> df2 = df[["A"]]
>>> df.iloc[0, 0] = 10
>>> df2.iloc[0, 0]  # what is this value?
```

Given that `df2` is _considered_ as a copy of df under this proposal (i.e. behaves as a
copy), also mutating the parent `df` will not mutate the subset `df2`.

### When do mutations propagate to other objects and when not?

This proposal basically means that mutations _never_ propagate to _other_ objects (as
would happen with views). The only way to modify a DataFrame or Series is to modify the
object itself directly.

But let's illustrate this in Python terms. Consider that we have a DataFrame `df1`, and we
assign that to another name `df2`:

```python
>>> df1 = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
>>> df2 = df1
```

Although we have now two variables (`df1` and `df2`), this assignment follows the standard
python semantics, and both names are pointing to the same object ("df1 and df2 are
_identical_"):

```python
>>> id(df1) == id(df2)  # or: df1 is df2
True
```

Thus, if you modify DataFrame `df2`, this is also reflected in the other variable `df1`, and
the other way around (since it's the same object):

```python
>>> df1.iloc[0, 0]
1
>>> df2.iloc[0, 0] = 10
>>> df1.iloc[0, 0]
10
```

In summary, modifications are only "propagated" between _identical_ objects (not just
equal (`==`), but identical (`is`) in python terms, see
[docs](https://docs.python.org/3/reference/expressions.html#is)). Propagation is not
really the proper term, since there is only one object that was modified.

However, when in some way creating a new object (even though it might be a DataFrame
with the same data, and thus be an "equal" DataFrame):

```python
>>> df1 = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
>>> df2 = df1[:]  # or df1.loc[...] with some indexer
```

Those objects are no longer identical:

```python
>>> id(df1) == id(df2)  # or df1 is df2
False
```

And thus modifications to one will not propagate to the other:

```python
>>> df1.iloc[0, 0]
1
>>> df2.iloc[0, 0] = 10
>>> df1.iloc[0, 0]  # not changed
1
```

Currently, any getitem indexing operation returns _new_ objects, and also almost all
DataFrame/Series methods return a _new_ object (except with `inplace=True` in some
cases), and thus follow the above logic of never modifying its parent/child DataFrame or
Series (using the lazy Copy-on-Write mechanism where possible).

## Copy / view behaviour in NumPy versus pandas

NumPy has the concept of "views" (an array that shares data with another array, viewing
the same memory, see e.g.
[this explanation](https://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html) for
more details). Typically you create views as a slice of another array. But other
indexing methods, often called "fancy indexing", do not return views but copies: using a
list of indices or a boolean mask.

Pandas, being built on NumPy, uses those concepts, and also exposes the behaviour
consequences to its users. This basically means that pandas users, to understand the
details of how indexing works, also need to understand those view / fancy indexing
concepts of numpy.

However, because DataFrames are not an array, the copy/view rules still differ from
NumPy's rules with current pandas. Slicing rows generally gives a view (following
NumPy), but slicing columns doesn't always give a view (this could be changed to match
NumPy however, see "Alternatives" 1b below). Fancy indexing rows (e.g. with a list of
(positional) labels) gives a copy, but fancy indexing columns _could_ give a view
(currently this gives a copy as well, but one of the "Alternatives" (1b) is to have this
always return a view).

The proposal in this document is to decouple the pandas user-facing behaviour from those
NumPy concepts. Creating a subset of a DataFrame with a slice or with a mask would
behave in a similar way for the user (both return a new object and behave as a copy of
the original). We still use the concept of views internally in pandas to optimize the
implementation, but this becomes hidden from the user.

## Alternatives

The [original document](https://docs.google.com/document/d/1csGE4qigPR2vzmU2--jwURn3sK5k5eVewinxd8OUPk0/edit) and GitHub issue ([Proposal for future copy / view semantics in indexing operations - #36195](https://github.com/pandas-dev/pandas/issues/36195)) discussed several options for making the copy/view situation more consistent and clear:

1. **Well-Defined copy/view rules:** ensure we have more consistent rules about which
   operations result in a copy and which in a view, and then views result in mutating
   the parent, copies not.
   a. A minimal change would be to officialize the current behaviour. This comes down to
      fixing some bugs and clearly documenting and testing which operations are views,
      and which are copies.
   b. An alternative would be to simplify the set of rules. For example: selecting
      columns is always a view, subsetting rows is always a copy. Or: selecting columns
      is always a view, subsetting rows as a slice is a view otherwise always a copy.

2. **Copy-on-Write**: The setitem operation would check if it's a view on another
   dataframe. If it is, then we would copy our data before mutating. (i.e. this
   proposal)

3. **Error-on-Write**: The setitem operation would check if it's a subset of another
   dataframe (both view of copy). Only rather than copying in case of a view we would
   raise an exception telling the user to either copy the data with
   ``.copy_if_needed()`` (name TBD) or mark the frame as "a mutable view" with
   ``.as_mutable_view()`` (name TBD).

This document basically proposes an extended version of option 2 (Copy-on-Write). Some
arguments in favor of Copy-on-Write compared to the other options:

* Copy-on-Write will improve the copy/view efficiency of _methods_ (e.g. rename,
  (re)set_index, drop columns, etc. See section above). This will result in
  lower memory usage and better performance.

* This proposal can also be seen as a clear "well-defined rule". Using Copy-on-Write
  under the hood is an implementation detail to delay the actual copy until it is
  needed. The rule of "always copy" is the simplest "well-defined rule" we can get.

  Other "well-defined rule" ideas above would always include some specific cases (and
  deviations from the NumPy rules). And even with clear rules a user still needs to know
  the details of those rules to understand that `df['a'][df['b'] < 0] = 0` or
  `df[df['b'] < 0]['a'] = 0` does something differently (switched order of column/row
  indexing: the first mutates df (if selecting a column is a view) and the second
  doesn't). While with the "always copy" rule with Copy-on-Write, neither of those
  examples will work to update `df`.

On the other hand, the proposal in this document does not give the user control over
whether a subset should be a view (when possible) that mutates the parent when being
mutated. The only way to modify the parent dataframe is with a direct indexing operation
on this dataframe itself.

See the GitHub comment with some more detailed argumentation:
[https://github.com/pandas-dev/pandas/issues/36195#issuecomment-786654449](https://github.com/pandas-dev/pandas/issues/36195#issuecomment-786654449)

## Disadvantages

Other than the fact that this proposal would result in a backwards incompatible,
breaking change in behaviour (see next section), there are some other potential
disadvantages:

* Deviation from NumPy: NumPy uses the copy and view concepts, while in this proposal
  views would basically not exist anymore in pandas (for the user, at least; we would
  still use it internally as an implementation detail)
  * But as a counter argument: many pandas users are probably not familiar with those
    concepts, and pandas already deviates from the exact rules in NumPy.
* Performance cost of indexing and methods becomes harder to predict: because the copy
  of the data doesn't happen at the moment when actually creating the new object, but
  can happen at a later stage when modifying either the parent or child object, it
  becomes less transparent about when pandas copies data (but in general we should copy
  less often). This is somewhat mitigated because Copy-on-Write will only copy the columns
  that are mutated. Unrelated columns won't get copied.
* Increased memory usage for some use cases: while the majority of use cases will
  see an improvement in memory usage with this proposal, there are a few use
  cases where this might not be the case. Specifically in cases where pandas currently
  does return a view (e.g. slicing rows) and in the case you are fine with (or don't care
  about) the current behaviour of it being a view when mutating that subset (i.e.
  mutating the sliced subset also mutates the parent dataframe), in such a case the
  proposal would introduce a new copy compared to the current behaviour. There is a
  workaround for this though: the copy is not needed if the previous object goes out
  of scope, e.g. the variable is reassigned to something else.

## Backward compatibility

The proposal in this document is clearly a backwards incompatible change that breaks
existing behaviour. Because of the current inconsistencies and subtleties around views
vs. copies and mutation, it would be difficult to change anything without breaking
changes. The current proposal is not the proposal with the minimal changes, though. A
change like this will in any case need to be accompanied with a major version bump (for
example pandas 3.0).

Doing a traditional deprecation cycle that lives in several minor feature releases will
be too noisy. Indexing is too common an operation to include a warning (even if we limit
it to just those operations that previously returned views). However, this proposal is
already implemented and thus available. Users can opt-in and test their code (this is
possible starting with version 1.5 with `pd.options.mode.copy_on_write = True`).

Further we will add a warning mode for pandas 2.2 that raises warnings for all cases that
will change behaviour under the Copy-on-Write proposal. We can
provide a clearly documented upgrade path to first enable the warnings, fix all
warnings, and then enable the Copy-on-Write mode and ensure your code is still working,
and then finally upgrade to the new major release.

## Implementation

The implementation is available since pandas 1.5 (and significantly improved starting
with pandas 2.0). It uses weakrefs to keep track of whether the
data of a Dataframe/Series are viewing the data of another (pandas) object or are being
viewed by another object. This way, whenever the series/dataframe gets modified, we can
check if its data first needs to be copied before mutating it
(see [here](https://pandas.pydata.org/docs/development/copy_on_write.html)).

To test the implementation and experiment with the new behaviour, you can
enable it with the following option:

```python
>>> pd.options.mode.copy_on_write = True
```

after importing pandas (or setting the `PANDAS_COPY_ON_WRITE=1` environment variable
before importing pandas).

## Concrete examples

### Chained assignment

Consider a "classic" case of chained indexing, which was the original motivation for the SettingWithCopy warning:

```python
>>> df[df['B'] > 3]['B'] = 10
```

That is roughly equivalent to

```python
>>> df2 = df[df['B'] > 3]  # Copy under NumPy's rules
>>> df2['B'] = 10  # Update (the copy) df2, df not changed
>>> del df2  # All references to df2 are lost, goes out of scope
```

And so `df` is not modified. For this reason, the SettingWithCopyWarning was introduced.

_With this proposal_, any result of an indexing operation behaves as a copy
(Copy-on-Write), and thus chained assignment will _never_ work. Given that there is then
no ambiguity, the idea is to drop the warning.

The above example is a case where chained assignment doesn't work with current pandas.
But there are of course also patterns with chained assignment that currently _do_ work
and are used. _With this proposal_, any chained assignment will not work, and so those
cases will stop working (e.g. the case above but switching the order):

```python
>>> df['B'][df['B'] > 3] = 10
# or
>>> df['B'][0:5] = 10
```

These cases will raise a warning ``ChainedAssignmentError``, because they can never
accomplish what the user intended. There will be false-positive cases when these
operations are triggered from Cython, because Cython uses a different reference counting
mechanism. These cases should be rare, since calling pandas code from Cython does not
have any performance benefits.

### Filtered dataframe

A typical example where the current SettingWithCopyWarning becomes annoying is when
filtering a DataFrame (which always already returns a copy):

```python
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
>>> df_filtered = df[df["A"] > 1]
>>> df_filtered["new_column"] = 1
SettingWithCopyWarning:
A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead
```

If you then modify your filtered dataframe (e.g. adding a column), you get the
unnecessary SettingWithCopyWarning (with confusing message). The only way to get rid of
the warning is by doing a defensive copy (`df_filtered = df[df["A"] > 1].copy()`, which
results in copying the data twice in the current implementation, Copy-on-Write would
not require ``.copy()`` anymore).

_With this proposal_, the filtered dataframe is never a view and the above
workflow would work as expected without warning (and thus without needing the extra
copy).

### Modifying a Series (from DataFrame column)

_Currently_, accessing a column of a DataFrame as a Series is one of the few cases that
is actually guaranteed to always be a view:

```python
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
>>> s = df["A"]
>>> s.loc[0] = 0   # will also modify df (but no longer with this proposal)
```

_With this proposal_, any indexing operation results in a copy, so also accessing a
column as a Series (in practice, it will still be a view of course, but behave as a copy
through Copy-on-Write). In the above example, mutating `s` will no longer modify the
parent `df`.

This situation is similar as the "chained assignment" case above, except with
an explicit intermediate variable. To actually change the original DataFrame,
the solution is the same: mutate directly the DataFrame in a single step.
For example:

```python
>>> df.loc[0, "A"] = 0
```

### "Shallow" copies

_Currently_, it is possible to create a "shallow" copy of a DataFrame with
`copy(deep=False)`. This creates a new DataFrame object but without copying the
underlying index and data. Any changes to the data of the original will be reflected in
the shallow copy (and vice versa). See the
[docs](https://pandas.pydata.org/pandas-docs/version/1.5/reference/api/pandas.DataFrame.copy.html).

```python
>>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
>>> df2 = df.copy(deep=False)
>>> df2.iloc[0, 0] = 0   # will also modify df (but no longer with this proposal)
```

_With this proposal_, this kind of shallow copy is no longer possible. Only "identical"
objects (in Python terms: `df2 is df`) can share data without triggering Copy-on-Write.
A shallow copy will rather become a "delayed" copy through Copy-on-Write.

See
[#36195 (comment)](https://github.com/pandas-dev/pandas/issues/36195#issuecomment-830579242)
for a more detailed comment on this.

### Methods returning a new DataFrame with the same data

This example is already shown above as well, but so _currently_ almost all methods on a
Series/DataFrame by default return a new object that is a copy of the original data:

```python
>>> df2 = df.rename(columns=str.lower)
>>> df3 = df2.set_index("a")
```

In the above example, df2 holds a copy of the data of df, and df3 holds a copy of the
data of df2. Mutating any of those DataFrames would not modify the parent dataframe.

_With this proposal_, those methods would continue to return new objects, but would use
the shallow copy mechanism with Copy-on-Write so that in practice, those methods don't
need to copy the data at each step, while preserving the current behaviour.

### Series and DataFrame constructors

_Currently_, the Series and DataFrame constructors don't always copy the input
(depending on the type of the input). For example:

```python
>>> s = pd.Series([1, 2, 3])
>>> s2 = pd.Series(s)
>>> s2.iloc[0] = 0   # will also modify the parent Series s
>>> s
0	0  # <-- modified
1	2
2	3
dtype: int64
```

_With this proposal_, we can also use the shallow copy with Copy-on-Write approach _by
default_ in the constructors. This would mean that by default, a new Series or DataFrame
(like `s2` in the above example) would not modify the data from which it is being
constructed (when being modified itself), honoring the proposed rules.

## More background: Current behaviour of views vs copy

To the best of our knowledge, indexing operations currently return views in the
following cases:

* Selecting a single column (as a Series) out of a DataFrame is always a view
  (``df['a']``)
* Slicing columns from a DataFrame creating a subset DataFrame (``df[['a':'b']]`` or
  ``df.loc[:, 'a': 'b']``) is a view _if_ the the original DataFrame consists of a
  single block (single dtype, consolidated) and _if_ you are slicing (so not a list
  selection). In all other cases, getting a subset is always a copy.
* Selecting rows _can_ return a view, when the row indexer is a `slice` object.

Remaining operations (subsetting rows with a list indexer or boolean mask) in practice
return a copy, and we will raise a SettingWithCopyWarning when the user tries to modify
the subset.

## More background: Previous attempts

We've discussed this general issue before. [https://github.com/pandas-dev/pandas/issues/10954](https://github.com/pandas-dev/pandas/issues/10954) and a few pull requests ([https://github.com/pandas-dev/pandas/pull/12036](https://github.com/pandas-dev/pandas/pull/12036), [https://github.com/pandas-dev/pandas/pull/11207](https://github.com/pandas-dev/pandas/pull/11207), [https://github.com/pandas-dev/pandas/pull/11500](https://github.com/pandas-dev/pandas/pull/11500)).

## Comparison with other languages / libraries

### R

For the user, R has somewhat similar behaviour. Most R objects can be considered
immutable, through "copy-on-modify"
([https://adv-r.hadley.nz/names-values.html#copy-on-modify](https://adv-r.hadley.nz/names-values.html#copy-on-modify)).
But in contrast to Python, in R this is a language feature, and any assignment (binding
a variable to a new name) or passing as function argument will essentially create a
"copy" (when mutating such an object, at that point the actual data get copied and
rebind to the name):

```r
x <- c(1, 2, 3)
y <- x
y[[1]] <- 10  # does not modify x
```

While if you would do the above example in Python with a list, x and y are "identical"
and mutating one will also mutate the other.

As a consequence of this language behaviour, modifying a `data.frame` will not modify
other data.frames that might share memory (before being copied with "copy-on-modify").

### Polars

Polars ([https://github.com/pola-rs/polars](https://github.com/pola-rs/polars)) is a
DataFrame library with a Python interface, mainly written in Rust on top of Arrow. It
explicitly
[mentions](https://pola-rs.github.io/polars-book/user-guide/introduction.html#current-status)
"Copy-on-Write" semantics as one its features.

Based on some experiments, the user-facing behaviour of Polars seems similar to the behaviour
described in this proposal (mutating a DataFrame/Series never mutates a parent/child
object, and so chained assignment also doesn't work)


## PDEP-7 History

- July 2021: Initial version
- February 2023: Converted into a PDEP

Note: this proposal has been discussed before it was turned into a PDEP. The main
discussion happened in [GH-36195](https://github.com/pandas-dev/pandas/issues/36195).
This document is modified from the original document discussing different options for
clear copy/view semantics started by Tom Augspurger
([google doc](https://docs.google.com/document/d/1csGE4qigPR2vzmU2--jwURn3sK5k5eVewinxd8OUPk0/edit)).

Related mailing list discussion: [https://mail.python.org/pipermail/pandas-dev/2021-July/001358.html](https://t.co/vT5dOMhNjV?amp=1)
