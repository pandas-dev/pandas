# PDEP-5: NoRowIndex

- Created: 14 November 2022
- Status: Withdrawn
- Discussion: [#49693](https://github.com/pandas-dev/pandas/pull/49693)
- Author: [Marco Gorelli](https://github.com/MarcoGorelli)
- Revision: 2

[TOC]

## Abstract

The suggestion is to add a ``NoRowIndex`` class. Internally, it would act a bit like
a ``RangeIndex``, but some methods would be stricter. This would be one
step towards enabling users who do not want to think about indices to not need to.

## Motivation

The Index can be a source of confusion and frustration for pandas users. For example, let's consider the inputs

```python
In[37]: ser1 = pd.Series([10, 15, 20, 25], index=[1, 2, 3, 5])

In[38]: ser2 = pd.Series([10, 15, 20, 25], index=[1, 2, 3, 4])
```

Then:

- it can be unexpected that adding `Series` with the same length (but different indices) produces `NaN`s in the result (https://stackoverflow.com/q/66094702/4451315):

  ```python
  In [41]: ser1 + ser2
  Out[41]:
  1    20.0
  2    30.0
  3    40.0
  4     NaN
  5     NaN
  dtype: float64
  ```

- concatenation, even with `ignore_index=True`, still aligns on the index (https://github.com/pandas-dev/pandas/issues/25349):

    ```python
    In [42]: pd.concat([ser1, ser2], axis=1, ignore_index=True)
    Out[42]:
          0     1
    1  10.0  10.0
    2  15.0  15.0
    3  20.0  20.0
    5  25.0   NaN
    4   NaN  25.0
    ```

- it can be frustrating to have to repeatedly call `.reset_index()` (https://twitter.com/chowthedog/status/1559946277315641345):

  ```python
  In [3]: ser1.reset_index(drop=True) + ser2.reset_index(drop=True)
  Out[3]:
  0    20
  1    30
  2    40
  3    50
  dtype: int64
  ```

If a user did not want to think about row labels (which they may have ended up after slicing / concatenating operations),
then ``NoRowIndex`` would enable the above to work in a more intuitive
manner (details and examples to follow below).

## Scope

This proposal deals exclusively with the ``NoRowIndex`` class. To allow users to fully "opt-out" of having to think
about row labels, the following could also be useful:
- a ``pd.set_option('mode.no_row_index', True)`` mode which would default to creating new ``DataFrame``s and
  ``Series`` with ``NoRowIndex`` instead of ``RangeIndex``;
- giving ``as_index`` options to methods which currently create an index
  (e.g. ``value_counts``, ``.sum()``, ``.pivot_table``) to just insert a new column instead of creating an
  ``Index``.

However, neither of the above will be discussed here.

## Detailed Description

The core pandas code would change as little as possible. The additional complexity should be handled
within the ``NoRowIndex`` object. It would act just like ``RangeIndex``, but would be a bit stricter
in some cases:
- `name` could only be `None`;
- `start` could only be `0`, `step` `1`;
- when appending one ``NoRowIndex`` to another ``NoRowIndex``, the result would still be ``NoRowIndex``.
  Appending a ``NoRowIndex`` to any other index (or vice-versa) would raise;
- the ``NoRowIndex`` class would be preserved under slicing;
- a ``NoRowIndex`` could only be aligned with another ``Index`` if it's also ``NoRowIndex`` and if it's of the same length;
- ``DataFrame`` columns cannot be `NoRowIndex` (so ``transpose`` would need some adjustments when called on a ``NoRowIndex`` ``DataFrame``);
- `insert` and `delete` should raise. As a consequence, if ``df`` is a ``DataFrame`` with a
  ``NoRowIndex``, then `df.drop` with `axis=0` would always raise;
- arithmetic operations (e.g. `NoRowIndex(3) + 2`) would always raise;
- when printing a ``DataFrame``/``Series`` with a ``NoRowIndex``, then the row labels would not be printed;
- a ``MultiIndex`` could not be created with a ``NoRowIndex`` as one of its levels.

Let's go into more detail for some of these. In the examples that follow, the ``NoRowIndex`` will be passed explicitly,
but this is not how users would be expected to use it (see "Usage and Impact" section for details).

### NoRowIndex.append

If one has two ``DataFrame``s with ``NoRowIndex``, then one would expect that concatenating them would
result in a ``DataFrame`` which still has ``NoRowIndex``. To do this, the following rule could be introduced:

> If appending a ``NoRowIndex`` of length ``y`` to a ``NoRowIndex`` of length ``x``, the result will be a
  ``NoRowIndex`` of length ``x + y``.

Example:

```python
In [6]: df1 = pd.DataFrame({'a': [1,  2], 'b': [4, 5]}, index=NoRowIndex(2))

In [7]: df2 = pd.DataFrame({'a': [4], 'b': [0]}, index=NoRowIndex(1))

In [8]: df1
Out[8]:
 a  b
 1  4
 2  5

In [9]: df2
Out[9]:
 a  b
 4  0

In [10]: pd.concat([df1, df2])
Out[10]:
 a  b
 1  4
 2  5
 4  0

In [11]: pd.concat([df1, df2]).index
Out[11]: NoRowIndex(len=3)
```

Appending anything other than another ``NoRowIndex`` would raise.

### Slicing a ``NoRowIndex``

If one has a ``DataFrame`` with ``NoRowIndex``, then one would expect that a slice of it would still have
a ``NoRowIndex``. This could be accomplished with:

> If a slice of length ``x`` is taken from a ``NoRowIndex`` of length ``y``, then one gets a
 ``NoRowIndex`` of length ``x``. Label-based slicing would not be allowed.

Example:

```python
In [12]: df = pd.DataFrame({'a': [1,  2, 3], 'b': [4, 5, 6]}, index=NoRowIndex(3))

In [13]: df.loc[df['a']>1, 'b']
Out[13]:
5
6
Name: b, dtype: int64

In [14]: df.loc[df['a']>1, 'b'].index
Out[14]: NoRowIndex(len=2)
```

Slicing by label, however, would be disallowed:
```python
In [15]: df.loc[0, 'b']
---------------------------------------------------------------------------
IndexError: Cannot use label-based indexing on NoRowIndex!
```

Note too that:
- other uses of ``.loc``, such as boolean masks, would still be allowed (see F.A.Q);
- ``.iloc`` and ``.iat`` would keep working as before;
- ``.at`` would raise.

### Aligning ``NoRowIndex``s

To minimise surprises, the rule would be:

> A ``NoRowIndex`` can only be aligned with another ``NoRowIndex`` of the same length.
> Attempting to align it with anything else would raise.

Example:
```python
In [1]: ser1 = pd.Series([1, 2, 3], index=NoRowIndex(3))

In [2]: ser2 = pd.Series([4, 5, 6], index=NoRowIndex(3))

In [3]: ser1 + ser2  # works!
Out[3]:
5
7
9
dtype: int64

In [4]: ser1 + ser2.iloc[1:]  # errors!
---------------------------------------------------------------------------
TypeError: Cannot join NoRowIndex of different lengths
```

### Columns cannot be NoRowIndex

This proposal deals exclusively with allowing users to not need to think about
row labels. There's no suggestion to remove the column labels.

In particular, calling ``transpose`` on a ``NoRowIndex`` ``DataFrame``
would error. The error would come with a helpful error message, informing
users that they should first set an index. E.g.:
```python
In [4]: df = pd.DataFrame({'a': [1,  2, 3], 'b': [4, 5, 6]}, index=NoRowIndex(3))

In [5]: df.transpose()
---------------------------------------------------------------------------
ValueError: Columns cannot be NoRowIndex.
If you got here via `transpose` or an `axis=1` operation, then you should first set an index, e.g.: `df.pipe(lambda _df: _df.set_axis(pd.RangeIndex(len(_df))))`
```

### DataFrameFormatter and SeriesFormatter changes

When printing an object with a ``NoRowIndex``, then the row labels would not be shown:

```python
In [15]: df = pd.DataFrame({'a': [1,  2, 3], 'b': [4, 5, 6]}, index=NoRowIndex(3))

In [16]: df
Out[16]:
 a  b
 1  4
 2  5
 3  6
```

Of the above changes, this may be the only one that would need implementing within
``DataFrameFormatter`` / ``SerieFormatter``, as opposed to within ``NoRowIndex``.

## Usage and Impact

Users would not be expected to work with the ``NoRowIndex`` class itself directly.
Usage would probably involve a mode which would change how the ``default_index``
function to return a ``NoRowIndex`` rather than a ``RangeIndex``.
Then, if a ``mode.no_row_index`` option was introduced and a user opted in to it with

```python
pd.set_option("mode.no_row_index", True)
```

then the following would all create a ``DataFrame`` with a ``NoRowIndex`` (as they
all call ``default_index``):

- ``df.reset_index()``;
- ``pd.concat([df1, df2], ignore_index=True)``
- ``df1.merge(df2, on=col)``;
- ``df = pd.DataFrame({'col_1': [1, 2, 3]})``

Further discussion of such a mode is out-of-scope for this proposal. A ``NoRowIndex`` would
just be a first step towards getting there.

## Implementation

Draft pull request showing proof of concept: https://github.com/pandas-dev/pandas/pull/49693.

Note that implementation details could well change even if this PDEP were
accepted. For example, ``NoRowIndex`` would not necessarily need to subclass
``RangeIndex``, and it would not necessarily need to be accessible to the user
(``df.index`` could well return ``None``)

## Likely FAQ

**Q: Could not users just use ``RangeIndex``? Why do we need a new class?**

**A**: ``RangeIndex`` is not preserved under slicing and appending, e.g.:
  ```python
  In[1]: ser = pd.Series([1, 2, 3])

  In[2]: ser[ser != 2].index
  Out[2]: Int64Index([0, 2], dtype="int64")
  ```
  If someone does not want to think about row labels and starts off
  with a ``RangeIndex``, they'll very quickly lose it.

**Q: Are indices not really powerful?**

**A:** Yes! And they're also confusing to many users, even experienced developers.
  Often users are using ``.reset_index`` to avoid issues with indices and alignment.
  Such users would benefit from being able to not think about indices
  and alignment. Indices would be here to stay, and ``NoRowIndex`` would not be the
  default.

**Q: How could one switch a ``NoRowIndex`` ``DataFrame`` back to one with an index?**

**A:** The simplest way would probably be:
  ```python
  df.set_axis(pd.RangeIndex(len(df)))
  ```
  There's probably no need to introduce a new method for this.

  Conversely, to get rid of the index, then if the ``mode.no_row_index`` option was introduced, then
  one could simply do ``df.reset_index(drop=True)``.

**Q: How would ``tz_localize`` and other methods which operate on the index work on a ``NoRowIndex`` ``DataFrame``?**

**A:** Same way they work on other ``NumericIndex``s, which would typically be to raise:

  ```python
  In [2]: ser.tz_localize('UTC')
  ---------------------------------------------------------------------------
  TypeError: index is not a valid DatetimeIndex or PeriodIndex
  ```

**Q: Why not let transpose switch ``NoRowIndex`` to ``RangeIndex`` under the hood before swapping index and columns?**

**A:** This is the kind of magic that can lead to surprising behaviour that's
  difficult to debug. For example, ``df.transpose().transpose()`` would not
  round-trip. It's easy enough to set an index after all, better to "force" users
  to be intentional about what they want and end up with fewer surprises later
  on.

**Q: What would df.sum(), and other methods which introduce an index, return?**

**A:** Such methods would still set an index and would work the same way they
  do now. There may be some way to change that (e.g. introducing ``as_index``
  arguments and introducing a mode to set its default) but that's out of scope
  for this particular PDEP.

**Q: How would a user opt-in to a ``NoRowIndex`` DataFrame?**

**A:** This PDEP would only allow it via the constructor, passing
  ``index=NoRowIndex(len(df))``. A mode could be introduced to toggle
  making that the default, but would be out-of-scope for the current PDEP.

**Q: Would ``.loc`` stop working?**

**A:** No. It would only raise if used for label-based selection. Other uses
  of ``.loc``, such as ``df.loc[:, col_1]`` or ``df.loc[boolean_mask, col_1]``, would
  continue working.

**Q: What's unintuitive about ``Series`` aligning indices when summing?**

**A:** Not sure, but I once asked a group of experienced developers what the
  output of
  ```python
  ser1 = pd.Series([1, 1, 1], index=[1, 2, 3])
  ser2 = pd.Series([1, 1, 1], index=[3, 4, 5])
  print(ser1 + ser2)
  ```
  would be, and _nobody_ got it right.

## Reasons for withdrawal

After some discussions, it has become clear there is not enough for support for the proposal in its current state.
In short, it would add too much complexity to justify the potential benefits. It would unacceptably increase
the maintenance burden, the testing requirements, and the benefits would be minimal.

Concretely:
- maintenance burden: it would not be possible to handle all the complexity within the ``NoRowIndex`` class itself, some
  extra logic would need to go into the pandas core codebase, which is already very complex and hard to maintain;
- the testing burden would be too high. Properly testing this would mean almost doubling the size of the test suite.
  Coverage for options already is not great: for example [this issue](https://github.com/pandas-dev/pandas/issues/49732)
  was caused by a PR which passed CI, but CI did not (and still does not) cover that option (plotting backends);
- it will not benefit most users, as users do not tend to use nor discover options which are not the default;
- it would be difficult to reconcile with some existing behaviours: for example, ``df.sum()`` returns a Series with the
  column names in the index.

In order to make no-index the pandas default and have a chance of benefiting users, a more comprehensive set of changes
would need to made at the same time. This would require a proposal much larger in scope, and would be a much more radical change.
It may be that this proposal will be revisited in the future, but in its current state (as an option) it cannot be accepted.

This has still been a useful exercise, though, as it has resulted in two related proposals (see below).

## Related proposals

- Deprecate automatic alignment, at least in some cases: https://github.com/pandas-dev/pandas/issues/49939;
- ``.value_counts`` behaviour change: https://github.com/pandas-dev/pandas/issues/49497

## PDEP History

- 14 November 2022: Initial draft
- 18 November 2022: First revision (limited the proposal to a new class, leaving a ``mode`` to a separate proposal)
- 14 December 2022: Withdrawal (difficulty reconciling with some existing methods, lack of strong support,
  maintenance burden increasing unjustifiably)
