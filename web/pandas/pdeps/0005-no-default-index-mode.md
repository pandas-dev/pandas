# PDEP-5: No-default-index mode

- Created: 14 November 2022
- Status: Draft
- Discussion: [#49693](https://github.com/pandas-dev/pandas/pull/49693)
- Author: [Marco Gorelli](https://github.com/MarcoGorelli)
- Revision: 1

## Abstract

The suggestion is to add a `mode.no_default_index` option which, if enabled,
would ensure:
- if a ``DataFrame`` / ``Series`` is created, then by default it won't have an ``Index``;
- nobody will get an ``Index`` unless they ask for one - this would affect the default behaviour of ``groupby``, ``value_counts``, ``pivot_table``, and more.

This option would not be the default. Users would need to explicitly opt-in to it, via ``pd.set_option('mode.no_default_index', True)``, via ``pd.option_context``, or via the ``PANDAS_NO_DEFAULT_INDEX`` environment variable.

## Motivation and Scope

The Index can be a source of confusion and frustration for pandas users. For example, let's consider the inputs

```python
In [37]: ser1 = df.groupby('sender')['amount'].sum()

In [38]: ser2 = df.groupby('receiver')['amount'].sum()

In [39]: ser1
Out[39]:
sender
1    10
2    15
3    20
5    25
Name: amount, dtype: int64

In [40]: ser2
Out[40]:
receiver
1    10
2    15
3    20
4    25
Name: amount, dtype: int64
```
. Then:

- it can be unexpected that summing `Series` with the same length (but different indices) produces `NaN`s in the result (https://stackoverflow.com/q/66094702/4451315):

  ```python
  In [41]: ser1 + ser2
  Out[41]:
  1    20.0
  2    30.0
  3    40.0
  4     NaN
  5     NaN
  Name: amount, dtype: float64
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
    In [45]: df.value_counts(['sender', 'receiver']).reset_index().rename(columns={0: 'count'})
    Out[45]:
       sender  receiver  count
    0       1         1      1
    1       2         2      1
    2       3         3      1
    3       5         4      1
    ```

With this option enabled, users who don't want to worry about indices wouldn't need to.

## Detailed Description

This would require 3 steps:
1. creation of a ``NoIndex`` object, which would be a subclass of ``RangeIndex`` on which
  some operations such as ``append`` would behave differently.
  The ``default_index`` function would then return ``NoIndex`` (rather than ``RangeIndex``) if this mode is enabled;
2. adjusting ``DataFrameFormatter`` and ``SeriesFormatter`` to not print row labels for objects with a ``NoIndex``;
3. adjusting methods which currently return an index to just insert a new column instead.

Let's expand on all three below.

### 1. NoIndex object

Most of the logic could be handled within the ``NoIndex`` object.
It would be like a ``RangeIndex``, but with the following differences:
- `name` could only be `None`;
- `start` could only be `0`, `step` `1`;
- when appending an extra element, the new `Index` would still be `NoIndex`;
- when slicing, one would still get a `NoIndex`;
- two ``NoIndex`` objects can't be aligned. Either they're the same length, or pandas raises;
- aligning a ``NoIndex`` object with one which has an index will raise, always;
- ``DataFrame`` columns can't be `NoIndex` (so ``transpose`` would need some adjustments when called on a ``NoIndex`` ``DataFrame``);
- `insert` and `delete` should raise. As a consequence, `.drop` with `axis=0` would always raise;
- arithmetic operations (e.g. `NoIndex(3) + 2`) would all raise.

### 2. DataFrameFormatter and SeriesFormatter changes

When printing an object with a ``NoIndex``, then the row labels wouldn't be shown:

```python
In [14]: pd.set_option('mode.no_default_index', True)

In [15]: df = pd.DataFrame({'a': [1,  2, 3], 'b': [4, 5, 6], 'c': [7, 8, 9]})

In [16]: df
Out[16]:
 a  b  c
 1  4  7
 2  5  8
 3  6  9
```

### 3. Nobody should get an index unless they ask for one

The following would work in the same way:
```python
pivot = (
    pd.pivot_table(df, values="D", index=["A", "B"], columns=["C"], aggfunc=np.sum)
).reset_index()

with pd.option_context('mode.no_default_index', True):
    pivot = (
        pd.pivot_table(df, values="D", index=["A", "B"], columns=["C"], aggfunc=np.sum)
    )
```

Likewise for ``value_counts``. In ``groupby``, the default would be ``as_index=False``.

## Usage and Impact

Users who like the power of the ``Index`` could continue using pandas exactly as it is,
without changing anything.

The addition of this mode would enable users who don't want to think about indices to
not have to.

The implementation would be quite simple: most of the logic would be handled within the
``NoIndex`` class, and only some minor adjustments (e.g. to the ``default_index`` function)
would be needed in core pandas.

## Implementation

Draft pull request showing proof of concept: https://github.com/pandas-dev/pandas/pull/49693.

## Likely FAQ

**Q: Aren't indices really powerful?**

**A:** Yes! And they're also confusing to many users, even experienced developers.
  It's fairly common to see pandas code with ``.reset_index`` scattered around every
  other line. Such users would benefit from a mode in which they wouldn't need to think
  about indices and alignment.

**Q: In this mode, could users still get an ``Index`` if they really wanted to?**

**A:** Yes! For example with
  ```python
  df.set_index(Index(range(len(df))))
  ```
  or, if they don't have a column named ``'index'``:
  ```python
  df.reset_index().set_index('index')
  ```

**Q: Why is it necessary to change the behaviour of ``value_counts``? Isn't the introduction of a ``NoIndex`` object enough?**

**A:** The objective of this mode is to enable users to not have to think about indices if they don't want to. If they have to call
  ``.reset_index`` after each ``value_counts`` / ``pivot_table`` call, or remember to pass ``as_index=False`` to each ``groupby``
  call, then this objective has arguably not quite been reached.

## PDEP History

- 14 November: Initial draft
