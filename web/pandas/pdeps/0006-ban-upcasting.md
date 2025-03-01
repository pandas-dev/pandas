# PDEP-6: Ban upcasting in setitem-like operations

- Created: 23 December 2022
- Status: Implemented
- Discussion: [#39584](https://github.com/pandas-dev/pandas/pull/50402)
- Author: [Marco Gorelli](https://github.com/MarcoGorelli) ([original issue](https://github.com/pandas-dev/pandas/issues/39584) by [Joris Van den Bossche](https://github.com/jorisvandenbossche))
- Revision: 1

[TOC]

## Abstract

The suggestion is that setitem-like operations would
not change a ``Series``' dtype (nor that of a ``DataFrame``'s column).

Current behaviour:
```python
In [1]: ser = pd.Series([1, 2, 3], dtype='int64')

In [2]: ser[2] = 'potage'

In [3]: ser  # dtype changed to 'object'!
Out[3]:
0         1
1         2
2    potage
dtype: object
```

Suggested behaviour:

```python
In [1]: ser = pd.Series([1, 2, 3])

In [2]: ser[2] = 'potage'  # raises!
---------------------------------------------------------------------------
ValueError: Invalid value 'potage' for dtype int64
```

## Motivation and Scope

Currently, pandas is extremely flexible in handling different dtypes.
However, this can potentially hide bugs, break user expectations, and copy data
in what looks like it should be an inplace operation.

An example of it hiding a bug is:
```python
In[9]: ser = pd.Series(pd.date_range("2000", periods=3))

In[10]: ser[2] = "2000-01-04"  # works, is converted to datetime64

In[11]: ser[2] = "2000-01-04x"  # typo - but pandas does not error, it upcasts to object
```

The scope of this PDEP is limited to setitem-like operations on Series (and DataFrame columns).
For example, starting with:
```python
df = DataFrame({"a": [1, 2, np.nan], "b": [4, 5, 6]})
ser = df["a"].copy()
```
then the following would all raise:

* setitem-like operations:

    - ``ser.fillna('foo', inplace=True)``
    - ``ser.where(ser.isna(), 'foo', inplace=True)``
    - ``ser.fillna('foo', inplace=False)``
    - ``ser.where(ser.isna(), 'foo', inplace=False)``

* setitem indexing operations (where ``indexer`` could be a slice, a mask,
  a single value, a list or array of values, or any other allowed indexer):

    - ``ser.iloc[indexer] = 'foo'``
    - ``ser.loc[indexer] = 'foo'``
    - ``df.iloc[indexer, 0] = 'foo'``
    - ``df.loc[indexer, 'a'] = 'foo'``
    - ``ser[indexer] = 'foo'``

It may be desirable to expand the top list to ``Series.replace`` and ``Series.update``,
but to keep the scope of the PDEP down, they are excluded for now.

Examples of operations which would not raise are:

- ``ser.diff()``
- ``pd.concat([ser, ser.astype(object)])``
- ``ser.mean()``
- ``ser[0] = 3``  # same dtype
- ``ser[0] = 3.``  # 3.0 is a 'round' float and so compatible with 'int64' dtype
- ``df['a'] = pd.date_range(datetime(2020, 1, 1), periods=3)``
- ``df.index.intersection(ser.index)``

## Detailed description

Concretely, the suggestion is:

- If a ``Series`` is of a given dtype, then a ``setitem``-like operation should not change its dtype.
- If a ``setitem``-like operation would previously have changed a ``Series``' dtype, it would now raise.

For a start, this would involve:

1. changing ``Block.setitem`` such that it does not have an ``except`` block in:

    <!-- language: python -->

        value = extract_array(value, extract_numpy=True)
        try:
            casted = np_can_hold_element(values.dtype, value)
        except LossSetitiemError:
            # current dtype cannot store value, coerce to common dtype
            nb = self.coerce_to_target_dtype(value)
            return nb.setitem(index, value)
        else:

2. making a similar change in:

    - ``Block.where``
    - ``Block.putmask``
    - ``EABackedBlock.setitem``
    - ``EABackedBlock.where``
    - ``EABackedBlock.putmask``

The above would already require several hundreds of tests to be adjusted. Note that once
implementation starts, the list of locations to change may turn out to be slightly
different.

### Ban upcasting altogether, or just upcasting to ``object``?

The trickiest part of this proposal concerns what to do when setting a float in an integer column:

```python
In[1]: ser = pd.Series([1, 2, 3])

In [2]: ser
Out[2]:
0    1
1    2
2    3
dtype: int64

In[3]: ser[0] = 1.5  # what should this do?
```

The current behaviour is to upcast to 'float64':
```python
In [4]: ser
Out[4]:
0    1.5
1    2.0
2    3.0
dtype: float64
```

This is not necessarily a sign of a bug, because the user might just be thinking of their ``Series`` as being
numeric (without much regard for ``int`` vs ``float``) - ``'int64'`` is just what pandas happened to infer
when constructing it.

Possible options could be:

1. Only accept round floats (e.g. ``1.0``) and raise on anything else (e.g. ``1.01``).
2. Convert the float value to ``int`` before setting it (i.e. silently round all float values).
3. Limit "banning upcasting" to when the upcasted dtype is ``object`` (i.e. preserve current behavior of upcasting the int64 Series to float64).

Let us compare with what other libraries do:

- ``numpy``: option 2
- ``cudf``: option 2
- ``polars``: option 2
- ``R data.frame``: just upcasts (like pandas does now for non-nullable dtypes);
- ``pandas`` (nullable dtypes): option 1
- ``datatable``: option 1
- ``DataFrames.jl``: option 1

Option ``2`` would be a breaking behaviour change in pandas. Further,
if the objective of this PDEP is to prevent bugs, then this is also not desirable:
someone might set ``1.5`` and later be surprised to learn that they actually set ``1``.

There are several downsides to option ``3``:

- It would be inconsistent with the nullable dtypes' behaviour.
- It would also add complexity to the codebase and to tests.
- It would be hard to teach, as instead of being able to teach a simple rule,
  There would be a rule with exceptions.
- There would be a risk of loss of precision and or overflow.
- It opens the door to other exceptions, such as not upcasting ``'int8'`` to ``'int16'``.

Option ``1`` is the maximally safe one in terms of protecting users from bugs, being
consistent with the current behaviour of nullable dtypes, and in being simple to teach.
Therefore, the option chosen by this PDEP is option 1.

## Usage and Impact

This would make pandas stricter, so there should not be any risk of introducing bugs. If anything, this would help prevent bugs.

Unfortunately, it would also risk annoying users who might have been intentionally upcasting.

Given that users could still get the current behaviour by first explicitly casting the Series
to float, it would be more beneficial to the community at large to err on the side
of strictness.

## Out of scope

Enlargement. For example:
```python
ser = pd.Series([1, 2, 3])
ser[len(ser)] = 4.5
```
There is arguably a larger conversation to be had about whether that should be allowed
at all. To keep this proposal focused, it is intentionally excluded from the scope.

## F.A.Q.

**Q: What happens if setting ``1.0`` in an ``int8`` Series?**

**A**: The current behavior would be to insert ``1.0`` as ``1`` and keep the dtype
  as ``int8``. So, this would not change.

**Q: What happens if setting ``1_000_000.0`` in an ``int8`` Series?**

**A**: The current behavior would be to upcast to ``int32``. So under this PDEP,
  it would instead raise.

**Q: What happens in setting ``16.000000000000001`` in an ``int8`` Series?**

**A**: As far as Python is concerned, ``16.000000000000001`` and ``16.0`` are the
  same number. So, it would be inserted as ``16`` and the dtype would not change
  (just like what happens now, there would be no change here).

**Q: What if I want ``1.0000000001`` to be inserted as ``1.0`` in an ``int8`` Series?**

**A**: You may want to define your own helper function, such as:

```python
def maybe_convert_to_int(x: int | float, tolerance: float):
    if np.abs(x - round(x)) < tolerance:
        return round(x)
    return x
```

which you could adapt according to your needs.


## Timeline

Deprecate sometime in the 2.x releases (after 2.0.0 has already been released), and enforce in 3.0.0.

### PDEP History

- 23 December 2022: Initial draft
- 4 July 2024: Change status to "implemented"
