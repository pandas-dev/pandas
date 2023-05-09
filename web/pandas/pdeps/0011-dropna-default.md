# PDEP-11: dropna default in pandas

- Created: 4 May 2023
- Status: Under discussion
- Discussion: [PR #53094](https://github.com/pandas-dev/pandas/pull/53094)
- Authors: [Richard Shadrach](https://github.com/rhshadrach)
- Revision: 1

## Abstract

Throughout pandas, almost all of the methods that have a `dropna` argument default
to `True`. Being the default, this can cause NA values to be silently dropped.
This PDEP proposes to deprecate the current default value of `True` and change it
to `False` in the next major release of pandas.

## Motivation and Scope

Upon seeing the output for a Series `ser`:

```python
print(ser.value_counts())

1    3
2    1
dtype: Int64
```

users may be surprised that the Series can contain NA values, as is argued in
[#21890](https://github.com/pandas-dev/pandas/issues/21890). By then operating
on data under the assumption NA values are not present, erroroneous results can
arise. The same issue can occur with `groupby`, which can also be used to produce
detailed summary statistics of data. We think it is not unreasonable that an
experienced pandas user seeing the code

    df[["a", "b"]].groupby("a").sum()

would describe this operation as something like the following.

> For each unique value in column `a`, compute the sum of corresponding values
> in column `b` and return the results in a DataFrame indexed by the unique
> values of `a`.

This is correct, except that NA values in the column `a` will be dropped from
the computation. That pandas is taking this additional step in the computation
is not apparent from the code, and can surprise users.

###

### Keeping the default `skipna=True`

Many reductions methods, such as `sum`, `mean`, and `var`, have a `skipna` argument.
In such operations, setting `skipna=False` would make the output of any operation
NA if a single NA value is encountered.

```python
df = pd.DataFrame({'a': [1, np.nan], 'b': [2, np.nan]})
print(df.sum(skipna=False))
# a   NaN
# b   NaN
# dtype: float64
```

This makes `skipna=False` an undesirable default. In the methods with `dropna`, this phenomena does not occur. By defaulting to `dropna=False` in these
methods, the results when NA values are encountered do not obscure the results of non-NA values.

### Possible deprecation of `dropna`

This PDEP takes no position on whether some methods with a `dropna` argument should have said argument deprecated.
However, if such a deprecation is to be pursued, then we believe that the final behavior should
be that of `dropna=False` across any of the methods listed below. With this, a necessary first step
in the deprecation process would be to change the default value to `dropna=False`.

## Detailed Description

The following methods have a dropna argument, those marked with a `*` already default to `False`.

```python
Series.groupby
Series.mode
Series.nunique
Series.to_hdf*
Series.value_counts
DataFrame.groupby
DataFrame.mode
DataFrame.nunique
DataFrame.pivot_table
DataFrame.stack
DataFrame.to_hdf*
DataFrame.value_counts
SeriesGroupBy.nunique
SeriesGroupBy.value_counts
DataFrameGroupBy.nunique
DataFrameGroupBy.value_counts
```

We propose to deprecate the current default of `dropna` and change it to
`False` across all methods listed above.

## Timeline

If accepted, the current `dropna` default would be deprecated as part of pandas
2.x and this deprecation would be enforced in pandas 3.0. In pandas 2.x, `FutureWarning` messages would
be emitted on any calls to these methods where the value of `dropna` is unspecified and
an NA value is present.

## PDEP History

- 4 May 2023: Initial draft
