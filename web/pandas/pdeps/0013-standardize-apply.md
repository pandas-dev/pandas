# PDEP-13: Make the Series.apply method operate Series-wise

- Created: 24 August 2023
- Status: Under discussion
- Discussion: [#52140](https://github.com/pandas-dev/pandas/issues/52509)
- Author: [Terji Petersen](https://github.com/topper-123)
- Revision: 1

## Abstract

Currently, giving an input to `Series.apply` is treated differently depending on the type of the input:

* if the input is a numpy `ufunc`, `series.apply(func)` is equivalent to `func(series)`, i.e. similar to `series.pipe(func)`.
* if the input is a callable, but not a numpy `ufunc`, `series.apply(func)` is similar to `Series([func(val) for val in series], index=series.index)`, i.e. similar to `series.map(func)`
* if the input is a list-like or dict-like, `series.apply(func)` is equivalent to `series.agg(func)` (which is subtly different than `series.apply`)

In contrast, `DataFrame.apply` has a consistent behavior:

* if the input is a callable, `df.apply(func)` always calls each columns in the DataFrame, so is similar to `func(col) for _, col in
df.items()` + wrapping functionality
* if the input is a list-like or dict-like, `df.apply` call each item in the list/dict and wraps the result as needed. So for example if the input is a list, `df.apply(func_list)` is equivalent to `[df.apply(func) for func in func_list]` + wrapping functionality

This PDEP proposes that:

- The current complex current behavior of `Series.apply` will be deprecated in Pandas 2.2.
- Single callables given to the `.apply` methods of `Series` will in Pandas 3.0 always be called on the whole `Series`, so `series.apply(func)` will become similar to `func(series)`,
- Lists or dicts of callables given to the `Series.apply` will in Pandas 3.0 always call `Series.apply` on each element of the list/dict

In short, this PDEP proposes changing `Series.apply` to be more similar to how `DataFrame.apply` works on single dataframe columns, i.e. operate on the whole series. If a user wants to map a callable to each element of a Series, they should be directed to use `Series.map` instead of using `Series.apply`.

## Motivation

`Series.apply` is currently a very complex method, whose behaviour will differ depending on the nature of its input.

`Series.apply` & `Series.map` currently often behave very similar, but differently enough for it to be confusing when it's a good idea to use one over the other and especially when `Series.apply` is a bad idea to use.

Also, calling `Series.apply` currently gives a different result than the per-column result from calling `DataFrame.apply`, which can be confusing for users who expect `Series.apply` to be the `Series` version of `DataFrame.apply`, similar to how `Series.agg` is the `Series` version of `DataFrame.agg`. For example, currently some functions may work fine with `DataFrame.apply`, but may fail, be very slow when given to `Series.apply` or give a different result than the per-column result from `DataFrame.apply`.

### Similarities and differences between `Series.apply` and `Series.map`

The similarity between the methods is especially that they both fall back to use `Series._map_values` and there use `algorithms.map_array` or `ExtensionArray.map` as relevant.

The differences are many, but each one is relative minor:

1. `Series.map` has a `na_action` parameter, which `Series.apply` doesn't
2. `Series.apply` can take advantage of numpy ufuncs, which `Series.map` can't
3. `Series.apply` can take `args` and `**kwargs`, which `Series.map` can't
4. `Series.apply` is more general and can take a string, e.g. `"sum"`, or lists or dicts of inputs which `Series.map` can't.
5. when given a numpy ufunc, the ufunc will be called on the whole Series, when given to `Series.apply` and on each element of the series, if given to `Series.map`.

In addition, `Series.apply` has some functionality, which `Series.map` does not, but which has already been deprecated:

6. `Series.apply` has a `convert_dtype` parameter, which has been deprecated (deprecated in pandas 2.1, see [GH52257](https://github.com/pandas-dev/pandas/pull/52257))
7. `Series.apply` will return a Dataframe, if its result is a list of Series (deprecated in pandas 2.1, see [GH52123]()https://github.com/pandas-dev/pandas/pull/52123)).

### Similarities and differences between `Series.apply` and `DataFrame.apply`

`Series.apply` and `DataFrame.apply` are similar when given numpy ufuncs as inputs, but when given non-ufuncs as inputs, `Series.apply` and `DataFrame.apply` will behave differently, because `series.apply(func)` will be similar to `series.map(func)` while `Dataframe.apply(func)` will call the input on each column series and combine the result.

If given a list-like or dict-like, `Series.apply` will behave similar to `Series.agg`, while `DataFrame.apply` will call each element in the list-like/dict-like on each column and combine the results.

Also `DataFrame.apply` has some parameters (`raw` and `result_type`) which are relevant for a 2D DataFrame, but may not be relevant for `Series.apply`, because `Series` is a 1D structure.

## Examples of problems with the current way `Series.apply` works

The above similarities and many minor differences makes for confusing and too complex rules for when its a good idea to use `Series.apply` over `Series.map` to do operations, and vica versa, and for when a callable will work well with `Series.apply` versus `DataFrame.apply`. Some examples will show some examples below.

First some setup:

```python
>>> import numpy as np
>>> import pandas as pd
>>>
>>> small_ser = pd.Series([1, 2, 3])
>>> large_ser = pd.Series(range(100_000))
```

### 1: string vs numpy funcs in `Series.apply`

```python
>>> small_ser.apply("sum")
6
>>> small_ser.apply(np.sum)
0    1
1    2
2    3
dtype: int64
```

It will surprise users that these two give different results. Also, anyone using the second pattern is probably making a mistake.

Note that giving `np.sum` to `DataFrame.apply` aggregates properly:

```python
>>> pd.DataFrame(small_ser).apply(np.sum)
0    6
dtype: int64
```

This PDEP proposes that callables will be applies to the whole `Series`, so we in the future will have:

```python
>>> small_ser.apply(np.sum)
6
```

### 2 Callables vs. list/dict of callables

Giving functions and lists/dicts of functions will give different results:

```python
>>> small_ser.apply(np.sum)
0    1
1    2
2    3
dtype: int64
>>> small_ser.apply([np.sum])
sum    6
dtype: int64
```

Also with non-numpy callables:

```python
>>> small_ser.apply(lambda x: x.sum())
AttributeError: 'int' object has no attribute 'sum'
>>> small_ser.apply([lambda x: x.sum()])
<lambda>    6
dtype: int64
```

In both cases above the difference is that `Series.apply` operates element-wise, if given a callable, but series-wise if given a list/dict of callables.

This PDEP proposes that callables will be applies to the whole `Series`, so we in the future will have:

```python
>>> small_ser.apply(lambda x: x.sum())
6
>>> small_ser.apply([lambda x: x.sum()])
<lambda>    6
dtype: int64
```

### 3. Functions in `Series.apply`

The `Series.apply` doc string have examples with using lambdas, but using lambdas in `Series.apply` is often a bad practices because of bad performance:

```python
>>> %timeit large_ser.apply(lambda x: x + 1)
24.1 ms ± 88.8 µs per loop
```

Currently, `Series` does not have a method that makes a callable operate on a series' data. Instead users need to use `Series.pipe` for that operation in order for the operation to be efficient:

```python
>>> %timeit large_ser.pipe(lambda x: x + 1)
44 µs ± 363 ns per loop
```

(The reason for the above performance differences is that apply gets called on each single element, while `pipe` calls `x.__add__(1)`, which operates on the whole array).

Note also that `.pipe` operates on the `Series` while `apply`currently operates on each element in the data, so there is some differences that may have some consequence in some cases.

This PDEP proposes that callables will be applies to the whole `Series`, so we in the future `Series.apply` will be as fast as `Series.pipe`.

### 4. ufuncs in `Series.apply` vs. normal functions

Performance-wise, ufuncs are fine in `Series.apply`, but non-ufunc functions are not:

```python
>>> %timeit large_ser.apply(np.sqrt)
71.6 µs ± 1.17 µs per loop
>>> %timeit large_ser.apply(lambda x:np.sqrt(x))
63.6 ms ± 1.03 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

It is difficult to understand why ufuncs are fast in `apply`, while other callables are slow in `apply` (answer: it's because ufuncs operate on the whole Series, while other callables operate elementwise).

This PDEP proposes that callables will be applies to the whole `Series`, so we in the future non-ufunc functions in `Series.apply` will be as fast as ufuncs.

### 5. callables in `Series.apply` is slow, callables in `DataFrame.apply` is fast

Above it was shown that using (non-ufunc) callables in `Series.apply` is bad performance-wise. OTOH using them in `DataFrame.apply` is fine:

```python
>>> %timeit large_ser.apply(lambda x: x + 1)
24.3 ms ± 24 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
>>> %timeit pd.DataFrame(large_ser).apply(lambda x: x + 1)
160 µs +- 1.17 µs per loop (mean ± std. dev. of 7 runs, 10,000 loops each)
```

Having callables being fast to use in the `DataFrame.apply` method, but slow in `Series.apply` is confusing for users.

This PDEP proposes that callables will be applies to the whole `Series`, so we in the future `Series.apply` will be as fast as `DataFrame.apply` already is.

### 6. callables in `Series.apply` may fail, while callables in `DataFrame.apply` do not and vica versa

```python
>>> ser.apply(lambda x: x.sum())
AttributeError: 'int' object has no attribute 'sum'
>>> pd.DataFrame(ser).apply(lambda x: x.sum())
0    6
dtype: int64
```

Having callables fail when used in `Series.apply`, but work in `DataFrame.Apply` or vica versa is confusing for users.

This PDEP proposes that callables will be applied to the whole `Series`, so callables given to `Series.apply` will work the same as when given to `DataFrame.apply`, so in the future we will have that:

```python
>>> ser.apply(lambda x: x.sum())
6
>>> pd.DataFrame(ser).apply(lambda x: x.sum())
0    6
dtype: int64
```

### 7.  `Series.apply` vs. `Series.agg`

The doc string for `Series.agg` says about the method's `func` parameter: "If a function, must ... work when passed ... to Series.apply". But compare these:

```python
>>> small_ser.agg(np.sum)
6
>>> small_ser.apply(np.sum)
0    1
1    2
2    3
dtype: int64
```

Users would expect these two to give the same result.

This PDEP proposes that callables will be applied to the whole `Series`, so in the future we will have:

```python
>>> small_ser.agg(np.sum)
6
>>> small_ser.apply(np.sum)
6
```

### 8. dictlikes vs. listlikes in `Series.apply`

Giving a *list* of transforming arguments to `Series.apply` returns a `DataFrame`:

```python
>>> small_ser.apply(["sqrt", np.abs])
       sqrt  absolute
0  1.000000         1
1  1.414214         2
2  1.732051         3
```

But giving a *dict* of transforming arguments returns a `Series` with a `MultiIndex`:

```python
>>> small_ser.apply({"sqrt" :"sqrt", "abs" : np.abs})
sqrt  0    1.000000
      1    1.414214
      2    1.732051
abs   0    1.000000
      1    2.000000
      2    3.000000
dtype: float64
```

These two should give same-shaped output for consistency. Using `Series.transform` instead of `Series.apply`, it returns a `DataFrame` in both cases and I think the dictlike example above should return a `DataFrame` similar to the listlike example.

Minor additional info: listlikes and dictlikes of aggregation arguments do behave the same, so this is only a problem with dictlikes of transforming arguments when using `apply`.

This PDEP proposes that the result from giving list-likes and dict-likes to `Series.apply` will have the same shape as when given list-likes currently:

```python
>>> small_ser.apply(["sqrt", np.abs])
       sqrt  absolute
0  1.000000         1
1  1.414214         2
2  1.732051         3
>>> small_ser.apply({"sqrt" :"sqrt", "abs" : np.abs})
       sqrt  absolute
0  1.000000         1
1  1.414214         2
2  1.732051         3
```

## Proposal

With the above in mind, it is proposed that:

1. When given a callable, `Series.apply` always operate on the series. I.e. let `series.apply(func)` be similar to `func(series)` + the needed additional functionality.
2. When given a list-like or dict-like, `Series.apply` will apply each element of the list-like/dict-like to the series. I.e. `series.apply(func_list)` will be similar to `[series.apply(func) for func in func_list]` + the needed additional functionality
3. The changes made to `Series.apply`will propagate to `Series.agg` and `Series.transform` as needed.

The difference between `Series.apply()` & `Series.map()` will then be that:

* `Series.apply()` makes the passed-in callable operate on the series, similarly to how `(DataFrame|SeriesGroupby|DataFrameGroupBy).apply` operate on series. This is very fast and can do almost anything,
* `Series.map()` makes the passed-in callable operate on each series data elements individually. This is very flexible, but can be very slow, so should only be used if `Series.apply` can't do it.

so, this API change will help make Pandas `Series.(apply|map)` API  clearer without losing functionality and let their functionality  be explainable in a simple manner, which would be a win for Pandas.

The result from the above change will be that `Series.apply` will operate similar to how `DataFrame.apply` works already per column, similar to how `Series.map` operates similar to how `DataFrame.map` works per column. This will give better coherence between same-named methods on `DataFrames` and `Series`.

## Deprecation process

To change the behavior to the current behavior will have to be deprecated. This can be done by adding a `by_row` parameter to `Series.apply`, which means, when `by_rows=False`, that `Series.apply` will not operate elementwise but Series-wise.

So we will have in pandas v2.2:

```python
>>> def apply(self, ..., by_row: bool | NoDefault=no_default, ...):
    if by_row is no_default:
        warn("The by_row parameter will be set to False in the future")
        by_row = True
    ...
```

In pandas v3.0 the signature will change to:

```python
>>> def apply(self, ..., by_row: NoDefault=no_default, ...):
    if by_row is not no_default:
        warn("Do not use the by_row parameter, it will be removed in the future")
    ...
```

I.e. the `by_row` parameter will be needed in the signature in v3.0 in order be backward compatible with v2.x, but will have no effect.

In Pandas v4.0, the `by_row` parameter will be removed.

## PDEP-13 History

- 24 august 2023: Initial version
