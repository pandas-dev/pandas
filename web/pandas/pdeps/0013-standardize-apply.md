# PDEP-13: Deprecate the apply method on Series and DataFrame and make the agg and transform methods operate on series data

- Created: 24 August 2023
- Status: Under discussion
- Discussion: [#52140](https://github.com/pandas-dev/pandas/issues/52509)
- Author: [Terji Petersen](https://github.com/topper-123)
- Revision: 2

## Abstract

The `apply`, `transform` and `agg` methods have very complex behavior when given callables because they in some cases operate on elements in series, in some cases on series and sometimes try one first, and it that fails, falls back to try the other. There is not a logical system how these behaviors are arranged and it can therefore be difficult for users to understand these methods.

It is proposed that `apply`, `transform` and `agg` in the future will work as follows:

1. the `agg` and `transform` methods of `Series`, `DataFrame` and `groupby` will always operate series-wise and never element-wise
2. `Series.apply` and `DataFrame.apply` will be deprecated.
3. The current behavior when supplying string to the methods will not be changed.
4. `groupby.apply` will not be deprecated (because it behaves differently than `Series.apply` and `DataFrame.apply`)

The above changes means that the future behavior, when users want to apply arbitrary callables in pandas, can be described as follows:

1. When users want to operate on single elements in a `Series` or `DataFrame`, they should use `Series.map` and `DataFrame.map` respectively.
2. When users want to aggregate a `Series`, columns/rows of a `DataFrame` or groups in `groupby` objects, they should use `Series.agg`, `DataFrame.agg` and `groupby.agg` respectively.
3. When users want to transform a `Series`, columns/rows of a `DataFrame` or groups in `groupby` objects, they should use `Series.transform`, `DataFrame.transform` and `groupby.transform` respectively.
4. Functions that are not applicable for  `map`, `agg` nor `transform` are considered relatively rare and in the future users should call these functions directly rather than use the `apply` method.

The use of `Series.apply` and  `DataFrame.apply` will after the proposed change in almost all cases be replaced by `map`, `agg` or `transform`. In the very few cases where `Series.apply` and  `DataFrame.apply` cannot be substituted by `map`, `agg` or `transform`, it is proposed that it will be accepted that users will have to find alternative ways to apply the functions, i.e. typically apply the functions manually and possibly concatenating the results.

## Motivation

The current behavior of `apply`, `agg` and `transform` is very complex and therefore difficult to understand for non-expert users. The difficulty is especially that the methods sometimes apply callables on elements of series/dataframes, sometimes on Series or columns/rows of Dataframes and sometimes try element-wise operation and if that fails, falls back to series-wise operations.

Below is an overview of the current behavior in table form when giving callables to `agg`, `transform` and `apply`. As an example on how to read the tables, when a non-ufunc callable is given to `Series.agg`, `Series.agg` will first try to apply the callable to each element in the series, and if that fails, will fall back to call the series using the callable.

(The description may not be 100 % accurate because of various special cases in the current implementation, but will give a good understanding of the current behavior).

### agg

|                                    | Series                           | DataFrame                        | groupby   |
|:-----------------------------------|:---------------------------------|:---------------------------------|:----------|
| ufunc or list/dict of ufuncs       | series                           | series                           | series    |
| other callables (non ufunc)        | Try elements, fallback to series | series                           | series    |
| list/dict of callables (non-ufunc) | Try elements, fallback to series | Try elements, fallback to series | series    |

### transform

|                                    | Series                           | DataFrame                        | groupby   |
|:-----------------------------------|:---------------------------------|:---------------------------------|:----------|
| ufunc or list/dict of ufuncs       | series                           | series                           | series    |
| other callables (non ufunc)        | Try elements, fallback to series | series                           | series    |
| list/dict of callables (non-ufunc) | Try elements, fallback to series | Try elements, fallback to series | series    |

### apply

|                                    | Series   | DataFrame   | groupby   |
|:-----------------------------------|:---------|:------------|:----------|
| ufunc or list/dict of ufuncs       | series   | series      | series    |
| other callables (non ufunc)        | elements | series      | series    |
| list/dict of callables (non-ufunc) | Try elements, fallback to series | series    | series    |

The 3 tables show that:

1. when given numpy ufuncs, callables given to `agg`/`transform`/`apply` operate on series data
2. when used on groupby objects, callables given to `agg`/`transform`/`apply` operate on series data
3. else, in some case it will try element-wise operation and fall back to series-wise operations if that fails, in some case will operate on series data and in some cases on element data.

The above differences result on some non-obvious differences in how the same callable given to `agg`/`transform`/`apply` will behave.

For example, calling `agg` using the same callable will give different results depending on context:

```python
>>> import pandas as pd
>>> df = pd.DataFrame({"A": range(3)})
>>>
>>> df.agg(lambda x: np.sum(x))  # ok
A    3
dtype: int64
>>> df.agg([lambda x: np.sum(x)])  # not ok
         A
  <lambda>
0        0
1        1
2        2
>>> df.A.agg(lambda x: np.sum(x))  # not ok
0    0
1    1
2    2
Name: A, dtype: int64
```

It can also have great effect on performance, even when the result is correct. For example:

```python
>>> df = pd.DataFrame({"A": range(1_000_000)})
>>> %tiemit df.transform(lambda x: x + 1)  # fast
1.43 ms ± 3.6 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)
 >>> %timeit df.transform([lambda x: x + 1])  # slow
163 ms ± 754 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
 >>> %timeit df.A.transform(lambda x: x + 1)  # slow
162 ms ± 980 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

The reason for the great performance difference is that `df.transform(func)` operates on series data, which is fast, while `df.transform(func_list)` will attempt elementwise operation first, and if that works (which is does here), will be much slower than series operations.

In addition to the above effects of the current implementation of `agg`/`transform` and `apply`, see [#52140](https://github.com/pandas-dev/pandas/issues/52140) for more examples of the unexpected effects of how `apply` is implemented.

It can also be noted that `Series.apply` and `DataFrame.apply` could almost always be replaced with calls to `agg`, `transform` or `map`, if `agg` and `transform` were to always operate on series data. For some examples, see the table below for alternatives using `apply(func)`:

| func                    | Series     | DataFrame   |
|:--------------------|:-----------|:------------|
| lambda x: str(x)    | .map       | .map        |
| lambda x: x + 1     | .transform | .transform  |
| [lambda x: x.sum()] | .agg       | .agg        |

So, for example, `ser.apply(lambda x: str(x))` can be replaced with `ser.map(lambda x: str(x))` while `df.apply([lambda x: x.sum()])` can be replaced with `df.agg([lambda x: x.sum()])`.

Overall, because of their flexibility, `Series.apply` and `DataFrame.apply` are considered unnecessarily complex, and it would be clearer for users to use `.map`, `.agg` or `.transform`, as appropriate in the given situation.

## Proposal

With the above in mind, it is proposed that in the future `apply`, `transform` and `agg` will work as follows:

1. the `agg` and `transform` methods of `Series`, `DataFrame` and `groupby` will always operate series-wise and never element-wise
2. `Series.apply` and `DataFrame.apply` will be deprecated.
3. `groupby.apply` will not be deprecated (because it behaves differently than `Series.apply` and `DataFrame.apply`)

The above changes means that the future behavior, when users want to apply arbitrary callables in pandas, can be described as follows:

1. When users want to operate on single elements in a `Series` or `DataFrame`, they should use `Series.map` and `DataFrame.map` respectively.
2. When users want to aggregate a `Series`, columns/rows of a `DataFrame` or groups in `groupby` objects, they should use `Series.agg`, `DataFrame.agg` and `groupby.agg` respectively.
3. When users want to transform a `Series`, columns/rows of a `DataFrame` or groups in `groupby` objects, they should use `Series.transform`, `DataFrame.transform` and `groupby.transform` respectively.
4. Functions that are not applicable for  `map`, `agg` nor `transform` are considered relatively rare and in the future users should call these functions directly rather than use the `apply` method.

The use of `Series.apply` and  `DataFrame.apply` will after the proposed change in almost all cases be replaced by `map`, `agg` or `transform`. In the very few cases where `Series.apply` and  `DataFrame.apply` cannot be substituted by `map`, `agg` or `transform`, it is proposed that it will be accepted that users will have to find alternative ways to apply the functions, i.e. typically apply the functions manually and possibly concatenating the results.

It can be noted that the behavior of `groupby.agg`, `groupby.transform` and `groupby.apply` are not proposed changed in this PDEP, because `groupby.agg`, `groupby.transform` already behave as desired and `groupby.apply` behaves differently than `Series.apply` and `DataFrame.apply`. Likewise, the behavior when given ufuncs (e.g. `np.sqrt`) and string input (e.g. `"sqrt"`) will remain unchanged, because the behavior is already as intended in all cases.

## Deprecation process

To change the current behavior, it will have to be deprecated. However, `Series.apply` and `DataFrame.apply` are very widely used methods, so will be deprecated very gradually:

This means that in v2.2:

1. Calls to `Series.apply` and `DataFrame.apply`will emit a `DeprecationWarning` with an appropriate deprecation message.
2. A `series_ops_only` with type `bool | lib.NoDefault` parameter will be added to the `agg` and `transform` methods of `Series` and `DataFrame`. When `series_ops_only` is set to False, `agg` and `transform` will behave as they do currently. When set to True, `agg` and `transform` will never operate on elements, but always on Series. When set to `no_default`, `agg` and `transform` will behave as `series_ops_only=False`, but will emit a FutureWarning the current behavior will be reoved in the future.

In Pandas v3.0:
1. Calls to `Series.apply` and `DataFrame.apply` will emit a `FutureWarning` and emit an appropriate deprecation message.
2. The `agg` and `transform` will always operate on series/columns/rows data and the `series_ops_only` parameter will have no effect and be deprecated.

In Pandas v4.0:
1. `Series.apply` and `DataFrame.apply` will be removed from the code base.
2. The `series_ops_only` parameter of agg` and `transform` will be removed from the code base.

## PDEP History

- 24 august 2023: Initial version (proposed to change `Series.apply` and `DataFrame.apply` to always operate on series/columns/rows)
- 17. september 2023: version 2 (renamed and proposing to deprecate `Series.apply` and `DataFrame.apply` and make `agg`/`transform` always operate on series/columns/rows)
