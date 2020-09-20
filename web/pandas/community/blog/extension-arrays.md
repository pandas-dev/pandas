Title: pandas extension arrays
Date: 2019-01-04

# pandas extension arrays

Extensibility was a major theme in pandas development over the last couple of
releases. This post introduces the pandas extension array interface: the
motivation behind it and how it might affect you as a pandas user. Finally, we
look at how extension arrays may shape the future of pandas.

Extension Arrays are just one of the changes in pandas 0.24.0. See the
[whatsnew][whatsnew] for a full changelog.

## The Motivation

Pandas is built on top of NumPy. You could roughly define a Series as a wrapper
around a NumPy array, and a DataFrame as a collection of Series with a shared
index. That's not entirely correct for several reasons, but I want to focus on
the "wrapper around a NumPy array" part. It'd be more correct to say "wrapper
around an array-like object".

Pandas mostly uses NumPy's builtin data representation; we've restricted it in
places and extended it in others. For example, pandas' early users cared greatly
about timezone-aware datetimes, which NumPy doesn't support. So pandas
internally defined a `DatetimeTZ` dtype (which mimics a NumPy dtype), and
allowed you to use that dtype in `Index`, `Series`, and as a column in a
`DataFrame`. That dtype carried around the tzinfo, but wasn't itself a valid
NumPy dtype.

As another example, consider `Categorical`. This actually composes *two* arrays:
one for the `categories` and one for the `codes`. But it can be stored in a
`DataFrame` like any other column.

Each of these extension types pandas added is useful on its own, but carries a
high maintenance cost. Large sections of the codebase need to be aware of how to
handle a NumPy array or one of these other kinds of special arrays. This made
adding new extension types to pandas very difficult.

Anaconda, Inc. had a client who regularly dealt with datasets with IP addresses.
They wondered if it made sense to add an [IPArray][IPArray] to pandas. In the
end, we didn't think it passed the cost-benefit test for inclusion in pandas
*itself*, but we were interested in defining an interface for third-party
extensions to pandas. Any object implementing this interface would be allowed in
pandas. I was able to write [cyberpandas][cyberpandas] outside of pandas, but it
feels like using any other dtype built into pandas.

## The Current State

As of pandas 0.24.0, all of pandas' internal extension arrays (Categorical,
Datetime with Timezone, Period, Interval, and Sparse) are now built on top of
the ExtensionArray interface. Users shouldn't notice many changes. The main
thing you'll notice is that things are cast to `object` dtype in fewer places,
meaning your code will run faster and your types will be more stable. This
includes storing `Period` and `Interval` data in `Series` (which were previously
cast to object dtype).

Additionally, we'll be able to add *new* extension arrays with relative ease.
For example, 0.24.0 (optionally) solved one of pandas longest-standing pain
points: missing values casting integer-dtype values to float.


```python
>>> int_ser = pd.Series([1, 2], index=[0, 2])
>>> int_ser
0    1
2    2
dtype: int64

>>> int_ser.reindex([0, 1, 2])
0    1.0
1    NaN
2    2.0
dtype: float64
```

With the new [IntegerArray][IntegerArray] and nullable integer dtypes, we can
natively represent integer data with missing values.

```python
>>> int_ser = pd.Series([1, 2], index=[0, 2], dtype=pd.Int64Dtype())
>>> int_ser
0    1
2    2
dtype: Int64

>>> int_ser.reindex([0, 1, 2])
0      1
1    NaN
2      2
dtype: Int64
```

One thing it does slightly change how you should access the raw (unlabeled)
arrays stored inside a Series or Index, which is occasionally useful. Perhaps
the method you're calling only works with NumPy arrays, or perhaps you want to
disable automatic alignment.

In the past, you'd hear things like "Use `.values` to extract the NumPy array
from a Series or DataFrame." If it were a good resource, they'd tell you that's
not *entirely* true, since there are some exceptions. I'd like to delve into
those exceptions.

The fundamental problem with `.values` is that it serves two purposes:

1. Extracting the array backing a Series, Index, or DataFrame
2. Converting the Series, Index, or DataFrame to a NumPy array

As we saw above, the "array" backing a Series or Index might not be a NumPy
array, it may instead be an extension array (from pandas or a third-party
library). For example, consider `Categorical`,

```python
>>> cat = pd.Categorical(['a', 'b', 'a'], categories=['a', 'b', 'c'])
>>> ser = pd.Series(cat)
>>> ser
0    a
1    b
2    a
dtype: category
Categories (3, object): ['a', 'b', 'c']

>>> ser.values
[a, b, a]
Categories (3, object): ['a', 'b', 'c']
```

In this case `.values` is a Categorical, not a NumPy array. For period-dtype
data, `.values` returns a NumPy array of `Period` objects, which is expensive to
create. For timezone-aware data, `.values` converts to UTC and *drops* the
timezone info. These kind of surprises (different types, or expensive or lossy
conversions) stem from trying to shoehorn these extension arrays into a NumPy
array. But the entire point of an extension array is for representing data NumPy
*can't* natively represent.

To solve the `.values` problem, we've split its roles into two dedicated methods:

1. Use `.array` to get a zero-copy reference to the underlying data
2. Use `.to_numpy()` to get a (potentially expensive, lossy) NumPy array of the
   data.

So with our Categorical example,

```python
>>> ser.array
[a, b, a]
Categories (3, object): ['a', 'b', 'c']

>>> ser.to_numpy()
array(['a', 'b', 'a'], dtype=object)
```

To summarize:

- `.array` will *always* be a an ExtensionArray, and is always a zero-copy
   reference back to the data.
- `.to_numpy()` is *always* a NumPy array, so you can reliably call
   ndarray-specific methods on it.

You shouldn't ever need `.values` anymore.

## Possible Future Paths

Extension Arrays open up quite a few exciting opportunities. Currently, pandas
represents string data using Python objects in a NumPy array, which is slow.
Libraries like [Apache Arrow][arrow] provide native support for variable-length
strings, and the [Fletcher][fletcher] library provides pandas extension arrays
for Arrow arrays. It will allow [GeoPandas][geopandas] to store geometry data
more efficiently. Pandas (or third-party libraries) will be able to support
nested data, data with units, geo data, GPU arrays. Keep an eye on the
[pandas ecosystem][eco] page, which will keep track of third-party extension
arrays. It's an exciting time for pandas development.

## Other Thoughts

I'd like to emphasize that this is an *interface*, and not a concrete array
implementation. We are *not* reimplementing NumPy here in pandas. Rather, this
is a way to take any array-like data structure (one or more NumPy arrays, an
Apache Arrow array, a CuPy array) and place it inside a DataFrame. I think
getting pandas out of the array business, and instead thinking about
higher-level tabular data things, is a healthy development for the project.

This works perfectly with NumPy's [`__array_ufunc__`][ufunc] protocol and
[NEP-18][nep18]. You'll be able to use the familiar NumPy API on objects that
aren't backed by NumPy memory.

## Upgrade

These new goodies are all available in the recently released pandas 0.24.

conda:

    conda install -c conda-forge pandas

pip:

    pip install --upgrade pandas

As always, we're happy to hear feedback on the [mailing list][ml],
[@pandas-dev][twitter], or [issue tracker][tracker].

Thanks to the many contributors, maintainers, and [institutional
partners][partners] involved in the pandas community.


[IPArray]: https://github.com/pandas-dev/pandas/issues/18767
[cyberpandas]: https://cyberpandas.readthedocs.io
[IntegerArray]: http://pandas.pydata.org/pandas-docs/version/0.24/reference/api/pandas.arrays.IntegerArray.html
[fletcher]: https://github.com/xhochy/fletcher
[arrow]: https://arrow.apache.org
[ufunc]: https://numpy.org/neps/nep-0013-ufunc-overrides.html
[nep18]: https://www.numpy.org/neps/nep-0018-array-function-protocol.html
[ml]: https://mail.python.org/mailman/listinfo/pandas-dev
[twitter]: https://twitter.com/pandas_dev
[tracker]: https://github.com/pandas-dev/pandas/issues
[partners]: https://github.com/pandas-dev/pandas-governance/blob/master/people.md
[eco]: http://pandas.pydata.org/pandas-docs/stable/ecosystem.html#extension-data-types
[whatsnew]: http://pandas.pydata.org/pandas-docs/version/0.24/whatsnew/v0.24.0.html
[geopandas]: https://github.com/geopandas/geopandas
