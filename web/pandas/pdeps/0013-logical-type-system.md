# PDEP-13: The pandas Logical Type System

- Created: 27 Apr 2024
- Status: Under discussion
- Discussion: [#58141](https://github.com/pandas-dev/pandas/issues/58141)
- Author: [Will Ayd](https://github.com/willayd),
- Revision: 3

## Abstract

This PDEP proposes a logical type system for pandas which will decouple user semantics (i.e. _this column should be an integer_) from the pandas implementation/internal (i.e. _we will use NumPy/pyarrow/X to store this array_).By decoupling these through a logical type system, the expectation is that this PDEP will:

  * Abstract underlying library differences from end users
  * More clearly define the types of data pandas supports
  * Allow pandas developers more flexibility to manage type implementations
  * Pave the way for continued adoption of Arrow in the pandas code base

## Background

When pandas was originally built, the data types that it exposed were a subset of the NumPy type system. Starting back in version 0.23.0, pandas introduced Extension types, which it also began to use internally as a way of creating arrays instead of exclusively relying upon NumPy. Over the course of the 1.x releases, pandas began using pyarrow for string storage and in version 1.5.0 introduced the high level ``pd.ArrowDtype`` wrapper.

While these new type systems have brought about many great features, they have surfaced three major problems.

### Problem 1: Inconsistent Type Naming / Behavior

There is no better example of our current type system being problematic than strings. Let's assess the number of string iterations a user could create (this is a non-exhaustive list):

```python
dtype=object
dtype=str
dtype="string"
dtype=pd.StringDtype()
dtype=pd.StringDtype("python", na_value=np.nan)
dtype=pd.StringDtype("python", na_value=pd.NA)
dtype=pd.StringDtype("pyarrow")
dtype="string[pyarrow]"
dtype="string[pyarrow_numpy]"  # added in 2.1, deprecated in 2.3
dtype=pd.ArrowDtype(pa.string())
dtype=pd.ArrowDtype(pa.large_string())
```

``dtype="string"`` was the first truly new string implementation starting back in pandas 0.23.0, and it is a common pitfall for new users not to understand that there is a huge difference between that and ``dtype=str``. The pyarrow strings have trickled in in more recent releases, but also are very difficult to reason about. The fact that ``dtype="string[pyarrow]"`` is not the same as ``dtype=pd.ArrowDtype(pa.string()`` or ``dtype=pd.ArrowDtype(pa.large_string())`` was a surprise [to the author of this PDEP](https://github.com/pandas-dev/pandas/issues/58321).

While some of these are aliases, the main reason why we have so many different string dtypes is because we have historically used NumPy and created custom missing value solutions around the ``np.nan`` marker, which are incompatible with ``pd.NA``. Our ``pd.StringDtype()`` uses the pd.NA sentinel, as do our pyarrow based solutions; bridging these into one unified solution has proven challenging.

To try and smooth over the different missing value semantics and how they affect the underlying type system, the status quo has always been to add another string dtype. With PDEP-14 we now have a "compatibility" string of ``pd.StringDtype("python|pyarrow", na_value=np.nan)`` that makes a best effort to move users towards all the benefits of PyArrow strings (assuming pyarrow is installed) while retaining backwards-compatible missing value handling with ``np.nan`` as the missing value marker. The usage of the ``pd.StringDtype`` in this manner is a good stepping stone towards the goals of this PDEP, although it is stuck in an "in-between" state without other types following suit.

For instance, if a user calls ``Series.value_counts()`` on the ``pd.StringDtype()``, the type of the returned Series can vary wildly, and in non-obvious ways:

```python
>>> pd.Series(["x"], dtype=pd.StringDtype("python", na_value=pd.NA)).value_counts().dtype
Int64Dtype()
>>> pd.Series(["x"], dtype=pd.StringDtype("pyarrow", na_value=pd.NA)).value_counts().dtype
int64[pyarrow]
>>> pd.Series(["x"], dtype=pd.StringDtype("python", na_value=np.nan)).value_counts().dtype
Int64Dtype()
>>> pd.Series(["x"], dtype=pd.StringDtype("pyarrow", na_value=np.nan)).value_counts().dtype
dtype('int64')
```

It is also worth noting that different methods will return different data types. For a pyarrow-backed string with pd.NA, ``Series.value_counts()`` returns a ``int64[pyarrow]`` but ``Series.str.len()`` returns a ``pd.Int64Dtype()``.

A logical type system can help us abstract all of these issues. At the end of the day, this PDEP assumes a user wants a string data type. If they call ``Series.str.value_counts()`` against a Series of that type with missing data, they should get back a Series with an integer data type.

### Problem 2: Inconsistent Constructors

The second problem is that the conventions for constructing types from the various _backends_ are inconsistent. Let's review string aliases used to construct certain types:

| logical type | NumPy   | pandas extension | pyarrow                  |
| int          | "int64" | "Int64"          | "int64[pyarrow]"         |
| string       | N/A     | "string"         | N/A                      |
| datetime     | N/A     | "datetime64[us]" | "timestamp[us][pyarrow]" |

"string[pyarrow]" is excluded from the above table because it is misleading; while "int64[pyarrow]" definitely gives you a pyarrow backed string, "string[pyarrow]" gives you a pandas extension array which itself then uses pyarrow. Subtleties like this then lead to behavior differences (see [issue 58321](https://github.com/pandas-dev/pandas/issues/58321)).

If you wanted to try and be more explicit about using pyarrow, you could use the ``pd.ArrowDtype`` wrapper. But this unfortunately exposes gaps when trying to use that pattern across all backends:

| logical type | NumPy    | pandas extension | pyarrow                           |
| int          | np.int64 | pd.Int64Dtype()  | pd.ArrowDtype(pa.int64())         |
| string       | N/A      | pd.StringDtype() | pd.ArrowDtype(pa.string())        |
| datetime     | N/A      | ???              | pd.ArrowDtype(pa.timestamp("us")) |

It would stand to reason in this approach that you could use a ``pd.DatetimeDtype()`` but no such type exists (there is a ``pd.DatetimeTZDtype`` which requires a timezone).

### Problem 3: Lack of Clarity on Type Support

The third issue is that the extent to which pandas may support any given type is unclear. Issue [#58307](https://github.com/pandas-dev/pandas/issues/58307) highlights one example. It would stand to reason that you could interchangeably use a pandas datetime64 and a pyarrow timestamp, but that is not always true. Another example is the use of NumPy fixed length strings, which users commonly try to use even though we claim no real support for them (see [#5764](https://github.com/pandas-dev/pandas/issues/57645)).

## Assessing the Current Type System(s)

A best effort at visualizing the current type system(s) with types that we currently "support" or reasonably may want to is shown [in this comment](https://github.com/pandas-dev/pandas/issues/58141#issuecomment-2047763186). Note that this does not include the ["pyarrow_numpy"](https://github.com/pandas-dev/pandas/pull/58451) string data type or the string data type that uses the NumPy 2.0 variable length string data type (see [comment](https://github.com/pandas-dev/pandas/issues/57073#issuecomment-2080798081)) as they are under active discussion.

## Proposal

### Proposed Logical Types

Derived from the hierarchical visual in the previous section, this PDEP proposes that pandas supports at least all of the following _logical_ types, excluding any type widths for brevity:

  - Signed Integer
  - Unsigned Integer
  - Floating Point
  - Decimal
  - Boolean
  - Date
  - Datetime
  - Duration
  - CalendarInterval
  - Period
  - Binary
  - String
  - Dict
  - List
  - Struct
  - Interval
  - Object
  - Null

One of the major problems this PDEP has tried to highlight is the historical tendency of our team to "create more types" to solve existing problems. To minimize the need for that, this PDEP proposes re-using our existing extension types where possible, and only adding new ones where they do not exist.

The existing extension types which will become our "logical types" are:

  - pd.StringDtype()
  - pd.IntXXDtype()
  - pd.UIntXXDtype()
  - pd.FloatXXDtype()
  - pd.BooleanDtype()
  - pd.PeriodDtype(freq)
  - pd.IntervalDtype()

To satisfy all of the types highlighted above, this would require the addition of:

  - pd.DecimalDtype()
  - pd.DateDtype()
  - pd.DatetimeDtype(unit, tz)
  - pd.Duration()
  - pd.CalendarInterval()
  - pd.BinaryDtype()
  - pd.DictDtype()
  - pd.ListDtype()
  - pd.StructDtype()
  - pd.ObjectDtype()
  - pd.NullDtype()

The storage / backend to each of these types is left as an implementation detail. The fact that ``pd.StringDtype()`` may be backed by Arrow while ``pd.PeriodDtype()`` continues to be a custom solution is of no concern to the end user. Over time this will allow us to adopt more Arrow behind the scenes without breaking the front end for our end users, but _still_ giving us the flexibility to produce data types that Arrow will not implement (e.g. ``pd.ObjectDtype()``).

The methods of each logical type are expected in turn to yield another logical type. This can enable us to smooth over differences between the NumPy and Arrow world, while also leveraging the best of both backends. To illustrate, let's look at some methods where the return type today deviates for end users depending on if they are using NumPy-backed data types or Arrow-backed data types. The equivalent PDEP-13 logical data type is presented as the last column:

| Method             | NumPy-backed result | Arrow Backed result type | PDEP-13 result type            |
|--------------------|---------------------|--------------------------|--------------------------------|
| Series.str.len()   | np.float64          | pa.int64()               | pd.Int64Dtype()                |
| Series.str.split() | object              | pa.list(pa.string())     | pd.ListDtype(pa.StringDtype()) |
| Series.dt.date     | object              | pa.date32()              | pd.DateDtype()                 |

The ``Series.dt.date`` example is worth an extra look - with a PDEP-13 logical type system in place we would theoretically have the ability to keep our default ``pd.DatetimeDtype()`` backed by our current NumPy-based array but leverage pyarrow for the ``Series.dt.date`` solution, rather than having to implement a DateArray ourselves.

To implement this PDEP, we expect all of the logical types to have at least the following metadata:

    * storage: Either "numpy" or "pyarrow". Describes the library used to create the data buffer
    * physical_type: Can expose the physical type being used. As an example, StringDtype could return pa.string_view
    * na_value: Either pd.NA, np.nan, or pd.NaT.

While these attributes are exposed as construction arguments to end users, users are highly discouraged from trying to control them directly. Put explicitly, this PDEP allows a user to request a ``pd.XXXDtype(storage="numpy")`` to request a NumPy-backed array, if possible. While pandas may respect that during construction, operations against that data make no guarantees that the storage backend will be persisted through, giving pandas the freedom to convert to whichever storage is internally optimal (Arrow will typically be preferred).

### Missing Value Handling

Missing value handling is a tricky area as developers are split between pd.NA semantics versus np.nan, and the transition path from one to the other is not always clear. This PDEP does not aim to "solve" that issue per se (for that discussion, please refer to PDEP-16), but aims to provide a go-forward path that strikes a reasonable balance between backwards compatibility and a consistent missing value approach in the future.

This PDEP proposes that the default missing value for logical types is ``pd.NA``. The reasoning is two-fold:

  1. We are in many cases re-using extension types as logical types, which mostly use pd.NA (StrDtype and datetimes are the exception)
  2. For new logical types that have nothing to do with NumPy, using np.nan as a missing value marker is an odd fit

However, to help with backwards compatibility for users that heavily rely on the semantics of ``np.nan`` or ``pd.NaT``, an option of ``pd.na_value = "legacy"`` can be set. This would mean that the missing value indicator for logical types would be:

| Logical Type      | Default Missing Value | Legacy Missing Value |
| pd.BooleanDtype() | pd.NA                 | np.nan               |
| pd.IntXXType()    | pd.NA                 | np.nan               |
| pd.FloatXXType()  | pd.NA                 | np.nan               |
| pd.StringDtype()  | pd.NA                 | np.nan               |
| pd.DatetimeType() | pd.NA                 | pd.NaT               |

However, all data types for which there is no legacy NumPy-backed equivalent will continue to use ``pd.NA``, even in "legacy" mode. Legacy is provided only for backwards compatibility, but ``pd.NA`` usage is encouraged going forward to give users one exclusive missing value indicator and better align with the goals of PDEP-16.

### Transitioning from Current Constructors

To maintain a consistent path forward, _all_ constructors with the implementation of this PDEP are expected to map to the logical types. This means that providing ``np.int64`` as the data type argument makes no guarantee that you actually get a NumPy managed storage buffer; pandas reserves the right to optimize as it sees fit and may decide instead to use PyArrow.

The theory behind this is that the majority of users are not expecting anything particular from NumPy to happen when they say ``dtype=np.int64``. The expectation is that a user just wants _integer_ data, and the ``np.int64`` specification owes to the legacy of pandas' evolution.

This PDEP makes no guarantee that we will stay that way forever; it is certainly reasonable that, in the future, we deprecate and fully stop support for backend-specifc constructors like ``np.int64`` or ``pd.ArrowDtype(pa.int64())``. However, for the execution of this PDEP, such an initiative is not in scope.

## PDEP-13 History

- 27 April 2024: Initial version
- 10 May 2024: First revision
- 01 Aug 2024: Revisions for PDEP-14 and PDEP-16
