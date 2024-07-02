# PDEP-13: The pandas Logical Type System

- Created: 27 Apr 2024
- Status: Under discussion
- Discussion: [#58141](https://github.com/pandas-dev/pandas/issues/58141)
- Author: [Will Ayd](https://github.com/willayd),
- Revision: 2

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
dtype=pd.StringDtype("pyarrow")
dtype="string[pyarrow]"
dtype="string[pyarrow_numpy]"
dtype=pd.ArrowDtype(pa.string())
dtype=pd.ArrowDtype(pa.large_string())
```

``dtype="string"`` was the first truly new string implementation starting back in pandas 0.23.0, and it is a common pitfall for new users not to understand that there is a huge difference between that and ``dtype=str``. The pyarrow strings have trickled in in more recent memory, but also are very difficult to reason about. The fact that ``dtype="string[pyarrow]"`` is not the same as ``dtype=pd.ArrowDtype(pa.string()`` or ``dtype=pd.ArrowDtype(pa.large_string())`` was a surprise [to the author of this PDEP](https://github.com/pandas-dev/pandas/issues/58321).

While some of these are aliases, the main reason why we have so many different string dtypes is because we have historically used NumPy and created custom missing value solutions around the ``np.nan`` marker, which are incompatible with the ``pd.NA`` sentinel introduced a few years back. Our ``pd.StringDtype()`` uses the pd.NA sentinel, as do our pyarrow based solutions; bridging these into one unified solution has proven challenging.

To try and smooth over the different missing value semantics and how they affect the underlying type system, the status quo has been to add another string dtype. ``string[pyarrow_numpy]`` was an attempt to use pyarrow strings but adhere to NumPy nullability semantics, under the assumption that the latter offers maximum backwards compatibility. However, being the exclusive data type that uses pyarrow for storage but NumPy for nullability handling, this data type just adds more inconsistency to how we handle missing data, a problem we have been attempting to solve back since discussions around pandas2. The name ``string[pyarrow_numpy]`` is not descriptive to end users, and unless it is inferred requires users to explicitly ``.astype("string[pyarrow_numpy]")``, again putting a burden on end users to know what ``pyarrow_numpy`` means and to understand the missing value semantics of both systems.

PDEP-14 has been proposed to smooth over that and change our ``pd.StringDtype()`` to be an alias for ``string[pyarrow_numpy]``. This would at least offer some abstraction to end users who just want strings, but on the flip side would be breaking behavior for users that have already opted into ``dtype="string"`` or ``dtype=pd.StringDtype()`` and the related pd.NA missing value marker for the prior 4 years of their existence.

A logical type system can help us abstract all of these issues. At the end of the day, this PDEP assumes a user wants a string data type. If they call ``Series.str.len()`` against a Series of that type with missing data, they should get back a Series with an integer data type.

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

The third issue is that the extent to which pandas may support any given type is unclear. Issue [#58307](https://github.com/pandas-dev/pandas/issues/58307) highlights one example. It would stand to reason that you could interchangeably use a pandas datetime64 and a pyarrow timestamp, but that is not always true. Another common example is the use of NumPy fixed length strings, which users commonly try to use even though we claim no real support for them (see [#5764](https://github.com/pandas-dev/pandas/issues/57645)).

## Assessing the Current Type System(s)

A best effort at visualizing the current type system(s) with types that we currently "support" or reasonably may want to is shown [in this comment](https://github.com/pandas-dev/pandas/issues/58141#issuecomment-2047763186). Note that this does not include the ["pyarrow_numpy"](https://github.com/pandas-dev/pandas/pull/58451) string data type or the string data type that uses the NumPy 2.0 variable length string data type (see [comment](https://github.com/pandas-dev/pandas/issues/57073#issuecomment-2080798081)) as they are under active discussion.

## Proposal

### Proposed Logical Types

Derived from the hierarchical visual in the previous section, this PDEP proposes that pandas supports at least all of the following _logical_ types, excluding any type widths for brevity:

  - Signed Integer
  - Unsigned Integer
  - Floating Point
  - Fixed Point
  - Boolean
  - Date
  - Datetime
  - Duration
  - CalendarInterval
  - Period
  - Binary
  - String
  - Map
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
  - pd.MapDtype()  # or pd.DictDtype()
  - pd.ListDtype()
  - pd.StructDtype()
  - pd.ObjectDtype()

The storage / backend to each of these types is left as an implementation detail. The fact that ``pd.StringDtype()`` may be backed by Arrow while ``pd.PeriodDtype()`` continues to be a custom solution is of no concern to the end user. Over time this will allow us to adopt more Arrow behind the scenes without breaking the front end for our end users, but _still_ giving us the flexibility to produce data types that Arrow will not implement (e.g. ``pd.ObjectDtype()``).

The methods of each logical type are expected in turn to yield another logical type. This can enable us to smooth over differences between the NumPy and Arrow world, while also leveraging the best of both backends. To illustrate, let's look at some methods where the return type today deviates for end users depending on if they are using NumPy-backed data types or Arrow-backed data types. The equivalent PDEP-13 logical data type is presented as the last column:

| Method             | NumPy-backed result | Arrow Backed result type | PDEP-13 result type            |
|--------------------|---------------------|--------------------------|--------------------------------|
| Series.str.len()   | np.float64          | pa.int64()               | pd.Int64Dtype()                |
| Series.str.split() | object              | pa.list(pa.string())     | pd.ListDtype(pa.StringDtype()) |
| Series.dt.date     | object              | pa.date32()              | pd.DateDtype()                 |

The ``Series.dt.date`` example is worth an extra look - with a PDEP-13 logical type system in place we would theoretically have the ability to keep our default ``pd.DatetimeDtype()`` backed by our current NumPy-based array but leverage pyarrow for the ``Series.dt.date`` solution, rather than having to implement a DateArray ourselves.

While this PDEP proposes reusing existing extension types, it also necessitates extending those types with extra metadata:

```python
class BaseType:

    @property
    def data_manager -> Literal["numpy", "pyarrow"]:
        """
        Who manages the data buffer - NumPy or pyarrow
        """
        ...

    @property
    def physical_type:
        """
        For logical types which may have different implementations, what is the
        actual implementation? For pyarrow strings this may mean pa.string() versus
        pa.large_string() versrus pa.string_view(); for NumPy this may mean object
        or their 2.0 string implementation.
        """
        ...

    @property
    def na_marker -> pd.NA|np.nan|pd.NaT:
        """
        Sentinel used to denote missing values
        """
        ...
```

``na_marker`` is expected to be read-only (see next section). For advanced users that have a particular need for a storage type, they may be able to construct the data type via ``pd.StringDtype(data_manager=np)`` to assert NumPy managed storage. While the PDEP allows constructing in this fashion, operations against that data make no guarantees that they will respect the storage backend and are free to convert to whichever storage the internals of pandas considers optimal (Arrow will typically be preferred).

### Missing Value Handling

Missing value handling is a tricky area as developers are split between pd.NA semantics versus np.nan, and the transition path from one to the other is not always clear.

Because this PDEP proposes reuse of the existing pandas extension type system, the default missing value marker will consistently be ``pd.NA``. However, to help with backwards compatibility for users that heavily rely on the equality semantics of np.nan, an option of ``pd.na_marker = "legacy"`` can be set. This would mean that the missing value indicator for logical types would be:

| Logical Type      | Default Missing Value | Legacy Missing Value |
| pd.BooleanDtype() | pd.NA                 | np.nan               |
| pd.IntXXType()    | pd.NA                 | np.nan               |
| pd.FloatXXType()  | pd.NA                 | np.nan               |
| pd.StringDtype()  | pd.NA                 | np.nan               |
| pd.DatetimeType() | pd.NA                 | pd.NaT               |

However, all data types for which there is no legacy NumPy-backed equivalent will continue to use ``pd.NA``, even in "legacy" mode. Legacy is provided only for backwards compatibility, but pd.NA usage is encouraged going forward to give users one exclusive missing value indicator.

### Transitioning from Current Constructors

To maintain a consistent path forward, _all_ constructors with the implementation of this PDEP are expected to map to the logical types. This means that providing ``np.int64`` as the data type argument makes no guarantee that you actually get a NumPy managed storage buffer; pandas reserves the right to optimize as it sees fit and may decide instead to just pyarrow.

The theory behind this is that the majority of users are not expecting anything particular from NumPy to happen when they say ``dtype=np.int64``. The expectation is that a user just wants _integer_ data, and the ``np.int64`` specification owes to the legacy of pandas' evolution.

This PDEP makes no guarantee that we will stay that way forever; it is certainly reasonable that a few years down the road we deprecate and fully stop support for backend-specifc constructors like ``np.int64`` or ``pd.ArrowDtype(pa.int64())``. However, for the execution of this PDEP, such an initiative is not in scope.

## PDEP-11 History

- 27 April 2024: Initial version
- 10 May 2024: First revision
