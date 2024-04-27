# PDEP-13: The pandas Logical Type System

- Created: 27 Apr 2024
- Status: Draft
- Discussion: [#58141](https://github.com/pandas-dev/pandas/issues/58141)
- Author: [Will Ayd](https://github.com/willayd),
- Revision: 1

## Abstract

This PDEP proposes a logical type system for pandas to abstract underlying library differences from end users, clarify the scope of pandas type support, and give pandas developers more flexibility to manage the implementation of types.

## Background

When pandas was originally built, the data types that it exposed were a subset of the NumPy type system. Starting back in version 0.23.0, pandas introduced Extension types, which it also began to use internally as a way of creating arrays instead of exclusively relying upon NumPy. Over the course of the 1.x releases, pandas began using pyarrow for string storage and in version 1.5.0 introduced the high level ``pd.ArrowDtype`` wrapper.

While these new type systems have brought about many great features, they have surfaced three major problems. The first is that we put the onus on users to understand the differences of the physical type implementations. Consider the many ways pandas allows you to create a "string":

```python
dtype=object
dtype=str
dtype="string"
dtype=pd.StringDtype()
dtype=pd.StringDtype("pyarrow")
dtype="string[pyarrow]"
dtype="string[pyarrow_numpy]"
dtype=pd.ArrowDtype(pa.string())
```

Keeping track of all of these iterations and their subtle differences is difficult even for [core maintainers](https://github.com/pandas-dev/pandas/issues/58321).

The second problem is that the conventions for constructing types from a given type backend are inconsistent. Let's review string aliases used to construct certain types:

| logical type | NumPy   | pandas extension | pyarrow                  |
| int          | "int64" | "Int64"          | "int64[pyarrow]"         |
| string       | N/A     | "string"         | N/A                      |
| datetime     | N/A     | "datetime64[us]" | "timestamp[us][pyarrow]" |

"string[pyarrow]" is excluded from the above table because it is misleading; while "int64[pyarrow]" definitely gives you a pyarrow backed string, "string[pyarrow]" gives you a pandas extension array which itself then uses pyarrow, which can introduce behavior differences (see [issue 58321](https://github.com/pandas-dev/pandas/issues/58321)).

If you wanted to try and be more explicit about using pyarrow, you could use the ``pd.ArrowDtype`` wrapper. But this unfortunately exposes gaps when trying to use that pattern across all backends:

| logical type | NumPy    | pandas extension | pyarrow                           |
| int          | np.int64 | pd.Int64Dtype()  | pd.ArrowDtype(pa.int64())         |
| string       | N/A      | pd.StringDtype() | pd.ArrowDtype(pa.string())        |
| datetime     | N/A      | ???              | pd.ArrowDtype(pa.timestamp("us")) |

It would stand to reason in this approach that you could use a ``pd.DatetimeDtype()`` but no such type exists (there is a ``pd.DatetimeTZDtype`` which requires a timezone).

The third issue is that the extent to which pandas may support any given type is unclear. Issue [#58307](https://github.com/pandas-dev/pandas/issues/58307) highlights one example. It would stand to reason that you could interchangeably use a pandas datetime64 and a pyarrow timestamp, but that is not always true. Another common example is the use of NumPy fixed length strings, which users commonly try to use even though we claim no real support for them (see [#5764](https://github.com/pandas-dev/pandas/issues/57645)).

## Assessing the Current Type System(s)

A best effort at visualizing the current type system(s) with types that we currently "support" or reasonably may want to is shown [in this comment](https://github.com/pandas-dev/pandas/issues/58141#issuecomment-2047763186). Note that this does not include the ["pyarrow_numpy"](https://github.com/pandas-dev/pandas/pull/58451) string data type or the string data type that uses the NumPy 2.0 variable length string data type (see [comment](https://github.com/pandas-dev/pandas/issues/57073#issuecomment-2080798081)) as they are under active discussion.

## Proposal

Derived from the hierarchical visual in the previous section, this PDEP proposes that pandas supports at least all of the following _logical_ types, excluding any type widths for brevity:

  - Signed Integer
  - Unsigned Integer
  - Floating Point
  - Fixed Point
  - Boolean
  - Date
  - Datetime
  - Duration
  - Interval
  - Period
  - Binary
  - String
  - Dictionary
  - List
  - Struct

To ensure we maintain all of the current functionality of our existing type system(s), a base type structure would need to look something like:

```python
class BaseType:

    @property
    def dtype_backend -> Literal["pandas", "numpy", "pyarrow"]:
        """
        Library is responsible for the array implementation
        """
        ...

    @property
    def physical_type:
        """
        How does the backend physically implement this logical type? i.e. our
        logical type may be a "string" and we are using pyarrow underneath -
        is it a pa.string(), pa.large_string(), pa.string_view() or something else?
        """
        ...

    @property
    def missing_value_marker -> pd.NA|np.nan:
        """
        Sentinel used to denote missing values
        """
        ...
```

The theory behind this PDEP is that most users /should not care/ about the physical type that is being used. But if the abstraction our logical type provides is too much, a user could at least inspect and potentially configure which physical type to use.

With regards to how we may expose such types to end users, there are two currently recognized proposals. The first would use factory functions to create the logical type, i.e. something like:

```python
pd.Series(["foo", "bar", "baz"], dtype=pd.string()  # assumed common case
pd.Series(["foo", "bar", "baz"], dtype=pd.string(missing_value_marker=np.nan)
pd.Series(["foo", "bar", "baz"], dtype=pd.string(physical_type=pa.string_view())
```

Another approach would be to classes:

```python
pd.Series(["foo", "bar", "baz"], dtype=pd.StringDtype()
pd.Series(["foo", "bar", "baz"], dtype=pd.StringDtype(missing_value_marker=np.nan)
pd.Series(["foo", "bar", "baz"], dtype=pd.StringDtype(physical_type=pa.string_view())
```
Note that the class-based approach would reuse existing classes like ``pd.StringDtype`` but change their purpose, whereas the factory function would more explicitly be a new approach. This is an area that requires more discussion amongst the team.

## String Type Arguments

This PDEP proposes that we maintain only a small subset of string arguments that can be used to construct logical types. Those string arguments are:

  - intXX
  - uintXX
  - floatXX
  - string
  - datetime64[unit]
  - datetime64[unit, tz]

However, new code should be encouraged to use the logical constructors outlined previously. Particularly for aggregate types, trying to encode all of the information into a string can become unwieldy. Instead, keyword argument use should be encouraged:

```python
pd.Series(dtype=pd.list(value_type=pd.string()))
```

## Bridging Type Systems

An interesting question arises when a user constructs two logical types with differing physical types. If one is backed by NumPy and the other is backed by pyarrow, what should happen?

This PDEP proposes the following backends should be prioritized in the following order (1. is the highest priority):

  1. Arrow
  2. pandas
  3. NumPy

One reason for this is that Arrow represents the most efficient and assumedly least-lossy physical representations. An obvious example comes when a pyarrow int64 array with missing data gets added to a NumPy int64 array; casting to the latter would lose data. Another reason is that Arrow represents the fastest growing ecosystem of tooling, and the PDEP author believes improving pandas's interoperability within that landscape is extremely important.

Aside from the backend, the C standard rules for [implicit conversion](https://en.cppreference.com/w/c/language/conversion) should apply to the data buffer, i.e. adding a pyarrow int8 array to a NumPy uint64 array should produce a pyarrow uint64 array.

For more expensive conversions, pandas retains the right to throw warnings or even error out when two of the same logical type with differing physical types is added. For example, attempting to do string concatenation of string arrays backed by pyarrow and Python objects may throw a ``PerformanceWarning``, or maybe even a ``MemoryError`` if such a conversion exhausts the available system resources.

The ``BaseType`` proposed above also has a property for the ``missing_value_marker``. Operations that use two logical types with different missing value markers should raise, as there is no clear way to prioritize between the various sentinels.

## PDEP-11 History

- 27 April 2024: Initial version
