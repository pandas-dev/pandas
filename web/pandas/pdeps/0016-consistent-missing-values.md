# PDEP-16: Consistent missing value handling (with a single NA scalar)

- Created: March 2024
- Status: Under discussion
- Discussion: [#32265](https://github.com/pandas-dev/pandas/issues/32265)
- Author: [Patrick Hoefler](https://github.com/phofl)
          [Joris Van den Bossche](https://github.com/jorisvandenbossche)
- Revision: 1

## Abstract

...

## Background

Currently, pandas handles missing data differently for different data types. We
use different types to indicate that a value is missing: ``np.nan`` for
floating-point data, ``np.nan`` or ``None`` for object-dtype data -- typically
strings or booleans -- with missing values, and ``pd.NaT`` for datetimelike
data. Some other data types, such as integer and bool, cannot store missing data
or are cast to float or object dtype. In addition, pandas 1.0 introduced a new
missing value sentinel, ``pd.NA``, which is being used for the experimental
nullable integer, float, boolean, and string data types, and more recently also
for the pyarrow-backed data types.

These different missing values also have different behaviors in user-facing
operations. Specifically, we introduced different semantics for the nullable
data types for certain operations (e.g. propagating in comparison operations
instead of comparing as False).

The nullable extension dtypes and the `pd.NA` scalar were originally designed to
solve these problems and to provide consistent missing value behavior between
different dtypes. Historically those are used as 1D arrays, which hinders usage
of those dtypes in certain scenarios that rely on the 2D block structure of the
pandas internals for fast operations (``axis=1`` operations, transposing, etc.).

Long term, we want to introduce consistent missing data handling for all data
types. This includes consistent behavior in all operations (indexing, arithmetic
operations, comparisons, etc.) and using a missing value scalar that behaves
consistently.

## Proposal

This proposal aims to unify the missing value handling across all dtypes. This
proposal is not meant to address implementation details, rather to provide a
high level way forward.

1. All data types support missing values and use `pd.NA` exclusively as the
   user-facing missing value indicator.

2. All data types implement consistent missing value "semantics" corresponding
   to the current nullable dtypes using `pd.NA` (i.e. regarding behaviour in
   comparisons, see below for details).

3. As a consequence, pandas will move to nullable extension arrays by default
   for all data types, instead of using the NumPy dtypes that are currently the
   default. To preserve the default 2D block structure of the DataFrame internals,
   the ExtensionArray interface will be extended to support 2D arrays.

4. For backwards compatibility, existing missing value indicators like `NaN` and
   `NaT` will be interpreted as `pd.NA` when introduced in user input, IO or
   through operations (to ensure it keeps being considered as missing).
   Specifically for floating dtypes, in practice this means a float column can
   for now only contain NA values. Potentially distinguishing NA and NaN is left
   for a separate discussion.

This will ensure that all dtypes have consistent missing value handling and there
is no need to upcast if a missing value is inserted into integers or booleans. Those
nullability semantics will be mostly consistent with how PyArrow treats nulls and thus
make switching between both set of dtypes easier. Additionally, it allows the usage of
other Arrow dtypes by default that use the same semantics (bytes, nested dtypes, ...).

In practice, this means solidifying the existing integer, float, boolean and
string nullable data types that already exist, and implementing (variants of)
the categorical, datetimelike and interval data types using `pd.NA`. The
proposal leaves the exact implementation details (e.g. whether to use a mask or
a sentinel (where the best strategy might vary by data type depending on
existing code), or whether to use byte masks vs bitmaps, or whether to use
PyArrow under the hood like the string dtype, etc) out of scope.

This PDEP also does not define the exact API for dtype constructors or
propose a new consistent interface; this is left for a separate discussion
(PDEP-13).

### The `NA` scalar

...

### Missing value semantics


...

## Backward compatibility

...

## Timeline

...

### PDEP History

- March 2024: Initial draft

Note: There is a very long discussion in [GH-32265](https://github.com/pandas-dev/pandas/issues/32265)
that concerns this topic.
