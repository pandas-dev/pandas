# PDEP-XX: Dedicated string data type for pandas 3.0

- Created: May 3, 2024
- Status: Under discussion
- Discussion: https://github.com/pandas-dev/pandas/pull/58551
- Author: [Joris Van den Bossche](https://github.com/jorisvandenbossche)
- Revision: 1

## Abstract

This PDEP proposes to introduce a dedicated string dtype that will be used by
default in pandas 3.0:

* In pandas 3.0, enable a "string" dtype by default, using PyArrow if available
  or otherwise the numpy object-dtype alternative.
* The default string dtype will use missing value semantics (using NaN) consistent
  with the other default data types.

This will give users a long-awaited proper string dtype for 3.0, while 1) not
(yet) making PyArrow a _hard_ dependency, but only a dependency used by default,
and 2) leaving room for future improvements (different missing value semantics,
using NumPy 2.0, etc).

# Dedicated string data type for pandas 3.0

## Background

Currently, pandas by default stores text data in an `object`-dtype NumPy array.
The current implementation has two primary drawbacks. First, `object` dtype is
not specific to strings: any Python object can be stored in an `object`-dtype
array, not just strings, and seeing `object` as the dtype for a column with
strings is confusing for users. Second: this is not efficient (all string
methods on a Series are eventually calling Python methods on the individual
string objects).

To solve the first issue, a dedicated extension dtype for string data has
already been
[added in pandas 1.0](https://pandas.pydata.org/docs/whatsnew/v1.0.0.html#dedicated-string-data-type).
This has always been opt-in for now, requiring users to explicitly request the
dtype (with `dtype="string"` or `dtype=pd.StringDtype()`). The array backing
this string dtype was initially almost the same as the default implementation,
i.e. an `object`-dtype NumPy array of Python strings.

To solve the second issue (performance), pandas contributed to the development
of string kernels in the PyArrow package, and a variant of the string dtype
backed by PyArrow was
[added in pandas 1.3](https://pandas.pydata.org/docs/whatsnew/v1.3.0.html#pyarrow-backed-string-data-type).
This could be specified with the `storage` keyword in the opt-in string dtype
(`pd.StringDtype(storage="pyarrow")`).

Since its introduction, the `StringDtype` has always been opt-in, and has used
the experimental `pd.NA` sentinel for missing values (which was also [introduced
in pandas 1.0](https://pandas.pydata.org/docs/whatsnew/v1.0.0.html#experimental-na-scalar-to-denote-missing-values)).
However, up to this date, pandas has not yet taken the step to use `pd.NA` by
default, and thus the `StringDtype` deviates in missing value behaviour compared
to the default data types.

In 2023, [PDEP-10](https://pandas.pydata.org/pdeps/0010-required-pyarrow-dependency.html)
proposed to start using a PyArrow-backed string dtype by default in pandas 3.0
(i.e. infer this type for string data instead of object dtype). To ensure we
could use the variant of `StringDtype` backed by PyArrow instead of Python
objects (for better performance), it proposed to make `pyarrow` a new required
runtime dependency of pandas.

In the meantime, NumPy has also been working on a native variable-width string
data type, which will be available [starting with NumPy
2.0](https://numpy.org/devdocs/release/2.0.0-notes.html#stringdtype-has-been-added-to-numpy).
This can provide a potential alternative to PyArrow for implementing a string
data type in pandas that is not backed by Python objects.

After acceptance of PDEP-10, two aspects of the proposal have been under
reconsideration:

- Based on user feedback (mostly around installation complexity and size), it
  has been considered to relax the new `pyarrow` requirement to not be a _hard_
  runtime dependency. In addition, NumPy 2.0 could in the future potentially
  reduce the need to make PyArrow a required dependency specifically for a
  dedicated pandas string dtype.
- The PDEP did not consider the usage of the experimental `pd.NA` as a
  consequence of adopting one of the existing implementations of the
  `StringDtype`.

For the second aspect, another variant of the `StringDtype` was
[introduced in pandas 2.1](https://pandas.pydata.org/docs/whatsnew/v2.1.0.html#whatsnew-210-enhancements-infer-strings)
that is still backed by PyArrow but follows the default missing values semantics
pandas uses for all other default data types (and using `NaN` as the missing
value sentinel) ([GH-54792](https://github.com/pandas-dev/pandas/issues/54792)).
At the time, the `storage` option for this new variant was called
`"pyarrow_numpy"` to disambiguate from the existing `"pyarrow"` option using `pd.NA`.

This last dtype variant is what you currently (pandas 2.2) get for string data
when enabling the ``future.infer_string`` option (to enable the behaviour which
is intended to become the default in pandas 3.0).

## Proposal

To be able to move forward with a string data type in pandas 3.0, this PDEP proposes:

1. For pandas 3.0, we enable a "string" dtype by default, which will use PyArrow
   if installed, and otherwise falls back to an in-house functionally-equivalent
   (but slower) version.
2. This default "string" dtype will follow the same behaviour for missing values
   as our other default data types, and use `NaN` as the missing value sentinel.
3. The version that is not backed by PyArrow can reuse the existing numpy
   object-dtype backed StringArray for its implementation.
4. We update installation guidelines to clearly encourage users to install
   pyarrow for the default user experience.

Those string dtypes enabled by default will then no longer be considered as
experimental.

### Default inference of a string dtype

By default, pandas will infer this new string dtype for string data (when
creating pandas objects, such as in constructors or IO functions).

The existing `future.infer_string` option can be used to opt-in to the future
default behaviour:

```python
>>> pd.options.future.infer_string = True
>>> pd.Series(["a", "b", None])
0      a
1      b
2    NaN
dtype: string
```

This option will be expanded to also work when PyArrow is not installed.

### Missing value semantics

Given that all other default data types use NaN semantics for missing values,
this proposal says that a new default string dtype should still use the same
default semantics. Further, it should result in default data types when doing
operations on the string column that result in a boolean or numeric data type
(e.g., methods like `.str.startswith(..)` or `.str.len(..)`, or comparison
operators like `==`, should result in default `int64` and `bool` data types).

Because the original `StringDtype` implementations already use `pd.NA` and
return masked integer and boolean arrays in operations, a new variant of the
existing dtypes that uses `NaN` and default data types is needed.

### Object-dtype "fallback" implementation

To avoid a hard dependency on PyArrow for pandas 3.0, this PDEP proposes to keep
a "fallback" option in case PyArrow is not installed. The original `StringDtype`
backed by a numpy object-dtype array of Python strings can be mostly reused for
this (adding a new variant of the dtype) and a new `StringArray` subclass only
needs minor changes to follow the above-mentioned missing value semantics
([GH-58451](https://github.com/pandas-dev/pandas/pull/58451)).

For pandas 3.0, this is the most realistic option given this implementation is
already available for a long time. Beyond 3.0, we can still explore further
improvements such as using NumPy 2.0 ([GH-58503](https://github.com/pandas-dev/pandas/issues/58503))
or nanoarrow ([GH-58552](https://github.com/pandas-dev/pandas/issues/58552)),
but at that point that is an implementation detail that should not have a
direct impact on users (except for performance).

### Naming

Given the long history of this topic, the naming of the dtypes is a difficult
topic.

In the first place, we need to acknowledge that most users should not need to
use storage-specific options. Users are expected to specify `pd.StringDtype()`
or `"string"`, and that will give them their default string dtype (which
depends on whether PyArrow is installed or not).

But for testing purposes and advanced use cases that want control over this, we
need some way to specify this and distinguish them from the other string dtypes.
Currently, the `StringDtype(storage="pyarrow_numpy")` is used, where
"pyarrow_numpy" is a rather confusing option.

TODO see if we can come up with a better naming scheme

## Alternatives

### Why not delay introducing a default string dtype?

To avoid introducing a new string dtype while other discussions and changes are
in flux (eventually making pyarrow a required dependency? adopting `pd.NA` as
the default missing value sentinel? using the new NumPy 2.0 capabilities?), we
could also delay introducing a default string dtype until there is more clarity
in those other discussions.

However:

1. Delaying has a cost: it further postpones introducing a dedicated string
   dtype that has massive benefits for our users, both in usability as (for the
   significant part of the user base that has PyArrow installed) in performance.
2. In case we eventually transition to use `pd.NA` as the default missing value
   sentinel, we will need a migration path for _all_ our data types, and thus
   the challenges around this will not be unique to the string dtype and
   therefore not a reason to delay this.

### Why not use the existing StringDtype with `pd.NA`?

Wouldn't adding even more variants of the string dtype will make things only more
confusing? Indeed, this proposal unfortunately introduces more variants of the
string dtype. However, the reason for this is to ensure the actual default user
experience is _less_ confusing, and the new string dtype fits better with the
other default data types.

If the new default string data type would use `pd.NA`, then after some
operations, a user can easily end up with a DataFrame that mixes columns using
`NaN` semantics and columns using `NA` semantics (and thus a DataFrame that
could have columns with two different int64, two different float64, two different
bool, etc dtypes). This would lead to a very confusing default experience.

With the proposed new variant of the StringDtype, this will ensure that for the
_default_ experience, a user will only see only 1 kind of integer dtype, only
kind of 1 bool dtype, etc. For now, a user should only get columns with an
`ArrowDtype` and/or using `pd.NA` when explicitly opting into this.

## Backward compatibility

The most visible backwards incompatible change will be that columns with string
data will no longer have an `object` dtype. Therefore, code that assumes
`object` dtype (such as `ser.dtype == object`) will need to be updated.

To allow testing your code in advance, the
`pd.options.future.infer_string = True` option is available.

Otherwise, the actual string-specific functionality (such as the `.str` accessor
methods) should all keep working as is. By preserving the current missing value
semantics, this proposal is also backwards compatible on this aspect.

One other backwards incompatible change is present for early adopters of the
existing `StringDtype`. In pandas 3.0, calling `pd.StringDtype()` will start
returning the new default string dtype, while up to now this returned the
experimental string dtype using `pd.NA` introduced in pandas 1.0. Those users
will need to start specifying a keyword in the dtype constructor if they want to
keep using `pd.NA` (but if they just want to have a dedicated string dtype, they
don't need to change their code).

## Timeline

The future PyArrow-backed string dtype was already made available behind a feature
flag in pandas 2.1 (by `pd.options.future.infer_string = True`).

Some small enhancements or fixes (or naming changes) might still be needed and
can be backported to pandas 2.2.x.

The variant using numpy object-dtype could potentially also be backported to
2.2.x to allow easier testing.

For pandas 3.0, this flag becomes enabled by default.


## PDEP-XX History

- 3 May 2024: Initial version
