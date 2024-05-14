# PDEP-14: Dedicated string data type for pandas 3.0

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
using NumPy 2.0 strings, etc).

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
default for any dtype, and thus the `StringDtype` deviates in missing value behaviour compared
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
- PDEP-10 did not consider the usage of the experimental `pd.NA` as a
  consequence of adopting one of the existing implementations of the
  `StringDtype`.

For the second aspect, another variant of the `StringDtype` was
[introduced in pandas 2.1](https://pandas.pydata.org/docs/whatsnew/v2.1.0.html#whatsnew-210-enhancements-infer-strings)
that is still backed by PyArrow but follows the default missing values semantics
pandas uses for all other default data types (and using `NaN` as the missing
value sentinel) ([GH-54792](https://github.com/pandas-dev/pandas/issues/54792)).
At the time, the `storage` option for this new variant was called
`"pyarrow_numpy"` to disambiguate from the existing `"pyarrow"` option using
`pd.NA` (but this PDEP proposes a better naming scheme, see the "Naming"
subsection below).

This last dtype variant is what users currently (pandas 2.2) get for string data
when enabling the ``future.infer_string`` option (to enable the behaviour which
is intended to become the default in pandas 3.0).

## Proposal

To be able to move forward with a string data type in pandas 3.0, this PDEP proposes:

1. For pandas 3.0,  a "string" dtype is enabled by default, which will use PyArrow
   if installed, and otherwise falls back to an in-house functionally-equivalent
   (but slower) version.
2. This default "string" dtype will follow the same behaviour for missing values
   as other default data types, and use `NaN` as the missing value sentinel.
3. The version that is not backed by PyArrow can reuse (with minor code
   additions) the existing numpy object-dtype backed StringArray for its
   implementation.
4. Installation guidelines are updated to clearly encourage users to install
   pyarrow for the default user experience.

Those string dtypes enabled by default will then no longer be considered as
experimental.

### Default inference of a string dtype

By default, pandas will infer this new string dtype instead of object dtype for
string data (when creating pandas objects, such as in constructors or IO
functions).

In pandas 2.2, the existing `future.infer_string` option can be used to opt-in to the future
default behaviour:

```python
>>> pd.options.future.infer_string = True
>>> pd.Series(["a", "b", None])
0      a
1      b
2    NaN
dtype: string
```

Right now (pandas 2.2), the existing option only enables the PyArrow-based
future dtype. For the remaining 2.x releases, this option will be expanded to
also work when PyArrow is not installed to enable the object-dtype fallback in
that case.

### Missing value semantics

As mentioned in the background section, the original `StringDtype` has used
the experimental `pd.NA` sentinel for missing values. In addition to using
`pd.NA` as the scalar for a missing value, this essentially means
that:

- String columns follow ["NA-semantics"](https://pandas.pydata.org/docs/user_guide/missing_data.html#na-semantics)
  for missing values, where `NA` propagates in boolean operations such as
  comparisons or predicates.
- Operations on the string column that give a numeric or boolean result use the
  nullable Integer/Float/Boolean data types (e.g. `ser.str.len()` returns the
  nullable `'Int64"` / `pd.Int64Dtype()` dtype instead of the numpy `int64`
  dtype (or `float64` in case of missing values)).

However, up to this date, all other default data types still use `NaN` semantics
for missing values. Therefore, this proposal says that a new default string
dtype should also still use the same default missing value semantics and return
default data types when doing operations on the string column, to be consistent
with the other default dtypes at this point.

In practice, this means that the default `"string"` dtype will use `NaN` as
the missing value sentinel, and:

- String columns will follow NaN-semantics for missing values, where `NaN` gives
  False in boolean operations such as comparisons or predicates.
- Operations on the string column that give a numeric or boolean result will use
  the default data types (i.e. numpy `int64`/`float64`/`bool`).

Because the original `StringDtype` implementations already use `pd.NA` and
return masked integer and boolean arrays in operations, a new variant of the
existing dtypes that uses `NaN` and default data types was needed. The original
variant of `StringDtype` using `pd.NA` will still be available for those who
want to keep using it (see below in the "Naming" subsection for how to specify
this).

### Object-dtype "fallback" implementation

To avoid a hard dependency on PyArrow for pandas 3.0, this PDEP proposes to keep
a "fallback" option in case PyArrow is not installed. The original `StringDtype`
backed by a numpy object-dtype array of Python strings can be mostly reused for
this (adding a new variant of the dtype) and a new `StringArray` subclass only
needs minor changes to follow the above-mentioned missing value semantics
([GH-58451](https://github.com/pandas-dev/pandas/pull/58451)).

For pandas 3.0, this is the most realistic option given this implementation has
already been available for a long time. Beyond 3.0, further improvements such as
using NumPy 2.0 ([GH-58503](https://github.com/pandas-dev/pandas/issues/58503))
or nanoarrow ([GH-58552](https://github.com/pandas-dev/pandas/issues/58552)) can
still be explored, but at that point that is an implementation detail that
should not have a direct impact on users (except for performance).

### Naming

Given the long history of this topic, the naming of the dtypes is a difficult
topic.

In the first place, it should be acknowledged that most users should not need to
use storage-specific options. Users are expected to specify `pd.StringDtype()`
or `"string"`, and that will give them their default string dtype (which
depends on whether PyArrow is installed or not).

But for testing purposes and advanced use cases that want control over this, we
need some way to specify this and distinguish them from the other string dtypes.
In addition, users that want to continue using the original NA-variant of the
dtype need a way to specify this.

Currently (pandas 2.2), `StringDtype(storage="pyarrow_numpy")` is used, where
the `"pyarrow_numpy"` storage was used to disambiguate from the existing
`"pyarrow"` option using `pd.NA`. However, `"pyarrow_numpy"` is a rather confusing
option and doesn't generalize well. Therefore, this PDEP proposes a new naming
scheme as outlined below, and `"pyarrow_numpy"` will be deprecated and removed
before pandas 3.0.

The `storage` keyword of `StringDtype` is kept to disambiguate the underlying
storage of the string data (using pyarrow or python objects), but an additional
`na_value` is introduced to disambiguate the the variants using NA semantics
and NaN semantics.

Overview of the different ways to specify a dtype and the resulting concrete
dtype of the data:

| User specification                       | Concrete dtype                                                | String alias                          | Note     |
|------------------------------------------|---------------------------------------------------------------|---------------------------------------|----------|
| Unspecified (inference)                  | `StringDtype(storage="pyarrow"\|"python", na_value=np.nan)`   | "string"                              | (1)      |
| `StringDtype()` or `"string"`            | `StringDtype(storage="pyarrow" \| "python", na_value=np.nan)` | "string"                              | (1), (2) |
| `StringDtype("pyarrow")`                 | `StringDtype(storage="pyarrow", na_value=np.nan)`             | "string"                              | (2)      |
| `StringDtype("python")`                  | `StringDtype(storage="python", na_value=np.nan)`              | "string"                              | (2)      |
| `StringDtype("pyarrow", na_value=pd.NA)` | `StringDtype(storage="pyarrow", na_value=pd.NA)`              | "string[pyarrow]"                     |          |
| `StringDtype("python", na_value=pd.NA)`  | `StringDtype(storage="python", na_value=pd.NA)`               | "string[python]"                      |          |
| `StringDtype(na_value=pd.NA)`            | `StringDtype(storage="pyarrow" \| "python", na_value=pd.NA)`  | "string[pyarrow]" or "string[python]" | (1)      |
| `StringDtype("pyarrow_numpy")`           | `StringDtype(storage="pyarrow", na_value=np.nan)`             | "string[pyarrow_numpy]"               | (3)      |

Notes:

- (1) You get "pyarrow" or "python" depending on pyarrow being installed.
- (2) Those three rows are backwards incompatible (i.e. they work now but give
  the NA-variant), see the "Backward compatibility" section below.
- (3) "pyarrow_numpy" is kept temporarily because this is already in a released
  version, but we can deprecate it in 2.x and have it removed for 3.0.

For the new default string dtype, only the `"string"` alias can be used to
specify the dtype as a string, i.e. a way would not be provided to make the
underlying storage (pyarrow or python) explicit through the string alias. This
string alias is only a convenience shortcut and for most users `"string"` is
sufficient (they don't need to specify the storage), and the explicit
`pd.StringDtype(...)` is still available for more fine-grained control.

## Alternatives

### Why not delay introducing a default string dtype?

To avoid introducing a new string dtype while other discussions and changes are
in flux (eventually making pyarrow a required dependency? adopting `pd.NA` as
the default missing value sentinel? using the new NumPy 2.0 capabilities?
overhauling all our dtypes to use a logical data type system?), introducing a
default string dtype could also be delayed until there is more clarity in those
other discussions.

However:

1. Delaying has a cost: it further postpones introducing a dedicated string
   dtype that has massive benefits for users, both in usability as (for the
   significant part of the user base that has PyArrow installed) in performance.
2. In case pandas eventually transitions to use `pd.NA` as the default missing value
   sentinel,  a migration path for _all_ pandas data types will be needed, and thus
   the challenges around this will not be unique to the string dtype and
   therefore not a reason to delay this.

Making this change now for 3.0 will benefit the majority of users, while coming
at a cost for a part of the users who already started using the `"string"` or
`pd.StringDtype()` dtype (they will have to update their code to continue to the
variant using `pd.NA`, see the "Backward compatibility" section below).

### Why not use the existing StringDtype with `pd.NA`?

Wouldn't adding even more variants of the string dtype make things only more
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
kind of 1 bool dtype, etc. For now, a user should only get columns using `pd.NA`
when explicitly opting into this.

### Naming alternatives

This PDEP now keeps the `pd.StringDtype` class constructor with the existing
`storage` keyword and with an additional `na_value` keyword.

During the discussion, several alternatives have been brought up. Both
alternative keyword names as using a different constructor. This PDEP opted to
keep using the existing `pd.StringDtype()` for now to keep the changes as
minimal as possible, leaving a larger overhaul of the dtype system (potentially
including different constructor functions or namespace) for a future discussion.
See [GH-58613](https://github.com/pandas-dev/pandas/issues/58613) for the full
discussion.

## Backward compatibility

The most visible backwards incompatible change will be that columns with string
data will no longer have an `object` dtype. Therefore, code that assumes
`object` dtype (such as `ser.dtype == object`) will need to be updated. This
change is done as a hard break in a major release, as warning in advance for the
changed inference is deemed too noisy.

To allow testing code in advance, the
`pd.options.future.infer_string = True` option is available for users.

Otherwise, the actual string-specific functionality (such as the `.str` accessor
methods) should generally all keep working as is. By preserving the current
missing value semantics, this proposal is also backwards compatible on this
aspect.

### For existing users of `StringDtype`

Users of the existing `StringDtype` will see more backwards incompatible
changes, though. In pandas 3.0, calling `pd.StringDtype()` (or specifying
`dtype="string"`) will start returning the new default string dtype using `NaN`,
while up to now this returned the string dtype using `pd.NA` introduced in
pandas 1.0.

For example, this code snippet returned the NA-variant of `StringDtype` with
pandas 1.x and 2.x:

```python
>>> pd.Series(["a", "b", None], dtype="string")
0      a
1      b
2   <NA>
dtype: string
```

but will start returning the new default NaN-variant of `StringDtype` with
pandas 3.0. This means that the missing value sentinel will change from `pd.NA`
to `NaN`, and that operations will no longer return nullable dtypes but default
numpy dtypes (see the "Missing value semantics" section above).

While this change will be transparent in many cases (e.g. checking for missing
values with `isna()`/`dropna()`/`fillna()` or filtering rows with the result of
a string predicate method keeps working regardless of the sentinel), this can be
a breaking change if users relied on the exact sentinel or resulting dtype. Since
pandas 1.0, the string dtype has been promoted quite a bit, and so we expect
that many users already have started using this dtype, even though officially
still labeled as "experimental".

To smooth the upgrade experience for those users, it is proposed to add a
deprecation warning before 3.0 when such dtype is created, giving them two
options:

- If the user just wants to have a dedicated "string" dtype (or the better
  performance when using pyarrow) but is fine with using the default NaN
  semantics, they can add `pd.options.future.infer_string = True` to their code
  to suppress the warning and already opt-in to the future behaviour of pandas
  3.0.
- If the user specifically wants the variant of the string dtype that uses
  `pd.NA` (and returns nullable numeric/boolean dtypes in operations), they will
  have to update their dtype specification from `"string"` / `pd.StringDtype()`
  to `pd.StringDtype(na_value=pd.NA)` to suppress the warning and further keep
  their code running as is.

## Timeline

The future PyArrow-backed string dtype was already made available behind a feature
flag in pandas 2.1 (enabled by `pd.options.future.infer_string = True`).

Some small enhancements or fixes might still be needed and can continue to be
backported to pandas 2.2.x.

The variant using numpy object-dtype can also be backported to the 2.2.x branch
to allow easier testing. It is proposed to release this as 2.3.0 (created from
the 2.2.x branch, given that the main branch already includes many other changes
targeted for 3.0), together with the deprecation warning when creating a dtype
from `"string"` / `pd.StringDtype()`.

The 2.3.0 release would then have all future string functionality available
(both the pyarrow and object-dtype based variants of the default string dtype),
and warn existing users of the `StringDtype` in advance of 3.0 about how to
update their code.

For pandas 3.0, this `future.infer_string` flag becomes enabled by default.

## PDEP-XX History

- 3 May 2024: Initial version
