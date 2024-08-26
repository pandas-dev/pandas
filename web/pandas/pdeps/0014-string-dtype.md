# PDEP-14: Dedicated string data type for pandas 3.0

- Created: May 3, 2024
- Status: Accepted
- Discussion: https://github.com/pandas-dev/pandas/pull/58551
- Author: [Joris Van den Bossche](https://github.com/jorisvandenbossche)
- Revision: 1

## Abstract

This PDEP proposes to introduce a dedicated string dtype that will be used by
default in pandas 3.0:

* In pandas 3.0, enable a string dtype (`"str"`) by default, using PyArrow if available
  or otherwise a string dtype using numpy object-dtype under the hood as fallback.
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
However, up to this date, pandas has not yet taken the step to use `pd.NA` for
for any default dtype, and thus the `StringDtype` deviates in missing value
behaviour compared to the default data types.

In 2023, [PDEP-10](https://pandas.pydata.org/pdeps/0010-required-pyarrow-dependency.html)
proposed to start using a PyArrow-backed string dtype by default in pandas 3.0
(i.e. infer this type for string data instead of object dtype). To ensure we
could use the variant of `StringDtype` backed by PyArrow instead of Python
objects (for better performance), it proposed to make `pyarrow` a new required
runtime dependency of pandas.

In the meantime, NumPy has also been working on a native variable-width string
data type, which was made available [starting with NumPy
2.0](https://numpy.org/devdocs/release/2.0.0-notes.html#stringdtype-has-been-added-to-numpy).
This can provide a potential alternative to PyArrow for implementing a string
data type in pandas that is not backed by Python objects.

After acceptance of PDEP-10, two aspects of the proposal have been under
reconsideration:

- Based on feedback from users and maintainers from other packages (mostly
  around installation complexity and size), it has been considered to relax the
  new `pyarrow` requirement to not be a _hard_ runtime dependency. In addition,
  NumPy 2.0 could in the future potentially reduce the need to make PyArrow a
  required dependency specifically for a dedicated pandas string dtype.
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

1. For pandas 3.0, a `"str"` string dtype is enabled by default, i.e. this
   string dtype will be used as the default dtype for text data when creating
   pandas objects (e.g. inference in constructors, I/O functions).
2. This default string dtype will follow the same behaviour for missing values
   as other default data types, and use `NaN` as the missing value sentinel.
3. The string dtype will use PyArrow if installed, and otherwise falls back to
   an in-house functionally-equivalent (but slower) version. This fallback can
   reuse (with minor code additions) the existing numpy object-dtype backed
   StringArray for its implementation.
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

As mentioned in the background section, the original `StringDtype` has always
used the experimental `pd.NA` sentinel for missing values. In addition to using
`pd.NA` as the scalar for a missing value, this essentially means that:

- String columns follow ["NA-semantics"](https://pandas.pydata.org/docs/user_guide/missing_data.html#na-semantics)
  for missing values, where `NA` propagates in boolean operations such as
  comparisons or predicates.
- Operations on the string column that give a numeric or boolean result use the
  nullable Integer/Float/Boolean data types (e.g. `ser.str.len()` returns the
  nullable `"Int64"` / `pd.Int64Dtype()` dtype instead of the numpy `int64`
  dtype (or `float64` in case of missing values)).

However, up to this date, all other default data types still use `NaN` semantics
for missing values. Therefore, this proposal says that a new default string
dtype should also still use the same default missing value semantics and return
default data types when doing operations on the string column, to be consistent
with the other default dtypes at this point.

In practice, this means that the default string dtype will use `NaN` as
the missing value sentinel, and:

- String columns will follow NaN-semantics for missing values, where `NaN` gives
  False in boolean operations such as comparisons or predicates.
- Operations on the string column that give a numeric or boolean result will use
  the default data types (i.e. numpy `int64`/`float64`/`bool`).

Because the original `StringDtype` implementations already use `pd.NA` and
return masked integer and boolean arrays in operations, a new variant of the
existing dtypes that uses `NaN` and default data types was needed. The original
variant of `StringDtype` using `pd.NA` will continue to be available for those
who were already using it.

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

For the original variant of `StringDtype` using `pd.NA`, currently the default
storage is `"python"` (the object-dtype based implementation). Also for this
variant, it is proposed to follow the same logic for determining the default
storage, i.e. default to `"pyarrow"` if available, and otherwise
fall back to `"python"`.

### Naming

Given the long history of this topic, the naming of the dtypes is a difficult
topic.

In the first place, it should be acknowledged that most users should not need to
use storage-specific options. Users are expected to specify a generic name (such
as `"str"` or `"string"`), and that will give them their default string dtype
(which depends on whether PyArrow is installed or not).

For the generic string alias to specify the dtype, `"string"` is already used
for the `StringDtype` using `pd.NA`. This PDEP proposes to use `"str"` for the
new default `StringDtype` using `NaN`. This ensures backwards compatibility for
code using `dtype="string"`, and was also chosen because `dtype="str"` or
`dtype=str` currently already works to ensure your data is converted to
strings (only using object dtype for the result).

But for testing purposes and advanced use cases that want control over the exact
variant of the `StringDtype`, we need some way to specify this and distinguish
them from the other string dtypes.

Currently (pandas 2.2), `StringDtype(storage="pyarrow_numpy")` is used for the new variant using `NaN`,
where the `"pyarrow_numpy"` storage was used to disambiguate from the existing
`"pyarrow"` option using `pd.NA`. However, `"pyarrow_numpy"` is a rather confusing
option and doesn't generalize well. Therefore, this PDEP proposes a new naming
scheme as outlined below, and `"pyarrow_numpy"` will be deprecated as an alias
in pandas 2.3 and removed in pandas 3.0.

The `storage` keyword of `StringDtype` is kept to disambiguate the underlying
storage of the string data (using pyarrow or python objects), but an additional
`na_value` is introduced to disambiguate the the variants using NA semantics
and NaN semantics.

Overview of the different ways to specify a dtype and the resulting concrete
dtype of the data:

| User specification                          | Concrete dtype                                                | String alias                          | Note     |
|---------------------------------------------|---------------------------------------------------------------|---------------------------------------|----------|
| Unspecified (inference)                     | `StringDtype(storage="pyarrow"\|"python", na_value=np.nan)`   | "str"                                 | (1)      |
| `"str"` or `StringDtype(na_value=np.nan)`   | `StringDtype(storage="pyarrow"\|"python", na_value=np.nan)`   | "str"                                 | (1)      |
| `StringDtype("pyarrow", na_value=np.nan)`   | `StringDtype(storage="pyarrow", na_value=np.nan)`             | "str"                                 |          |
| `StringDtype("python", na_value=np.nan)`    | `StringDtype(storage="python", na_value=np.nan)`              | "str"                                 |          |
| `StringDtype("pyarrow")`                    | `StringDtype(storage="pyarrow", na_value=pd.NA)`              | "string[pyarrow]"                     |          |
| `StringDtype("python")`                     | `StringDtype(storage="python", na_value=pd.NA)`               | "string[python]"                      |          |
| `"string"` or `StringDtype()`               | `StringDtype(storage="pyarrow"\|"python", na_value=pd.NA)`    | "string[pyarrow]" or "string[python]" | (1)      |
| `StringDtype("pyarrow_numpy")`              | `StringDtype(storage="pyarrow", na_value=np.nan)`             | "string[pyarrow_numpy]"               | (2)      |

Notes:

- (1) You get "pyarrow" or "python" depending on pyarrow being installed.
- (2) "pyarrow_numpy" is kept temporarily because this is already in a released
  version, but it will be deprecated in 2.x and removed for 3.0.

For the new default string dtype, only the `"str"` alias can be used to
specify the dtype as a string, i.e. pandas would not provide a way to make the
underlying storage (pyarrow or python) explicit through the string alias. This
string alias is only a convenience shortcut and for most users `"str"` is
sufficient (they don't need to specify the storage), and the explicit
`pd.StringDtype(storage=..., na_value=np.nan)` is still available for more
fine-grained control.

Also for the existing variant using `pd.NA`, specifying the storage through the
string alias could be deprecated, but that is left for a separate decision.

## Alternatives

### Why not delay introducing a default string dtype?

To avoid introducing a new string dtype while other discussions and changes are
in flux (eventually making pyarrow a required dependency? adopting `pd.NA` as
the default missing value sentinel? using the new NumPy 2.0 capabilities?
overhauling all our dtypes to use a logical data type system?), introducing a
default string dtype could also be delayed until there is more clarity in those
other discussions. Specifically, it would avoid temporarily switching to use
`NaN` for the string dtype, while in a future version we might switch back
to `pd.NA` by default.

However:

1. Delaying has a cost: it further postpones introducing a dedicated string
   dtype that has significant benefits for users, both in usability as (for the
   part of the user base that has PyArrow installed) in performance.
2. In case pandas eventually transitions to use `pd.NA` as the default missing value
   sentinel, a migration path for _all_ pandas data types will be needed, and thus
   the challenges around this will not be unique to the string dtype and
   therefore not a reason to delay this.

Making this change now for 3.0 will benefit the majority of users, and the PDEP
author believes this is worth the cost of the added complexity around "yet
another dtype" (also for other data types we already have multiple variants).

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

An initial version of this PDEP proposed to use the `"string"` alias and the
default `pd.StringDtype()` class constructor for the new default dtype.
However, that caused a lot of discussion around backwards compatibility for
existing users of `dtype=pd.StringDtype()` and `dtype="string"`, that uses
`pd.NA` to represent missing values.

During the discussion, several alternatives have been brought up. Both
alternative keyword names as using a different constructor. In the end,
this PDEP proposes to use a different string alias (`"str"`) but to keep
using the existing `pd.StringDtype` (with the existing `storage` keyword but
with an additional `na_value` keyword) for now to keep the changes as
minimal as possible, leaving a larger overhaul of the dtype system (potentially
including different constructor functions or namespace) for a future discussion.
See [GH-58613](https://github.com/pandas-dev/pandas/issues/58613) for the full
discussion.

One consequence is that when using the class constructor for the default dtype,
it has to be used with non-default arguments, i.e. a user needs to specify
`pd.StringDtype(na_value=np.nan)` to get the default dtype using `NaN`.
Therefore, the pandas documentation will focus on the usage of `dtype="str"`.

## Backward compatibility

The most visible backwards incompatible change will be that columns with string
data will no longer have an `object` dtype. Therefore, code that assumes
`object` dtype (such as `ser.dtype == object`) will need to be updated. This
change is done as a hard break in a major release, as warning in advance for the
changed inference is deemed too noisy.

To allow testing code in advance, the
`pd.options.future.infer_string = True` option is available for users.

Otherwise, the actual string-specific functionality (such as the `.str` accessor
methods) should generally all keep working as is.

By preserving the current missing value semantics, this proposal is also mostly
backwards compatible on this aspect. When storing strings in object dtype, pandas
however did allow using `None` as the missing value indicator as well (and in
certain cases such as the `shift` method, pandas even introduced this itself).
For all the cases where currently `None` was used as the missing value sentinel,
this will change to consistently use `NaN`.

### For existing users of `StringDtype`

Existing code that already opted in to use the `StringDtype` using `pd.NA`
should generally keep working as is. The latest version of this PDEP preserves
the behaviour of `dtype="string"` or `dtype=pd.StringDtype()` to mean the
`pd.NA` variant of the dtype.

It does propose the change the default storage to `"pyarrow"` (if available) for
the opt-in `pd.NA` variant as well, but this should have limited, if any,
user-visible impact.

## Timeline

The future PyArrow-backed string dtype was already made available behind a feature
flag in pandas 2.1 (enabled by `pd.options.future.infer_string = True`).

The variant using numpy object-dtype can also be backported to the 2.2.x branch
to allow easier testing. It is proposed to release this as 2.3.0 (created from
the 2.2.x branch, given that the main branch already includes many other changes
targeted for 3.0), together with the changes to the naming scheme.

The 2.3.0 release would then have all future string functionality available
(both the pyarrow and object-dtype based variants of the default string dtype).

For pandas 3.0, this `future.infer_string` flag becomes enabled by default.

## PDEP-14 History

- 3 May 2024: Initial version
