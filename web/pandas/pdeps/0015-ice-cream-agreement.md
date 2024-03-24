# PDEP-15: Ice Cream Agreement

- Created: March 2024
- Status: Under discussion
- Discussion: [#32265](https://github.com/pandas-dev/pandas/issues/32265)
- Author: [Patrick Hoefler](https://github.com/phofl)
          [Joris Van den Bossche](https://github.com/jorisvandenbossche)
- Revision: 1

## Abstract

Short summary of the proposal:

1. The pandas Extension Array interface will fully support 2D arrays. Currently, pandas publicly only
   supports 1D Extension Arrays. Additionally, pandas will make all internal NumPy-based
   Extension Arrays 2D. This specifically includes our nullable Extension Arrays and
   __excludes__ Arrow-based extension arrays. Consequently, pandas will move to the nullable
   extension dtypes by default to provide consistent missing value handling across all dtypes.

2. The NumPy based Extension Arrays will exclusively use ``pd.NA`` as a missing value indicator.
   ``np.nan`` will not allowed to be present, which removes the need to distinguish between
   ``pd.NA`` and ``np.nan``. The ``FloatingArray`` will thus only use ``NA`` and not ``nan``.

This addresses several issues that have been open for years:

1) Clear and consistent missing value handling across all dtypes.
2) A resolution of the discussion how to treat ``NA`` and ``NaN`` in FloatingArrays.
3) Making the NA-scalar easier to use through no longer raising on ``bool(pd.NA)``.
4) The ExtensionArray interface will be a first class citizen, which simplifies 3rd-party
   extensions.

## Background

pandas currently maintains three different sets of dtypes next to each other:

- NumPy dtypes that use NumPy arrays to store the data
- Arrow dtypes that use PyArrow Arrays to store the data
- Nullable extension dtypes that use pandas Extension Arrays to store the data. These
  arrays add a layer on top of NumPy to modify the behavior.

The NumPy dtypes are currently default and the most widely used. They use NaN as the missing
value indicator, which is a float and can't be stored in an integer or boolean array. Consequently,
these dtypes are cast to float/object if a missing value is inserted into them.

The nullable extension dtypes were originally designed to solve these problems and to provide
consistent missing value behavior between different dtypes. These arrays use a strict 1D layout
and store missing values through an accompanying mask. The integer and boolean dtypes are
supported well across the pandas API, but the float dtypes still have many inconsistencies
with respect to missing value handling and the behavior of ``pd.NA`` and ``np.nan``. The
nullable arrays generally are hindered in some scenarios because of the 1D layout (``axis=1``
operations, transposing, etc.).

The Arrow dtypes are the most recent addition to pandas. They are currently separate from the
other two sets of dtypes since they user a different data model under the hood and are strictly
1D.

## Proposal

This proposal aims to unify the missing value handling across all dtypes and to resolve
outstanding issues for the FloatingArray implementation. This proposal is not meant to
address implementation details, rather to provide a high level way forward.

1. The ``FloatingArray`` implementation will exclusively use ``pd.NA`` was missing value
   indicator. ``np.nan`` will not be allowed to be present in the array. The missing value
   behavior will follow the semantics of the other nullable extension dtypes.

2. The ExtensionArray interface will be extended to support 2D arrays. This will allow
   us to make our internal nullable ExtensionArrays 2D and also make this option available
   to 3rd party arrays.

3. pandas will move to nullable extension arrays by default instead of using the NumPy
   dtypes that are currently the default. Every constructor and IO method will infer
   extension dtypes by default if not explicitly specified by the user. This is
   similar to the current ``dtype_backend="numpy_nullable"`` keyword in IO methods,
   but will be made the new default and extended to the constructors.

We will obey the following dtype mapping:

- int*/uint* -> Int*/UInt*
- float* -> Float*
- bool -> boolean
- object dtype will be mapped to string, but this is covered by PDEP10
- object dtype will be used for values that aren't strings

This will ensure that all dtypes have consistent missing value handling and there
is no need to upcast if a missing value is inserted into integers or booleans. Those
nullability semantics will be mostly consistent with how PyArrow treats nulls and thus
make switching between both set of dtypes easier.  Additionally, it allows the usage of
other Arrow dtypes by default that user the same semantics (bytes, nested dtypes, ...).

This proposal formalizes the results of the pandas core sprint in 2023.

## Backward compatibility

Making Extension Arrays 2D can be considered an implementation detail and shouldn't
impact users negatively.

The ``FloatingArray`` implementation is still experimental and currently riddled with
bugs with respect to handling of ``pd.NA`` and ``np.nan``. It's experimental status allows
us to change this without worrying to much about backwards compatibility. Additionally,
because of the bugs related to NA handling makes it unlikely that it is used in serious
applications.

Switching to nullable dtypes by default will be a huge change for pandas. It will deviate
from the current NumPy dtypes and change nullability semantics for users. This will require
care when implementing this change to make the change in behavior as small as possible and
to ensure that the new implementation is well tested and easy to opt in for users before
we make this switch.

## Considerations

### 2D Extension Arrays

The current restriction of 1D Extension Arrays only has a number of limitations internally.
``axis=1`` operations and more generally operations that transpose the data in some way
tend to fall back to object. Additionally, the 1D limitation requires copies when converting
between NumPy and pandas in all cases for DataFrames. Our internal algorithms are more
performant for 2D arrays like groupby aggregations. There are currently 35
TODOs across the code base that are related to 2D extension arrays.

I am not aware of any drawbacks compared to the current default dtypes at the point
of writing.

### FloatingArray

The FloatingArray implementation is currently experimental and has a number of bugs.
The main source of issues stems from the fact that both ``np.nan`` and ``pd.NA``
are allowed and not properly handled.

**Status quo**

When constructing a FloatingArray from a NumPy array, a Series with a
NumPy dtype or another list-like the constructor converts ``np.nan`` to ``pd.NA``.

```python

In [3]: pd.array(np.array([1.5, np.nan]), dtype="Float64")
Out[3]:
<FloatingArray>
[1.5, <NA>]
Length: 2, dtype: Float64
```

This is done because NumPy doesn't have a missing value sentinel and pandas
considers ``np.nan`` to be missing.

Inserting ``np.nan`` into a FloatingArray will also coerce to ``pd.NA``.

```python
In [4]: arr = pd.array([1.5, np.nan], dtype="Float64")
In [5]: arr[0] = np.nan

In [6]: arr
Out[6]:
<FloatingArray>
[<NA>, <NA>]
Length: 2, dtype: Float64
```

You can introduce NaN values through ``0/0`` for example, but having NaN
causes other issues. None of our na-detection methods (fillna, isna, ...)
will match NaN values, they only match ``pd.NA``. A non exhaustive list of
issues this behavior causes can be found on the
[pandas issue tracker](https://github.com/pandas-dev/pandas/issues?q=is%3Aopen+is%3Aissue+label%3A%22Ice+Cream+Agreement%22).

**Solution**

The current state makes the FloatingArray unusable if you rely on missing values
in any way. We solve this problem through disallowing ``np.nan`` in the FloatingArray.
Only ``NA`` will be allowed to be present.

- This solution makes the implementation of all methods that interact with NA
  simpler and more consistent. This includes methods like ``fillna`` but also
  sorting operations.
- Users are used to only having ``np.nan`` as a missing value indicator in pandas.
  Staying with one missing value indicator in ``pd.NA`` will make the behavior
  less confusing for users.
- In- and Output to and from NumPy is unambiguous. Every NaN is converted to NA
  and back.


**Drawbacks**

- There is no option to distinguish between missing and invalid values. This is currently not
  possible either and generally would require increasing the API surface to handle both cases.
  Methods interacting with missing values would need to be configurable. There was never much
  demand for this feature, so the additional complexity does not seem justified.

Distinguishing NA and NaN adds a lot of complexities:

- Roundtripping through NumPy is not really possible. Currently, we are converting to NA and then
  converting back. This is potentially very confusing if users end up with the same value after
  going to NumPy.
- Increases the API surface and complexity for all methods that interact with NA values. Additionally,
  it also concerns all methods that have the ``skipna`` argument.
- Having both values in pandas is confusing for users. Historically, pandas used only NaN.
  Differentiating between NA and NaN would make behavior less intuitive for non-expert users.
- It adds maintenance burden.
- NA and NaN have different semantics in comparison operations, which adds further mental complexity.

## Timeline

Make Extension Arrays 2D and fix all inconsistencies in the FloatingArray. Ensure that
this is done by the time that pandas 4.0 is released and then prepare the migration
to nullable dtypes by default in the next major release.

### PDEP History

- March 2024: Initial draft

Note: There is a very long discussion in [GH-32265](https://github.com/pandas-dev/pandas/issues/32265)
that concerns this topic.
