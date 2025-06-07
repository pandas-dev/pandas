# PDEP-18: Nullable Object Dtype for Pandas

- Created: 07 June 2025
- Status: Draft
- Discussion: [#32931](https://github.com/pandas-dev/pandas/issues/32931)
- Author: [Simon Hawkins](https://github.com/simonjayhawkins)
- Revision: 1

## Abstract

This proposal outlines the introduction of a nullable object dtype to the pandas library. The goal is to provide a dedicated dtype for handling arbitrary Python objects with consistent missing value semantics using `pd.NA`. Unlike the traditional `object` dtype which lacks robust missing data handling, this new nullable dtype will add clarity and consistency in representing missing or undefined values within object arrays.

## Motivation

Currently, the `object` dtype in pandas is a catch-all for heterogeneous Python objects, but it does not enforce any particular missing-value semantics. As pandas has evolved to include extension types (like `string[python]`, `Int64`, or `boolean`), there is a clear benefit in extending these improvements to the object datatype. A nullable object dtype would help:
- **Consistency**: Enforce a uniform approach to managing missing values with `pd.NA` across all dtypes.
- **Interoperability**: Enable cleaner and more predictable behavior when performing operations on data previously stored as generic objects.
- **Clarity**: Help users distinguish between truly “object” data and data that is better represented by a nullable container supporting missing values.

This proposal is driven by frequent community discussions and development efforts that aim to unify missing value handling across pandas data types.

## Detailed Proposal

### Definition

The proposal introduces a new extension type, tentatively named `"object_nullable"`, that stores an underlying array of Python objects alongside a boolean mask that indicates missing (i.e., `pd.NA`) values. The API should mimic that of existing extension arrays, ensuring that missing value propagation, casting, and arithmetic comparisons (where applicable) behave consistently with other nullable types.

### Key Features
1. **Consistent Missing Value Semantics**:
    - Missing entries will be represented by `pd.NA`, ensuring compatibility with pandas nullable dtypes that use `pd.NA` as the missing value indicator as well as the experimental `ArrowDType`.
    - Operations that encounter missing values will handle `pd.NA` uniformly consistent with other pandas nullable dtypes that use `pd.NA` as the missing value indicator.
2. **Underlying Data Storage**:
    - The core data structure will consist of a NumPy array of Python objects and an associated boolean mask. (not so different from the current `object` backed nullable string array variant that uses `pd.NA` as the missing value.)
    - Consideration should be given to performance, ensuring that operations remain as vectorized as possible despite the inherent overhead of handling Python objects.
3. **API Integration**:
    - The new dtype will implement the ExtensionArray interface.
    - Methods such as `astype`, `isna`, `fillna`, and element-wise operations are already defined to respect missing values in the other pandas nullable dtypes.
    - All operations on a nullable object array will return a pandas nullable array except where requested, such as `astype`. Methods like `fillna` would still return a nullable object array even though there are no missing values to avoid introducing mixed-propagation behavior.
    - Ensure compatibility with pandas functions, like groupby, concatenation, and merging, where the semantics of missing values are critical.
4. **Transition and Interoperability**:
    - Users should be able to convert from the legacy object dtype to object_nullable using a constructor or an explicit method (e.g., `pd.array(old_array, dtype="object_nullable")`) using the exisiting api.
    - Operations on existing pandas nullable dtypes that would normally produce an object dtype should be updated (or made configurable as a transition path) to yield "object_nullable" in all cases even when missing values are not present to avoid introducing mixed-propagation behavior.
    - `ArrowDType` does not offer an `object` dtype for hetrogeneous Python objects and therefore a user requesting arrow dtypes could be given "object_nullable" arrays where appropriate to avoid mixed `pd.NA`/`np.nan` semantics when using `dtype_backend="pyarrow"`.


### Implementation Considerations
1. **Performance**:
    - Handling arbitrary Python objects is inherently slower than operations on native numerical types.
    - Expanding the EA interface to 2D is outside the scope of this PDEP.

2. **Backward Compatibility**:
    - Existing code that uses the traditional object dtype should not break. (Making the pandas nullable object dtype the default is not part of this proposal and would be discussed in conjunction with moving the other pandas nullable dtypes to be default.)
    - Existing code that uses the pandas nullable dtypes should not break without warnings, even though they are considered experimental, as these dtypes have been available to users for a long time. The new dtype can be offered as an opt-in feature initially.
3. **Testing and Documentation**:
    - Extensive tests will be required to validate behavior against edge cases.
    - Updated documentation should explain differences between the legacy object dtype and object_nullable, including examples and migration tips.
4. **Community Feedback**:
    - Continuous discussions on GitHub, mailing lists, and related channels will inform refinements. The nullable object dtype should be available as opt-in for at least 2 minor versions to allow sufficient time for feedback before the return types of the existing pandas nullable dtypes are changed.

## Alternatives Considered
- Continuing with the Legacy Object Dtype:
    - Retaining the ambiguous missing value semantics of the legacy object dtype does not provide a robust and consistent solution, aligning with the design of other extension arrays.
    - Not having a nullable object dtype could potentially be a blocker for a potential future nullable by default policy.

## Drawbacks and Future Directions
1. **Overhead Cost**:
The additional memory required for a boolean mask and possible performance penalties in highly heterogeneous arrays are acknowledged trade-offs.
2. **Integration Complexity**:
Ensuring seamless integration with the full suite of pandas functionality may reveal edge cases that require careful handling.
3. **Incompatibility**:
The existing object array can hold any python object, even `pd.NA` itself. The proposed nullable object array will be unable to hold `np.nan`, `None` or `pd.NaT` as these will be considered missing in the constructors and other conversions when following the existing API for the other nullable types. Users will not be able to round-trip between the legacy and nullable object dtypes.

## Conclusion
Introducing a nullable object dtype in pandas will offer a clearer semantic for missing values and align the behavior of object arrays with other nullable types. This proposal is aimed at fostering discussion and soliciting community feedback to refine the design and implementation roadmap.



## PDEP-18 History

- 07 June 2025: Initial version.
