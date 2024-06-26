# PDEP-15: Do not require PyArrow as a required dependency (for pandas 3.0)

- Created: 8 May 2024
- Status: Under Discussion
- Discussion:  [#58623](https://github.com/pandas-dev/pandas/pull/58623)
               [#52711](https://github.com/pandas-dev/pandas/pull/52711)
               [#52509](https://github.com/pandas-dev/pandas/issues/52509)
               [#54466](https://github.com/pandas-dev/pandas/issues/54466)
- Author: [Thomas Li](https://github.com/lithomas1)
- Revision: 1

## Abstract

This PDEP was supersedes PDEP-10, which stipulated that PyArrow should become a required dependency
for pandas 3.0. After reviewing feedback posted
on the feedback issue [#54466](https://github.com/pandas-dev/pandas/issues/54466), we, the members of
the core team, have decided against moving forward with this PDEP for pandas 3.0.

The primary reasons for rejecting this PDEP are twofold:

1) Requiring pyarrow as a dependency causes installation problems.
   - Pyarrow does not fit or has a hard time fitting in space-constrained environments
such as AWS Lambda and WASM, due to its large size of around ~40 MB for a compiled wheel
(which is larger than pandas' own wheel sizes)
   - Installation of pyarrow is not possible on some platforms. We provide support for some
less widely used platforms such as Alpine Linux (and there is third party support for pandas in
pyodide, a WASM distribution of pandas), both of which pyarrow does not provide wheels for.

   While both of these reasons are mentioned in the drawbacks section of this PDEP, at the time of the writing
of the PDEP, we underestimated the impact this would have on users, and also downstream developers.

2) Many of the benefits presented in PDEP-10 can be materialized for users that have pyarrow installed, without
   forcing a pyarrow requirement on other users.

   In PDEP-10, there are three primary benefits listed:

    - First class support for strings.

      - This is covered by PDEP-14, which will enable the usage of a pyarrow backed string dtype by default,
        (for users who have pyarrow instaleld) and the use of a Python object based fallback in the case .

    - Support for dtypes not present in pandas (e.g. nested dtypes, decimals)
      - Users can already create arrays with these dtypes if they have pyarrow installed, but we cannot infer
      arrays to those dtypes by default without pyarrow installed (as there is no Python/numpy equivalent).

    - Interoperability
      - The Arrow C Data Interface would allow us to import/export pandas DataFrames to and from other libraries
        that support Arrow in a zero-copy manner.

        Support for the Arrow C Data interface in pandas and other libraries in the ecosystem is still very new, though,
        (support in pandas itself was only added as of pandas 2.2), and the dataframe interchange protocol, which allows
        for dataframe interchange between Python dataframe implementations is currently better supported in downstream
        libraries.

Although this PR recommends not adopting pyarrow as a required dependency in pandas 3.0, this does not mean that we are
abandoning pyarrow support and integration in pandas. Adopting support for pyarrow arrays
and data types in more of pandas will lead to greater interoperability with the
ecosystem and better performance for users. Furthermore, a lot of the drawbacks, such as the large installation size of pyarrow
and the lack of support for certain platforms, can be solved, and potential solutions have been proposed for them, allowing us
to potentially revisit this decision in the future.

However, at this point in time, it is clear that we are not ready to require pyarrow
as a dependency in pandas.
