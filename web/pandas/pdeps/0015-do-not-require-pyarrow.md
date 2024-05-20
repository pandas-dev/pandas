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

2) Many of the benefits presented in this PDEP can be materialized even with payrrow as an optional dependency.

   For example, as detailed in PDEP-14, it is possible to create a new string data type with the same semantics
   as our current default object string data type, but that allows users to experience faster performance and memory savings
   compared to the object strings (if pyarrow is installed).

While we've decided to not move forward with requiring pyarrow in pandas 3.0, the rejection of this PDEP
does not mean that we are abandoning pyarrow support and integration in pandas. We, as the core team, still believe
that adopting support for pyarrow arrays and data types in more of pandas will lead to greater interoperability with the
ecosystem and better performance for users. Furthermore, a lot of the drawbacks, such as the large installation size of pyarrow
and the lack of support for certain platforms, can be solved, and potential solutions have been proposed for them, allowing us
to potentially revisit this decision in the future.

However, at this point in time, it is clear that we are not ready to require pyarrow
as a dependency in pandas.
