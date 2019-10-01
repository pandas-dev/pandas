# Roadmap

This page provides an overview of the major themes in pandas'
development. Each of these items requires a relatively large amount of
effort to implement. These may be achieved more quickly with dedicated
funding or interest from contributors.

An item being on the roadmap does not mean that it will *necessarily*
happen, even with unlimited funding. During the implementation period we
may discover issues preventing the adoption of the feature.

Additionally, an item *not* being on the roadmap does not exclude it
from inclusion in pandas. The roadmap is intended for larger,
fundamental changes to the project that are likely to take months or
years of developer time. Smaller-scoped items will continue to be
tracked on our [issue tracker](https://github.com/pandas-dev/pandas/issues).

See [Roadmap evolution](#roadmap-evolution) for proposing
changes to this document.

## Extensibility

Pandas `extending.extension-types` allow
for extending NumPy types with custom data types and array storage.
Pandas uses extension types internally, and provides an interface for
3rd-party libraries to define their own custom data types.

Many parts of pandas still unintentionally convert data to a NumPy
array. These problems are especially pronounced for nested data.

We'd like to improve the handling of extension arrays throughout the
library, making their behavior more consistent with the handling of
NumPy arrays. We'll do this by cleaning up pandas' internals and
adding new methods to the extension array interface.

## String data type

Currently, pandas stores text data in an `object` -dtype NumPy array.
The current implementation has two primary drawbacks: First, `object`
-dtype is not specific to strings: any Python object can be stored in an
`object` -dtype array, not just strings. Second: this is not efficient.
The NumPy memory model isn't especially well-suited to variable width
text data.

To solve the first issue, we propose a new extension type for string
data. This will initially be opt-in, with users explicitly requesting
`dtype="string"`. The array backing this string dtype may initially be
the current implementation: an `object` -dtype NumPy array of Python
strings.

To solve the second issue (performance), we'll explore alternative
in-memory array libraries (for example, Apache Arrow). As part of the
work, we may need to implement certain operations expected by pandas
users (for example the algorithm used in, `Series.str.upper`). That work
may be done outside of pandas.

## Apache Arrow interoperability

[Apache Arrow](https://arrow.apache.org) is a cross-language development
platform for in-memory data. The Arrow logical types are closely aligned
with typical pandas use cases.

We'd like to provide better-integrated support for Arrow memory and
data types within pandas. This will let us take advantage of its I/O
capabilities and provide for better interoperability with other
languages and libraries using Arrow.

## Block manager rewrite

We'd like to replace pandas current internal data structures (a
collection of 1 or 2-D arrays) with a simpler collection of 1-D arrays.

Pandas internal data model is quite complex. A DataFrame is made up of
one or more 2-dimensional "blocks", with one or more blocks per dtype.
This collection of 2-D arrays is managed by the BlockManager.

The primary benefit of the BlockManager is improved performance on
certain operations (construction from a 2D array, binary operations,
reductions across the columns), especially for wide DataFrames. However,
the BlockManager substantially increases the complexity and maintenance
burden of pandas.

By replacing the BlockManager we hope to achieve

-   Substantially simpler code
-   Easier extensibility with new logical types
-   Better user control over memory use and layout
-   Improved micro-performance
-   Option to provide a C / Cython API to pandas' internals

See [these design
documents](https://dev.pandas.io/pandas2/internal-architecture.html#removal-of-blockmanager-new-dataframe-internals)
for more.

## Decoupling of indexing and internals

The code for getting and setting values in pandas' data structures
needs refactoring. In particular, we must clearly separate code that
converts keys (e.g., the argument to `DataFrame.loc`) to positions from
code that uses these positions to get or set values. This is related to
the proposed BlockManager rewrite. Currently, the BlockManager sometimes
uses label-based, rather than position-based, indexing. We propose that
it should only work with positional indexing, and the translation of
keys to positions should be entirely done at a higher level.

Indexing is a complicated API with many subtleties. This refactor will
require care and attention. More details are discussed at
<https://github.com/pandas-dev/pandas/wiki/(Tentative)-rules-for-restructuring-indexing-code>

## Numba-accelerated operations

[Numba](https://numba.pydata.org) is a JIT compiler for Python code.
We'd like to provide ways for users to apply their own Numba-jitted
functions where pandas accepts user-defined functions (for example,
`Series.apply`,
`DataFrame.apply`,
`DataFrame.applymap`, and in groupby and
window contexts). This will improve the performance of
user-defined-functions in these operations by staying within compiled
code.

## Documentation improvements

We'd like to improve the content, structure, and presentation of the
pandas documentation. Some specific goals include

-   Overhaul the HTML theme with a modern, responsive design
    (`15556`)
-   Improve the "Getting Started" documentation, designing and writing
    learning paths for users different backgrounds (e.g. brand new to
    programming, familiar with other languages like R, already familiar
    with Python).
-   Improve the overall organization of the documentation and specific
    subsections of the documentation to make navigation and finding
    content easier.

## Package docstring validation

To improve the quality and consistency of pandas docstrings, we've
developed tooling to check docstrings in a variety of ways.
<https://github.com/pandas-dev/pandas/blob/master/scripts/validate_docstrings.py>
contains the checks.

Like many other projects, pandas uses the
[numpydoc](https://numpydoc.readthedocs.io/en/latest/) style for writing
docstrings. With the collaboration of the numpydoc maintainers, we'd
like to move the checks to a package other than pandas so that other
projects can easily use them as well.

## Performance monitoring

Pandas uses [airspeed velocity](https://asv.readthedocs.io/en/stable/)
to monitor for performance regressions. ASV itself is a fabulous tool,
but requires some additional work to be integrated into an open source
project's workflow.

The [asv-runner](https://github.com/asv-runner) organization, currently
made up of pandas maintainers, provides tools built on top of ASV. We
have a physical machine for running a number of project's benchmarks,
and tools managing the benchmark runs and reporting on results.

We'd like to fund improvements and maintenance of these tools to

-   Be more stable. Currently, they're maintained on the nights and
    weekends when a maintainer has free time.
-   Tune the system for benchmarks to improve stability, following
    <https://pyperf.readthedocs.io/en/latest/system.html>
-   Build a GitHub bot to request ASV runs *before* a PR is merged.
    Currently, the benchmarks are only run nightly.

## Roadmap Evolution

Pandas continues to evolve. The direction is primarily determined by
community interest. Everyone is welcome to review existing items on the
roadmap and to propose a new item.

Each item on the roadmap should be a short summary of a larger design
proposal. The proposal should include

1.  Short summary of the changes, which would be appropriate for
    inclusion in the roadmap if accepted.
2.  Motivation for the changes.
3.  An explanation of why the change is in scope for pandas.
4.  Detailed design: Preferably with example-usage (even if not
    implemented yet) and API documentation
5.  API Change: Any API changes that may result from the proposal.

That proposal may then be submitted as a GitHub issue, where the pandas
maintainers can review and comment on the design. The [pandas mailing
list](https://mail.python.org/mailman/listinfo/pandas-dev) should be
notified of the proposal.

When there's agreement that an implementation would be welcome, the
roadmap should be updated to include the summary and a link to the
discussion issue.
