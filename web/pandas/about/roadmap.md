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

The roadmap is defined as a set of major enhancement proposals named PDEPs.
For more information about PDEPs, and how to submit one, please refer to
[PEDP-1]({{ base_url }}pdeps/0001-purpose-and-guidelines.html).

## PDEPs

{% for pdep_type in ["Under discussion", "Accepted", "Implemented", "Rejected"] %}

<h3 id="pdeps-{{pdep_type}}">{{ pdep_type.replace("_", " ").capitalize() }}</h3>

<ul>
{% for pdep in pdeps[pdep_type] %}
    <li><a href="{% if not pdep.url.startswith("http") %}{{ base_url }}{% endif %}{{ pdep.url }}">{{ pdep.title }}</a></li>
{% else %}
    <li>There are currently no PDEPs with this status</li>
{% endfor %}
</ul>

{% endfor %}

## Roadmap points pending a PDEP

<div class="alert alert-warning" role="alert">
  pandas is in the process of moving roadmap points to PDEPs (implemented in
  August 2022). During the transition, some roadmap points will exist as PDEPs,
  while others will exist as sections below.
</div>

### Extensibility

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

### String data type

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

### Apache Arrow interoperability

[Apache Arrow](https://arrow.apache.org) is a cross-language development
platform for in-memory data. The Arrow logical types are closely aligned
with typical pandas use cases.

We'd like to provide better-integrated support for Arrow memory and
data types within pandas. This will let us take advantage of its I/O
capabilities and provide for better interoperability with other
languages and libraries using Arrow.

### Block manager rewrite

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

### Decoupling of indexing and internals

The code for getting and setting values in pandas' data structures
needs refactoring. In particular, we must clearly separate code that
converts keys (e.g., the argument to `DataFrame.loc`) to positions from
code that uses these positions to get or set values. This is related to
the proposed BlockManager rewrite. Currently, the BlockManager sometimes
uses label-based, rather than position-based, indexing. We propose that
it should only work with positional indexing, and the translation of
keys to positions should be entirely done at a higher level.

Indexing is a complicated API with many subtleties. This refactor will require care
and attention. The following principles should inspire refactoring of indexing code and
should result on cleaner, simpler, and more performant code.

1. Label indexing must never involve looking in an axis twice for the same label(s).
This implies that any validation step must either:

  * limit validation to general features (e.g. dtype/structure of the key/index), or
  * reuse the result for the actual indexing.

2. Indexers must never rely on an explicit call to other indexers.
For instance, it is OK to have some internal method of `.loc` call some
internal method of `__getitem__` (or of their common base class),
but never in the code flow of `.loc` should `the_obj[something]` appear.

3. Execution of positional indexing must never involve labels (as currently, sadly, happens).
That is, the code flow of a getter call (or a setter call in which the right hand side is non-indexed)
to `.iloc` should never involve the axes of the object in any way.

4. Indexing must never involve accessing/modifying values (i.e., act on `._data` or `.values`) more than once.
The following steps must hence be clearly decoupled:

  * find positions we need to access/modify on each axis
  * (if we are accessing) derive the type of object we need to return (dimensionality)
  * actually access/modify the values
  * (if we are accessing) construct the return object

5. As a corollary to the decoupling between 4.i and 4.iii, any code which deals on how data is stored
(including any combination of handling multiple dtypes, and sparse storage, categoricals, third-party types)
must be independent from code that deals with identifying affected rows/columns,
and take place only once step 4.i is completed.

  * In particular, such code should most probably not live in `pandas/core/indexing.py`
  * ... and must not depend in any way on the type(s) of axes (e.g. no `MultiIndex` special cases)

6. As a corollary to point 1.i, `Index` (sub)classes must provide separate methods for any desired validity check of label(s) which does not involve actual lookup,
on the one side, and for any required conversion/adaptation/lookup of label(s), on the other.

7. Use of trial and error should be limited, and anyway restricted to catch only exceptions
which are actually expected (typically `KeyError`).

  * In particular, code should never (intentionally) raise new exceptions in the `except` portion of a `try... exception`

8. Any code portion which is not specific to setters and getters must be shared,
and when small differences in behavior are expected (e.g. getting with `.loc` raises for
missing labels, setting still doesn't), they can be managed with a specific parameter.

### Numba-accelerated operations

[Numba](https://numba.pydata.org) is a JIT compiler for Python code.
We'd like to provide ways for users to apply their own Numba-jitted
functions where pandas accepts user-defined functions (for example,
`Series.apply`,
`DataFrame.apply`,
`DataFrame.applymap`, and in groupby and
window contexts). This will improve the performance of
user-defined-functions in these operations by staying within compiled
code.

### Documentation improvements

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

### Performance monitoring

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
