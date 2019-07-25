.. _roadmap:

=======
Roadmap
=======

This page provides an overview of the major themes pandas' development. Implementation
of these goals may be hastened with dedicated funding.

Extensibility
-------------

Pandas :ref:`extending.extension-types` provide 3rd-party libraries the ability to
extend pandas' supported data types. Many parts of pandas still unintentionally
convert data to a NumPy array. These problems are especially pronounced for nested data.

We'd like to improve the handling of extension arrays throughout the library,
making their behavior more consistent with the handling of NumPy arrays. We'll do this
by cleaning up pandas' internals and adding new methods do the extension array interface.

String Dtype
------------

Currently, pandas stores text data in an ``object`` -dtype NumPy array.
Each array stores Python strings. While pragmatic, since we rely on NumPy
for storage and Python for string operations, this is memory inefficient
and slow. We'd like to provide a native string type for pandas.

The most obvious alternative is Apache Arrow. Currently, Arrow provides
storage for string data. We can work with the Arrow developers to implement
the operations needed for pandas users (for example, ``Series.str.upper``).
These operations could be implemented in Numba (
as prototyped in `Fletcher <https://github.com/xhochy/fletcher>`__)
or in the Apache Arrow C++ library.

Apache Arrow Interoperability
-----------------------------

`Apache Arrow <https://arrow.apache.org>`__ is a cross-language development
platform for in-memory data. The Arrow logical types are closely aligned with
typical pandas use cases.

We'd like to provide better-integrated support for Arrow memory and data types
within pandas. This will let us take advantage of its I/O capabilities and
provides for better interoperability with other languages and libraries
using Arrow.

Block Manager Rewrite
---------------------

We'd like to replace pandas current internal data structures (a collection of
1 or 2-D arrays) with a simpler collection of 1-D arrays.

Pandas internal data model is quite complex. A DataFrame is made up of
one or more 2-dimension "blocks", with one or more blocks per dtype. This
collection of 2-D arrays is managed by the BlockManager.

The primary benefit of the BlockManager is improved performance on certain
operations (construction from a 2D array, binary operations, reductions across the columns),
especially for wide DataFrames. However, the BlockManager substantially increases the
complexity and maintenance burden of pandas'.

By replacing the BlockManager we hope to achieve

* Substantially simpler code
* Easier extensibility with new logical types
* Better user control over memory use and layout
* Improved microperformance
* Option to provide a C / Cython API to pandas' internals

See `these design documents <https://dev.pandas.io/pandas2/internal-architecture.html#removal-of-blockmanager-new-dataframe-internals>`__
for more.

Decoupling of Indexing and Internals
------------------------------------

The code for getting and setting values in pandas' data structures needs refactoring.
In particular, a clear separation must be made between code that
converts keys (e.g., the argument to ``DataFrame.loc``) to positions from code that uses
uses these positions to get or set values. This is related to the proposed BlockManager rewrite.
Currently, the BlockManager sometimes uses label-based, rather than position-based, indexing.
We propose that it should only work with positional indexing, and the translation of keys
to positions should be entirely done at a higher level.

Indexing is a complicated API with many subtleties. This refactor will require care
and attention. More details are discussed at
https://github.com/pandas-dev/pandas/wiki/(Tentative)-rules-for-restructuring-indexing-code

Numba-Accelerated Operations
----------------------------

[Numba](https://numba.pydata.org) is a JIT compiler for Python code. We'd like to provide
ways for users to apply their own Numba-jitted functions within pandas' groupby and window
contexts. This will improve the performance of user-defined-functions in these operations
by staying within compiled code.


Weighted Operations
-------------------

In many fields, sample weights are necessary to correctly estimate population
statistics. We'd like to support weighted operations (like `mean`, `sum`, `std`,
etc.), possibly with an API similar to `DataFrame.groupby`.

See https://github.com/pandas-dev/pandas/issues/10030 for more.

Documentation Improvements
--------------------------

We'd like to improve the content, structure, and presentation of pandas documentation.
Some specific goals include

* Overhaul the HTML theme with a modern, responsive design.
* Improve the "Getting Started" documentation, designing and writing learning paths
  for users different backgrounds (e.g. brand new to programming, familiar with
  other languages like R, already familiar with Python).
* Improve the overall organization of the documentation and specific subsections
  of the documentation to make navigation and finding content easier.

Package Docstring Validation
----------------------------

To improve the quality and consistency of pandas docstrings, we've developed
tooling to check docstrings in a variety of ways.
https://github.com/pandas-dev/pandas/blob/master/scripts/validate_docstrings.py
contains the checks.

Like many other projects, pandas uses the
`numpydoc <https://numpydoc.readthedocs.io/en/latest/>`__ style for writing
docstrings. With the collaboration of the numpydoc maintainers, we'd like to
move the checks to a package other than pandas so that other projects can easily
use them as well.

Performance Monitoring
----------------------

Pandas uses `airspeed velocity <https://asv.readthedocs.io/en/stable/>`__ to
monitor for performance regressions. ASV itself is a fabulous tool, but requires
some additional work to be integrated into an open source project's workflow.

The `asv-runner <https://github.com/asv-runner>`__ organization, currently made up
of pandas maintainers, provides tools built on top of ASV. We have a physical
machine for running a number of project's benchmarks, and tools managing the
benchmark runs and reporting on results.

We'd like to fund improvements and maintenance of these tools to

* Be more stable. Currently, they're maintained on the nights and weekends when
  a maintainer has free time.
* Tune the system for benchmarks to improve stability, following
  https://pyperf.readthedocs.io/en/latest/system.html
* Build a GitHub bot to request ASV runs *before* a PR is merged. Currently, the
  benchmarks are only run nightly.
