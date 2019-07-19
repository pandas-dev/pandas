
.. _roadmap:

=======
Roadmap
=======

This page provides an overview of the major themes pandas' development.

Extensibility
-------------

Pandas Extension Arrays provide 3rd party libraries the ability to
extend pandas' supported types. In theory, these 3rd party types can do
everything one of pandas. In practice, many places in pandas will unintentionally
convert the ExtensionArray to a NumPy array of Python objects, causing
performance issues and the loss of type information.

Additionally, there are known issues when the scalar type of an ExtensionArray
is actual a sequence (nested data). Developing APIs for working with nested data,
and ensuring that pandas' internal routines can handled it, will be a major effort.

String Dtype
------------

Currently, pandas stores text data in an ``object`` -dtype NumPy array.
Each array stores Python strings. While pragmatic, since we rely on NumPy
for storage and Python for string operations, this is memory inefficient
and slow. We'd like to provide a native string type for pandas.

The most obvious candidate is Apache Arrow. Currently, Arrow provides
storage for string data. We can work with the Arrow developers to implement
the operations needed for pandas users (for example, ``Series.str.upper``).

Apache Arrow Interoperability
-----------------------------

`Apache Arrow <https://arrow.apache.org>`__ is a cross-language development
platform for in-memory data. The Arrow logical types are closely aligned with
typical pandas use cases (for example, support for nullable integers).

We'd like have a pandas DataFrame be backed by a collection of Apache Arrow
arrays. This should simplify pandas internals and ensure more consistent
handling of data types through operations.

Block Manager Rewrite
---------------------

Pandas internal data model is quite complex. A DataFrame is made up of
one or more 2-dimension "blocks", with one or more blocks per dtype. This
collection of 2-D arrays is managed by the BlockManager.

The primary benefit of the BlockManager is improved performance on certain
operations, especially for wide DataFrames. Consider summing a table with ``P``
columns. When stored as a 2-D array, this results in a single call to
``numpy.sum``. If this were stored as ``P`` arrays, we'd have a Python for loop
going calling ``numpy.sum`` P times.

By replacing the BlockManager we hope to achieve

* Substantially simpler code
* Easier extensibility with new logical types
* Better user control over memory use and layout
* Improved microperformance

See `these design documents <https://dev.pandas.io/pandas2/internal-architecture.html#removal-of-blockmanager-new-dataframe-internals>`__
for more.

Weighted Operations
-------------------

In many fields, sample weights are necessary to correctly estimate population
statistics. We'd like to support weighted operations (like `mean`, `sum`, `std`,
etc.), possibly with an API similar to `DataFrame.groupby`.

See https://github.com/pandas-dev/pandas/issues/10030 for more.

Package Docstring Validation
----------------------------

To improve the quality and consistency of pandas docstrings, we've developed
tooling to check docstrings in a variety of ways.
https://github.com/pandas-dev/pandas/blob/master/scripts/validate_docstrings.py
contains the checks.

Like many other projects, pandas uses the
[numpydoc](https://numpydoc.readthedocs.io/en/latest/) style for writing
docstrings. With the collaboration of the numpydoc maintainers, we'd like to
move the checks to a package other than pandas so that other projects can easily
use them as well.

Performance Monitoring
----------------------

Pandas uses [airspeed velocity](https://asv.readthedocs.io/en/stable/) to
monitor for performance regressions. ASV itself is a fabulous tool, but requires
some additional work to be integrated into an open source project's workflow.

The [asv-runner](https://github.com/asv-runner) organization, currently made up
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

