# PDEP-2: Build System Overhaul


- Created: 5 August 2022
- Status: Under discussion
- Discussion: [#47380](https://github.com/pandas-dev/pandas/pull/47380)
- Author: [Will Ayd](https://github.com/willayd)
- Revision: 1

# Abstract

This PDEP proposes replacing pandas build system away from its legacy setuptools-based implementation to either a CMake or meson build system.

# Detailed Description

Pandas has long used setuptools as its build system, likely due to this historically being the only feasible option. Within setuptools, much of the build system logic ends up getting placed in ``setup.py`` with limited options for extension or customization. Pandas in turn has a monolithic setup.py that has its own quirks. For instance, parallel compilation "works" but only during development (yielding very long source distribution installs) and not without its own [bugs](see https://github.com/pandas-dev/pandas/issues/30873)).

As Python packaging has evolved, so too has its ability to leverage other build, packaging and installation systems. [PEP-517](https://peps.python.org/pep-0517/) and [PEP-518](https://peps.python.org/pep-0518/) (which provides more history on setuptools) have standardized how libraries can implement and extend their own build systems.

Given the opportunity to decouple a hard setuptools dependency, this PEP proposes pandas explores other popular build tools, with the goals of providing a better developer and end-user experience. The two main tools discussed so far are [CMake](https://cmake.org/) and [Meson](https://mesonbuild.com/Python-module.html)

# CMake

CMake is a popular build system, particularly in the C++ world. It is the build tool of choice for [Apache Arrow](https://arrow.apache.org/). CMake hearkens back to the early 2000s and has a very large ecosystem of documentation, books and training resources. [Why CMake?](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Why%20CMake.html) within their own documentation will provide users with history and motivation for use.

A reference implementation of CMake for pandas can be found at [#47380](https://github.com/pandas-dev/pandas/pull/47380).

# Meson

Meson is a more recent entry into the build system world. Comparisons to other build systems can be found in the [Meson documentation](https://mesonbuild.com/Comparisons.html). One of the more attractive points to Meson is that its syntax is much more Pythonic in nature than the DSL offered by CMake. Meson is used by SciPy and NumPy, with scikit-learn [likely moving](https://github.com/pandas-dev/pandas/pull/47380#issuecomment-1162817318) to Meson in the future.

A reference implementation of Meson for pandas can be found [here](https://github.com/lithomas1/pandas/pull/19)

# Consideration Points

## Syntax Differences

Meson has a language heavily influenced by Python. CMake by comparison has its own DSL, with relatively limited support for lists and little to no dict support. For a Python developer, Meson would seem more natural at first glance.

## Popularity

When looking at [StackOverflow data,](https://data.stackexchange.com/stackoverflow/query/1280378/number-of-stackoverflow-quesions-total-and-answered-on-bazel-cmake-meson-and), CMake is far and away the most tagged build system. As of August 5, 2022, the number of tagged answers for popular build systems is:

|Build System|Answered Questions|
|---|---|
|CMake|18363|
|MSBuild|12351|
|Bazel|2224|
|Meson|295|
