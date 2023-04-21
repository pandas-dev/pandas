# PDEP-10: PyArrow as a required dependency

- Created: 17 April 2023
- Status: Under discussion
- Discussion: [#52711](https://github.com/pandas-dev/pandas/pull/52711)
              [#52509](https://github.com/pandas-dev/pandas/issues/52509)
- Author: [Matthew Roeschke](https://github.com/mroeschke)
- Revision: 1

## Abstract

This PDEP proposes that:

- PyArrow becomes a runtime dependency starting with pandas 2.1
- The minimum version of PyArrow supported starting with pandas 2.1 is version 7 of PyArrow.
- The minimum version of PyArrow will be bumped every major pandas release to the highest
  PyArrow version that has been released for at least 2 years, and the minimum PyArrow version will be
  maintained for every minor version in the major version series.

## Background

PyArrow is an optional dependency of pandas that provides a wide range of supplemental features to pandas:

- Since pandas 0.21.0, PyArrow provided I/O reading functionality for Parquet
- Since pandas 1.2.0, pandas integrated PyArrow into the `ExtensionArray` interface to provide an optional string data type backed by PyArrow
- Since pandas 1.4.0, PyArrow provided I/0 reading functionality for CSV
- Since pandas 1.5.0, pandas provided an `ArrowExtensionArray` and `ArrowDtype` to support all PyArrow data types within the `ExtensionArray` interface
- Since pandas 2.0.0, All I/O readers have the option to return PyArrow-backed data types, and many methods now utilize PyArrow compute functions to
accelerate PyArrow-backed data in pandas, notibly string and datetime types.

As of pandas 2.0, one can feasibly utilize PyArrow as an alternative data representation to NumPy with advantages such as:

1. Consistent ``NA`` support for all data types
2. Broader support of data types such as ``decimal``, ``date`` and nested types

## Motivation

While all the functionality described in the previous paragraph is currently optional, PyArrow has significant integration into many areas
of pandas. With our roadmap noting that pandas strives for better Apache Arrow interoperability [^1] and many projects [^2], within or beyond the Python ecosystem, adopting or interacting with the Arrow format, making PyArrow a required dependency provides an additional signal of confidence in the Arrow
ecosystem to pandas users.

Additionally, requiring PyArrow would simplify the related development within pandas and potentially improve NumPy functionality that would be better suited
by PyArrow including:

- Avoiding runtime checking if PyArrow is available to perform PyArrow object inference during constructor or indexing operations
- Avoiding NumPy object data types more by default for analogous types that have native PyArrow support such as decimal, binary, and nested types

## Drawbacks

Including PyArrow would naturally increase the installation size of pandas. For example, installing pandas and PyArrow using pip from wheels, numpy and pandas
are about `70MB`, and PyArrow is around `120MB`. An increase of installation size would have negative impliciation using pandas in space-constrained development
or deployment environments such as AWS Lambda.

Additionally, if a user is installing pandas in an environment where wheels are not available and needs to build from source, the user will need to build Arrow C++ and related dependencies. These environments include

- Alpine linux (commonly used as a base for Docker containers)
- WASM (pyodide and pyscript)
- Python development versions

Lastly, pandas development and releases will need to be mindful of PyArrow's development and release cadance. For example when supporting a newly released Python version, pandas will also need to be mindful of PyArrow's wheel support for that Python version before releasing a new pandas version.

### PDEP-1 History

- 17 April 2023: Initial version

[^1] <https://pandas.pydata.org/docs/development/roadmap.html#apache-arrow-interoperability>
[^2] <https://arrow.apache.org/powered_by/>
