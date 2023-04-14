# PDEP-10: PyArrow as a required dependency

- Created: 13 April 2023
- Status: Under discussion
- Discussion: [#?](https://github.com/pandas-dev/pandas/pull/?)
              [#52509](https://github.com/pandas-dev/pandas/issues/52509)
- Author: [Matthew Roeschke](https://github.com/mroeschke)
- Revision: 1

## Abstract

This PDEP proposes that:

- PyArrow becomes a runtime dependency starting pandas 2.1
- The minimum version of PyArrow supported starting pandas 2.1 is version 6.
- The minimum version of PyArrow will be bumped every major pandas release to the highest
  PyArrow version that has been released for at least 2 years.

## Background

PyArrow has been an optional dependency of pandas since version 0.21.0. PyArrow
initially provided I/O reading functionality for formats such as Parquet and CSV. In pandas version 1.2,
pandas integrated PyArrow into the ExtensionArray interface to provide an optional string data type backed by PyArrow.
In pandas version 1.5 this functionality was expanded to support all data types that PyArrow supports. As of pandas version 2.0,
all I/O readers have the option to return PyArrow-backed data types, and a lot of methods now utilize PyArrow compute functions to
accelerate PyArrow-backed data in pandas, notibly string and datetime types.

## Motivation

While all the functionality described in the previous paragraph is currently optional, PyArrow has significant integration into many areas
of pandas. With our roadmap noting that pandas strives for better Apache Arrow interoperability [^1] and many projects [^2], within or beyond the Python ecosystem, adopting or interacting with the Arrow format, making PyArrow a required dependency provides an additional signal of confidence in the Arrow
ecosystem to pandas users.

Additionally, requiring PyArrow would simplify the related development within pandas and potentially improve NumPy functionality that would be better suited
by PyArrow including:

- Avoiding runtime checking if PyArrow is available to perform PyArrow object inference during constructor or indexing operations
- Improve NumPy object data type support by default for analogous types that have native PyArrow support such as decimal, binary, and nested types

## Drawbacks

Including PyArrow would naturally increase the installation size of pandas.

### PDEP-1 History

- 13 April 2023: Initial version

[^1] <https://pandas.pydata.org/docs/development/roadmap.html#apache-arrow-interoperability>
[^2] <https://arrow.apache.org/powered_by/>
