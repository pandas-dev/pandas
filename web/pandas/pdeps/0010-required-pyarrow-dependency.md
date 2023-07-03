# PDEP-10: PyArrow as a required dependency for default string inference implementation

- Created: 17 April 2023
- Status: Under discussion
- Discussion: [#52711](https://github.com/pandas-dev/pandas/pull/52711)
              [#52509](https://github.com/pandas-dev/pandas/issues/52509)
- Author: [Matthew Roeschke](https://github.com/mroeschke)
          [Patrick Hoefler](https://github.com/phofl)
- Revision: 1

## Abstract

This PDEP proposes that:

- PyArrow becomes a required runtime dependency starting with pandas 3.0
- The minimum version of PyArrow supported starting with pandas 3.0 is version 7 of PyArrow.
- When the minimum version of PyArrow is bumped, PyArrow will be bumped to the highest version that has
  been released for at least 2 years.
- The pandas 2.1 release notes will have a big warning that PyArrow will become a required dependency starting
  with pandas 3.0.
- Starting in pandas 2.2, pandas raises a ``FutureWarning`` when PyArrow is not installed in the users
  environment when pandas is imported. This will ensure that only one warning is raised and users can
  easily silence it if necessary.
- Starting in pandas 3.0, the default type inferred for string data will be `ArrowDtype` with `pyarrow.string`
  instead of `object`

## Background

PyArrow is an optional dependency of pandas that provides a wide range of supplemental features to pandas:

- Since pandas 0.21.0, PyArrow provided I/O reading functionality for Parquet
- Since pandas 1.2.0, pandas integrated PyArrow into the `ExtensionArray` interface to provide an
  optional string data type backed by PyArrow
- Since pandas 1.4.0, PyArrow provided I/0 reading functionality for CSV
- Since pandas 1.5.0, pandas provided an `ArrowExtensionArray` and `ArrowDtype` to support all PyArrow
  data types within the `ExtensionArray` interface
- Since pandas 2.0.0, all I/O readers have the option to return PyArrow-backed data types, and many methods
  now utilize PyArrow compute functions to
accelerate PyArrow-backed data in pandas, notibly string and datetime types.

As of pandas 2.0, one can feasibly utilize PyArrow as an alternative data representation to NumPy with advantages such as:

1. Consistent `NA` support for all data types
2. Broader support of data types such as `decimal`, `date` and nested types

Currently, when users pass string data into pandas constructors without specifying a data type, the resulting data type
is `object`. With pyarrow string support available since 1.2.0, requiring pyarrow for 3.0 will allow pandas to default
the inferred type to the more efficient pyarrow string type.

```python
In [1]: import pandas as pd

In [2]: pd.Series(["a"]).dtype
# Current behavior
Out[2]: dtype('O')

# Future behavior in 3.0
Out[2]: string[pyarrow]
```

## Motivation

While all the functionality described in the previous paragraph is currently optional, PyArrow has significant
integration into many areas of pandas. With our roadmap noting that pandas strives for better Apache Arrow
interoperability [^1] and many projects [^2], within or beyond the Python ecosystem, adopting or interacting with
the Arrow format, making PyArrow a required dependency provides an additional signal of confidence in the Arrow
ecosystem to pandas users.

Additionally, requiring PyArrow would simplify the related development within pandas and potentially improve NumPy
functionality that would be better suited by PyArrow including:

- Avoiding runtime checking if PyArrow is available to perform PyArrow object inference during constructor or indexing operations

- Removing redundant functionality:
  - fastparquet engine in `read_parquet`
  - potentially simplifying the `read_csv` logic (needs more investigation)

- Avoiding NumPy object data types more by default for analogous types that have native PyArrow support such as:
  - decimal
  - binary
  - nested types (list or dict data)
  - strings

Out of this group, strings offer the most advantages for users. They use significantly less memory and are faster:

**Performance:**

```python
import string
import random

import pandas as pd


def random_string() -> str:
    return "".join(random.choices(string.printable, k=random.randint(10, 100)))


ser_object = pd.Series([random_string() for _ in range(1_000_000)])
ser_string = ser_object.astype("string[pyarrow]")\
```

PyArrow backed strings are significantly faster than NumPy object strings:

*str.len*

```python
In[1]: %timeit ser_object.str.len()
118 ms ± 260 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

In[2]: %timeit ser_string.str.len()
24.2 ms ± 187 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

*str.startswith*

```python
In[3]: %timeit ser_object.str.startswith("a")
136 ms ± 300 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

In[4]: %timeit ser_string.str.startswith("a")
11 ms ± 19.8 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

Another advantage is I/O. PyArrow engines in pandas can provide a significant speedup. Currently, the data
are cast to NumPy dtypes, which requires roundtripping when converting back to PyArrow strings explicitly, which
hinders performance.

**Memory**

PyArrow backed strings use significantly less memory. Dask developers investigated this [here](https://www.coiled.io/blog/pyarrow-strings-in-dask-dataframes).

Short summary: PyArrow strings required 1/3 of the original memory.


## Drawbacks

Including PyArrow would naturally increase the installation size of pandas. For example, installing pandas and PyArrow
using pip from wheels, numpy and pandas requires about `70MB`, and including PyArrow requires around `120MB`. An increase
of installation size would have negative implication using pandas in space-constrained development or deployment environments
such as AWS Lambda.

Additionally, if a user is installing pandas in an environment where wheels are not available through a `pip install` or `conda install`,
the user will need to also build Arrow C++ and related dependencies when installing from source. These environments include

- Alpine linux (commonly used as a base for Docker containers)
- WASM (pyodide and pyscript)
- Python development versions

Lastly, pandas development and releases will need to be mindful of PyArrow's development and release cadance. For example when
supporting a newly released Python version, pandas will also need to be mindful of PyArrow's wheel support for that Python version
before releasing a new pandas version.

### PDEP-10 History

- 17 April 2023: Initial version
- 8 May 2023: Changed proposal to make pyarrow required in pandas 3.0 instead of 2.1

[^1] <https://pandas.pydata.org/docs/development/roadmap.html#apache-arrow-interoperability>
[^2] <https://arrow.apache.org/powered_by/>
