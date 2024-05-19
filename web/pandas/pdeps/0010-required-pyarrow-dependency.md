# PDEP-10: PyArrow as a required dependency for default string inference implementation

- Created: 17 April 2023
- Status: Accepted
- Discussion: [#52711](https://github.com/pandas-dev/pandas/pull/52711)
              [#52509](https://github.com/pandas-dev/pandas/issues/52509)
- Author: [Matthew Roeschke](https://github.com/mroeschke)
          [Patrick Hoefler](https://github.com/phofl)
- Revision: 1

[TOC]

## Abstract

This PDEP proposes that:

- PyArrow becomes a required runtime dependency starting with pandas 3.0
- The minimum version of PyArrow supported starting with pandas 3.0 is version 7 of PyArrow.
- When the minimum version of PyArrow is bumped, PyArrow will be bumped to the highest version that has
  been released for at least 2 years.
- The pandas 2.1 release notes will have a big warning that PyArrow will become a required dependency starting
  with pandas 3.0. We will pin a feedback issue on the pandas issue tracker. The note in the release notes will point
  to that issue.
- Starting in pandas 2.2, pandas raises a ``FutureWarning`` when PyArrow is not installed in the users
  environment when pandas is imported. This will ensure that only one warning is raised and users can
  easily silence it if necessary. This warning will point to the feedback issue.
- Starting in pandas 3.0, the default type inferred for string data will be `ArrowDtype` with `pyarrow.string`
  instead of `object`. Additionally, we will infer all dtypes that are listed below as well instead of storing as object.

This will bring **immediate benefits to users**, as well as opening up the door for significant further
benefits in the future.

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

1. Consistent `NA` support for all data types;
2. Broader support of data types such as `decimal`, `date` and nested types;
3. Better interoperability with other dataframe libraries based on Arrow.

## Motivation

While all the functionality described in the previous paragraph is currently optional, PyArrow has significant
integration into many areas of pandas. With our roadmap noting that pandas strives for better Apache Arrow
interoperability [^1] and many projects [^2], within or beyond the Python ecosystem, adopting or interacting with
the Arrow format, making PyArrow a required dependency provides an additional signal of confidence in the Arrow
ecosystem (as well as improving interoperability with it).

### Immediate User Benefit 1: pyarrow strings

Currently, when users pass string data into pandas constructors without specifying a data type, the resulting data type
is `object`, which has significantly much worse memory usage and performance as compared to pyarrow strings.
With pyarrow string support available since 1.2.0, requiring pyarrow for 3.0 will allow pandas to default
the inferred type to the more efficient pyarrow string type.

```python
In [1]: import pandas as pd

In [2]: pd.Series(["a"]).dtype
# Current behavior
Out[2]: dtype('O')

# Future behavior in 3.0
Out[2]: string[pyarrow]
```

Dask developers investigated performance and memory of pyarrow strings [here](https://www.coiled.io/blog/pyarrow-strings-in-dask-dataframes),
and found them to be a significant improvement over the current `object` dtype.

Little demo:
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

### Immediate User Benefit 2: Nested Datatypes

Currently, if you try storing `dict`s in a pandas `Series`, you will again get the horrendeous `object` dtype:
```python
In [6]: pd.Series([{'a': 1, 'b': 2}, {'a': 2, 'b': 99}])
Out[6]:
0     {'a': 1, 'b': 2}
1    {'a': 2, 'b': 99}
dtype: object
```

If `pyarrow` were required, this could have been auto-inferred to be `pyarrow.struct`, which again
would come with memory and performance improvements.

### Immediate User Benefit 3: Interoperability

Other Arrow-backed dataframe libraries are growing in popularity. Having the same memory representation
would improve interoperability with them, as operations such as:
```python
import pandas as pd
import polars as pl

df = pd.DataFrame(
  {
    'a': ['one', 'two'],
    'b': [{'name': 'Billy', 'age': 3}, {'name': 'Bob', 'age': 4}],
  }
)
pl.from_pandas(df)
```
could be zero-copy. Users making use of multiple dataframe libraries would more easily be able to
switch between them.

### Future User Benefits:

Requiring PyArrow would simplify the related development within pandas and potentially improve NumPy
functionality that would be better suited by PyArrow including:

- Avoiding runtime checking if PyArrow is available to perform PyArrow object inference during constructor or indexing operations

- NumPy object dtype will be avoided as much as possible. This means that every dtype that has a PyArrow equivalent is inferred automatically as such. This includes:
  - decimal
  - binary
  - nested types (list or dict data)
  - strings
  - time
  - date

#### Developer benefits

First, this would simplify development of pyarrow-backed datatypes, as it would avoid
optional dependency checks.

Second, it could potentially remove redundant functionality:
- fastparquet engine in `read_parquet`;
- potentially simplifying the `read_csv` logic (needs more investigation);
- factorization;
- datetime/timezone ops.

## Drawbacks

Including PyArrow would naturally increase the installation size of pandas. For example, installing pandas and PyArrow
using pip from wheels, numpy and pandas requires about `70MB`, and including PyArrow requires an additional `120MB`.
An increase of installation size would have negative implication using pandas in space-constrained development or deployment environments
such as AWS Lambda.

Additionally, if a user is installing pandas in an environment where wheels are not available through a `pip install` or `conda install`,
the user will need to also build Arrow C++ and related dependencies when installing from source. These environments include

- Alpine linux (commonly used as a base for Docker containers)
- WASM (pyodide and pyscript)
- Python development versions

Lastly, pandas development and releases will need to be mindful of PyArrow's development and release cadance. For example when
supporting a newly released Python version, pandas will also need to be mindful of PyArrow's wheel support for that Python version
before releasing a new pandas version.

## F.A.Q.

**Q: Why can't pandas just use numpy string and numpy void datatypes instead of pyarrow string and pyarrow struct?**

**A**: NumPy strings aren't yet available, whereas pyarrow strings are. NumPy void datatype would be different to pyarrow struct,
  not bringing the same interoperabitlity benefit with other arrow-based dataframe libraries.

**Q: Are all pyarrow dtypes ready? Isn't it too soon to make them the default?**

**A**: They will likely be ready by 3.0 - however, we're not making them the default (yet).
  For example, `pd.Series([1, 2, 3])` will continue to be auto-inferred to be
  `np.int64`.  We will only change the default for dtypes which currently have no `numpy`-backed equivalent and which are
  stored as `object` dtype, such as strings and nested datatypes.

### PDEP-10 History

- 17 April 2023: Initial version
- 8 May 2023: Changed proposal to make pyarrow required in pandas 3.0 instead of 2.1

[^1] <https://pandas.pydata.org/docs/development/roadmap.html#apache-arrow-interoperability>
[^2] <https://arrow.apache.org/powered_by/>
