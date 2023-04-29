# PDEP-9: Allow third-party projects to register pandas connectors with a standard API

- Created: 5 March 2023
- Status: Rejected
- Discussion: [#51799](https://github.com/pandas-dev/pandas/pull/51799)
- Author: [Marc Garcia](https://github.com/datapythonista)
- Revision: 1

## PDEP Summary

This document proposes that third-party projects implementing I/O or memory
connectors can register them using Python's entrypoint system, and make them
available to pandas users with the existing I/O interface. For example:

```python
import pandas

df = pandas.DataFrame.read_duckdb("SELECT * FROM 'my_dataset.parquet';")

df.to_deltalake('/delta/my_dataset')
```

This would allow to easily extend the existing number of connectors, adding
support to new formats and database engines, data lake technologies,
out-of-core connectors, the new ADBC interface and at the same time reduce the
maintenance cost of the pandas core.

## Current state

pandas supports importing and exporting data from different formats using
I/O connectors, currently implemented in `pandas/io`, as well as connectors
to in-memory structures, like Python structures or other library formats.
In many cases, those connectors wrap an existing Python library, while in
some others, pandas implements the logic to read and write to a particular
format.

In some cases, different engines exist for the same format. The API to use
those connectors is `pandas.read_<format>(engine='<engine-name>', ...)` to
import data, and `DataFrame.to_<format>(engine='<engine-name>', ...)` to
export data.

For objects exported to memory (like a Python dict) the API is the same as
for I/O, `DataFrame.to_<format>(...)`. For formats imported from objects in
memory, the API is different using the `from_` prefix instead of `read_`,
`DataFrame.from_<format>(...)`.

In some cases, the pandas API provides `DataFrame.to_*` methods that are not
used to export the data to a disk or memory object, but instead to transform
the index of a `DataFrame`: `DataFrame.to_period` and `DataFrame.to_timestamp`.

Dependencies of the connectors are not loaded by default, and will be
imported when the connector is used. If the dependencies are not installed,
an `ImportError` is raised.

```python
>>> pandas.read_gbq(query)
Traceback (most recent call last):
  ...
ImportError: Missing optional dependency 'pandas-gbq'.
pandas-gbq is required to load data from Google BigQuery.
See the docs: https://pandas-gbq.readthedocs.io.
Use pip or conda to install pandas-gbq.
```

### Supported formats

The list of formats can be found in the
[IO guide](https://pandas.pydata.org/docs/dev/user_guide/io.html).
A more detailed table, including in memory objects, and I/O connectors in the
DataFrame styler is presented next:

| Format       | Reader | Writer | Engines                                                                           |
|--------------|--------|--------|-----------------------------------------------------------------------------------|
| CSV          | X      | X      | `c`, `python`, `pyarrow`                                                          |
| FWF          | X      |        | `c`, `python`, `pyarrow`                                                          |
| JSON         | X      | X      | `ujson`, `pyarrow`                                                                |
| HTML         | X      | X      | `lxml`, `bs4/html5lib` (parameter `flavor`)                                       |
| LaTeX        |        | X      |                                                                                   |
| XML          | X      | X      | `lxml`, `etree` (parameter `parser`)                                              |
| Clipboard    | X      | X      |                                                                                   |
| Excel        | X      | X      | `xlrd`, `openpyxl`, `odf`, `pyxlsb` (each engine supports different file formats) |
| HDF5         | X      | X      |                                                                                   |
| Feather      | X      | X      |                                                                                   |
| Parquet      | X      | X      | `pyarrow`, `fastparquet`                                                          |
| ORC          | X      | X      |                                                                                   |
| Stata        | X      | X      |                                                                                   |
| SAS          | X      |        |                                                                                   |
| SPSS         | X      |        |                                                                                   |
| Pickle       | X      | X      |                                                                                   |
| SQL          | X      | X      | `sqlalchemy`, `dbapi2` (inferred from the type of the `con` parameter)            |
| BigQuery     | X      | X      |                                                                                   |
| dict         | X      | X      |                                                                                   |
| records      | X      | X      |                                                                                   |
| string       |        | X      |                                                                                   |
| markdown     |        | X      |                                                                                   |
| xarray       |        | X      |                                                                                   |

At the time of writing this document, the `io/` module contains
close to 100,000 lines of Python, C and Cython code.

There is no objective criteria for when a format is included
in pandas, and the list above is mostly the result of a developer
being interested in implementing the connectors for a certain
format in pandas.

The number of existing formats available for data that can be processed with
pandas is constantly increasing, and its difficult for pandas to keep up to
date even with popular formats. It could possibly make sense to have connectors
to PyArrow, PySpark, Iceberg, DuckDB, Hive, Polars, and many others.

At the same time, some of the formats are not frequently used as shown in the
[2019 user survey](https://pandas.pydata.org//community/blog/2019-user-survey.html).
Those less popular formats include SPSS, SAS, Google BigQuery and
Stata. Note that only I/O formats (and not memory formats like records or xarray)
were included in the survey.

The maintenance cost of supporting all formats is not only in maintaining the
code and reviewing pull requests, but also it has a significant cost in time
spent on CI systems installing dependencies, compiling code, running tests, etc.

In some cases, the main maintainers of some of the connectors are not part of
the pandas core development team, but people specialized in one of the formats
without commit rights.

## Proposal

While the current pandas approach has worked reasonably well, it is difficult
to find a stable solution where the maintenance incurred in pandas is not
too big, while at the same time users can interact with all different formats
and representations they are interested in, in an easy and intuitive way.

Third-party packages are already able to implement connectors to pandas, but
there are some limitations to it:

- Given the large number of formats supported by pandas itself, third-party
  connectors are likely seen as second class citizens, not important enough
  to be used, or not well supported.
- There is no standard API for external I/O connectors, and users need
  to learn each of them individually.
- Method chaining is not possible with third-party I/O connectors to export
  data, unless authors monkey patch the `DataFrame` class, which should not
  be encouraged.

This document proposes to open the development of pandas I/O connectors to
third-party libraries in a standard way that overcomes those limitations.

### Proposal implementation

Implementing this proposal would not require major changes to pandas, and
the API defined next would be used.

#### User API

Users will be able to install third-party packages implementing pandas
connectors using the standard packaging tools (pip, conda, etc.). These
connectors should implement an entrypoint that pandas will use to
automatically create the corresponding methods `pandas.read_*` and
`pandas.DataFrame.to_*`. Arbitrary function or method names will not
be created by this interface, only the `read_*` and `to_*` pattern will
be allowed. By simply installing the appropriate packages users will
be able to use code like this:

```python
import pandas

df = pandas.read_duckdb("SELECT * FROM 'dataset.parquet';")

df.to_hive(hive_conn, "hive_table")
```

This API allows for method chaining:

```python
(pandas.read_duckdb("SELECT * FROM 'dataset.parquet';")
       .to_hive(hive_conn, "hive_table"))
```

The total number of I/O functions and methods is expected to be small, as users
in general use only a small subset of formats. The number could actually be
reduced from the current state if the less popular formats (such as SAS, SPSS,
BigQuery, etc.) are removed from the pandas core into third-party packages.
Moving these connectors is not part of this proposal, and will be discussed
later in a separate proposal.

#### Plugin registration

Third-party packages would implement an
[entrypoint](https://setuptools.pypa.io/en/latest/userguide/entry_point.html#entry-points-for-plugins)
to define the connectors that they implement, under a group `dataframe.io`.

For example, a hypothetical project `pandas_duckdb` implementing a `read_duckdb`
function, could use `pyproject.toml` to define the next entry point:

```toml
[project.entry-points."dataframe.io"]
from_duckdb = "pandas_duckdb:read_duckdb"
```

On import of the pandas module, it would read the entrypoint registry for the
`dataframe.io` group, and would dynamically create methods in the `pandas` and
`pandas.DataFrame` namespace for them. Method not starting with `read_` or `to_`
would make pandas raise an exception. This would guarantee a reasonably
consistent API among third-party connectors.

#### Internal API

Connectors would use one of two different interface options: an Apache Arrow table
or a pandas `DataFrame`.

The Apache Arrow format would allow that connectors do not need to use pandas to
create the data, making them more robust and less likely to break when changes to
pandas internals happen. It would also allow other possible consumers of the
connectors to not have pandas as a dependency. Testing also becomes simpler by
using Apache Arrow, since connectors can be tested independently, and pandas does
not need to be tested for each connector. If the Apache Arrow specification is
respected in both sides, the communication between connectors and pandas is
guaranteed to work. If pandas eventually has a hard dependency on an Apache
Arrow implementation, this should be the preferred interface.

Allowing connectors to use pandas dataframes directly makes users not have to
depend on PyArrow for connectors that do not use an Apache Arrow object. It also
helps move existing connectors to this new API, since they are using pandas
dataframes as an exchange object now. It has the disadvantages stated in the
previous paragraph, and a future proposal may be created to discuss deprecating
pandas dataframes as a possible interface for connectors.

In case a `read` method returned something different than a PyArrow table,
pandas would raise an exception. pandas would expect all `to_` methods to have
`table: pyarrow.Table | pandas.DataFrame` as the first parameter, and it would
raise an exception otherwise. The `table` parameter would be exposed as the
`self` parameter of the `to_*` method in pandas.

In case the Apache Arrow interface is used, metadata not supported by Apache
Arrow may be provided by users. For example the column to use for row indices
or the data type backend to use in the object being created. This would be
managed independently from the connectors. Given the previous example, a new
argument `index_col` could be added directly into pandas to the function or
method automatically generated from the entrypoint. Since this would apply to
all functions and methods automatically generated, it would also improve the
consistency of pandas connectors. For example:

```python
def read_myformat(index_col=None, *args, **kwargs):
    # The third-party connector doesn't need to know about functionality
    # specific to pandas like the row index
    arrow_table = third_party_connector.from_myformat(*args, **kwargs)

    df = convert_arrow_table_to_pandas_dataframe(arrow_table)

    # Transformations to the dataframe with custom parameters is possible
    if index_col is not None:
        df = df.set_index(index_col)

    return df
```

#### Connector guidelines

In order to provide a better and more consistent experience to users, guidelines
will be created to unify terminology and behavior. Some of the topics to unify are
defined next.

**Existence and naming of columns**, since many connectors are likely to provide
similar features, like loading only a subset of columns in the data, or dealing
with file names. Examples of recommendations to connector developers:

- `columns`: Use this argument to let the user load a subset of columns. Allow a
  list or tuple.
- `path`: Use this argument if the dataset is a file in the file disk. Allow a string,
  a `pathlib.Path` object, or a file descriptor. For a string object, allow URLs that
  will be automatically download, compressed files that will be automatically
  uncompressed, etc. A library can be provided to deal with those in an easier and
  more consistent way.
- `schema`: For datasets that don't have a schema (e.g. `csv`), allow providing an
  Apache Arrow schema instance, and automatically infer types if not provided.

Note that the above are only examples of guidelines for illustration, and not
a proposal of the guidelines.

**Guidelines to avoid name conflicts**. Since it is expected that more than one
implementation exists for certain formats, as it already happens, guidelines on
how to name connectors would be created. The easiest approach is probably to use
the format `from_<format>_<implementation-id>` / `to_<format>_<implementation-id>`.

For example a `csv` loader based on PyArrow could be named as `from_csv_pyarrow`,
and an implementation that does not infer types and raises an exception in case
of mistake or ambiguity it could be named `from_csv_strict`. Exact guidelines
would be developed independently from this proposal.

**Connector registry and documentation**. To simplify the discovery of connectors
and its documentation, connector developers can be encourage to register their
projects in a central location, and to use a standard structure for documentation.
This would allow the creation of a unified website to find the available
connectors, and their documentation. It would also allow to customize the
documentation for specific implementations, and include their final API, and
include arguments specific to the implementation. In the case of pandas, it
would allow to add arguments such as `index_col` to all loader methods, and
to potentially build the API reference of certain third-party connectors as part
as the pandas own documentation. That may or may not be a good idea, but
standardizing the documentation of connectors would allow it.

### Connector examples

This section lists specific examples of connectors that could immediately
benefit from this proposal.

**PyArrow** currently provides `Table.from_pandas` and `Table.to_pandas`.
With the new interface, it could also register `DataFrame.from_pyarrow`
and `DataFrame.to_pyarrow`, so pandas users can use the converters with
the interface they are used to, when PyArrow is installed in the environment.
Better integration with PyArrow tables was discussed in
[#51760](https://github.com/pandas-dev/pandas/issues/51760).

_Current API_:

```python
pyarrow.Table.from_pandas(table.to_pandas()
                               .query('my_col > 0'))
```

_Proposed API_:

```python
(pandas.DataFrame.io.from_pyarrow(table)
                 .query('my_col > 0')
                 .io.to_pyarrow())
```

**Polars**, **Vaex** and other dataframe frameworks could benefit from
third-party projects that make the interoperability with pandas use a
more explicitly API. Integration with Polars was requested in
[#47368](https://github.com/pandas-dev/pandas/issues/47368).

_Current API_:

```python
polars.DataFrame(df.to_pandas()
                   .query('my_col > 0'))
```

_Proposed API_:

```python
(pandas.DataFrame.io.from_polars(df)
                 .query('my_col > 0')
                 .io.to_polars())
```

**DuckDB** provides an out-of-core engine able to push predicates before
the data is loaded, making much better use of memory and significantly
decreasing loading time. pandas, because of its eager nature is not able
to easily implement this itself, but could benefit from a DuckDB loader.
The loader can already be implemented inside pandas (it has already been
proposed in [#45678](https://github.com/pandas-dev/pandas/issues/45678),
or as a third-party extension with an arbitrary API. But this proposal would
let the creation of a third-party extension with a standard and intuitive API:

```python
pandas.DataFrame.io.from_duckdb("SELECT *
                                 FROM 'dataset.parquet'
                                 WHERE my_col > 0")
```

**Out-of-core algorithms** push some operations like filtering or grouping
to the loading of the data. While this is not currently possible, connectors
implementing out-of-core algorithms could be developed using this interface.

**Big data** systems such as Hive, Iceberg, Presto, etc. could benefit
from a standard way to load data to pandas. Also regular **SQL databases**
that can return their query results as Arrow, would benefit from better
and faster connectors than the existing ones based on SQL Alchemy and
Python structures.

Any other format, including **domain-specific formats** could easily
implement pandas connectors with a clear and intuitive API.

## Proposal extensions

The scope of the current proposal is limited to the registration of functions
defined by third-party projects, if an entrypoint is defined.

Any changes to the current connectors of pandas (e.g. `read_csv`,
`from_records`, etc.) or their migration to the new system are out of scope for
this proposal, but the next tasks can be considered for future work and proposals:

- Move out of the pandas repository and into their own third-party projects
  some of the existing I/O connectors. This would require a transition period
  to let users know that future versions of pandas will require a dependency
  installed for a particular connector to exist.
- Implement with the new interface some of the data structures that the
  `DataFrame` constructor accepts.

## PDEP-9 History

- 5 March 2023: Initial version
