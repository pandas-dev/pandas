# PDEP-9: Allow third-party projects to register pandas connectors with a standard API

- Created: 5 March 2023
- Status: Rejected
- Discussion: [#51799](https://github.com/pandas-dev/pandas/pull/51799)
              [#53005](https://github.com/pandas-dev/pandas/pull/53005)
- Author: [Marc Garcia](https://github.com/datapythonista)
- Revision: 1

[TOC]

## PDEP Summary

This document proposes that third-party projects implementing I/O or memory
connectors to pandas can register them using Python's entrypoint system,
and make them available to pandas users with the usual pandas I/O interface.
For example, packages independent from pandas could implement readers from
DuckDB and writers to Delta Lake, and when installed in the user environment
the user would be able to use them as if they were implemented in pandas.
For example:

```python
import pandas

pandas.load_io_plugins()

df = pandas.DataFrame.read_duckdb("SELECT * FROM 'my_dataset.parquet';")

df.to_deltalake('/delta/my_dataset')
```

This would allow to easily extend the existing number of connectors, adding
support to new formats and database engines, data lake technologies,
out-of-core connectors, the new ADBC interface, and others, and at the
same time reduce the maintenance cost of the pandas code base.

## Current state

pandas supports importing and exporting data from different formats using
I/O connectors, currently implemented in `pandas/io`, as well as connectors
to in-memory structures like Python structures or other library formats.
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

Dependencies of the connectors are not loaded by default, and are
imported when the connector is used. If the dependencies are not installed
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
date even with popular formats. It possibly makes sense to have connectors
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
the pandas core development team, but people specialized in one of the formats.

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
  to learn each of them individually. Since the pandas I/O API is inconsistent
  by using read/to instead of read/write or from/to, developers in many cases
  ignore the convention. Also, even if developers follow the pandas convention
  the namespaces would be different, since developers of connectors will rarely
  monkeypatch their functions into the `pandas` or `DataFrame` namespace.
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
connectors should implement entrypoints that pandas will use to
automatically create the corresponding methods `pandas.read_*`,
`pandas.DataFrame.to_*` and `pandas.Series.to_*`. Arbitrary function or
method names will not be created by this interface, only the `read_*`
and `to_*` pattern will be allowed.

By simply installing the appropriate packages and calling the function
`pandas.load_io_plugins()` users will be able to use code like this:

```python
import pandas

pandas.load_io_plugins()

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
Moving these connectors is not part of this proposal, and could be discussed
later in a separate proposal.

#### Plugin registration

Third-party packages would implement
[entrypoints](https://setuptools.pypa.io/en/latest/userguide/entry_point.html#entry-points-for-plugins)
to define the connectors that they implement, under a group `dataframe.io`.

For example, a hypothetical project `pandas_duckdb` implementing a `read_duckdb`
function, could use `pyproject.toml` to define the next entry point:

```toml
[project.entry-points."dataframe.io"]
reader_duckdb = "pandas_duckdb:read_duckdb"
```

When the user calls `pandas.load_io_plugins()`, it would read the entrypoint registry for the
`dataframe.io` group, and would dynamically create methods in the `pandas`,
`pandas.DataFrame` and `pandas.Series` namespaces for them. Only entrypoints with
name starting by `reader_` or `writer_` would be processed by pandas, and the functions
registered in the entrypoint would be made available to pandas users in the corresponding
pandas namespaces. The text after the keywords `reader_` and `writer_` would be used
for the name of the function. In the example above, the entrypoint name `reader_duckdb`
would create `pandas.read_duckdb`. An entrypoint with name `writer_hive` would create
the methods `DataFrame.to_hive` and `Series.to_hive`.

Entrypoints not starting with `reader_` or `writer_` would be ignored by this interface,
but will not raise an exception since they can be used for future extensions of this
API, or other related dataframe I/O interfaces.

#### Internal API

Connectors will use the dataframe interchange API to provide data to pandas. When
data is read from a connector, and before returning it to the user as a response
to `pandas.read_<format>`, data will be parsed from the data interchange interface
and converted to a pandas DataFrame. In practice, connectors are likely to return
a pandas DataFrame or a PyArrow Table, but the interface will support any object
implementing the dataframe interchange API.

#### Connector guidelines

In order to provide a better and more consistent experience to users, guidelines
will be created to unify terminology and behavior. Some of the topics to unify are
defined next.

**Guidelines to avoid name conflicts**. Since it is expected that more than one
implementation exists for certain formats, as it already happens, guidelines on
how to name connectors would be created. The easiest approach is probably to use
as the format a string of the type `to_<format>_<implementation-id>` if it is
expected that more than one connector can exist. For example, for LanceDB it is likely
that only one connector exist, and the name `lance` can be used (which would create
`pandas.read_lance` or `DataFrame.to_lance`. But if a new `csv` reader based in the
Arrow2 Rust implementation, the guidelines can recommend to use `csv_arrow2` to
create `pandas.read_csv_arrow2`, etc.

**Existence and naming of parameters**, since many connectors are likely to provide
similar features, like loading only a subset of columns in the data, or dealing
with paths. Examples of recommendations to connector developers could be:

- `columns`: Use this argument to let the user load a subset of columns. Allow a
  list or tuple.
- `path`: Use this argument if the dataset is a file in the file disk. Allow a string,
  a `pathlib.Path` object, or a file descriptor. For a string object, allow URLs that
  will be automatically download, compressed files that will be automatically
  uncompressed, etc. Specific libraries can be recommended to deal with those in an
  easier and more consistent way.
- `schema`: For datasets that don't have a schema (e.g. `csv`), allow providing an
  Apache Arrow schema instance, and automatically infer types if not provided.

Note that the above are only examples of guidelines for illustration, and not
a proposal of the guidelines, which would be developed independently after this
PDEP is approved.

**Connector registry and documentation**. To simplify the discovery of connectors
and its documentation, connector developers can be encourage to register their
projects in a central location, and to use a standard structure for documentation.
This would allow the creation of a unified website to find the available
connectors, and their documentation. It would also allow to customize the
documentation for specific implementations, and include their final API.

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
(pandas.read_pyarrow(table)
       .query('my_col > 0')
       .to_pyarrow())
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
(pandas.read_polars(df)
       .query('my_col > 0')
       .to_polars())
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
pandas.read_duckdb("SELECT *
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

### Limitations

The implementation of this proposal has some limitations discussed here:

- **Lack of support for multiple engines.** The current pandas I/O API
  supports multiple engines for the same format (for the same function or
  method name). For example `read_csv(engine='pyarrow', ...)`. Supporting
  engines requires that all engines for a particular format use the same
  signature (the same parameters), which is not ideal. Different connectors
  are likely to have different parameters and using `*args` and `**kwargs`
  provides users with a more complex and difficult experience. For this
  reason this proposal prefers that function and method names are unique
  instead of supporting an option for engines.
- **Lack of support for type checking of connectors.** This PDEP proposes
  creating functions and methods dynamically, and those are not supported
  for type checking using stubs. This is already the case for other
  dynamically created components of pandas, such as custom accessors.
- **No improvements to the current I/O API**. In the discussions of this
  proposal it has been considered to improve the current pandas I/O API to
  fix the inconsistency of using `read` / `to` (instead of for example
  `read` / `write`), avoid using `to_` prefixed methods for non-I/O
  operations, or using a dedicated namespace (e.g. `DataFrame.io`) for
  the connectors. All of these changes are out of scope for this PDEP.

## Future plans

This PDEP is exclusively to support a better API for existing of future
connectors. It is out of scope for this PDEP to implement changes to any
connectors existing in the pandas code base.

Some ideas for future discussion related to this PDEP include:

- Automatically loading of I/O plugins when pandas is imported.

- Removing from the pandas code base some of the least frequently used connectors,
such as SAS, SPSS or Google BigQuery, and move them to third-party connectors
registered with this interface.

- Discussing a better API for pandas connectors. For example, using `read_*`
methods instead of `from_*` methods, renaming `to_*` methods not used as I/O
connectors, using a consistent terminology like from/to, read/write, load/dump, etc.
or using a dedicated namespace for connectors (e.g. `pandas.io` instead of the
general `pandas` namespace).

- Implement as I/O connectors some of the formats supported by the `DataFrame`
constructor.

## PDEP-9 History

- 5 March 2023: Initial version
- 30 May 2023: Major refactoring to use the pandas existing API,
  the dataframe interchange API and to make the user be explicit to load
  the plugins
- 13 June 2023: The PDEP did not get any support after several iterations,
  and its been closed as rejected by the author
