# PDEP-9: Allow third-party projects to register pandas connectors with a standard API

- Created: 5 March 2023
- Status: Draft
- Discussion: [#XXXX](https://github.com/pandas-dev/pandas/pull/XXXX)
- Author: [Marc Garcia](https://github.com/datapythonista)
- Revision: 1

## PDEP Summary

This document proposes that third-party projects implementing I/O or memory
connectors, can register them using Python's entrypoint system, and make them
available to pandas users with a standard interface in a dedicated namespace
`DataFrame.io`. For example:

```python
import pandas

df = pandas.DataFrame.io.from_duckdb("SELECT * FROM 'dataset.parquet';")

df.io.to_hive(hive_conn, "hive_table")
```

## Current state

pandas supports importing and exporting data from different formats using
I/O connectors, currently implemented in `pandas/io`, as well as connectors
to in-memory structure, like Python structures or other library formats.
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

| Format       | Reader | Writer |
|--------------|--------|--------|
| CSV          | X      | X      |
| FWF          | X      |        |
| JSON         | X      | X      |
| HTML         | X      | X      |
| LaTeX        |        | X      |
| XML          | X      | X      |
| Clipboard    | X      | X      |
| Excel        | X      | X      |
| HDF5         | X      | X      |
| Feather      | X      | X      |
| Parquet      | X      | X      |
| ORC          | X      | X      |
| Stata        | X      | X      |
| SAS          | X      |        |
| SPSS         | X      |        |
| Pickle       | X      | X      |
| SQL          | X      | X      |
| BigQuery     |        |        |
| dict         | X      | X      |
| records      | X      | X      |
| string       |        | X      |
| markdown     |        | X      |
| xarray       |        | X      |

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
where included in the survey.

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
- There is no standard API for I/O connectors, and users of them need to learn
  each of them individually.
- Method chaining, is not possible with third-party I/O connectors to export
  data, unless authors monkey patch the `DataFrame` class, which should not
  be encouraged.

This document proposes to open the development of pandas I/O connectors to
third-party libraries in a standard way that overcomes those limitations.

### Proposal implementation

Implementing this proposal would not require major changes to pandas, and
the API defined next would be used.

A new `.io` accessor would be created for the `DataFrame` class, where all
I/O connector methods from third-parties would be loaded. Nothing else would
live under that namespace.

Third-party packages would implement a
[setuptools entrypoint](https://setuptools.pypa.io/en/latest/userguide/entry_point.html#entry-points-for-plugins)
to define the connectors that they implement, under a group `dataframe.io`.

For example, a hypothetical project `pandas_duckdb` implementing a `from_duckdb`
function, could use `pyproject.toml` to define the next entry point:

```toml
[project.entry-points."dataframe.io"]
from_duckdb = "pandas_duckdb:from_duckdb"
```

On import of the pandas module, it would read the entrypoint registry for the
`dataframe.io` group, and would dynamically create methods in the `DataFrame.io`
namespace for them. Method names would only be allowed to start with `from_` or
`to_`, and any other prefix would make pandas raise an exception. This would
guarantee a reasonably consistent API among third-party I/O connectors.

Connectors would use Apache Arrow as the only interface to load data from and
to pandas. This would prevent that changes to the pandas API affecting
connectors in any way. This would simplify the development of connectors,
make testing of them much more reliable and resistant to changes in pandas, and
allow connectors to be reused by other projects of the ecosystem.

In case a `from_` method returned something different than a PyArrow table,
pandas would raise an exception. pandas would expect all `to_` methods to have
`table: pyarrow.Table` as the first parameter, and it would raise an exception
otherwise. The `table` parameter would be exposed as the `self` parameter in
pandas, when the original function is registered as a method of the `.io`
accessor.

### Connector examples

This section lists specific examples of connectors that could immediately
benefit from this proposal.

**PyArrow** currently provides `Table.from_pandas` and `Table.to_pandas`.
With the new interface, it could also register `DataFrame.from_pyarrow`
and `DataFrame.to_pyarrow`, so pandas users can use the converters with
the interface they are used to, when PyArrow is installed in the environment.

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
more explicitly API.

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
The loader can already be implemented inside pandas, or as a third-party
extension with an arbitrary API. But this proposal would let the creation
of a third-party extension with a standard and intuitive API:

```python
pandas.DataFrame.io.from_duckdb("SELECT *
                                 FROM 'dataset.parquet'
                                 WHERE my_col > 0")
```

**Big data** systems such as Hive, Iceberg, Presto, etc. could benefit
from a standard way to load data to pandas. Also regular **SQL databases**
that can return their query results as Arrow, would benefit from better
and faster connectors than the existing ones based on SQL Alchemy and
Python structures.

Any other format, including **domain-specific formats** could easily
implement pandas connectors with a clear an intuitive API.

## Proposal extensions

The scope of the current proposal is limited to the addition of the
`DataFrame.io` namespace, and the automatic registration of functions defined
by third-party projects, if an entrypoint is defined.

Any changes to the current I/O of pandas are out of scope for this proposal,
but the next tasks can be considered for future work and proposals:

- Migrate I/O connectors currently implemented in pandas to the new interface.
  This would require a transition period where users would be warned that
  existing `DataFrame.read_*` may have been moved to `DataFrame.io.from_*`,
  and that the old API will stop working in a future version.
- Move out of the pandas repository and into their own third-party projects
  some of the existing I/O connectors.
- Implement with the new interface some of the data structures that the
  `DataFrame` constructor accepts.

## PDEP-9 History

- 5 March 2023: Initial version
