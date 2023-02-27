# PDEP-9: Implement pandas I/O connectors as extensions

- Created: 26 February 2023
- Status: Draft
- Discussion: [#XXXX](https://github.com/pandas-dev/pandas/pull/XXXX)
- Author: [Marc Garcia](https://github.com/datapythonista)
- Revision: 1

## Introduction

pandas supports importing and exporting data from different formats using
connectors, currently implemented in `pandas/io`. In many cases, those
connectors wrap an existing Python library, while in some others, pandas
implements the format logic.

In some cases, different engines exist for the same format. The API to use
those connectors is `pandas.read_<format>(engine='<engine-name>', ...)` to
import data, and `DataFrame.to_<format>(engine='<engine-name>', ...)` to
export data.

For objects exported to memory (like a Python dict) the API is the same as
for I/O, `DataFrame.to_<format>(...)`. For formats imported from objects in
memory, the API is different, `DataFrame.from_<format>(...)`.

In some cases, the pandas API provides `DataFrame.to_*` methods that are not
used to export the data to a disk or memory object, but instead to transform
the index of a `DataFrame`: `DataFrame.to_period` and `DataFrame.to_timestamp`.

Dependencies of the I/O connectors are not loaded by default, and will be
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

The list of formats can be found in the [IO guide](https://pandas.pydata.org/docs/dev/user_guide/io.html).
A more detailed table, including in memory objects, with
the engines and dependencies is presented next.

| Format       | Reader | Writer | Engines                  | Dependencies |
|--------------|--------|--------|--------------------------|--------------|
| CSV          | X      | X      | `c`, `python`, `pyarrow` | `pyarrow`    |
| FWF          | X      |        |                          |              |
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
| Pickle       | X      | X      |                          |               |
| SQL          | X      | X      |
| BigQuery     |        |        |
| dict         | X      | X      |
| records      | X      | X      |
| string       |        | X      |
| markdown     |        | X      |
| xarray       |        | X      |

### Inclusion criteria

There is no objective criteria for when a format is included
in pandas, and the list above is mostly the result of developers
being interested in implementing the connectors for a certain
format in pandas.

The number of existing formats is constantly increasing, and its
difficult for pandas to keep up to date even with popular formats.
It could probably make sense to have connectors to pyarrow,
pyspark, Iceberg, DuckDB, Polars, and others.

At the same time, some of the formats are not frequently used as
shown in the [2019 user survey](https://pandas.pydata.org//community/blog/2019-user-survey.html).
Those less popular formats include SPSS, SAS, Google BigQuery and
Stata. Note that only I/O formats (and not memory formats like
records or xarray) where included in the survey.

## Proposal

The main proposal in this PDEP is to open the development of pandas
connectors to third-parties. This would not only allow the development
of new connectors in a faster and easier way, without the intervention of
the pandas team, but also remove from the pandas code base a number of the
existing connectors, simplifying the code, the CI and the builds.
While a limited set of core connectors could live in the pandas code base,
most of the existing connectors would be moved to third-party projects.

The user experience would remain similar to the existing one, but making
better use of namespaces, and adding consistency. Any pandas connector
(regardless of being implemented as a third-party module or not) would define
a Python entrypoint specifying the format they connect to, the operations
they support (read and/or write) and the name of the engine to be used.
On load, pandas would access this registry of connectors, and would create
the corresponding import and export methods.

To use the connectors for the format, users would install the third-party
connector package, instead of installing the required dependencies as they
need to do now.

### Python API

The Python API can be improved from the current one to make better use
of namespaces, and avoid inconsistencies. The proposed API is:

```python
import pandas

df = pandas.DataFrame.io.read_<format>(engine='<engine>', ...)

df.io.write_<format>(engine='<engine>', ...)
```
The `engine` parameter would only be required when more than an engine
is available for a format. This is similar to the the current API, that
would use the default engine if not specified.

For example:

```python
import pandas

df = pandas.DataFrame.io.read_hdf5('input.hdf5')

df.io.write_parquet('output.parquet')
```

All the I/O connectors would be accessed via `DataFrame.io`, significantly
reducing the number of items in the namespace of the `pandas` module, and
the `DataFrame` class. Introspection would make it fast and simple to
list the existing connectors `dir(pandas.DataFrame.io)`.

The API is more intuitive than the current one, as it would be used for
both in memory formats and disk formats, and does not mix read/to (users
in general would expect read/write, from/to, import/export, input/output,
and not a mix of those pairs).

### Ecosystem of connectors

In the same way Python can be extended with third-party modules, pandas
would be extendable with I/O plugins. This has some advantages:

- **Supression of the pandas maintainers bottleneck.** Everybody would be
  able to develop and promote their own I/O connectors, without the
  approval or intervention of pandas maintainers.
- **Lower the entry barrier to pandas code.** Since pandas is a huge and
  mature project, writing code in pandas itself is complex. Several
  linters and autoformatters are required, policies like adding release
  notes need to be followed. Proper testing must be implemented.
  CI is slow and takes hours to complete. pandas needs to be compiled
  due to its C extensions. All those would not be necessary, and
  creating new I/O connectors would be faster and simpler.
- **CI and packaging simplification.** pandas has currently around 20
  dependencies required by connectors. And a significant number of
  tests, some of them requiring a high level of customization (such as
  an available database server to test `read_sql`, or a virtual
  clipboard to test `read_clipboard`). Moving connectors out of
  pandas would make the CI faster, and the number of problems caused
  by updates in dependencies smaller.
- **Competition and alternatives for I/O operations.** Some of the
  supported formats allow for different approaches in terms of
  implementation. For example, `csv` connectors can be optimized
  for performance and reliability, or for easiness of use. When
  building a production pipeline, users would often appreciate a
  loader that requires an expected schema, loads faster because of
  it, and fails if the file contains errors. While Jupyter users
  may prefer inference and magic that helps them write code faster.
- **Reusability with other projects.** In some cases, it can make
  sense to load a format into for example Apache Arrow, and then
  convert it to a pandas `DataFrame` in the connector. It could
  also be quite simple when that is implemented to return a Vaex
  or a Polars object. Having connectors as third-party packages
  would allow to implement this, as opposed as our current
  connectors. This reusability would not only benefit other
  dataframe projects, but it would also have better maintained
  connectors, as they will be shared by a larger ecosystem.

## Disadvantages

The main disadvantages to implementing this PDEP are:

- **Backward compatibility**.
- **More verbose API.**
- **Fragmented documentation.**

## Transition period

This proposal involves some important changes regarding user
facing code.

The implementation of connectors as third-party packages is quite
small for users, who would just need to install `pandas-xarray`
instead of `xarray` to be able to use `DataFrame.to_xarray`. Also,
the `ImportError` message users would get in case it was not
properly installed, can provide the required information for users
to install the right package without issues.

The part that requires more careful management and a long transition
period is the change to the Python API proposed here. The
new API does not overlap with the old one (everything would be in
the new `DataFrame.io` accessor). This allows to easily implement
both the new and old API in parallel, raising `FutureWarning`
warnings in the old API, so users can slowly adapt their code,
and get used to the new API. Since the changes affect all pandas
users, keeping the old behavior until at least pandas 4.0 seems
a reasonable transition period.

## PDEP-9 History

- 26 February 2023: Initial version
