"""
Load I/O plugins from third-party libraries into the pandas namespace.

Third-party libraries defining I/O plugins register an entrypoint in
the `dataframe.io` group. For example:

```
[project.entry-points."dataframe.io"]
repr = "pandas_repr:ReprDataFrameIO"
```

The class `ReprDataFrameIO` will implement at least one of a reader
and a writer that supports the dataframe interchange protocol:

https://data-apis.org/dataframe-protocol/latest/API.html

For example:

```python
class ReprDataFrameIO:
    @staticmethod
    def reader(self, fname):
        with open(fname) as f:
            # for simplicity this assumes eval will create a DataFrame object
            return eval(f.read())

    def writer(self, fname, mode='w'):
        with open(fname, mode) as f:
            f.write(repr(self))
```

pandas will create wrapper functions or methods to call the reader or
writer from the pandas standard I/O namespaces. For example, for the
entrypoint above with name `repr` and both methods `reader` and
`writer` implemented, pandas will create the next functions and methods:

- `pandas.read_repr(...)`
- `pandas.Series.to_repr(...)`
- `pandas.DataFrame.to_repr(...)`

The reader wrappers make sure that the returned object is a pandas
DataFrame, since the user always expects the return of `read_*()`
to be a pandas DataFrame, not matter what the connector returns.
In few cases, the return can be a list or dict of dataframes, which
is supported.

If more than one reader or writer with the same name is loaded, pandas
raises an exception. For example, if two connectors use the name
`arrow` pandas will raise when `load_io_plugins()` is called, since
only one `pandas.read_arrow` function can exist, and pandas should not
make an arbitrary decision on which to use.
"""
import functools
from importlib.metadata import entry_points

import pandas as pd


def _create_reader_function(io_plugin):
    """
    Create and return a wrapper function for the original I/O reader.

    We can't directly call the original reader implemented in
    the connector, since the return of third-party connectors is not necessarily
    a pandas DataFrame but any object supporting the dataframe interchange
    protocol. We make sure here that `read_<whatever>` returns a pandas DataFrame.
    """

    # TODO: Create this function dynamically so the resulting signature contains
    # the original parameters and not `*args` and `**kwargs`
    @functools.wraps(io_plugin.reader)
    def reader_wrapper(*args, **kwargs):
        result = io_plugin.reader(*args, **kwargs)

        if isinstance(result, list):
            result = [pd.api.interchange.from_dataframe(df) for df in result]
        elif isinstance(result, dict):
            result = {k: pd.api.interchange.from_dataframe(df)
                      for k, df in result.items()}
        else:
            result = pd.api.interchange.from_dataframe(result)

        return result

    # TODO `function.wraps` changes the name of the wrapped function to the
    # original `pandas_reader`, change it to the function exposed in pandas.
    return reader_wrapper


def _create_series_writer_function(format_name):
    """
    When calling `Series.to_<whatever>` we call the dataframe writer, so
    we need to convert the Series to a one column dataframe.
    """
    def series_writer_wrapper(self, *args, **kwargs):
        dataframe_writer = getattr(self.to_frame(), f"to_{format_name}")
        dataframe_writer(*args, **kwargs)

    return series_writer_wrapper


def load_io_plugins():
    """
    Looks for entrypoints in the `dataframe.io` group and creates the
    corresponding pandas I/O methods.
    """
    for dataframe_io_entry_point in entry_points().get("dataframe.io", []):
        format_name = dataframe_io_entry_point.name
        io_plugin = dataframe_io_entry_point.load()

        if hasattr(io_plugin, "reader"):
            if hasattr(pd, f"read_{format_name}"):
                raise RuntimeError(
                    "More than one installed library provides the "
                    "`read_{format_name}` reader. Please uninstall one of "
                    "the I/O plugins providing connectors for this format."
                )
            setattr(
                pd,
                f"read_{format_name}",
                _create_reader_function(io_plugin),
            )

        if hasattr(io_plugin, "writer"):
            if hasattr(pd.DataFrame, f"to_{format_name}"):
                raise RuntimeError(
                    "More than one installed library provides the "
                    "`to_{format_name}` reader. Please uninstall one of "
                    "the I/O plugins providing connectors for this format."
                )
            setattr(
                pd.DataFrame,
                f"to_{format_name}",
                getattr(io_plugin, "writer"),
            )
            setattr(
                pd.Series,
                f"to_{format_name}",
                _create_series_writer_function(format_name),
            )
