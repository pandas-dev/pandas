"""
Load I/O plugins from third-party libraries into the pandas namespace.

Third-party libraries defining I/O plugins register an entrypoint in
the `dataframe.io` group. For example:

```
[project.entry-points."dataframe.io"]
repr = "pandas_repr:ReprDataFrameIO"
```

The class `ReprDataFrameIO` will implement readers and writers in at
least one of different exchange formats supported by the protocol.
For now a pandas DataFrame or a PyArrow table, in the future probably
nanoarrow, or a Python wrapper to the Arrow Rust implementations.
For example:

```python
class FancyFormatDataFrameIO:
    @staticmethod
    def pandas_reader(self, fname):
        with open(fname) as f:
            return eval(f.read())

    def pandas_writer(self, fname, mode='w'):
        with open(fname, mode) as f:
            f.write(repr(self))
```

If the I/O plugin implements a reader or writer supported by pandas,
pandas will create a wrapper function or method to call the reader or
writer from the pandas standard I/O namespaces. For example, for the
entrypoint above with name `repr` and methods `pandas_reader` and
`pandas_writer` pandas will create the next functions and methods:

- `pandas.read_repr(...)`
- `pandas.Series.to_repr(...)`
- `pandas.DataFrame.to_repr(...)`

The reader wrappers validates that the returned object is a pandas
DataFrame when the exchange format is `pandas`, and will convert the
other supported objects (e.g. a PyArrow Table) to a pandas DataFrame,
so the result of `pandas.read_repr` is a pandas DataFrame, as the
user would expect.

If more than one reader or writer with the same name is loaded, pandas
raises an exception.
"""
import functools
from importlib.metadata import entry_points

import pandas as pd

supported_exchange_formats = ["pandas"]

try:
    import pyarrow as pa
except ImportError:
    pa = None
else:
    supported_exchange_formats.append("pyarrow")


def _create_reader_function(io_plugin, exchange_format):
    original_reader = getattr(io_plugin, f"{exchange_format}_reader")

    # TODO: Create this function dynamically so the resulting signature contains
    # the original parameters and not `*args` and `**kwargs`
    @functools.wraps(original_reader)
    def reader_wrapper(*args, **kwargs):
        result = original_reader(*args, **kwargs)
        if exchange_format == "pyarrow":
            result = result.to_pandas()

        # validate output type
        if isinstance(result, list):
            assert all((isinstance(item, pd.DataFrame) for item in result))
        elif isinstance(result, dict):
            assert all(
                (
                    (isinstance(k, str) and isinstance(v, pd.DataFrame))
                    for k, v in result.items()
                )
            )
        elif not isinstance(result, pd.DataFrame):
            raise AssertionError("Returned object is not a DataFrame")
        return result

    return reader_wrapper


def _create_series_writer_function(format_name):
    def series_writer_wrapper(self, *args, **kwargs):
        dataframe_writer = getattr(self.to_frame(), f"to_{format_name}")
        dataframe_writer(*args, **kwargs)

    return series_writer_wrapper


def load_io_plugins():
    for dataframe_io_entry_point in entry_points().get("dataframe.io", []):
        format_name = dataframe_io_entry_point.name
        io_plugin = dataframe_io_entry_point.load()

        for exchange_format in supported_exchange_formats:
            if hasattr(io_plugin, f"{exchange_format}_reader"):
                if hasattr(pd, f"read_{format_name}"):
                    raise RuntimeError(
                        "More than one installed library provides the "
                        "`read_{format_name}` reader. Please uninstall one of "
                        "the I/O plugins to be able to load the pandas I/O plugins."
                    )
                setattr(
                    pd,
                    f"read_{format_name}",
                    _create_reader_function(io_plugin, exchange_format),
                )

            if hasattr(io_plugin, f"{exchange_format}_writer"):
                setattr(
                    pd.DataFrame,
                    f"to_{format_name}",
                    getattr(io_plugin, f"{exchange_format}_writer"),
                )
                setattr(
                    pd.Series,
                    f"to_{format_name}",
                    _create_series_writer_function(format_name),
                )
