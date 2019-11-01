import io
from typing import Callable
from collections.abc import Mapping
from functools import wraps
import pandas as pd


def read_input_wrapper(read_func=None, **kwargs):
    """A decorator to make the `pandas.io.parser.read` functions
    take `bytes` as input.

    Parameters
    ----------
    `read_func` : `Callable[..., pd.DataFrame]`
        The `pandas.io.parsers` function to decorate.
        If not set `read_input_wrapper` will return a decorator.
    **`kwargs` : `dict`
        `kwargs` to pass on to `read_func`.

    Returns
    -------
    function : `Callable[input: bytes, pd.DataFrame]` |
               `Callable[[Callable[..., pd.DataFrame]],
                         Callable[input: bytes, pd.DataFrame]]`
        either return a decorator which will wrap a pandas parser function
        or a wrapped parser function:

    Examples
    -------
    Examples should be written in doctest format, and
    should illustrate how to use the function/class.
    >>> read_csv = read_input_wrapper(pd.read_csv)
    >>> read_tsv = read_input_wrapper(pd.read_csv, sep='\t')

    or as a decorator

        @read_input_wrapper
        def read_foo(file, **kwargs) -> pd.DataFrame:
            # some custom foo
            return pd.DataFrame()

    or

        @read_input_wrapper(sep='\t')
        def read_bar(file, **kwargs) -> pd.DataFrame:
            # some custom bar
            return pd.DataFrame()
    """

    def wrapper(func: Callable[..., pd.DataFrame]):

        @wraps(func)
        def reader(input: bytes) -> pd.DataFrame:
            return func(io.BytesIO(input), **kwargs)
        return reader

    if read_func is None:
        return wrapper
    return wrapper(read_func)


read_csv = read_input_wrapper(pd.read_csv)
read_tsv = read_input_wrapper(pd.read_csv, sep='\t')
read_excel = read_input_wrapper(pd.read_excel, sheet_name=None)


class ApplicationMapper(Mapping):
    """A `Mapping` class to map MIME application types to a pandas reader."""

    application_mapping = {
        "text/plain": read_tsv,
        "text/csv": read_csv,
        "application/vnd.ms-excel": read_excel
    }

    aplication_prefixed = (
      (
        'application/vnd.ms-excel.sheet',
        read_excel
      )
      (
        'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        read_excel
      )
    )

    def __getitem__(self, key):
        func = self.application_mapping.get(key)
        if func is not None:
            return func
        for prefix, func in self.aplication_prefixed:
            if key.startswith(prefix):
                return read_excel

    def __iter__(self):
        return iter(self.application_mapping)

    def __len__(self):
        return len(self.application_mapping)


application_mapping = ApplicationMapper()
