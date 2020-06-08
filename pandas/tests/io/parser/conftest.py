import os
from typing import List, Optional

import pytest

from pandas import read_csv, read_table


class BaseParser:
    engine: Optional[str] = None
    low_memory = True
    float_precision_choices: List[Optional[str]] = []

    def update_kwargs(self, kwargs):
        kwargs = kwargs.copy()
        kwargs.update(dict(engine=self.engine, low_memory=self.low_memory))

        return kwargs

    def read_csv(self, *args, **kwargs):
        kwargs = self.update_kwargs(kwargs)
        return read_csv(*args, **kwargs)

    def read_table(self, *args, **kwargs):
        kwargs = self.update_kwargs(kwargs)
        return read_table(*args, **kwargs)


class CParser(BaseParser):
    engine = "c"
    float_precision_choices = [None, "high", "round_trip"]


class CParserHighMemory(CParser):
    low_memory = False


class CParserLowMemory(CParser):
    low_memory = True


class PythonParser(BaseParser):
    engine = "python"
    float_precision_choices = [None]


@pytest.fixture
def csv_dir_path(datapath):
    """
    The directory path to the data files needed for parser tests.
    """
    return datapath("io", "parser", "data")


@pytest.fixture
def csv1(datapath):
    """
    The path to the data file "test1.csv" needed for parser tests.
    """
    return os.path.join(datapath("io", "data", "csv"), "test1.csv")


_cParserHighMemory = CParserHighMemory()
_cParserLowMemory = CParserLowMemory()
_pythonParser = PythonParser()

_py_parsers_only = [_pythonParser]
_c_parsers_only = [_cParserHighMemory, _cParserLowMemory]
_all_parsers = [*_c_parsers_only, *_py_parsers_only]

_py_parser_ids = ["python"]
_c_parser_ids = ["c_high", "c_low"]
_all_parser_ids = [*_c_parser_ids, *_py_parser_ids]


@pytest.fixture(params=_all_parsers, ids=_all_parser_ids)
def all_parsers(request):
    """
    Fixture all of the CSV parsers.
    """
    return request.param


@pytest.fixture(params=_c_parsers_only, ids=_c_parser_ids)
def c_parser_only(request):
    """
    Fixture all of the CSV parsers using the C engine.
    """
    return request.param


@pytest.fixture(params=_py_parsers_only, ids=_py_parser_ids)
def python_parser_only(request):
    """
    Fixture all of the CSV parsers using the Python engine.
    """
    return request.param


_utf_values = [8, 16, 32]

_encoding_seps = ["", "-", "_"]
_encoding_prefixes = ["utf", "UTF"]

_encoding_fmts = [
    f"{prefix}{sep}" + "{0}" for sep in _encoding_seps for prefix in _encoding_prefixes
]


@pytest.fixture(params=_utf_values)
def utf_value(request):
    """
    Fixture for all possible integer values for a UTF encoding.
    """
    return request.param


@pytest.fixture(params=_encoding_fmts)
def encoding_fmt(request):
    """
    Fixture for all possible string formats of a UTF encoding.
    """
    return request.param
