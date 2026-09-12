import importlib.util
import os
import shutil
import warnings

import pytest

import pandas.util._exceptions as pue

_SRC = (
    "import warnings\n"
    "from pandas.util._exceptions import find_stack_level\n"
    "def caller():\n"
    "    warnings.warn('boom', UserWarning, stacklevel=find_stack_level())\n"
)


@pytest.mark.parametrize(
    "suffix, tail",
    [
        # third-party distributions whose install directory merely starts
        # with the pandas package directory, e.g. site-packages/pandas_ta
        ("_ta", "strategy.py"),
        ("_gbq", "io.py"),
        # a pandas-prefixed sibling with its own tests directory
        ("2/tests", "mod.py"),
    ],
)
def test_find_stack_level_pandas_prefixed_siblings(tmp_path, suffix, tail):
    # GH#67992: find_stack_level must not treat pandas-prefixed sibling
    # directories as pandas-internal, else the reported stacklevel points
    # past the caller and user code can end up with no warning at all.
    pue.find_stack_level()  # populate the _pkg_dir/_test_dir cache

    pkg_dir = os.path.normpath(os.fspath(pue._pkg_dir))
    parent = os.path.dirname(pkg_dir)
    pkg_name = os.path.basename(pkg_dir)
    module_dir = os.path.join(parent, pkg_name + suffix)
    module_path = os.path.join(module_dir, *tail.split("/"))

    try:
        os.makedirs(os.path.dirname(module_path), exist_ok=True)
        with open(module_path, "w", encoding="utf-8") as fh:
            fh.write(_SRC)
        spec = importlib.util.spec_from_file_location(f"gh67992{suffix}", module_path)
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            module.caller()
    except OSError:
        pytest.skip(f"could not create {module_dir}")
    finally:
        shutil.rmtree(module_dir, ignore_errors=True)

    assert len(record) == 1
    assert record[0].filename == module_path
