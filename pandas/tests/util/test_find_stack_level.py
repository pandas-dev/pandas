import os
import types
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
def test_find_stack_level_pandas_prefixed_siblings(suffix, tail):
    # GH#67992: find_stack_level must not treat pandas-prefixed sibling
    # directories as pandas-internal, else the reported stacklevel points
    # past the caller and user code can end up with no warning at all.
    pue.find_stack_level()  # populate the _pkg_dir/_test_dir cache

    pkg_dir = os.path.normpath(os.fspath(pue._pkg_dir))
    module_path = os.path.normpath(
        os.path.join(
            os.path.dirname(pkg_dir),
            os.path.basename(pkg_dir) + suffix,
            *tail.split("/"),
        )
    )

    # An installed pandas_ta/pandas_gbq has exactly module_path as its file,
    # so instead of creating that file (which would clobber a real pandas_*
    # install next to pandas) compile the source with module_path as its
    # filename: the frame then reports the same co_filename it would have
    # when imported from such an install.
    module = types.ModuleType(f"gh67992{suffix}")
    module.__file__ = module_path
    exec(compile(_SRC, module_path, "exec"), module.__dict__)  # noqa: S102

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        module.caller()

    assert len(record) == 1
    assert record[0].filename == module_path
