import importlib
import sys

import pytest

import pandas._testing as tm

depr_msg_template = (
    "pandas.core.{0} is deprecated and will be moved to pandas._core.{0} "
    "in a future version. Import from the public API instead."
)


def clear_sys_modules(name):
    keys = list(sys.modules.keys())
    for key in keys:
        if name in key:
            del sys.modules[key]


@pytest.mark.parametrize(
    "import_statement",
    [
        "import pandas.core.groupby",
        "from pandas.core.groupby import generic",
        "import pandas.core.groupby.generic",
        "from pandas.core.groupby import GroupBy",
        "from pandas.core.groupby.groupby import GroupBy",
    ],
)
def test_core_depr(import_statement):
    msg = depr_msg_template.format("groupby")
    clear_sys_modules("pandas.core.groupby")
    with tm.assert_produces_warning(
        DeprecationWarning, match=msg, check_stacklevel="<string>"
    ):
        # Use of exec makes warning raised from str
        exec(import_statement)  # noqa: PDF003

    # importing twice should not raise a second warning
    with tm.assert_produces_warning(None):
        exec(import_statement)  # noqa: PDF003


@pytest.mark.parametrize("submodule", ["groupby"])
def test_core_depr_importlib(submodule):
    msg = depr_msg_template.format(submodule)
    clear_sys_modules("pandas.core")
    # Anything following `from pandas.core import ...` will be accessed
    # via importlib in pandas.core.__init__.
    with tm.assert_produces_warning(
        DeprecationWarning, match=msg, check_stacklevel=importlib.__file__
    ):
        # Use of exec makes warning raised from importlib
        exec(f"from pandas.core import {submodule}")  # noqa: PDF003

    # importing twice should not raise a second warning
    with tm.assert_produces_warning(None):
        exec(f"from pandas.core import {submodule}")  # noqa: PDF003
