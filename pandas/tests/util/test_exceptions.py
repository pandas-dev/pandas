import os
import types
import warnings

import pandas as pd

import pandas.util._exceptions as pue


def test_find_stack_level_stops_at_sibling_pandas_distributions() -> None:
    # GH#67992: a sibling distribution whose install directory merely starts
    # with "pandas" (e.g. pandas_ta) must not be treated as pandas-internal,
    # otherwise the warning is reported past the real caller.
    pue.find_stack_level()  # populate the cached _pkg_dir/_test_dir
    pkg_dir = os.path.dirname(pd.__file__)

    src = (
        "def caller():\n"
        "    warnings.warn('boom', UserWarning, stacklevel=find_stack_level())\n"
    )

    for path in (pkg_dir + "_ta/strategy.py", "/somewhere/else/strategy.py"):
        code = compile(src, path, "exec")
        caller_code = next(
            const for const in code.co_consts if isinstance(const, types.CodeType)
        )
        caller = types.FunctionType(
            caller_code,
            {
                "__builtins__": __builtins__,
                "warnings": warnings,
                "find_stack_level": pue.find_stack_level,
            },
        )
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            caller()
        assert rec[0].filename == path
        assert rec[0].lineno == 2
