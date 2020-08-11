import pytest

import pandas._testing as tm

from pandas.io.excel import ExcelWriter

odf = pytest.importorskip("odf")

pytestmark = pytest.mark.parametrize("ext", [".ods"])


def setup_module(module):
    import pandas.util._test_decorators as td

    yield from td.check_file_leaks()


def test_write_append_mode_raises(ext):
    msg = "Append mode is not supported with odf!"

    with tm.ensure_clean(ext) as f:
        with pytest.raises(ValueError, match=msg):
            ExcelWriter(f, engine="odf", mode="a")
