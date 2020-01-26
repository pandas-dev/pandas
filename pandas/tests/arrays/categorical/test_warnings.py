import pytest

from pandas.util._test_decorators import async_mark

import pandas._testing as tm


class TestCategoricalWarnings:
    @async_mark()
    async def test_tab_complete_warning(self, ip):
        # https://github.com/pandas-dev/pandas/issues/16409
        pytest.importorskip("IPython", minversion="6.0.0")
        from IPython.core.completer import provisionalcompleter

        code = "import pandas as pd; c = Categorical([])"
        await ip.run_code(code)

        # deprecation warning is raised due to jedi upgrades
        with tm.assert_produces_warning(DeprecationWarning):
            with provisionalcompleter("ignore"):
                list(ip.Completer.completions("c.", 1))
