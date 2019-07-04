import pytest

import pandas as pd
import pandas.util.testing as tm


class TestCategoricalWarnings:
    def test_tab_complete_warning(self, ip):
        # https://github.com/pandas-dev/pandas/issues/16409
        pytest.importorskip("IPython", minversion="6.0.0")
        from IPython.core.completer import provisionalcompleter

        code = "import pandas as pd; c = Categorical([])"
        ip.run_code(code)
        with tm.assert_produces_warning(None):
            with provisionalcompleter("ignore"):
                list(ip.Completer.completions("c.", 1))

    def test_CategoricalAccessor_categorical_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            pd.Series(["a", "b"], dtype="category").cat.categorical

    def test_CategoricalAccessor_name_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            pd.Series(["a", "b"], dtype="category").cat.name

    def test_CategoricalAccessor_index_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            pd.Series(["a", "b"], dtype="category").cat.index
