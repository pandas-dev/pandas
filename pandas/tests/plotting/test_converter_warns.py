# This test ensure that we warn for implicitly registered converters.
# We isolate it completely, because its behavior depends on some global state
# set at import time, which is tricky to get right.
# We use unittest instead of pytest, since pytest will import pandas when
# loading our conftest.py
# https://github.com/pandas-dev/pandas/issues/24963
# https://github.com/pandas-dev/pandas/pull/18307
# TODO: Remove the special runner in ci/run_tests.sh
import sys
import unittest


class TestConverter(unittest.TestCase):
    @unittest.skipIf(sys.version_info[0] == 2, "CI Failure")
    @unittest.skipIf("pandas" in sys.modules, "pandas musn't be imported.")
    def test_converter_warning(self):
        import pandas as pd
        import pandas.util.testing as tm
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise unittest.SkipTest("No matplotlib")

        fig, ax = plt.subplots()
        ser = pd.Series(range(12), index=pd.date_range('2000', periods=12))
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            ax.plot(ser)


if __name__ == '__main__':
    unittest.main()
