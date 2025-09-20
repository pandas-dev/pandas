import importlib
import sys

import unittest


class TestLlvmVersion(unittest.TestCase):

    def test_llvmlite_version(self):
        # test the system it's running on
        import llvmlite
        import numba
        self.assertTrue(numba.__version__)

        llvmlite_version = llvmlite.__version__
        def cleanup():
            llvmlite.__version__ = llvmlite_version
        self.addCleanup(cleanup)

        # explicitly test all 3 cases of version string
        ver = numba._min_llvmlite_version
        version_pass = '%d.%d.%d' % ver
        git_version_pass = '%d.%d.%d-10-g92584ed' % ver
        rc_version_pass = '%d.%d.%drc1' % (ver[0], ver[1], ver[2] + 1)
        version_fail = '%d.%d.0' % (ver[0], ver[1] - 1)
        git_version_fail = '%d.%d.9-10-g92584ed' % (ver[0], ver[1] - 1)

        ver_pass = (version_pass, git_version_pass, rc_version_pass)
        ver_fail = (version_fail, git_version_fail)
        for v in ver_pass:
            llvmlite.__version__ = v
            importlib.reload(numba)
            self.assertTrue(numba.__version__)

        for v in ver_fail:
            with self.assertRaises(ImportError):
                llvmlite.__version__ = v
                importlib.reload(numba)


if __name__ == '__main__':
    unittest.main()
