import unittest
import os
import tempfile
from nose2.tools import params
from vortexa_utils.versioning.versioner import Versioner

specs = [
    ((0, 0, 0), (0, 0, 1)),
    ((0, 0, 1), (0, 0, 2)),
    ((0, 1, 0), (0, 1, 0)),
    ((0, 1, 1), (0, 1, 0)),
    ((1, 0, 0), (1, 0, 0)),
    ((1, 0, 1), (1, 0, 0)),
    ((1, 1, 0), (1, 0, 0)),
    ((1, 1, 1), (1, 0, 0))
]


class TestVersioner(unittest.TestCase):
    def setUp(self):
        fh, filename = tempfile.mkstemp()
        os.fdopen(fh).close()
        self.version: Versioner = Versioner(filename)

    def tearDown(self):
        os.remove(self.version.VERSION_FILE)

    def test_version_none(self):
        self.assertEqual(self.version.__version__, None)

    def test_version_init(self):
        self.assertEqual(
            self.version.version,
            self.version.SemanticVersion(0, 0, 1)
        )
        self.assertTrue(os.path.isfile(self.version.VERSION_FILE))
        with open(self.version.VERSION_FILE, "r") as f:
            self.assertEqual(f.readline(), "0.0.1")

    @params(*specs)
    def test_version_incriment(self, flags, output):
        self.test_version_init()
        self.version.update_version(flags)
        self.assertEqual(
            self.version.version,
            self.version.SemanticVersion(*output)
        )
