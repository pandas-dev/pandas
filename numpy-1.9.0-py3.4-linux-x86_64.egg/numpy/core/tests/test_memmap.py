from __future__ import division, absolute_import, print_function

import sys
from tempfile import NamedTemporaryFile, TemporaryFile, mktemp, mkdtemp
import os
import shutil

from numpy import memmap
from numpy import arange, allclose, asarray
from numpy.testing import *

class TestMemmap(TestCase):
    def setUp(self):
        self.tmpfp = NamedTemporaryFile(prefix='mmap')
        self.tempdir = mkdtemp()
        self.shape = (3, 4)
        self.dtype = 'float32'
        self.data = arange(12, dtype=self.dtype)
        self.data.resize(self.shape)

    def tearDown(self):
        self.tmpfp.close()
        shutil.rmtree(self.tempdir)

    def test_roundtrip(self):
        # Write data to file
        fp = memmap(self.tmpfp, dtype=self.dtype, mode='w+',
                    shape=self.shape)
        fp[:] = self.data[:]
        del fp # Test __del__ machinery, which handles cleanup

        # Read data back from file
        newfp = memmap(self.tmpfp, dtype=self.dtype, mode='r',
                       shape=self.shape)
        assert_(allclose(self.data, newfp))
        assert_array_equal(self.data, newfp)

    def test_open_with_filename(self):
        tmpname = mktemp('', 'mmap', dir=self.tempdir)
        fp = memmap(tmpname, dtype=self.dtype, mode='w+',
                       shape=self.shape)
        fp[:] = self.data[:]
        del fp

    def test_unnamed_file(self):
        with TemporaryFile() as f:
            fp = memmap(f, dtype=self.dtype, shape=self.shape)
            del fp

    def test_attributes(self):
        offset = 1
        mode = "w+"
        fp = memmap(self.tmpfp, dtype=self.dtype, mode=mode,
                    shape=self.shape, offset=offset)
        self.assertEqual(offset, fp.offset)
        self.assertEqual(mode, fp.mode)
        del fp

    def test_filename(self):
        tmpname = mktemp('', 'mmap', dir=self.tempdir)
        fp = memmap(tmpname, dtype=self.dtype, mode='w+',
                       shape=self.shape)
        abspath = os.path.abspath(tmpname)
        fp[:] = self.data[:]
        self.assertEqual(abspath, fp.filename)
        b = fp[:1]
        self.assertEqual(abspath, b.filename)
        del b
        del fp

    def test_filename_fileobj(self):
        fp = memmap(self.tmpfp, dtype=self.dtype, mode="w+",
                    shape=self.shape)
        self.assertEqual(fp.filename, self.tmpfp.name)

    @dec.knownfailureif(sys.platform=='gnu0', "This test is known to fail on hurd")
    def test_flush(self):
        fp = memmap(self.tmpfp, dtype=self.dtype, mode='w+',
                    shape=self.shape)
        fp[:] = self.data[:]
        assert_equal(fp[0], self.data[0])
        fp.flush()

    def test_del(self):
        # Make sure a view does not delete the underlying mmap
        fp_base = memmap(self.tmpfp, dtype=self.dtype, mode='w+',
                    shape=self.shape)
        fp_base[0] = 5
        fp_view = fp_base[0:1]
        assert_equal(fp_view[0], 5)
        del fp_view
        # Should still be able to access and assign values after
        # deleting the view
        assert_equal(fp_base[0], 5)
        fp_base[0] = 6
        assert_equal(fp_base[0], 6)

    def test_arithmetic_drops_references(self):
        fp = memmap(self.tmpfp, dtype=self.dtype, mode='w+',
                    shape=self.shape)
        tmp = (fp + 10)
        if isinstance(tmp, memmap):
            assert tmp._mmap is not fp._mmap

    def test_indexing_drops_references(self):
        fp = memmap(self.tmpfp, dtype=self.dtype, mode='w+',
                    shape=self.shape)
        tmp = fp[[(1, 2), (2, 3)]]
        if isinstance(tmp, memmap):
            assert tmp._mmap is not fp._mmap

    def test_slicing_keeps_references(self):
        fp = memmap(self.tmpfp, dtype=self.dtype, mode='w+',
                    shape=self.shape)
        assert fp[:2, :2]._mmap is fp._mmap

    def test_view(self):
        fp = memmap(self.tmpfp, dtype=self.dtype, shape=self.shape)
        new1 = fp.view()
        new2 = new1.view()
        assert(new1.base is fp)
        assert(new2.base is fp)
        new_array = asarray(fp)
        assert(new_array.base is fp)

if __name__ == "__main__":
    run_module_suite()
