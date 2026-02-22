import sys
import zlib
import itertools

import numpy as np

import tables as tb
from tables.tests import common


class ArrayDirectChunkingTestCase(
    common.TempFileMixin, common.PyTablesTestCase
):
    obj = np.arange(25, dtype="uint8")

    def setUp(self):
        super().setUp()
        self.array = self.h5file.create_array("/", "array", self.obj)

    def test_chunk_info(self):
        self.assertRaises(
            tb.NotChunkedError, self.array.chunk_info, (0,) * self.array.ndim
        )

    def test_read_chunk(self):
        self.assertRaises(
            tb.NotChunkedError, self.array.read_chunk, (0,) * self.array.ndim
        )

    def test_read_chunk_out(self):
        arr = np.zeros(self.obj.shape, dtype=self.obj.dtype)
        self.assertRaises(
            tb.NotChunkedError,
            self.array.read_chunk,
            (0,) * self.array.ndim,
            out=memoryview(arr),
        )

    def test_write_chunk(self):
        arr = self.obj // 2
        self.assertRaises(
            tb.NotChunkedError,
            self.array.write_chunk,
            (0,) * self.array.ndim,
            arr,
        )


# For enlargeable and non-enlargeable datasets.
class DirectChunkingTestCase(common.TempFileMixin, common.PyTablesTestCase):
    # Class attributes:
    shape: tuple[int, ...]
    chunkshape: tuple[int, ...]
    shuffle: bool
    obj: np.ndarray

    # Instance attributes:
    array: tb.Leaf  # set by ``setUp()`` and ``_reopen()``
    filters: tb.Filters

    def setUp(self):
        super().setUp()
        self.filters = tb.Filters(
            complib="zlib", complevel=1, shuffle=self.shuffle
        )

    def modified(self, obj):
        # Return altered copy with same dtype and shape.
        raise NotImplementedError

    def iter_chunks(self):
        chunk_ranges = list(
            range(0, s, cs) for (s, cs) in zip(self.shape, self.chunkshape)
        )
        yield from itertools.product(*chunk_ranges)

    def test_chunk_info_aligned(self):
        for chunk_start in self.iter_chunks():
            chunk_info = self.array.chunk_info(chunk_start)
            self.assertEqual(chunk_info.start, chunk_start)
            self.assertIsNotNone(chunk_info.filter_mask)
            self.assertIsNotNone(chunk_info.offset)
            self.assertIsNotNone(chunk_info.size)

    def test_chunk_info_unaligned(self):
        chunk_info_a = self.array.chunk_info((0,) * self.array.ndim)
        chunk_info_u = self.array.chunk_info((1,) * self.array.ndim)
        self.assertIsNotNone(chunk_info_a.start)
        self.assertEqual(chunk_info_a, chunk_info_u)

    def test_chunk_info_aligned_beyond(self):
        beyond = tuple(
            (1 + s // cs) * cs for (s, cs) in zip(self.shape, self.chunkshape)
        )
        self.assertRaises(IndexError, self.array.chunk_info, beyond)

    def test_chunk_info_unaligned_beyond(self):
        beyond = tuple(
            1 + (1 + s // cs) * cs
            for (s, cs) in zip(self.shape, self.chunkshape)
        )
        self.assertRaises(IndexError, self.array.chunk_info, beyond)

    def shuffled(self, bytes_):
        itemsize = self.obj.dtype.itemsize
        return b"".join(bytes_[d::itemsize] for d in range(itemsize))

    def filter_chunk(self, bytes_, shuffle=None):
        assert self.filters.complib == "zlib"
        if shuffle is None:
            shuffle = self.shuffle
        maybe_shuffled = self.shuffled(bytes_) if shuffle else bytes_
        return zlib.compress(maybe_shuffled, self.filters.complevel)

    def test_read_chunk(self):
        # Extended to fit chunk boundaries.
        ext_obj = np.pad(
            self.obj,
            [(0, s % cs) for (s, cs) in zip(self.shape, self.chunkshape)],
        )
        for chunk_start in self.iter_chunks():
            chunk = self.array.read_chunk(chunk_start)
            self.assertIsInstance(chunk, bytes)
            obj_slice = tuple(
                slice(s, s + cs)
                for (s, cs) in zip(chunk_start, self.chunkshape)
            )
            obj_bytes = self.filter_chunk(ext_obj[obj_slice].tobytes())
            self.assertEqual(chunk, obj_bytes)

    def test_read_chunk_out(self):
        # Extended to fit chunk boundaries.
        ext_obj = np.pad(
            self.obj,
            [(0, s % cs) for (s, cs) in zip(self.shape, self.chunkshape)],
        )
        chunk_start = (0,) * self.obj.ndim
        obj_slice = tuple(
            slice(s, s + cs) for (s, cs) in zip(chunk_start, self.chunkshape)
        )
        obj_bytes = self.filter_chunk(ext_obj[obj_slice].tobytes())
        chunk_size = len(obj_bytes)

        chunk_out = bytearray(chunk_size - 1)  # too short
        self.assertRaises(
            ValueError, self.array.read_chunk, chunk_start, out=chunk_out
        )

        chunk_out = bytearray(chunk_size)
        chunk = self.array.read_chunk(chunk_start, out=chunk_out)
        self.assertIsInstance(chunk, memoryview)
        self.assertEqual(chunk, obj_bytes)
        self.assertEqual(chunk_out, obj_bytes)

    def test_read_chunk_unaligned(self):
        self.assertRaises(
            tb.NotChunkAlignedError,
            self.array.read_chunk,
            (1,) * self.array.ndim,
        )

    def test_read_chunk_beyond(self):
        beyond = tuple(
            (1 + s // cs) * cs for (s, cs) in zip(self.shape, self.chunkshape)
        )
        self.assertRaises(IndexError, self.array.read_chunk, beyond)

    def test_write_chunk(self):
        new_obj = self.modified(self.obj)
        # Extended to fit chunk boundaries.
        ext_obj = np.pad(
            new_obj,
            [(0, s % cs) for (s, cs) in zip(self.shape, self.chunkshape)],
        )
        for chunk_start in self.iter_chunks():
            obj_slice = tuple(
                slice(s, s + cs)
                for (s, cs) in zip(chunk_start, self.chunkshape)
            )
            obj_bytes = self.filter_chunk(ext_obj[obj_slice].tobytes())
            self.array.write_chunk(chunk_start, obj_bytes)

        self._reopen()
        self.assertTrue(common.areArraysEqual(self.array[:], new_obj))

    def test_write_chunk_filtermask(self):
        no_shuffle_mask = 0x00000004  # to turn shuffle off

        chunk_start = (0,) * self.obj.ndim
        obj_slice = tuple(
            slice(s, s + cs) for (s, cs) in zip(chunk_start, self.chunkshape)
        )
        new_obj = self.obj.copy()
        new_obj[obj_slice] = self.modified(new_obj[obj_slice])
        obj_bytes = self.filter_chunk(
            new_obj[obj_slice].tobytes(), shuffle=False
        )
        self.array.write_chunk(
            chunk_start, obj_bytes, filter_mask=no_shuffle_mask
        )

        self._reopen()
        arr_obj = self.array[:]  # first chunk is shuffled, fix it
        fixed_bytes = self.shuffled(arr_obj[obj_slice].tobytes())
        fixed_chunk = np.ndarray(
            self.chunkshape, dtype=self.obj.dtype, buffer=fixed_bytes
        )
        arr_obj[obj_slice] = fixed_chunk
        self.assertTrue(common.areArraysEqual(arr_obj, new_obj))

        chunk_info = self.array.chunk_info(chunk_start)
        self.assertEqual(chunk_info.filter_mask, no_shuffle_mask)

    def test_write_chunk_unaligned(self):
        self.assertRaises(
            tb.NotChunkAlignedError,
            self.array.write_chunk,
            (1,) * self.array.ndim,
            b"foobar",
        )

    def test_write_chunk_beyond(self):
        beyond = tuple(
            (1 + s // cs) * cs for (s, cs) in zip(self.shape, self.chunkshape)
        )
        self.assertRaises(
            IndexError, self.array.write_chunk, beyond, b"foobar"
        )


# For enlargeable datasets only.
class XDirectChunkingTestCase(DirectChunkingTestCase):
    def test_chunk_info_miss_extdim(self):
        # Next chunk in the enlargeable dimension.
        assert self.array.extdim == 0
        chunk_start = (
            ((1 + self.shape[0] // self.chunkshape[0]) * self.chunkshape[0]),
            *((0,) * (self.array.ndim - 1)),
        )
        self.assertRaises(IndexError, self.array.chunk_info, chunk_start)

        # Enlarge the array to put the (missing) chunk within the shape.
        self.array.truncate(chunk_start[0] + self.chunkshape[0])
        chunk_info = self.array.chunk_info(chunk_start)
        self.assertIsNone(chunk_info.filter_mask)
        self.assertIsNone(chunk_info.offset)
        self.assertIsNone(chunk_info.size)

    def test_chunk_info_miss_noextdim(self):
        if self.array.ndim < 2:
            raise common.unittest.SkipTest(
                "missing chunk always within enlargeable dimension"
            )

        # Next chunk in the first non-enlargeable dimension.
        assert self.array.extdim != 1
        chunk_start = (
            0,
            ((1 + self.shape[1] // self.chunkshape[1]) * self.chunkshape[1]),
            *((0,) * (self.array.ndim - 2)),
        )
        self.assertRaises(IndexError, self.array.chunk_info, chunk_start)

    def test_read_chunk_miss_extdim(self):
        # Next chunk in the enlargeable dimension.
        assert self.array.extdim == 0
        chunk_start = (
            ((1 + self.shape[0] // self.chunkshape[0]) * self.chunkshape[0]),
            *((0,) * (self.array.ndim - 1)),
        )
        self.assertRaises(IndexError, self.array.read_chunk, chunk_start)

        # Enlarge the array to put the (missing) chunk within the shape.
        self.array.truncate(chunk_start[0] + self.chunkshape[0])
        self.assertRaises(
            tb.NoSuchChunkError, self.array.read_chunk, chunk_start
        )

    def _test_write_chunk_missing(self, shrink_after):
        # Enlarge array by two chunk rows,
        # copy first old chunk in first chunk of new last chunk row.
        assert self.array.extdim == 0
        chunk_start = (
            (
                (1 + self.shape[0] // self.chunkshape[0]) * self.chunkshape[0]
                + self.chunkshape[0]
            ),
            *((0,) * (self.array.ndim - 1)),
        )
        chunk = self.array.read_chunk((0,) * self.array.ndim)
        self.array.truncate(chunk_start[0] + self.chunkshape[0])
        self.array.write_chunk(chunk_start, chunk)
        if shrink_after:
            self.array.truncate(self.shape[0] + 1)
            self.array.truncate(self.shape[0] - 1)

        new_obj = self.obj.copy()
        new_obj.resize(self.array.shape, refcheck=False)
        obj_slice = tuple(
            slice(s, s + cs) for (s, cs) in zip(chunk_start, self.chunkshape)
        )
        if not shrink_after:
            new_obj[obj_slice] = new_obj[
                tuple(slice(0, cs) for cs in self.chunkshape)
            ]

        self._reopen()
        self.assertTrue(common.areArraysEqual(self.array[:], new_obj))

    def test_write_chunk_missing1(self):
        return self._test_write_chunk_missing(shrink_after=False)

    def test_write_chunk_missing2(self):
        return self._test_write_chunk_missing(shrink_after=True)


class CArrayDirectChunkingTestCase(DirectChunkingTestCase):
    shape = (5, 5)
    chunkshape = (2, 2)  # 3 x 3 chunks, incomplete at right/bottom boundaries
    shuffle = True
    obj = np.arange(np.prod(shape), dtype="u2").reshape(shape)

    def setUp(self):
        super().setUp()
        self.array = self.h5file.create_carray(
            "/",
            "carray",
            chunkshape=self.chunkshape,
            obj=self.obj,
            filters=self.filters,
        )

    def _reopen(self):
        super()._reopen()
        self.array = self.h5file.root.carray

    def modified(self, obj):
        return obj * 2


class EArrayDirectChunkingTestCase(XDirectChunkingTestCase):
    shape = (5, 5)  # enlargeable along first dimension
    chunkshape = (2, 2)  # 3 x 3 chunks, incomplete at right/bottom boundaries
    shuffle = True
    obj = np.arange(np.prod(shape), dtype="u2").reshape(shape)

    def setUp(self):
        super().setUp()
        atom = tb.Atom.from_dtype(self.obj.dtype)
        shape = (0, *self.shape[1:])
        self.array = self.h5file.create_earray(
            "/",
            "earray",
            atom,
            shape,
            chunkshape=self.chunkshape,
            filters=self.filters,
        )
        self.array.append(self.obj)

    def _reopen(self):
        super()._reopen()
        self.array = self.h5file.root.earray

    def modified(self, obj):
        return obj * 2


class TableDirectChunkingTestCase(XDirectChunkingTestCase):
    shape = (5,)  # enlargeable along first dimension
    chunkshape = (2,)  # 3 chunks, incomplete at bottom boundary
    shuffle = True
    obj = np.array(
        [(i, float(i)) for i in range(np.prod(shape))], dtype="u4,f4"
    )

    def setUp(self):
        super().setUp()
        desc, _ = tb.descr_from_dtype(self.obj.dtype)
        self.array = self.h5file.create_table(
            "/",
            "table",
            desc,
            chunkshape=self.chunkshape,
            filters=self.filters,
        )
        self.array.append(self.obj)

    def _reopen(self):
        super()._reopen()
        self.array = self.h5file.root.table

    def modified(self, obj):
        flat = obj.copy().reshape((np.prod(obj.shape),))
        fnames = flat.dtype.names
        for i in range(len(flat)):
            for f in fnames:
                flat[i][f] *= 2
        return flat.reshape(obj.shape)


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1

    for i in range(niter):
        theSuite.addTest(common.make_suite(ArrayDirectChunkingTestCase))
        theSuite.addTest(common.make_suite(CArrayDirectChunkingTestCase))
        theSuite.addTest(common.make_suite(EArrayDirectChunkingTestCase))
        theSuite.addTest(common.make_suite(TableDirectChunkingTestCase))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
