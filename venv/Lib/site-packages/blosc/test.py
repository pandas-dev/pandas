import gc
import os
from ._version import LooseVersion
import ctypes
import blosc
import unittest

try:
    import numpy
except ImportError:
    has_numpy = False
else:
    has_numpy = True

try:
    import psutil
except ImportError:
    psutil = None


class TestCodec(unittest.TestCase):
    def setUp(self):
        self.PY_27_INPUT = b'\x02\x01\x03\x02\x85\x00\x00\x00\x84\x00\x00\x00\x95\x00\x00\x00\x80\x02cnumpy.core.multiarray\n_reconstruct\nq\x01cnumpy\nndarray\nq\x02K\x00\x85U\x01b\x87Rq\x03(K\x01K\x05\x85cnumpy\ndtype\nq\x04U\x02S2K\x00K\x01\x87Rq\x05(K\x03U\x01|NNNK\x02K\x01K\x00tb\x89U\n\xc3\xa5\xc3\xa7\xc3\xb8\xcf\x80\xcb\x9atb.'

    def test_basic_codec(self):
        s = b'0123456789'
        c = blosc.compress(s, typesize=1)
        d = blosc.decompress(c)
        self.assertEqual(s, d)

    def test_get_clib(self):
        s = b'0123456789'
        for cname in blosc.compressor_list():
            c = blosc.compress(s, typesize=1, cname=cname)
            clib = blosc.get_clib(c)
            self.assertEqual(clib, blosc.cname2clib[cname])

    def test_all_compressors(self):
        s = b'0123456789'*100
        for cname in blosc.compressor_list():
            c = blosc.compress(s, typesize=1, cname=cname)
            d = blosc.decompress(c)
            self.assertEqual(s, d)

    def test_all_filters(self):
        s = b'0123456789'*100
        filters = [blosc.NOSHUFFLE, blosc.SHUFFLE]
        # BITFILTER only works properly from 1.8.0 on
        if LooseVersion(blosc.blosclib_version) >= LooseVersion("1.8.0"):
            filters.append(blosc.BITSHUFFLE)
        for filter_ in filters:
            c = blosc.compress(s, typesize=1, shuffle=filter_)
            d = blosc.decompress(c)
            self.assertEqual(s, d)

    def test_set_nthreads_exceptions(self):
        self.assertRaises(ValueError, blosc.set_nthreads,
                          blosc.MAX_THREADS + 1)

    def test_compress_input_types(self):
        import numpy as np
        # assume the expected answer was compressed from bytes
        expected = blosc.compress(b'0123456789', typesize=1)

        # now for all the things that support the buffer interface
        self.assertEqual(expected,
                         blosc.compress(memoryview(b'0123456789'),
                                        typesize=1))

        self.assertEqual(expected, blosc.compress(
            bytearray(b'0123456789'), typesize=1))
        self.assertEqual(expected, blosc.compress(
            np.array([b'0123456789']), typesize=1))

    def test_decompress_input_types(self):
        import numpy as np
        # assume the expected answer was compressed from bytes
        expected = b'0123456789'
        compressed = blosc.compress(expected, typesize=1)

        # now for all the things that support the buffer interface
        self.assertEqual(expected, blosc.decompress(compressed))
        self.assertEqual(expected, blosc.decompress(memoryview(compressed)))

        self.assertEqual(expected, blosc.decompress(bytearray(compressed)))
        self.assertEqual(expected, blosc.decompress(np.array([compressed])))

    def test_decompress_ptr_input_types(self):
        import numpy as np
        # assume the expected answer was compressed from bytes
        expected = b'0123456789'
        out = np.zeros(len(expected), dtype=np.byte)
        compressed = blosc.compress(expected, typesize=1)

        # now for all the things that support the buffer interface
        out[:] = 0  # reset the output array
        nout = blosc.decompress_ptr(compressed, out.ctypes.data)
        self.assertEqual(expected, out.tobytes())
        self.assertEqual(len(expected), nout)  # check that we didn't write too many bytes

        out[:] = 0
        nout = blosc.decompress_ptr(memoryview(compressed), out.ctypes.data)
        self.assertEqual(expected, out.tobytes())
        self.assertEqual(len(expected), nout)

        out[:] = 0
        nout = blosc.decompress_ptr(bytearray(compressed), out.ctypes.data)
        self.assertEqual(expected, out.tobytes())
        self.assertEqual(len(expected), nout)

        out[:] = 0
        nout = blosc.decompress_ptr(np.frombuffer(compressed, dtype=np.byte), out.ctypes.data)
        self.assertEqual(expected, out.tobytes())
        self.assertEqual(len(expected), nout)

    def test_decompress_releasegil(self):
        import numpy as np
        # assume the expected answer was compressed from bytes
        blosc.set_releasegil(True)
        expected = b'0123456789'
        compressed = blosc.compress(expected, typesize=1)

        # now for all the things that support the buffer interface
        self.assertEqual(expected, blosc.decompress(compressed))
        self.assertEqual(expected, blosc.decompress(memoryview(compressed)))

        self.assertEqual(expected, blosc.decompress(bytearray(compressed)))
        self.assertEqual(expected, blosc.decompress(np.array([compressed])))
        blosc.set_releasegil(False)

    def test_decompress_input_types_as_bytearray(self):
        import numpy as np
        # assume the expected answer was compressed from bytes
        expected = bytearray(b'0123456789')
        compressed = blosc.compress(expected, typesize=1)

        # now for all the things that support the buffer interface
        self.assertEqual(expected, blosc.decompress(compressed,
                                                    as_bytearray=True))
        self.assertEqual(expected,
                         blosc.decompress(memoryview(compressed),
                                          as_bytearray=True))

        self.assertEqual(expected, blosc.decompress(bytearray(compressed),
                                                    as_bytearray=True))
        self.assertEqual(expected, blosc.decompress(np.array([compressed]),
                                                    as_bytearray=True))

    def test_compress_exceptions(self):
        s = b'0123456789'

        self.assertRaises(ValueError, blosc.compress, s, typesize=0)
        self.assertRaises(ValueError, blosc.compress, s,
                          typesize=blosc.MAX_TYPESIZE+1)

        self.assertRaises(ValueError, blosc.compress, s, typesize=1, clevel=-1)
        self.assertRaises(ValueError, blosc.compress, s, typesize=1, clevel=10)

        self.assertRaises(TypeError, blosc.compress, 1.0, 1)
        self.assertRaises(TypeError, blosc.compress, ['abc'], 1)

        self.assertRaises(ValueError, blosc.compress, 'abc',
                          typesize=1, cname='foo')

        # Create a simple mock to avoid having to create a buffer of 2 GB
        class LenMock:
            def __len__(self):
                return blosc.MAX_BUFFERSIZE+1
        self.assertRaises(ValueError, blosc.compress, LenMock(), typesize=1)

    def test_compress_ptr_exceptions(self):
        # Make sure we do have a valid address, to reduce the chance of a
        # segfault if we do actually start compressing because the exceptions
        # aren't raised.
        typesize, items = 8, 8
        data = [float(i) for i in range(items)]
        Array = ctypes.c_double * items
        array = Array(*data)
        address = ctypes.addressof(array)

        self.assertRaises(ValueError, blosc.compress_ptr, address, items,
                          typesize=-1)
        self.assertRaises(ValueError, blosc.compress_ptr, address, items,
                          typesize=blosc.MAX_TYPESIZE+1)

        self.assertRaises(ValueError, blosc.compress_ptr, address, items,
                          typesize=typesize, clevel=-1)
        self.assertRaises(ValueError, blosc.compress_ptr, address, items,
                          typesize=typesize, clevel=10)

        self.assertRaises(TypeError, blosc.compress_ptr, 1.0, items,
                          typesize=typesize)
        self.assertRaises(TypeError, blosc.compress_ptr, ['abc'], items,
                          typesize=typesize)

        self.assertRaises(ValueError, blosc.compress_ptr, address, -1,
                          typesize=typesize)
        self.assertRaises(ValueError, blosc.compress_ptr, address,
                          blosc.MAX_BUFFERSIZE+1, typesize=typesize)

    def test_decompress_exceptions(self):
        self.assertRaises(TypeError, blosc.decompress, 1.0)
        self.assertRaises(TypeError, blosc.decompress, ['abc'])

    def test_decompress_ptr_exceptions(self):
        # make sure we do have a valid address
        typesize, items = 8, 8
        data = [float(i) for i in range(items)]
        Array = ctypes.c_double * items
        in_array = Array(*data)
        c = blosc.compress_ptr(ctypes.addressof(in_array), items, typesize)
        out_array = ctypes.create_string_buffer(items*typesize)

        self.assertRaises(TypeError, blosc.decompress_ptr, 1.0,
                          ctypes.addressof(out_array))
        self.assertRaises(TypeError, blosc.decompress_ptr, ['abc'],
                          ctypes.addressof(out_array))

        self.assertRaises(TypeError, blosc.decompress_ptr, c, 1.0)
        self.assertRaises(TypeError, blosc.decompress_ptr, c, ['abc'])

    @unittest.skipIf(not has_numpy, "Numpy not available")
    def test_pack_array_exceptions(self):

        self.assertRaises(TypeError, blosc.pack_array, 'abc')
        self.assertRaises(TypeError, blosc.pack_array, 1.0)

        items = (blosc.MAX_BUFFERSIZE // 8) + 1
        one = numpy.ones(1, dtype=numpy.int64)
        self.assertRaises(ValueError, blosc.pack_array, one, clevel=-1)
        self.assertRaises(ValueError, blosc.pack_array, one, clevel=10)

        # use stride trick to make an array that looks like a huge one
        ones = numpy.lib.stride_tricks.as_strided(one, shape=(1, items),
                                                  strides=(8, 0))[0]

        # This should always raise an error
        self.assertRaises(ValueError, blosc.pack_array, ones)

    def test_unpack_array_with_unicode_characters(self):
        import numpy as np
        input_array = np.array(['å', 'ç', 'ø', 'π', '˚'])
        packed_array = blosc.pack_array(input_array)
        np.testing.assert_array_equal(input_array, blosc.unpack_array(packed_array, encoding='UTF-8'))

    def test_unpack_array_with_from_py27_exceptions(self):
        self.assertRaises(UnicodeDecodeError, blosc.unpack_array, self.PY_27_INPUT)

    def test_unpack_array_with_unicode_characters_from_py27(self):
        import numpy as np
        out_array = np.array(['å', 'ç', 'ø', 'π', '˚'])
        np.testing.assert_array_equal(out_array, blosc.unpack_array(self.PY_27_INPUT, encoding='bytes'))

    def test_unpack_array_exceptions(self):
        self.assertRaises(TypeError, blosc.unpack_array, 1.0)

    @unittest.skipIf(not psutil, "psutil not available, cannot test for leaks")
    def test_no_leaks(self):

        num_elements = 10000000
        typesize = 8
        data = [float(i) for i in range(num_elements)]  # ~76MB
        Array = ctypes.c_double * num_elements
        array = Array(*data)
        address = ctypes.addressof(array)

        def leaks(operation, repeats=3):
            gc.collect()
            used_mem_before = psutil.Process(os.getpid()).memory_info()[0]
            for _ in range(repeats):
                operation()
            gc.collect()
            used_mem_after = psutil.Process(os.getpid()).memory_info()[0]
            # We multiply by an additional factor of .01 to account for
            # storage overhead of Python classes
            return (used_mem_after - used_mem_before) >= num_elements * 8.01

        def compress():
            blosc.compress(array, typesize, clevel=1)

        def compress_ptr():
            blosc.compress_ptr(address, num_elements, typesize, clevel=0)

        def decompress():
            cx = blosc.compress(array, typesize, clevel=1)
            blosc.decompress(cx)

        def decompress_ptr():
            cx = blosc.compress_ptr(address, num_elements, typesize, clevel=0)
            blosc.decompress_ptr(cx, address)

        self.assertFalse(leaks(compress), msg='compress leaks memory')
        self.assertFalse(leaks(compress_ptr), msg='compress_ptr leaks memory')
        self.assertFalse(leaks(decompress), msg='decompress leaks memory')
        self.assertFalse(leaks(decompress_ptr), msg='decompress_ptr leaks memory')

    def test_get_blocksize(self):
        s = b'0123456789' * 1000
        blosc.set_blocksize(2**14)
        blosc.compress(s, typesize=1)
        d = blosc.get_blocksize()
        self.assertEqual(d, 2**14)

    def test_get_cbuffer_sizes(self):
        s = b'0123456789' * 100000
        blosc.set_blocksize(2**16)
        c = blosc.compress(s, typesize=1)
        t = blosc.get_cbuffer_sizes(c)
        self.assertEqual(t[0], 1000000)
        # One cannot be sure of the exact compressed bytes, so round to KB
        self.assertEqual(t[1] // 2**10, 4354 // 2**10)
        self.assertEqual(t[2], 2**16)

    def test_bitshuffle_not_multiple(self):
        # Check the fix for #133
        x = numpy.ones(27266, dtype='uint8')
        xx = x.tobytes()
        zxx = blosc.compress(xx, typesize=8, shuffle=blosc.BITSHUFFLE)
        last_xx = blosc.decompress(zxx)[-3:]
        self.assertEqual(last_xx, b'\x01\x01\x01')

    def test_cbuffer_validate(self):
        import numpy as np
        expected = b'0123456789' * 1000
        compressed = blosc.compress(expected)

        # now for all the things that support the buffer interface
        self.assertTrue(blosc.cbuffer_validate(compressed))
        self.assertTrue(blosc.cbuffer_validate(memoryview(compressed)))

        self.assertTrue(blosc.cbuffer_validate(bytearray(compressed)))
        self.assertTrue(blosc.cbuffer_validate(np.array([compressed])))

    def test_cbuffer_validate_failures(self):
        import numpy as np
        compressed = b'this_is_total_garbage'

        # now for all the things that support the buffer interface
        self.assertFalse(blosc.cbuffer_validate(compressed))
        self.assertFalse(blosc.cbuffer_validate(memoryview(compressed)))

        self.assertFalse(blosc.cbuffer_validate(bytearray(compressed)))
        self.assertFalse(blosc.cbuffer_validate(np.array([compressed])))

    def test_bithuffle_leftovers(self):
        # Test for https://github.com/Blosc/c-blosc2/pull/100
        import numpy as np
        buffer = b" " * 641091  # a buffer that is not divisible by 8
        cbuffer = blosc.compress(buffer, typesize=8, shuffle=blosc.BITSHUFFLE, clevel=1)
        dbuffer = blosc.decompress(cbuffer)
        self.assertTrue(buffer == dbuffer)


def run(verbosity=2):
    import blosc
    import blosc.toplevel
    blosc.print_versions()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCodec)
    # If in the future we split this test file in several, the auto-discover
    # might be interesting

    # suite = unittest.TestLoader().discover(start_dir='.', pattern='test*.py')
    suite.addTests(unittest.TestLoader().loadTestsFromModule(blosc.toplevel))
    assert unittest.TextTestRunner(verbosity=verbosity).\
        run(suite).wasSuccessful()


if __name__ == '__main__':
    run()
