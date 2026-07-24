"""Tests for gdb interacting with the DWARF numba generates"""
from numba.tests.support import TestCase, linux_only
from numba.tests.gdb_support import needs_gdb, skip_unless_pexpect, GdbMIDriver
from unittest.mock import patch, Mock
from numba.core import datamodel
import numpy as np
from numba import typeof
import ctypes as ct
import unittest


@linux_only
@needs_gdb
@skip_unless_pexpect
class TestGDBDwarf(TestCase):
    # This runs the tests in numba.tests.gdb, each submodule must contain one
    # test class called "Test" and it must contain one test called "test".
    # Variation is provided by the module name. The reason this convention exits
    # is because gdb tests tend to be line number sensitive (breakpoints etc
    # care about this) and doing this prevents constant churn and permits the
    # reuse of the existing subprocess_test_runner harness.
    _NUMBA_OPT_0_ENV = {'NUMBA_OPT': '0'}

    def _gdb_has_numpy(self):
        """Returns True if gdb has NumPy support, False otherwise"""
        driver = GdbMIDriver(__file__, debug=False,)
        has_numpy = driver.supports_numpy()
        driver.quit()
        return has_numpy

    def _subprocess_test_runner(self, test_mod):
        themod = f'numba.tests.gdb.{test_mod}'
        self.subprocess_test_runner(test_module=themod,
                                    test_class='Test',
                                    test_name='test',
                                    envvars=self._NUMBA_OPT_0_ENV)

    def test_basic(self):
        self._subprocess_test_runner('test_basic')

    def test_array_arg(self):
        self._subprocess_test_runner('test_array_arg')

    def test_conditional_breakpoint(self):
        self._subprocess_test_runner('test_conditional_breakpoint')

    def test_break_on_symbol(self):
        self._subprocess_test_runner('test_break_on_symbol')

    def test_break_on_symbol_version(self):
        self._subprocess_test_runner('test_break_on_symbol_version')

    def test_pretty_print(self):
        if not self._gdb_has_numpy():
            _msg = "Cannot find gdb with NumPy support"
            self.skipTest(_msg)

        self._subprocess_test_runner('test_pretty_print')


class TestGDBPrettyPrinterLogic(TestCase):
    # Tests the logic in numba.misc.gdb_print_extension.NumbaArrayPrinter
    # it's quite involved and susceptible to changes to the string
    # representation of Numba array and dtypes as it parses these
    # representations and recreates NumPy array/dtypes based on them!

    def setUp(self):
        # Patch sys.modules with mock gdb modules such that the
        # numba.misc.gdb_print_extension can import ok, the rest of the gdb
        # classes etc are implemented later

        mock_modules = {'gdb': Mock(),
                        'gdb.printing': Mock()}
        self.patched_sys = patch.dict('sys.modules', mock_modules)
        self.patched_sys.start()

        # Now sys.modules has a gdb in it, patch the gdb.selected_inferior.
        # This function should return a process wrapping object that has a
        # read_memory method that can read a memory region from a given address
        # in the process' address space.

        import gdb

        class SelectedInferior():

            def read_memory(self, data, extent):
                buf = (ct.c_char * extent).from_address(data)
                return buf.raw # this is bytes

        si = SelectedInferior()
        gdb.configure_mock(**{'selected_inferior': lambda :si})

    def tearDown(self):
        # drop the sys.modules patch
        self.patched_sys.stop()

    def get_gdb_repr(self, array):
        # Returns the gdb repr of an array as reconstructed via the
        # gdb_print_extension (should be the same as NumPy!).

        # This is the module being tested, it uses gdb and gdb.printing, both
        # of which are mocked in self.setUp()
        from numba.misc import gdb_print_extension

        # The following classes are ducks for the gdb classes (which are not
        # easily/guaranteed importable from the test suite). They implement the
        # absolute bare minimum necessary to test the gdb_print_extension.

        class DISubrange():
            def __init__(self, lo, hi):
                self._lo = lo
                self._hi = hi

            @property
            def type(self):
                return self

            def range(self):
                return self._lo, self._hi

        class DW_TAG_array_type():
            def __init__(self, lo, hi):
                self._lo, self._hi = lo, hi

            def fields(self):
                return [DISubrange(self._lo, self._hi),]

        class DIDerivedType_tuple():
            def __init__(self, the_tuple):
                self._type = DW_TAG_array_type(0, len(the_tuple) - 1)
                self._tuple = the_tuple

            @property
            def type(self):
                return self._type

            def __getitem__(self, item):
                return self._tuple[item]

        class DICompositeType_Array():
            def __init__(self, arr, type_str):
                self._arr = arr
                self._type_str = type_str

            def __getitem__(self, item):
                return getattr(self, item)

            @property
            def data(self):
                return self._arr.ctypes.data

            @property
            def itemsize(self):
                return self._arr.itemsize

            @property
            def shape(self):
                return DIDerivedType_tuple(self._arr.shape)

            @property
            def strides(self):
                return DIDerivedType_tuple(self._arr.strides)

            @property
            def type(self):
                return self._type_str

        # The type string encoded into the DWARF is the string repr of the Numba
        # type followed by the LLVM repr of the data model in brackets.
        dmm = datamodel.default_manager
        array_model = datamodel.models.ArrayModel(dmm, typeof(array))
        data_type = array_model.get_data_type()
        type_str = f"{typeof(array)} ({data_type.structure_repr()})"
        fake_gdb_arr = DICompositeType_Array(array, type_str)

        printer = gdb_print_extension.NumbaArrayPrinter(fake_gdb_arr)

        return printer.to_string().strip() # strip, there's new lines

    def check(self, array):
        gdb_printed = self.get_gdb_repr(array)
        self.assertEqual(str(gdb_printed), str(array))

    def test_np_array_printer_simple_numeric_types(self):
        # Tests printer over a selection of basic types
        n = 4
        m = 3

        for dt in (np.int8, np.uint16, np.int64, np.float32, np.complex128):
            arr = np.arange(m * n, dtype=dt).reshape(m, n)
            self.check(arr)

    def test_np_array_printer_simple_numeric_types_strided(self):
        # Tests printer over randomized strided arrays
        n_tests = 30
        np.random.seed(0)

        for _ in range(n_tests):

            shape = np.random.randint(1, high=12, size=np.random.randint(1, 5))
            tmp = np.arange(np.prod(shape)).reshape(shape)

            slices = []
            for x in shape:
                start = np.random.randint(0, x)
                # x + 3 is to ensure that sometimes the stop is beyond the
                # end of the size in a given dimension
                stop = np.random.randint(start + 1, max(start + 1, x + 3))
                step = np.random.randint(1, 3) # step as 1, 2
                strd = slice(start, stop, step)
                slices.append(strd)

            arr = tmp[tuple(slices)]
            self.check(arr)

    def test_np_array_printer_simple_structured_dtype(self):
        # Tests printer over a selection of basic types
        n = 4
        m = 3

        aligned = np.dtype([("x", np.int16), ("y", np.float64)], align=True)
        unaligned = np.dtype([("x", np.int16), ("y", np.float64)], align=False)

        for dt in (aligned, unaligned):
            arr = np.empty(m * n, dtype=dt).reshape(m, n)
            arr['x'] = np.arange(m * n, dtype=dt['x']).reshape(m, n)
            arr['y'] = 100 * np.arange(m * n, dtype=dt['y']).reshape(m, n)
            self.check(arr)

    def test_np_array_printer_chr_array(self):
        # Test unichr array
        arr = np.array(['abcde'])
        self.check(arr)

    def test_np_array_printer_unichr_structured_dtype(self):
        # Not supported yet
        n = 4
        m = 3

        dt = np.dtype([("x", '<U5'), ("y", np.float64)], align=True)
        arr = np.zeros(m * n, dtype=dt).reshape(m, n)
        rep = self.get_gdb_repr(arr)
        self.assertIn("array[Exception:", rep)
        self.assertIn("Unsupported sub-type", rep)
        self.assertIn("[unichr x 5]", rep)

    def test_np_array_printer_nested_array_structured_dtype(self):
        # Not supported yet
        n = 4
        m = 3

        dt = np.dtype([("x", np.int16, (2,)), ("y", np.float64)], align=True)
        arr = np.zeros(m * n, dtype=dt).reshape(m, n)
        rep = self.get_gdb_repr(arr)
        self.assertIn("array[Exception:", rep)
        self.assertIn("Unsupported sub-type", rep)
        self.assertIn("nestedarray(int16", rep)


if __name__ == '__main__':
    unittest.main()
