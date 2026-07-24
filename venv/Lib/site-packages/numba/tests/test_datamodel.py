from llvmlite import ir, binding as ll

from numba.core import types, datamodel
from numba.core.datamodel.testing import test_factory
from numba.core.datamodel.manager import DataModelManager
from numba.core.datamodel.models import OpaqueModel
import unittest


class TestBool(test_factory()):
    fe_type = types.boolean


class TestPyObject(test_factory()):
    fe_type = types.pyobject


class TestInt8(test_factory()):
    fe_type = types.int8


class TestInt16(test_factory()):
    fe_type = types.int16


class TestInt32(test_factory()):
    fe_type = types.int32


class TestInt64(test_factory()):
    fe_type = types.int64


class TestUInt8(test_factory()):
    fe_type = types.uint8


class TestUInt16(test_factory()):
    fe_type = types.uint16


class TestUInt32(test_factory()):
    fe_type = types.uint32


class TestUInt64(test_factory()):
    fe_type = types.uint64


class TestFloat(test_factory()):
    fe_type = types.float32


class TestDouble(test_factory()):
    fe_type = types.float64


class TestComplex(test_factory()):
    fe_type = types.complex64


class TestDoubleComplex(test_factory()):
    fe_type = types.complex128


class TestPointerOfInt32(test_factory()):
    fe_type = types.CPointer(types.int32)


class TestUniTupleOf2xInt32(test_factory()):
    fe_type = types.UniTuple(types.int32, 2)


class TestUniTupleEmpty(test_factory()):
    fe_type = types.UniTuple(types.int32, 0)


class TestTupleInt32Float32(test_factory()):
    fe_type = types.Tuple([types.int32, types.float32])


class TestTupleEmpty(test_factory()):
    fe_type = types.Tuple([])


class Test1DArrayOfInt32(test_factory()):
    fe_type = types.Array(types.int32, 1, 'C')


class Test2DArrayOfComplex128(test_factory()):
    fe_type = types.Array(types.complex128, 2, 'C')


class Test0DArrayOfInt32(test_factory()):
    fe_type = types.Array(types.int32, 0, 'C')


class TestArgInfo(unittest.TestCase):

    def _test_as_arguments(self, fe_args):
        """
        Test round-tripping types *fe_args* through the default data model's
        argument conversion and unpacking logic.
        """
        dmm = datamodel.default_manager
        fi = datamodel.ArgPacker(dmm, fe_args)

        module = ir.Module()
        fnty = ir.FunctionType(ir.VoidType(), [])
        function = ir.Function(module, fnty, name="test_arguments")
        builder = ir.IRBuilder()
        builder.position_at_end(function.append_basic_block())

        args = [ir.Constant(dmm.lookup(t).get_value_type(), None)
                for t in fe_args]

        # Roundtrip
        values = fi.as_arguments(builder, args)
        asargs = fi.from_arguments(builder, values)

        self.assertEqual(len(asargs), len(fe_args))
        valtys = tuple([v.type for v in values])
        self.assertEqual(valtys, fi.argument_types)

        expect_types = [a.type for a in args]
        got_types = [a.type for a in asargs]

        self.assertEqual(expect_types, got_types)

        # Assign names (check this doesn't raise)
        fi.assign_names(values, ["arg%i" for i in range(len(fe_args))])

        builder.ret_void()

        ll.parse_assembly(str(module))

    def test_int32_array_complex(self):
        fe_args = [types.int32,
                   types.Array(types.int32, 1, 'C'),
                   types.complex64]
        self._test_as_arguments(fe_args)

    def test_two_arrays(self):
        fe_args = [types.Array(types.int32, 1, 'C')] * 2
        self._test_as_arguments(fe_args)

    def test_two_0d_arrays(self):
        fe_args = [types.Array(types.int32, 0, 'C')] * 2
        self._test_as_arguments(fe_args)

    def test_tuples(self):
        fe_args = [types.UniTuple(types.int32, 2),
                   types.UniTuple(types.int32, 3)]
        self._test_as_arguments(fe_args)
        # Tuple of struct-likes
        arrty = types.Array(types.int32, 1, 'C')
        fe_args = [types.UniTuple(arrty, 2),
                   types.UniTuple(arrty, 3)]
        self._test_as_arguments(fe_args)
        # Nested tuple
        fe_args = [types.UniTuple(types.UniTuple(types.int32, 2), 3)]
        self._test_as_arguments(fe_args)

    def test_empty_tuples(self):
        # Empty tuple
        fe_args = [types.UniTuple(types.int16, 0),
                   types.Tuple(()),
                   types.int32]
        self._test_as_arguments(fe_args)

    def test_nested_empty_tuples(self):
        fe_args = [types.int32,
                   types.UniTuple(types.Tuple(()), 2),
                   types.int64]
        self._test_as_arguments(fe_args)


class TestMemInfo(unittest.TestCase):
    def setUp(self):
        self.dmm = datamodel.default_manager

    def test_number(self):
        ty = types.int32
        dm = self.dmm[ty]
        self.assertFalse(dm.contains_nrt_meminfo())

    def test_array(self):
        ty = types.int32[:]
        dm = self.dmm[ty]
        self.assertTrue(dm.contains_nrt_meminfo())

    def test_tuple_of_number(self):
        ty = types.UniTuple(dtype=types.int32, count=2)
        dm = self.dmm[ty]
        self.assertFalse(dm.contains_nrt_meminfo())

    def test_tuple_of_array(self):
        ty = types.UniTuple(dtype=types.int32[:], count=2)
        dm = self.dmm[ty]
        self.assertTrue(dm.contains_nrt_meminfo())


class TestMisc(unittest.TestCase):

    def test_issue2921(self):
        import numpy as np
        from numba import njit

        @njit
        def copy(a, b):
            for i in range(a.shape[0]):
                a[i] = b[i]

        b = np.arange(5, dtype=np.uint8).view(np.bool_)
        a = np.zeros_like(b)
        copy(a, b)
        np.testing.assert_equal(a, np.array((False,) + (True,) * 4))


class TestDMMChaining(unittest.TestCase):
    def test_basic(self):
        dmm = DataModelManager()

        class int_handler(OpaqueModel):
            pass

        class float_handler(OpaqueModel):
            pass

        dmm.register(types.Integer, int_handler)
        dmm.register(types.Float, float_handler)

        inter_dmm = DataModelManager()

        class new_int_handler(OpaqueModel):
            pass

        inter_dmm.register(types.Integer, new_int_handler)
        chained_dmm = inter_dmm.chain(dmm)

        # Check that the chained DMM has the new handler
        self.assertIsInstance(chained_dmm.lookup(types.intp), new_int_handler)
        # and not the old handler
        self.assertNotIsInstance(chained_dmm.lookup(types.intp), int_handler)
        # Check that the base DMM has the old handler
        self.assertIsInstance(dmm.lookup(types.intp), int_handler)
        # Check that float goes to the float_handler
        self.assertIsInstance(chained_dmm.lookup(types.float32), float_handler)
        self.assertIsInstance(dmm.lookup(types.float32), float_handler)
        # Check the intermediate DMM
        self.assertIsInstance(inter_dmm.lookup(types.intp), new_int_handler)
        with self.assertRaises(KeyError):
            inter_dmm.lookup(types.float32)


if __name__ == '__main__':
    unittest.main()
