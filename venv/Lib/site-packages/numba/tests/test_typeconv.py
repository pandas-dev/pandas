import itertools

from numba.core import types
from numba.core.typeconv.typeconv import TypeManager, TypeCastingRules
from numba.core.typeconv import rules
from numba.core.typeconv import castgraph, Conversion
import unittest


class CompatibilityTestMixin(unittest.TestCase):

    def check_number_compatibility(self, check_compatible):
        b = types.boolean
        i8 = types.int8
        i16 = types.int16
        i32 = types.int32
        i64 = types.int64
        u8 = types.uint8
        u16 = types.uint16
        u32 = types.uint32
        u64 = types.uint64
        f16 = types.float16
        f32 = types.float32
        f64 = types.float64
        c64 = types.complex64
        c128 = types.complex128

        self.assertEqual(check_compatible(i32, i32), Conversion.exact)

        self.assertEqual(check_compatible(b, i8), Conversion.safe)
        self.assertEqual(check_compatible(b, u8), Conversion.safe)
        self.assertEqual(check_compatible(i8, b), Conversion.unsafe)
        self.assertEqual(check_compatible(u8, b), Conversion.unsafe)

        self.assertEqual(check_compatible(i32, i64), Conversion.promote)
        self.assertEqual(check_compatible(i32, u32), Conversion.unsafe)
        self.assertEqual(check_compatible(u32, i32), Conversion.unsafe)
        self.assertEqual(check_compatible(u32, i64), Conversion.safe)

        self.assertEqual(check_compatible(i16, f16), Conversion.unsafe)
        self.assertEqual(check_compatible(i32, f32), Conversion.unsafe)
        self.assertEqual(check_compatible(u32, f32), Conversion.unsafe)
        self.assertEqual(check_compatible(i32, f64), Conversion.safe)
        self.assertEqual(check_compatible(u32, f64), Conversion.safe)
        # Note this is inconsistent with i32 -> f32...
        self.assertEqual(check_compatible(i64, f64), Conversion.safe)
        self.assertEqual(check_compatible(u64, f64), Conversion.safe)

        self.assertEqual(check_compatible(f32, c64), Conversion.safe)
        self.assertEqual(check_compatible(f64, c128), Conversion.safe)
        self.assertEqual(check_compatible(f64, c64), Conversion.unsafe)

        # Propagated compatibility relationships
        self.assertEqual(check_compatible(i16, f64), Conversion.safe)
        self.assertEqual(check_compatible(i16, i64), Conversion.promote)
        self.assertEqual(check_compatible(i32, c64), Conversion.unsafe)
        self.assertEqual(check_compatible(i32, c128), Conversion.safe)
        self.assertEqual(check_compatible(i32, u64), Conversion.unsafe)

        for ta, tb in itertools.product(types.number_domain,
                                        types.number_domain):
            if ta in types.complex_domain and tb not in types.complex_domain:
                continue
            self.assertTrue(check_compatible(ta, tb) is not None,
                            msg="No cast from %s to %s" % (ta, tb))


class TestTypeConv(CompatibilityTestMixin, unittest.TestCase):

    def test_typeconv(self):
        tm = TypeManager()

        i32 = types.int32
        i64 = types.int64
        f32 = types.float32

        tm.set_promote(i32, i64)
        tm.set_unsafe_convert(i32, f32)

        sig = (i32, f32)
        ovs = [
            (i32, i32),
            (f32, f32),
            (i64, i64),
        ]

        # allow_unsafe = True => a conversion from i32 to f32 is chosen
        sel = tm.select_overload(sig, ovs, True, False)
        self.assertEqual(sel, 1)
        # allow_unsafe = False => no overload available
        with self.assertRaises(TypeError):
            sel = tm.select_overload(sig, ovs, False, False)

    def test_default_rules(self):
        tm = rules.default_type_manager
        self.check_number_compatibility(tm.check_compatible)

    def test_overload1(self):
        tm = rules.default_type_manager

        i32 = types.int32
        i64 = types.int64

        sig = (i64, i32, i32)
        ovs = [
            (i32, i32, i32),
            (i64, i64, i64),
        ]
        # The first overload is unsafe, the second is safe => the second
        # is always chosen, regardless of allow_unsafe.
        self.assertEqual(tm.select_overload(sig, ovs, True, False), 1)
        self.assertEqual(tm.select_overload(sig, ovs, False, False), 1)

    def test_overload2(self):
        tm = rules.default_type_manager

        i16 = types.int16
        i32 = types.int32
        i64 = types.int64

        sig = (i32, i16, i32)
        ovs = [
            # Three promotes
            (i64, i64, i64),
            # One promotes, two exact types
            (i32, i32, i32),
            # Two unsafe converts, one exact type
            (i16, i16, i16),
        ]
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=False,
                                            exact_match_required=False), 1)
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=True,
                                            exact_match_required=False), 1)

        # The same in reverse order
        ovs.reverse()
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=False,
                                            exact_match_required=False), 1)
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=True,
                                            exact_match_required=False), 1)

    def test_overload3(self):
        # Promotes should be preferred over safe converts
        tm = rules.default_type_manager

        i32 = types.int32
        i64 = types.int64
        f64 = types.float64

        sig = (i32, i32)
        ovs = [
            # Two promotes
            (i64, i64),
            # Two safe converts
            (f64, f64),
        ]
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=False,
                                            exact_match_required=False), 0)
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=True,
                                            exact_match_required=False), 0)

        # The same in reverse order
        ovs.reverse()
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=False,
                                            exact_match_required=False), 1)
        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=True,
                                            exact_match_required=False), 1)

    def test_overload4(self):
        tm = rules.default_type_manager

        i16 = types.int16
        i32 = types.int32
        i64 = types.int64
        f16 = types.float16
        f32 = types.float32

        sig = (i16, f16, f16)
        ovs = [
            # One unsafe, one promote, one exact
            (f16, f32, f16),
            # Two unsafe, one exact types
            (f32, i32, f16),
        ]

        self.assertEqual(tm.select_overload(sig, ovs, allow_unsafe=True,
                                            exact_match_required=False), 0)

    def test_type_casting_rules(self):
        tm = TypeManager()
        tcr = TypeCastingRules(tm)

        i16 = types.int16
        i32 = types.int32
        i64 = types.int64
        f64 = types.float64
        f32 = types.float32
        f16 = types.float16
        made_up = types.Dummy("made_up")

        tcr.promote_unsafe(i32, i64)
        tcr.safe_unsafe(i32, f64)
        tcr.promote_unsafe(f32, f64)
        tcr.promote_unsafe(f16, f32)
        tcr.unsafe_unsafe(i16, f16)

        def base_test():
            # As declared
            self.assertEqual(tm.check_compatible(i32, i64), Conversion.promote)
            self.assertEqual(tm.check_compatible(i32, f64), Conversion.safe)
            self.assertEqual(tm.check_compatible(f16, f32), Conversion.promote)
            self.assertEqual(tm.check_compatible(f32, f64), Conversion.promote)
            self.assertEqual(tm.check_compatible(i64, i32), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(f64, i32), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(f64, f32), Conversion.unsafe)

            # Propagated
            self.assertEqual(tm.check_compatible(i64, f64), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(f64, i64), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(i64, f32), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(i32, f32), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(f32, i32), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(i16, f16), Conversion.unsafe)
            self.assertEqual(tm.check_compatible(f16, i16), Conversion.unsafe)

        # Test base graph
        base_test()

        self.assertIsNone(tm.check_compatible(i64, made_up))
        self.assertIsNone(tm.check_compatible(i32, made_up))
        self.assertIsNone(tm.check_compatible(f32, made_up))
        self.assertIsNone(tm.check_compatible(made_up, f64))
        self.assertIsNone(tm.check_compatible(made_up, i64))

        # Add new test
        tcr.promote(f64, made_up)
        tcr.unsafe(made_up, i32)

        # Ensure the graph did not change by adding the new type
        base_test()

        # To "made up" type
        self.assertEqual(tm.check_compatible(i64, made_up), Conversion.unsafe)
        self.assertEqual(tm.check_compatible(i32, made_up), Conversion.safe)
        self.assertEqual(tm.check_compatible(f32, made_up), Conversion.promote)
        self.assertEqual(tm.check_compatible(made_up, f64), Conversion.unsafe)
        self.assertEqual(tm.check_compatible(made_up, i64), Conversion.unsafe)

    def test_castgraph_propagate(self):
        saved = []

        def callback(src, dst, rel):
            saved.append((src, dst, rel))

        tg = castgraph.TypeGraph(callback)

        i32 = types.int32
        i64 = types.int64
        f64 = types.float64
        f32 = types.float32

        tg.insert_rule(i32, i64, Conversion.promote)
        tg.insert_rule(i64, i32, Conversion.unsafe)

        saved.append(None)

        tg.insert_rule(i32, f64, Conversion.safe)
        tg.insert_rule(f64, i32, Conversion.unsafe)

        saved.append(None)

        tg.insert_rule(f32, f64, Conversion.promote)
        tg.insert_rule(f64, f32, Conversion.unsafe)

        self.assertIn((i32, i64, Conversion.promote), saved[0:2])
        self.assertIn((i64, i32, Conversion.unsafe), saved[0:2])
        self.assertIs(saved[2], None)

        self.assertIn((i32, f64, Conversion.safe), saved[3:7])
        self.assertIn((f64, i32, Conversion.unsafe), saved[3:7])
        self.assertIn((i64, f64, Conversion.unsafe), saved[3:7])
        self.assertIn((i64, f64, Conversion.unsafe), saved[3:7])
        self.assertIs(saved[7], None)

        self.assertIn((f32, f64, Conversion.promote), saved[8:14])
        self.assertIn((f64, f32, Conversion.unsafe), saved[8:14])
        self.assertIn((f32, i32, Conversion.unsafe), saved[8:14])
        self.assertIn((i32, f32, Conversion.unsafe), saved[8:14])
        self.assertIn((f32, i64, Conversion.unsafe), saved[8:14])
        self.assertIn((i64, f32, Conversion.unsafe), saved[8:14])
        self.assertEqual(len(saved[14:]), 0)


if __name__ == '__main__':
    unittest.main()
