"""
Test mutable struct, aka, structref
"""
import warnings

import numpy as np

from numba import typed, njit, errors, typeof
from numba.core import types
from numba.experimental import structref
from numba.extending import overload_method, overload_attribute
from numba.tests.support import (
    MemoryLeakMixin, TestCase, temp_directory, override_config,
)


@structref.register
class MySimplerStructType(types.StructRef):
    """
    Test associated with this type represent the lowest level uses of structref.
    """
    pass


my_struct_ty = MySimplerStructType(
    fields=[("values", types.intp[:]), ("counter", types.intp)]
)

structref.define_boxing(MySimplerStructType, structref.StructRefProxy)


class MyStruct(structref.StructRefProxy):

    def __new__(cls, values, counter):
        # Define this method to customize the constructor.
        # The default takes `*args`. Customizing allow the use of keyword-arg.
        # The impl of the method calls `StructRefProxy.__new__`
        return structref.StructRefProxy.__new__(cls, values, counter)

    # The below defines wrappers for attributes and methods manually

    @property
    def values(self):
        return get_values(self)

    @values.setter
    def values(self, val):
        return set_values(self, val)

    @property
    def counter(self):
        return get_counter(self)

    def testme(self, arg):
        return self.values * arg + self.counter

    @property
    def prop(self):
        return self.values, self.counter

    def __hash__(self):
        return compute_fields(self)


@structref.register
class MyStructType(types.StructRef):
    """Test associated with this type represent the higher-level uses of
    structef.
    """
    pass


# Call to define_proxy is needed to register the use of `MyStruct` as a
# PyObject proxy for creating a Numba-allocated structref.
# The `MyStruct` class can then be used in both jit-code and interpreted-code.
structref.define_proxy(
    MyStruct,
    MyStructType,
    ['values', 'counter'],
)


@njit
def my_struct(values, counter):
    st = structref.new(my_struct_ty)
    my_struct_init(st, values, counter)
    return st


@njit
def my_struct_init(self, values, counter):
    self.values = values
    self.counter = counter


@njit
def ctor_by_intrinsic(vs, ctr):
    st = my_struct(vs, counter=ctr)
    st.values += st.values
    st.counter *= ctr
    return st


@njit
def ctor_by_class(vs, ctr):
    return MyStruct(values=vs, counter=ctr)


@njit
def get_values(st):
    return st.values


@njit
def set_values(st, val):
    st.values = val


@njit
def get_counter(st):
    return st.counter


@njit
def compute_fields(st):
    return st.values + st.counter


@njit
def test_structref_is():
    c = MyStruct(3, 4)
    d = MyStruct(3, 4)
    return (c is c, c is d)


class TestStructRefIs(TestCase):
    def test_is_identity(self):
        res = test_structref_is()
        self.assertTrue(res[0])
        self.assertFalse(res[1])


class TestStructRefBasic(MemoryLeakMixin, TestCase):
    def test_structref_type(self):
        sr = types.StructRef([('a', types.int64)])
        self.assertEqual(sr.field_dict['a'], types.int64)
        sr = types.StructRef([('a', types.int64), ('b', types.float64)])
        self.assertEqual(sr.field_dict['a'], types.int64)
        self.assertEqual(sr.field_dict['b'], types.float64)
        # bad case
        with self.assertRaisesRegex(ValueError,
                                    "expecting a str for field name"):
            types.StructRef([(1, types.int64)])
        with self.assertRaisesRegex(ValueError,
                                    "expecting a Numba Type for field type"):
            types.StructRef([('a', 123)])

    def test_invalid_uses(self):
        with self.assertRaisesRegex(ValueError, "cannot register"):
            structref.register(types.StructRef)
        with self.assertRaisesRegex(ValueError, "cannot register"):
            structref.define_boxing(types.StructRef, MyStruct)

    def test_MySimplerStructType(self):
        vs = np.arange(10, dtype=np.intp)
        ctr = 13

        first_expected = vs + vs
        first_got = ctor_by_intrinsic(vs, ctr)
        # the returned instance is a structref.StructRefProxy
        # but not a MyStruct
        self.assertNotIsInstance(first_got, MyStruct)
        self.assertPreciseEqual(first_expected, get_values(first_got))

        second_expected = first_expected + (ctr * ctr)
        second_got = compute_fields(first_got)
        self.assertPreciseEqual(second_expected, second_got)

    def test_MySimplerStructType_wrapper_has_no_attrs(self):
        vs = np.arange(10, dtype=np.intp)
        ctr = 13
        wrapper = ctor_by_intrinsic(vs, ctr)
        self.assertIsInstance(wrapper, structref.StructRefProxy)
        with self.assertRaisesRegex(AttributeError, 'values'):
            wrapper.values
        with self.assertRaisesRegex(AttributeError, 'counter'):
            wrapper.counter

    def test_MyStructType(self):
        vs = np.arange(10, dtype=np.float64)
        ctr = 11

        first_expected_arr = vs.copy()
        first_got = ctor_by_class(vs, ctr)
        self.assertIsInstance(first_got, MyStruct)
        self.assertPreciseEqual(first_expected_arr, first_got.values)

        second_expected = first_expected_arr + ctr
        second_got = compute_fields(first_got)
        self.assertPreciseEqual(second_expected, second_got)
        self.assertEqual(first_got.counter, ctr)

    def test_MyStructType_mixed_types(self):
        # structref constructor is generic
        @njit
        def mixed_type(x, y, m, n):
            return MyStruct(x, y), MyStruct(m, n)

        a, b = mixed_type(1, 2.3, 3.4j, (4,))
        self.assertEqual(a.values, 1)
        self.assertEqual(a.counter, 2.3)
        self.assertEqual(b.values, 3.4j)
        self.assertEqual(b.counter, (4,))

    def test_MyStructType_in_dict(self):
        td = typed.Dict()
        td['a'] = MyStruct(1, 2.3)
        self.assertEqual(td['a'].values, 1)
        self.assertEqual(td['a'].counter, 2.3)
        # overwrite
        td['a'] = MyStruct(2, 3.3)
        self.assertEqual(td['a'].values, 2)
        self.assertEqual(td['a'].counter, 3.3)
        # mutate
        td['a'].values += 10
        self.assertEqual(td['a'].values, 12)    # changed
        self.assertEqual(td['a'].counter, 3.3)  # unchanged
        # insert
        td['b'] = MyStruct(4, 5.6)

    def test_MyStructType_in_dict_mixed_type_error(self):
        self.disable_leak_check()

        td = typed.Dict()
        td['a'] = MyStruct(1, 2.3)
        self.assertEqual(td['a'].values, 1)
        self.assertEqual(td['a'].counter, 2.3)

        # ERROR: store different types
        with self.assertRaisesRegex(errors.TypingError,
                                    r"Cannot cast numba.MyStructType"):
            # because first field is not a float;
            # the second field is now an integer.
            td['b'] = MyStruct(2.3, 1)

    def test_MyStructType_hash_no_typeof_recursion(self):
        # Tests that __hash__ is not called prematurely in typeof
        # causing infinite recursion (see #8241).
        st = MyStruct(1, 2)
        typeof(st)

        self.assertEqual(hash(st), 3)


@overload_method(MyStructType, "testme")
def _ol_mystructtype_testme(self, arg):
    def impl(self, arg):
        return self.values * arg + self.counter
    return impl


@overload_attribute(MyStructType, "prop")
def _ol_mystructtype_prop(self):
    def get(self):
        return self.values, self.counter
    return get


class TestStructRefExtending(MemoryLeakMixin, TestCase):
    def test_overload_method(self):
        @njit
        def check(x):
            vs = np.arange(10, dtype=np.float64)
            ctr = 11
            obj = MyStruct(vs, ctr)
            return obj.testme(x)

        x = 3
        got = check(x)
        expect = check.py_func(x)
        self.assertPreciseEqual(got, expect)

    def test_overload_attribute(self):
        @njit
        def check():
            vs = np.arange(10, dtype=np.float64)
            ctr = 11
            obj = MyStruct(vs, ctr)
            return obj.prop

        got = check()
        expect = check.py_func()
        self.assertPreciseEqual(got, expect)


def caching_test_make(x, y):
    struct = MyStruct(values=x, counter=y)
    return struct


def caching_test_use(struct, z):
    return struct.testme(z)


class TestStructRefCaching(MemoryLeakMixin, TestCase):
    def setUp(self):
        self._cache_dir = temp_directory(TestStructRefCaching.__name__)
        self._cache_override = override_config('CACHE_DIR', self._cache_dir)
        self._cache_override.__enter__()
        warnings.simplefilter("error")
        warnings.filterwarnings(action="ignore", module="typeguard")

    def tearDown(self):
        self._cache_override.__exit__(None, None, None)
        warnings.resetwarnings()

    def test_structref_caching(self):
        def assert_cached(stats):
            self.assertEqual(len(stats.cache_hits), 1)
            self.assertEqual(len(stats.cache_misses), 0)

        def assert_not_cached(stats):
            self.assertEqual(len(stats.cache_hits), 0)
            self.assertEqual(len(stats.cache_misses), 1)

        def check(cached):
            check_make = njit(cache=True)(caching_test_make)
            check_use = njit(cache=True)(caching_test_use)
            vs = np.random.random(3)
            ctr = 17
            factor = 3
            st = check_make(vs, ctr)
            got = check_use(st, factor)
            expect = vs * factor + ctr
            self.assertPreciseEqual(got, expect)

            if cached:
                assert_cached(check_make.stats)
                assert_cached(check_use.stats)
            else:
                assert_not_cached(check_make.stats)
                assert_not_cached(check_use.stats)

        check(cached=False)
        check(cached=True)


@structref.register
class PolygonStructType(types.StructRef):

    def preprocess_fields(self, fields):
        # temp name to allow Optional instantiation
        self.name = f"numba.PolygonStructType#{id(self)}"
        fields = tuple([
            ('value', types.Optional(types.int64)),
            ('parent', types.Optional(self)),
        ])

        return fields


polygon_struct_type = PolygonStructType(fields=(
    ('value', types.Any),
    ('parent', types.Any)
))


class PolygonStruct(structref.StructRefProxy):
    def __new__(cls, value, parent):
        return structref.StructRefProxy.__new__(cls, value, parent)

    @property
    def value(self):
        return PolygonStruct_get_value(self)

    @property
    def parent(self):
        return PolygonStruct_get_parent(self)


@njit
def PolygonStruct_get_value(self):
    return self.value


@njit
def PolygonStruct_get_parent(self):
    return self.parent


structref.define_proxy(
    PolygonStruct,
    PolygonStructType,
    ["value", "parent"]
)


@overload_method(PolygonStructType, "flip")
def _ol_polygon_struct_flip(self):
    def impl(self):
        if self.value is not None:
            self.value = -self.value
    return impl


@overload_attribute(PolygonStructType, "prop")
def _ol_polygon_struct_prop(self):
    def get(self):
        return self.value, self.parent
    return get


class TestStructRefForwardTyping(MemoryLeakMixin, TestCase):
    def test_same_type_assignment(self):
        @njit
        def check(x):
            poly = PolygonStruct(None, None)
            p_poly = PolygonStruct(None, None)
            poly.value = x
            poly.parent = p_poly
            p_poly.value = x
            return poly.parent.value

        x = 11
        got = check(x)
        expect = x
        self.assertPreciseEqual(got, expect)

    def test_overload_method(self):
        @njit
        def check(x):
            poly = PolygonStruct(None, None)
            p_poly = PolygonStruct(None, None)
            poly.value = x
            poly.parent = p_poly
            p_poly.value = x
            poly.flip()
            poly.parent.flip()
            return poly.parent.value

        x = 3
        got = check(x)
        expect = -x
        self.assertPreciseEqual(got, expect)

    def test_overload_attribute(self):
        @njit
        def check():
            obj = PolygonStruct(5, None)
            return obj.prop[0]

        got = check()
        expect = 5
        self.assertPreciseEqual(got, expect)
