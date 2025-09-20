# "magictoken" is used for markers as beginning and ending of example text.

import unittest

# magictoken.ex_structref_type_definition.begin
import numpy as np

from numba import njit
from numba.core import types
from numba.experimental import structref

from numba.tests.support import skip_unless_scipy


# Define a StructRef.
# `structref.register` associates the type with the default data model.
# This will also install getters and setters to the fields of
# the StructRef.
@structref.register
class MyStructType(types.StructRef):
    def preprocess_fields(self, fields):
        # This method is called by the type constructor for additional
        # preprocessing on the fields.
        # Here, we don't want the struct to take Literal types.
        return tuple((name, types.unliteral(typ)) for name, typ in fields)


# Define a Python type that can be use as a proxy to the StructRef
# allocated inside Numba. Users can construct the StructRef via
# the constructor for this type in python code and jit-code.
class MyStruct(structref.StructRefProxy):
    def __new__(cls, name, vector):
        # Overriding the __new__ method is optional, doing so
        # allows Python code to use keyword arguments,
        # or add other customized behavior.
        # The default __new__ takes `*args`.
        # IMPORTANT: Users should not override __init__.
        return structref.StructRefProxy.__new__(cls, name, vector)

    # By default, the proxy type does not reflect the attributes or
    # methods to the Python side. It is up to users to define
    # these. (This may be automated in the future.)

    @property
    def name(self):
        # To access a field, we can define a function that simply
        # return the field in jit-code.
        # The definition of MyStruct_get_name is shown later.
        return MyStruct_get_name(self)

    @property
    def vector(self):
        # The definition of MyStruct_get_vector is shown later.
        return MyStruct_get_vector(self)


@njit
def MyStruct_get_name(self):
    # In jit-code, the StructRef's attribute is exposed via
    # structref.register
    return self.name


@njit
def MyStruct_get_vector(self):
    return self.vector


# This associates the proxy with MyStructType for the given set of
# fields. Notice how we are not constraining the type of each field.
# Field types remain generic.
structref.define_proxy(MyStruct, MyStructType, ["name", "vector"])
# magictoken.ex_structref_type_definition.end


@skip_unless_scipy
class TestStructRefUsage(unittest.TestCase):
    def test_type_definition(self):
        np.random.seed(0)
        # Redirect print
        buf = []

        def print(*args):
            buf.append(args)

        # magictoken.ex_structref_type_definition_test.begin
        # Let's test our new StructRef.

        # Define one in Python
        alice = MyStruct("Alice", vector=np.random.random(3))

        # Define one in jit-code
        @njit
        def make_bob():
            bob = MyStruct("unnamed", vector=np.zeros(3))
            # Mutate the attributes
            bob.name = "Bob"
            bob.vector = np.random.random(3)
            return bob

        bob = make_bob()

        # Out: Alice: [0.5488135  0.71518937 0.60276338]
        print(f"{alice.name}: {alice.vector}")
        # Out: Bob: [0.88325739 0.73527629 0.87746707]
        print(f"{bob.name}: {bob.vector}")

        # Define a jit function to operate on the structs.
        @njit
        def distance(a, b):
            return np.linalg.norm(a.vector - b.vector)

        # Out: 0.4332647200356598
        print(distance(alice, bob))
        # magictoken.ex_structref_type_definition_test.end

        self.assertEqual(len(buf), 3)

    def test_overload_method(self):
        # magictoken.ex_structref_method.begin
        from numba.core.extending import overload_method
        from numba.core.errors import TypingError

        # Use @overload_method to add a method for
        # MyStructType.distance(other)
        # where *other* is an instance of MyStructType.
        @overload_method(MyStructType, "distance")
        def ol_distance(self, other):
            # Guard that *other* is an instance of MyStructType
            if not isinstance(other, MyStructType):
                raise TypingError(
                    f"*other* must be a {MyStructType}; got {other}"
                )

            def impl(self, other):
                return np.linalg.norm(self.vector - other.vector)

            return impl

        # Test
        @njit
        def test():
            alice = MyStruct("Alice", vector=np.random.random(3))
            bob = MyStruct("Bob", vector=np.random.random(3))
            # Use the method
            return alice.distance(bob)
        # magictoken.ex_structref_method.end

        self.assertIsInstance(test(), float)
