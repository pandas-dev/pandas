# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import contextlib
import os
import shutil
import subprocess
import weakref
from uuid import uuid4, UUID
import sys

import pytest
try:
    import numpy as np
except ImportError:
    np = None

import pyarrow as pa
from pyarrow.vendored.version import Version


@contextlib.contextmanager
def registered_extension_type(ext_type):
    pa.register_extension_type(ext_type)
    try:
        yield
    finally:
        pa.unregister_extension_type(ext_type.extension_name)


class TinyIntType(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.int8(), 'pyarrow.tests.TinyIntType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        assert storage_type == pa.int8()
        return cls()


class IntegerType(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.int64(), 'pyarrow.tests.IntegerType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        assert storage_type == pa.int64()
        return cls()


class IntegerEmbeddedType(pa.ExtensionType):

    def __init__(self):
        super().__init__(IntegerType(), 'pyarrow.tests.IntegerType')

    def __arrow_ext_serialize__(self):
        # XXX pa.BaseExtensionType should expose C++ serialization method
        return self.storage_type.__arrow_ext_serialize__()

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        deserialized_storage_type = storage_type.__arrow_ext_deserialize__(
            serialized)
        assert deserialized_storage_type == storage_type
        return cls()


class ExampleUuidScalarType(pa.ExtensionScalar):
    def as_py(self, *, maps_as_pydicts=None):
        return None if self.value is None else UUID(bytes=self.value.as_py())


class ExampleUuidType(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.binary(16), 'pyarrow.tests.ExampleUuidType')

    def __reduce__(self):
        return ExampleUuidType, ()

    def __arrow_ext_scalar_class__(self):
        return ExampleUuidScalarType

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        return cls()


class ExampleUuidType2(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.binary(16), 'pyarrow.tests.ExampleUuidType2')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        return cls()


class LabelType(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.string(), 'pyarrow.tests.LabelType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        return cls()


class ParamExtType(pa.ExtensionType):

    def __init__(self, width):
        self._width = width
        super().__init__(pa.binary(width), 'pyarrow.tests.ParamExtType')

    @property
    def width(self):
        return self._width

    def __arrow_ext_serialize__(self):
        return str(self._width).encode()

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        width = int(serialized.decode())
        assert storage_type == pa.binary(width)
        return cls(width)


class MyStructType(pa.ExtensionType):
    storage_type = pa.struct([('left', pa.int64()),
                              ('right', pa.int64())])

    def __init__(self):
        super().__init__(self.storage_type, 'pyarrow.tests.MyStructType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        assert storage_type == cls.storage_type
        return cls()


class MyListType(pa.ExtensionType):

    def __init__(self, storage_type):
        assert isinstance(storage_type, pa.ListType)
        super().__init__(storage_type, 'pyarrow.tests.MyListType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        return cls(storage_type)


class MyFixedListType(pa.ExtensionType):

    def __init__(self, storage_type):
        assert isinstance(storage_type, pa.FixedSizeListType)
        super().__init__(storage_type, 'pyarrow.tests.MyFixedListType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        return cls(storage_type)


class AnnotatedType(pa.ExtensionType):
    """
    Generic extension type that can store any storage type.
    """

    def __init__(self, storage_type, annotation):
        self.annotation = annotation
        super().__init__(storage_type, 'pyarrow.tests.AnnotatedType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        return cls(storage_type)


def ipc_write_batch(batch):
    stream = pa.BufferOutputStream()
    writer = pa.RecordBatchStreamWriter(stream, batch.schema)
    writer.write_batch(batch)
    writer.close()
    return stream.getvalue()


def ipc_read_batch(buf):
    reader = pa.RecordBatchStreamReader(buf)
    return reader.read_next_batch()


def test_ext_type_basics():
    ty = ExampleUuidType()
    assert ty.extension_name == "pyarrow.tests.ExampleUuidType"


def test_ext_type_str():
    ty = IntegerType()
    expected = "extension<pyarrow.tests.IntegerType<IntegerType>>"
    assert str(ty) == expected
    assert pa.DataType.__str__(ty) == expected


def test_ext_type_repr():
    ty = IntegerType()
    assert repr(ty) == "IntegerType(DataType(int64))"


def test_ext_type_lifetime():
    ty = ExampleUuidType()
    wr = weakref.ref(ty)
    del ty
    assert wr() is None


def test_ext_type_storage_type():
    ty = ExampleUuidType()
    assert ty.storage_type == pa.binary(16)
    assert ty.__class__ is ExampleUuidType
    ty = ParamExtType(5)
    assert ty.storage_type == pa.binary(5)
    assert ty.__class__ is ParamExtType


def test_ext_type_byte_width():
    # Test for fixed-size binary types
    ty = pa.uuid()
    assert ty.byte_width == 16
    ty = ParamExtType(5)
    assert ty.byte_width == 5

    # Test for non fixed-size binary types
    ty = LabelType()
    with pytest.raises(ValueError, match="Non-fixed width type"):
        _ = ty.byte_width


def test_ext_type_bit_width():
    # Test for fixed-size binary types
    ty = pa.uuid()
    assert ty.bit_width == 128
    ty = ParamExtType(5)
    assert ty.bit_width == 40

    # Test for non fixed-size binary types
    ty = LabelType()
    with pytest.raises(ValueError, match="Non-fixed width type"):
        _ = ty.bit_width


def test_ext_type_as_py():
    ty = ExampleUuidType()
    expected = uuid4()
    scalar = pa.ExtensionScalar.from_storage(ty, expected.bytes)
    assert scalar.as_py() == expected

    # test array
    uuids = [uuid4() for _ in range(3)]
    storage = pa.array([uuid.bytes for uuid in uuids], type=pa.binary(16))
    arr = pa.ExtensionArray.from_storage(ty, storage)

    # Works for __get_item__
    for i, expected in enumerate(uuids):
        assert arr[i].as_py() == expected

    # Works for __iter__
    for result, expected in zip(arr, uuids):
        assert result.as_py() == expected

    # test chunked array
    data = [
        pa.ExtensionArray.from_storage(ty, storage),
        pa.ExtensionArray.from_storage(ty, storage)
    ]
    carr = pa.chunked_array(data)
    for i, expected in enumerate(uuids + uuids):
        assert carr[i].as_py() == expected

    for result, expected in zip(carr, uuids + uuids):
        assert result.as_py() == expected


def test_uuid_type_pickle(pickle_module):
    for proto in range(0, pickle_module.HIGHEST_PROTOCOL + 1):
        ty = ExampleUuidType()
        ser = pickle_module.dumps(ty, protocol=proto)
        del ty
        ty = pickle_module.loads(ser)
        wr = weakref.ref(ty)
        assert ty.extension_name == "pyarrow.tests.ExampleUuidType"
        del ty
        assert wr() is None

    for proto in range(0, pickle_module.HIGHEST_PROTOCOL + 1):
        ty = pa.uuid()
        ser = pickle_module.dumps(ty, protocol=proto)
        del ty
        ty = pickle_module.loads(ser)
        wr = weakref.ref(ty)
        assert ty.extension_name == "arrow.uuid"
        del ty
        assert wr() is None


def test_ext_type_equality():
    a = ParamExtType(5)
    b = ParamExtType(6)
    c = ParamExtType(6)
    assert a != b
    assert b == c
    d = ExampleUuidType()
    e = ExampleUuidType()
    assert a != d
    assert d == e


def test_ext_array_basics():
    ty = ParamExtType(3)
    storage = pa.array([b"foo", b"bar"], type=pa.binary(3))
    arr = pa.ExtensionArray.from_storage(ty, storage)
    arr.validate()
    assert arr.type is ty
    assert arr.storage.equals(storage)


def test_ext_array_lifetime():
    ty = ParamExtType(3)
    storage = pa.array([b"foo", b"bar"], type=pa.binary(3))
    arr = pa.ExtensionArray.from_storage(ty, storage)

    refs = [weakref.ref(ty), weakref.ref(arr), weakref.ref(storage)]
    del ty, storage, arr
    for ref in refs:
        assert ref() is None


def test_ext_array_to_pylist():
    ty = ParamExtType(3)
    storage = pa.array([b"foo", b"bar", None], type=pa.binary(3))
    arr = pa.ExtensionArray.from_storage(ty, storage)

    assert arr.to_pylist() == [b"foo", b"bar", None]


def test_ext_array_errors():
    ty = ParamExtType(4)
    storage = pa.array([b"foo", b"bar"], type=pa.binary(3))
    with pytest.raises(TypeError, match="Incompatible storage type"):
        pa.ExtensionArray.from_storage(ty, storage)


def test_ext_array_equality():
    storage1 = pa.array([b"0123456789abcdef"], type=pa.binary(16))
    storage2 = pa.array([b"0123456789abcdef"], type=pa.binary(16))
    storage3 = pa.array([], type=pa.binary(16))
    ty1 = ExampleUuidType()
    ty2 = ParamExtType(16)

    a = pa.ExtensionArray.from_storage(ty1, storage1)
    b = pa.ExtensionArray.from_storage(ty1, storage2)
    assert a.equals(b)
    c = pa.ExtensionArray.from_storage(ty1, storage3)
    assert not a.equals(c)
    d = pa.ExtensionArray.from_storage(ty2, storage1)
    assert not a.equals(d)
    e = pa.ExtensionArray.from_storage(ty2, storage2)
    assert d.equals(e)
    f = pa.ExtensionArray.from_storage(ty2, storage3)
    assert not d.equals(f)


def test_ext_array_wrap_array():
    ty = ParamExtType(3)
    storage = pa.array([b"foo", b"bar", None], type=pa.binary(3))
    arr = ty.wrap_array(storage)
    arr.validate(full=True)
    assert isinstance(arr, pa.ExtensionArray)
    assert arr.type == ty
    assert arr.storage == storage

    storage = pa.chunked_array([[b"abc", b"def"], [b"ghi"]],
                               type=pa.binary(3))
    arr = ty.wrap_array(storage)
    arr.validate(full=True)
    assert isinstance(arr, pa.ChunkedArray)
    assert arr.type == ty
    assert arr.chunk(0).storage == storage.chunk(0)
    assert arr.chunk(1).storage == storage.chunk(1)

    # Wrong storage type
    storage = pa.array([b"foo", b"bar", None])
    with pytest.raises(TypeError, match="Incompatible storage type"):
        ty.wrap_array(storage)

    # Not an array or chunked array
    with pytest.raises(TypeError, match="Expected array or chunked array"):
        ty.wrap_array(None)


def test_ext_scalar_from_array():
    data = [b"0123456789abcdef", b"0123456789abcdef",
            b"zyxwvutsrqponmlk", None]
    storage = pa.array(data, type=pa.binary(16))
    ty1 = ExampleUuidType()
    ty2 = ParamExtType(16)
    ty3 = ExampleUuidType2()

    a = pa.ExtensionArray.from_storage(ty1, storage)
    b = pa.ExtensionArray.from_storage(ty2, storage)
    c = pa.ExtensionArray.from_storage(ty3, storage)

    scalars_a = list(a)
    assert len(scalars_a) == 4

    assert ty1.__arrow_ext_scalar_class__() == ExampleUuidScalarType
    assert isinstance(a[0], ExampleUuidScalarType)
    assert isinstance(scalars_a[0], ExampleUuidScalarType)

    for s, val in zip(scalars_a, data):
        assert isinstance(s, pa.ExtensionScalar)
        assert s.is_valid == (val is not None)
        assert s.type == ty1
        if val is not None:
            assert s.value == pa.scalar(val, storage.type)
            assert s.as_py() == UUID(bytes=val)
        else:
            assert s.value is None

    scalars_b = list(b)
    assert len(scalars_b) == 4

    for sa, sb in zip(scalars_a, scalars_b):
        assert isinstance(sb, pa.ExtensionScalar)
        assert sa.is_valid == sb.is_valid
        if sa.as_py() is None:
            assert sa.as_py() == sb.as_py()
        else:
            assert sa.as_py().bytes == sb.as_py()
        assert sa != sb

    scalars_c = list(c)
    assert len(scalars_c) == 4

    for s, val in zip(scalars_c, data):
        assert isinstance(s, pa.ExtensionScalar)
        assert s.is_valid == (val is not None)
        assert s.type == ty3
        if val is not None:
            assert s.value == pa.scalar(val, storage.type)
            assert s.as_py() == val
        else:
            assert s.value is None

    assert a.to_pylist() == [UUID(bytes=x) if x else None for x in data]


def test_ext_scalar_from_storage():
    ty = ExampleUuidType()

    s = pa.ExtensionScalar.from_storage(ty, None)
    assert isinstance(s, pa.ExtensionScalar)
    assert s.type == ty
    assert s.is_valid is False
    assert s.value is None

    s = pa.ExtensionScalar.from_storage(ty, b"0123456789abcdef")
    assert isinstance(s, pa.ExtensionScalar)
    assert s.type == ty
    assert s.is_valid is True
    assert s.value == pa.scalar(b"0123456789abcdef", ty.storage_type)

    s = pa.ExtensionScalar.from_storage(ty, pa.scalar(None, ty.storage_type))
    assert isinstance(s, pa.ExtensionScalar)
    assert s.type == ty
    assert s.is_valid is False
    assert s.value is None

    s = pa.ExtensionScalar.from_storage(
        ty, pa.scalar(b"0123456789abcdef", ty.storage_type))
    assert isinstance(s, pa.ExtensionScalar)
    assert s.type == ty
    assert s.is_valid is True
    assert s.value == pa.scalar(b"0123456789abcdef", ty.storage_type)


def test_ext_array_pickling(pickle_module):
    for proto in range(0, pickle_module.HIGHEST_PROTOCOL + 1):
        ty = ParamExtType(3)
        storage = pa.array([b"foo", b"bar"], type=pa.binary(3))
        arr = pa.ExtensionArray.from_storage(ty, storage)
        ser = pickle_module.dumps(arr, protocol=proto)
        del ty, storage, arr
        arr = pickle_module.loads(ser)
        arr.validate()
        assert isinstance(arr, pa.ExtensionArray)
        assert arr.type == ParamExtType(3)
        assert arr.type.storage_type == pa.binary(3)
        assert arr.storage.type == pa.binary(3)
        assert arr.storage.to_pylist() == [b"foo", b"bar"]


@pytest.mark.numpy
def test_ext_array_conversion_to_numpy():
    storage1 = pa.array([1, 2, 3], type=pa.int64())
    storage2 = pa.array([b"123", b"456", b"789"], type=pa.binary(3))
    ty1 = IntegerType()
    ty2 = ParamExtType(3)

    arr1 = pa.ExtensionArray.from_storage(ty1, storage1)
    arr2 = pa.ExtensionArray.from_storage(ty2, storage2)

    result = arr1.to_numpy()
    expected = np.array([1, 2, 3], dtype="int64")
    np.testing.assert_array_equal(result, expected)

    with pytest.raises(ValueError, match="zero_copy_only was True"):
        arr2.to_numpy()
    result = arr2.to_numpy(zero_copy_only=False)
    expected = np.array([b"123", b"456", b"789"])
    np.testing.assert_array_equal(result, expected)


@pytest.mark.pandas
def test_ext_array_conversion_to_pandas():
    import pandas as pd

    storage1 = pa.array([1, 2, 3], type=pa.int64())
    storage2 = pa.array([b"123", b"456", b"789"], type=pa.binary(3))
    ty1 = IntegerType()
    ty2 = ParamExtType(3)

    arr1 = pa.ExtensionArray.from_storage(ty1, storage1)
    arr2 = pa.ExtensionArray.from_storage(ty2, storage2)

    result = arr1.to_pandas()
    expected = pd.Series([1, 2, 3], dtype="int64")
    pd.testing.assert_series_equal(result, expected)

    result = arr2.to_pandas()
    expected = pd.Series([b"123", b"456", b"789"], dtype=object)
    pd.testing.assert_series_equal(result, expected)


@pytest.fixture
def struct_w_ext_data():
    storage1 = pa.array([1, 2, 3], type=pa.int64())
    storage2 = pa.array([b"123", b"456", b"789"], type=pa.binary(3))
    ty1 = IntegerType()
    ty2 = ParamExtType(3)

    arr1 = pa.ExtensionArray.from_storage(ty1, storage1)
    arr2 = pa.ExtensionArray.from_storage(ty2, storage2)

    sarr1 = pa.StructArray.from_arrays([arr1], ["f0"])
    sarr2 = pa.StructArray.from_arrays([arr2], ["f1"])

    return [sarr1, sarr2]


@pytest.mark.numpy
def test_struct_w_ext_array_to_numpy(struct_w_ext_data):
    # ARROW-15291
    # Check that we don't segfault when trying to build
    # a numpy array from a StructArray with a field being
    # an ExtensionArray

    result = struct_w_ext_data[0].to_numpy(zero_copy_only=False)
    expected = np.array([{'f0': 1}, {'f0': 2},
                         {'f0': 3}], dtype=object)
    np.testing.assert_array_equal(result, expected)

    result = struct_w_ext_data[1].to_numpy(zero_copy_only=False)
    expected = np.array([{'f1': b'123'}, {'f1': b'456'},
                         {'f1': b'789'}], dtype=object)
    np.testing.assert_array_equal(result, expected)


@pytest.mark.pandas
def test_struct_w_ext_array_to_pandas(struct_w_ext_data):
    # ARROW-15291
    # Check that we don't segfault when trying to build
    # a Pandas dataframe from a StructArray with a field
    # being an ExtensionArray
    import pandas as pd

    result = struct_w_ext_data[0].to_pandas()
    expected = pd.Series([{'f0': 1}, {'f0': 2},
                         {'f0': 3}], dtype=object)
    pd.testing.assert_series_equal(result, expected)

    result = struct_w_ext_data[1].to_pandas()
    expected = pd.Series([{'f1': b'123'}, {'f1': b'456'},
                         {'f1': b'789'}], dtype=object)
    pd.testing.assert_series_equal(result, expected)


def test_cast_kernel_on_extension_arrays():
    # test array casting
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(IntegerType(), storage)

    # test that no allocation happens during identity cast
    allocated_before_cast = pa.total_allocated_bytes()
    casted = arr.cast(pa.int64())
    assert pa.total_allocated_bytes() == allocated_before_cast

    cases = [
        (pa.int64(), pa.Int64Array),
        (pa.int32(), pa.Int32Array),
        (pa.int16(), pa.Int16Array),
        (pa.uint64(), pa.UInt64Array),
        (pa.uint32(), pa.UInt32Array),
        (pa.uint16(), pa.UInt16Array)
    ]
    for typ, klass in cases:
        casted = arr.cast(typ)
        assert casted.type == typ
        assert isinstance(casted, klass)

    # test chunked array casting
    arr = pa.chunked_array([arr, arr])
    casted = arr.cast(pa.int16())
    assert casted.type == pa.int16()
    assert isinstance(casted, pa.ChunkedArray)


@pytest.mark.parametrize("data,ty", (
    ([1, 2], pa.int32),
    ([1, 2], pa.int64),
    (["1", "2"], pa.string),
    ([b"1", b"2"], pa.binary),
    ([1.0, 2.0], pa.float32),
    ([1.0, 2.0], pa.float64)
))
def test_casting_to_extension_type(data, ty):
    arr = pa.array(data, ty())
    out = arr.cast(IntegerType())
    assert isinstance(out, pa.ExtensionArray)
    assert out.type == IntegerType()
    assert out.to_pylist() == [1, 2]


def test_cast_between_extension_types():
    array = pa.array([1, 2, 3], pa.int8())

    tiny_int_arr = array.cast(TinyIntType())
    assert tiny_int_arr.type == TinyIntType()

    # Casting between extension types w/ different storage types not okay.
    msg = ("Casting from 'extension<.*?<TinyIntType>>' "
           "to different extension type "
           "'extension<.*?<IntegerType>>' not permitted. "
           "One can first cast to the storage type, "
           "then to the extension type."
           )
    with pytest.raises(TypeError, match=msg):
        tiny_int_arr.cast(IntegerType())
    tiny_int_arr.cast(pa.int64()).cast(IntegerType())

    # Between the same extension types is okay
    array = pa.array([b'1' * 16, b'2' * 16], pa.binary(16)).cast(ExampleUuidType())
    out = array.cast(ExampleUuidType())
    assert out.type == ExampleUuidType()

    # Will still fail casting between extensions who share storage type,
    # can only cast between exactly the same extension types.
    with pytest.raises(TypeError, match='Casting from *'):
        array.cast(ExampleUuidType2())


def test_cast_to_extension_with_extension_storage():
    # Test casting directly, and IntegerType -> IntegerEmbeddedType
    array = pa.array([1, 2, 3], pa.int64())
    array.cast(IntegerEmbeddedType())
    array.cast(IntegerType()).cast(IntegerEmbeddedType())


@pytest.mark.parametrize("data,type_factory", (
    # list<extension>
    ([[1, 2, 3]], lambda: pa.list_(IntegerType())),
    # struct<extension>
    ([{"foo": 1}], lambda: pa.struct([("foo", IntegerType())])),
    # list<struct<extension>>
    ([[{"foo": 1}]], lambda: pa.list_(pa.struct([("foo", IntegerType())]))),
    # struct<list<extension>>
    ([{"foo": [1, 2, 3]}], lambda: pa.struct(
        [("foo", pa.list_(IntegerType()))])),
))
def test_cast_nested_extension_types(data, type_factory):
    ty = type_factory()
    a = pa.array(data)
    b = a.cast(ty)
    assert b.type == ty  # casted to target extension
    assert b.cast(a.type)  # and can cast back


def test_casting_dict_array_to_extension_type():
    storage = pa.array([b"0123456789abcdef"], type=pa.binary(16))
    arr = pa.ExtensionArray.from_storage(ExampleUuidType(), storage)
    dict_arr = pa.DictionaryArray.from_arrays(pa.array([0, 0], pa.int32()),
                                              arr)
    out = dict_arr.cast(ExampleUuidType())
    assert isinstance(out, pa.ExtensionArray)
    assert out.to_pylist() == [UUID('30313233-3435-3637-3839-616263646566'),
                               UUID('30313233-3435-3637-3839-616263646566')]


def test_cast_to_extension_with_nested_storage():
    # https://github.com/apache/arrow/issues/37669

    # With fixed-size list
    array = pa.array([[1, 2], [3, 4], [5, 6]], pa.list_(pa.float64(), 2))
    result = array.cast(MyFixedListType(pa.list_(pa.float64(), 2)))
    expected = pa.ExtensionArray.from_storage(MyFixedListType(array.type), array)
    assert result.equals(expected)

    ext_type = MyFixedListType(pa.list_(pa.float32(), 2))
    result = array.cast(ext_type)
    expected = pa.ExtensionArray.from_storage(
        ext_type, array.cast(ext_type.storage_type)
    )
    assert result.equals(expected)

    # With variable-size list
    array = pa.array([[1, 2], [3], [4, 5, 6]], pa.list_(pa.float64()))
    result = array.cast(MyListType(pa.list_(pa.float64())))
    expected = pa.ExtensionArray.from_storage(MyListType(array.type), array)
    assert result.equals(expected)

    ext_type = MyListType(pa.list_(pa.float32()))
    result = array.cast(ext_type)
    expected = pa.ExtensionArray.from_storage(
        ext_type, array.cast(ext_type.storage_type)
    )
    assert result.equals(expected)


def test_concat():
    arr1 = pa.array([1, 2, 3], IntegerType())
    arr2 = pa.array([4, 5, 6], IntegerType())

    result = pa.concat_arrays([arr1, arr2])
    expected = pa.array([1, 2, 3, 4, 5, 6], IntegerType())
    assert result.equals(expected)

    # nested in a struct
    struct_arr1 = pa.StructArray.from_arrays([arr1], names=["a"])
    struct_arr2 = pa.StructArray.from_arrays([arr2], names=["a"])
    result = pa.concat_arrays([struct_arr1, struct_arr2])
    expected = pa.StructArray.from_arrays([expected], names=["a"])
    assert result.equals(expected)


def test_null_storage_type():
    ext_type = AnnotatedType(pa.null(), {"key": "value"})
    storage = pa.array([None] * 10, pa.null())
    arr = pa.ExtensionArray.from_storage(ext_type, storage)
    assert arr.null_count == 10
    arr.validate(full=True)


def example_batch():
    ty = ParamExtType(3)
    storage = pa.array([b"foo", b"bar"], type=pa.binary(3))
    arr = pa.ExtensionArray.from_storage(ty, storage)
    return pa.RecordBatch.from_arrays([arr], ["exts"])


def check_example_batch(batch, *, expect_extension):
    arr = batch.column(0)
    if expect_extension:
        assert isinstance(arr, pa.ExtensionArray)
        assert arr.type.storage_type == pa.binary(3)
        assert arr.storage.to_pylist() == [b"foo", b"bar"]
    else:
        assert arr.type == pa.binary(3)
        assert arr.to_pylist() == [b"foo", b"bar"]
    return arr


def test_ipc_unregistered():
    batch = example_batch()
    buf = ipc_write_batch(batch)
    del batch

    batch = ipc_read_batch(buf)
    batch.validate(full=True)
    check_example_batch(batch, expect_extension=False)


def test_ipc_registered():
    with registered_extension_type(ParamExtType(1)):
        batch = example_batch()
        buf = ipc_write_batch(batch)
        del batch

        batch = ipc_read_batch(buf)
        batch.validate(full=True)
        arr = check_example_batch(batch, expect_extension=True)
        assert arr.type == ParamExtType(3)


class PeriodArray(pa.ExtensionArray):
    pass


class PeriodType(pa.ExtensionType):
    def __init__(self, freq):
        # attributes need to be set first before calling
        # super init (as that calls serialize)
        self._freq = freq
        pa.ExtensionType.__init__(self, pa.int64(), 'test.period')

    @property
    def freq(self):
        return self._freq

    def __arrow_ext_serialize__(self):
        return f"freq={self.freq}".encode()

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        serialized = serialized.decode()
        assert serialized.startswith("freq=")
        freq = serialized.split('=')[1]
        return PeriodType(freq)

    def __eq__(self, other):
        if isinstance(other, pa.BaseExtensionType):
            return (isinstance(self, type(other)) and
                    self.freq == other.freq)
        else:
            return NotImplemented


class PeriodTypeWithClass(PeriodType):
    def __init__(self, freq):
        PeriodType.__init__(self, freq)

    def __arrow_ext_class__(self):
        return PeriodArray

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        freq = PeriodType.__arrow_ext_deserialize__(
            storage_type, serialized).freq
        return PeriodTypeWithClass(freq)


class PeriodTypeWithToPandasDtype(PeriodType):
    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        freq = PeriodType.__arrow_ext_deserialize__(
            storage_type, serialized).freq
        return PeriodTypeWithToPandasDtype(freq)

    def to_pandas_dtype(self):
        import pandas as pd
        return pd.PeriodDtype(freq=self.freq)


@pytest.fixture(params=[PeriodType('D'),
                        PeriodTypeWithClass('D'),
                        PeriodTypeWithToPandasDtype('D')])
def registered_period_type(request):
    # setup
    period_type = request.param
    period_class = period_type.__arrow_ext_class__()
    pa.register_extension_type(period_type)
    yield period_type, period_class
    # teardown
    try:
        pa.unregister_extension_type('test.period')
    except KeyError:
        pass


def test_generic_ext_type():
    period_type = PeriodType('D')
    assert period_type.extension_name == "test.period"
    assert period_type.storage_type == pa.int64()
    # default ext_class expected.
    assert period_type.__arrow_ext_class__() == pa.ExtensionArray


def test_generic_ext_type_ipc(registered_period_type):
    period_type, period_class = registered_period_type
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(period_type, storage)
    batch = pa.RecordBatch.from_arrays([arr], ["ext"])
    # check the built array has exactly the expected clss
    assert isinstance(arr, period_class)

    buf = ipc_write_batch(batch)
    del batch
    batch = ipc_read_batch(buf)

    result = batch.column(0)
    # check the deserialized array class is the expected one
    assert isinstance(result, period_class)
    assert result.type.extension_name == "test.period"
    assert arr.storage.to_pylist() == [1, 2, 3, 4]

    # we get back an actual PeriodType
    assert isinstance(result.type, PeriodType)
    assert result.type.freq == 'D'
    assert result.type == period_type

    # using different parametrization as how it was registered
    period_type_H = period_type.__class__('H')
    assert period_type_H.extension_name == "test.period"
    assert period_type_H.freq == 'H'

    arr = pa.ExtensionArray.from_storage(period_type_H, storage)
    batch = pa.RecordBatch.from_arrays([arr], ["ext"])

    buf = ipc_write_batch(batch)
    del batch
    batch = ipc_read_batch(buf)
    result = batch.column(0)
    assert isinstance(result.type, PeriodType)
    assert result.type.freq == 'H'
    assert isinstance(result, period_class)


def test_generic_ext_type_ipc_unknown(registered_period_type):
    period_type, _ = registered_period_type
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(period_type, storage)
    batch = pa.RecordBatch.from_arrays([arr], ["ext"])

    buf = ipc_write_batch(batch)
    del batch

    # unregister type before loading again => reading unknown extension type
    # as plain array (but metadata in schema's field are preserved)
    pa.unregister_extension_type('test.period')

    batch = ipc_read_batch(buf)
    result = batch.column(0)

    assert isinstance(result, pa.Int64Array)
    ext_field = batch.schema.field('ext')
    assert ext_field.metadata == {
        b'ARROW:extension:metadata': b'freq=D',
        b'ARROW:extension:name': b'test.period'
    }


def test_generic_ext_type_equality():
    period_type = PeriodType('D')
    assert period_type.extension_name == "test.period"

    period_type2 = PeriodType('D')
    period_type3 = PeriodType('H')
    assert period_type == period_type2
    assert not period_type == period_type3


def test_generic_ext_type_pickling(registered_period_type, pickle_module):
    # GH-36038
    for proto in range(0, pickle_module.HIGHEST_PROTOCOL + 1):
        period_type, _ = registered_period_type
        ser = pickle_module.dumps(period_type, protocol=proto)
        period_type_pickled = pickle_module.loads(ser)
        assert period_type == period_type_pickled


def test_generic_ext_array_pickling(registered_period_type, pickle_module):
    for proto in range(0, pickle_module.HIGHEST_PROTOCOL + 1):
        period_type, _ = registered_period_type
        storage = pa.array([1, 2, 3, 4], pa.int64())
        arr = pa.ExtensionArray.from_storage(period_type, storage)
        ser = pickle_module.dumps(arr, protocol=proto)
        del storage, arr
        arr = pickle_module.loads(ser)
        arr.validate()
        assert isinstance(arr, pa.ExtensionArray)
        assert arr.type == period_type
        assert arr.type.storage_type == pa.int64()
        assert arr.storage.type == pa.int64()
        assert arr.storage.to_pylist() == [1, 2, 3, 4]


def test_generic_ext_type_register(registered_period_type):
    # test that trying to register other type does not segfault
    with pytest.raises(TypeError):
        pa.register_extension_type(pa.string())

    # register second time raises KeyError
    period_type = PeriodType('D')
    with pytest.raises(KeyError):
        pa.register_extension_type(period_type)


@pytest.mark.parquet
def test_parquet_period(tmpdir, registered_period_type):
    # Parquet support for primitive extension types
    period_type, period_class = registered_period_type
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(period_type, storage)
    table = pa.table([arr], names=["ext"])

    import pyarrow.parquet as pq

    filename = tmpdir / 'period_extension_type.parquet'
    pq.write_table(table, filename)

    # Stored in parquet as storage type but with extension metadata saved
    # in the serialized arrow schema
    meta = pq.read_metadata(filename)
    assert meta.schema.column(0).physical_type == "INT64"
    assert b"ARROW:schema" in meta.metadata

    import base64
    decoded_schema = base64.b64decode(meta.metadata[b"ARROW:schema"])
    schema = pa.ipc.read_schema(pa.BufferReader(decoded_schema))
    # Since the type could be reconstructed, the extension type metadata is
    # absent.
    assert schema.field("ext").metadata == {}

    # When reading in, properly create extension type if it is registered
    result = pq.read_table(filename)
    result.validate(full=True)
    assert result.schema.field("ext").type == period_type
    assert result.schema.field("ext").metadata == {}
    # Get the exact array class defined by the registered type.
    result_array = result.column("ext").chunk(0)
    assert type(result_array) is period_class

    # When the type is not registered, read in as storage type
    pa.unregister_extension_type(period_type.extension_name)
    result = pq.read_table(filename)
    result.validate(full=True)
    assert result.schema.field("ext").type == pa.int64()
    # The extension metadata is present for roundtripping.
    assert result.schema.field("ext").metadata == {
        b'ARROW:extension:metadata': b'freq=D',
        b'ARROW:extension:name': b'test.period'
    }


@pytest.mark.parquet
def test_parquet_extension_with_nested_storage(tmpdir):
    # Parquet support for extension types with nested storage type
    import pyarrow.parquet as pq

    struct_array = pa.StructArray.from_arrays(
        [pa.array([0, 1], type="int64"), pa.array([4, 5], type="int64")],
        names=["left", "right"])
    list_array = pa.array([[1, 2, 3], [4, 5]], type=pa.list_(pa.int32()))

    mystruct_array = pa.ExtensionArray.from_storage(MyStructType(),
                                                    struct_array)
    mylist_array = pa.ExtensionArray.from_storage(
        MyListType(list_array.type), list_array)

    orig_table = pa.table({'structs': mystruct_array,
                           'lists': mylist_array})
    filename = tmpdir / 'nested_extension_storage.parquet'
    pq.write_table(orig_table, filename)

    # Unregistered
    table = pq.read_table(filename)
    table.validate(full=True)
    assert table.column('structs').type == struct_array.type
    assert table.column('structs').combine_chunks() == struct_array
    assert table.column('lists').type == list_array.type
    assert table.column('lists').combine_chunks() == list_array

    # Registered
    with registered_extension_type(mystruct_array.type):
        with registered_extension_type(mylist_array.type):
            table = pq.read_table(filename)
            table.validate(full=True)
            assert table.column('structs').type == mystruct_array.type
            assert table.column('lists').type == mylist_array.type
            assert table == orig_table

            # Cannot select a subfield of an extension type with
            # a struct storage type.
            with pytest.raises(pa.ArrowInvalid,
                               match='without all of its fields'):
                pq.ParquetFile(filename).read(columns=['structs.left'])


@pytest.mark.parquet
def test_parquet_nested_extension(tmpdir):
    # Parquet support for extension types nested in struct or list
    import pyarrow.parquet as pq

    ext_type = IntegerType()
    storage = pa.array([4, 5, 6, 7], type=pa.int64())
    ext_array = pa.ExtensionArray.from_storage(ext_type, storage)

    # Struct of extensions
    struct_array = pa.StructArray.from_arrays(
        [storage, ext_array],
        names=['ints', 'exts'])

    orig_table = pa.table({'structs': struct_array})
    filename = tmpdir / 'struct_of_ext.parquet'
    pq.write_table(orig_table, filename)

    table = pq.read_table(filename)
    table.validate(full=True)
    assert table.column(0).type == pa.struct({'ints': pa.int64(),
                                              'exts': pa.int64()})
    with registered_extension_type(ext_type):
        table = pq.read_table(filename)
        table.validate(full=True)
        assert table.column(0).type == struct_array.type
        assert table == orig_table

    # List of extensions
    list_array = pa.ListArray.from_arrays([0, 1, None, 3], ext_array)

    orig_table = pa.table({'lists': list_array})
    filename = tmpdir / 'list_of_ext.parquet'
    pq.write_table(orig_table, filename)

    table = pq.read_table(filename)
    table.validate(full=True)
    assert table.column(0).type == pa.list_(pa.int64())
    with registered_extension_type(ext_type):
        table = pq.read_table(filename)
        table.validate(full=True)
        assert table.column(0).type == list_array.type
        assert table == orig_table

    # Large list of extensions
    list_array = pa.LargeListArray.from_arrays([0, 1, None, 3], ext_array)

    orig_table = pa.table({'lists': list_array})
    filename = tmpdir / 'list_of_ext.parquet'
    pq.write_table(orig_table, filename)

    table = pq.read_table(filename)
    table.validate(full=True)
    assert table.column(0).type == pa.large_list(pa.int64())
    with registered_extension_type(ext_type):
        table = pq.read_table(filename)
        table.validate(full=True)
        assert table.column(0).type == list_array.type
        assert table == orig_table


@pytest.mark.parquet
def test_parquet_extension_nested_in_extension(tmpdir):
    # Parquet support for extension<list<extension>>
    import pyarrow.parquet as pq

    inner_ext_type = IntegerType()
    inner_storage = pa.array([4, 5, 6, 7], type=pa.int64())
    inner_ext_array = pa.ExtensionArray.from_storage(inner_ext_type,
                                                     inner_storage)

    list_array = pa.ListArray.from_arrays([0, 1, None, 3], inner_ext_array)
    mylist_array = pa.ExtensionArray.from_storage(
        MyListType(list_array.type), list_array)

    orig_table = pa.table({'lists': mylist_array})
    filename = tmpdir / 'ext_of_list_of_ext.parquet'
    pq.write_table(orig_table, filename)

    table = pq.read_table(filename)
    assert table.column(0).type == pa.list_(pa.int64())
    with registered_extension_type(mylist_array.type):
        with registered_extension_type(inner_ext_array.type):
            table = pq.read_table(filename)
            assert table.column(0).type == mylist_array.type
            assert table == orig_table


@pytest.mark.numpy
def test_to_numpy():
    period_type = PeriodType('D')
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(period_type, storage)

    expected = storage.to_numpy()
    result = arr.to_numpy()
    np.testing.assert_array_equal(result, expected)

    result = np.asarray(arr)
    np.testing.assert_array_equal(result, expected)

    # chunked array
    a1 = pa.chunked_array([arr, arr])
    a2 = pa.chunked_array([arr, arr], type=period_type)
    expected = np.hstack([expected, expected])

    for charr in [a1, a2]:
        assert charr.type == period_type
        for result in [np.asarray(charr), charr.to_numpy()]:
            assert result.dtype == np.int64
            np.testing.assert_array_equal(result, expected)

    # zero chunks
    charr = pa.chunked_array([], type=period_type)
    assert charr.type == period_type

    for result in [np.asarray(charr), charr.to_numpy()]:
        assert result.dtype == np.int64
        np.testing.assert_array_equal(result, np.array([], dtype='int64'))


def test_empty_take():
    # https://issues.apache.org/jira/browse/ARROW-13474
    ext_type = IntegerType()
    storage = pa.array([], type=pa.int64())
    empty_arr = pa.ExtensionArray.from_storage(ext_type, storage)

    result = empty_arr.filter(pa.array([], pa.bool_()))
    assert len(result) == 0
    assert result.equals(empty_arr)

    result = empty_arr.take(pa.array([], pa.int32()))
    assert len(result) == 0
    assert result.equals(empty_arr)


@pytest.mark.parametrize("data,ty", (
    ([1, 2, 3], IntegerType),
    (["cat", "dog", "horse"], LabelType)
))
@pytest.mark.parametrize(
    "into", [
        pytest.param("to_numpy", marks=pytest.mark.numpy),
        pytest.param("to_pandas", marks=pytest.mark.pandas)
    ]
)
def test_extension_array_to_numpy_pandas(data, ty, into):
    storage = pa.array(data)
    ext_arr = pa.ExtensionArray.from_storage(ty(), storage)
    offsets = pa.array([0, 1, 2, 3])
    list_arr = pa.ListArray.from_arrays(offsets, ext_arr)
    result = getattr(list_arr, into)(zero_copy_only=False)

    list_arr_storage_type = list_arr.cast(pa.list_(ext_arr.type.storage_type))
    expected = getattr(list_arr_storage_type, into)(zero_copy_only=False)
    if into == "to_pandas":
        assert result.equals(expected)
    else:
        assert np.array_equal(result, expected)


@pytest.mark.numpy
def test_array_constructor():
    ext_type = IntegerType()
    storage = pa.array([1, 2, 3], type=pa.int64())
    expected = pa.ExtensionArray.from_storage(ext_type, storage)

    result = pa.array([1, 2, 3], type=IntegerType())
    assert result.equals(expected)

    result = pa.array(np.array([1, 2, 3]), type=IntegerType())
    assert result.equals(expected)

    result = pa.array(np.array([1.0, 2.0, 3.0]), type=IntegerType())
    assert result.equals(expected)


@pytest.mark.pandas
def test_array_constructor_from_pandas():
    import pandas as pd

    ext_type = IntegerType()
    storage = pa.array([1, 2, 3], type=pa.int64())
    expected = pa.ExtensionArray.from_storage(ext_type, storage)

    result = pa.array(pd.Series([1, 2, 3]), type=IntegerType())
    assert result.equals(expected)

    result = pa.array(
        pd.Series([1, 2, 3], dtype="category"), type=IntegerType()
    )
    assert result.equals(expected)


@pytest.mark.numpy
@pytest.mark.cython
def test_cpp_extension_in_python(tmpdir):
    from .test_cython import (
        setup_template, compiler_opts, test_ld_path, test_util, here)
    with tmpdir.as_cwd():
        # Set up temporary workspace
        pyx_file = 'extensions.pyx'
        shutil.copyfile(os.path.join(here, pyx_file),
                        os.path.join(str(tmpdir), pyx_file))
        # Create setup.py file
        setup_code = setup_template.format(pyx_file=pyx_file,
                                           compiler_opts=compiler_opts,
                                           test_ld_path=test_ld_path)
        with open('setup.py', 'w') as f:
            f.write(setup_code)

        subprocess_env = test_util.get_modified_env_with_pythonpath()

        # Compile extension module
        subprocess.check_call([sys.executable, 'setup.py',
                               'build_ext', '--inplace'],
                              env=subprocess_env)

    sys.path.insert(0, str(tmpdir))
    mod = __import__('extensions')

    uuid_type = mod._make_uuid_type()
    assert uuid_type.extension_name == "example-uuid"
    assert uuid_type.storage_type == pa.binary(16)

    array = mod._make_uuid_array()
    assert array.type == uuid_type
    assert array.to_pylist() == [b'abcdefghijklmno0', b'0onmlkjihgfedcba']
    assert array[0].as_py() == b'abcdefghijklmno0'
    assert array[1].as_py() == b'0onmlkjihgfedcba'

    buf = ipc_write_batch(pa.RecordBatch.from_arrays([array], ["example-uuid"]))

    batch = ipc_read_batch(buf)
    reconstructed_array = batch.column(0)
    assert reconstructed_array.type == uuid_type
    assert reconstructed_array == array


def test_uuid_extension():
    data = [b"0123456789abcdef", b"0123456789abcdef",
            b"zyxwvutsrqponmlk", None]

    uuid_type = pa.uuid()
    assert uuid_type.extension_name == "arrow.uuid"
    assert uuid_type.storage_type == pa.binary(16)
    assert uuid_type.__class__ is pa.UuidType

    storage = pa.array(data, pa.binary(16))
    array = pa.ExtensionArray.from_storage(uuid_type, storage)
    assert array.type == uuid_type

    assert array.to_pylist() == [x if x is None else UUID(bytes=x) for x in data]
    assert array[0].as_py() == UUID(bytes=data[0])
    assert array[3].as_py() is None

    buf = ipc_write_batch(pa.RecordBatch.from_arrays([array], ["uuid"]))

    batch = ipc_read_batch(buf)
    reconstructed_array = batch.column(0)
    assert reconstructed_array.type == uuid_type
    assert reconstructed_array == array

    assert uuid_type.__arrow_ext_scalar_class__() == pa.UuidScalar
    assert isinstance(array[0], pa.UuidScalar)


def test_tensor_type():
    tensor_type = pa.fixed_shape_tensor(pa.int8(), [2, 3])
    assert tensor_type.extension_name == "arrow.fixed_shape_tensor"
    assert tensor_type.storage_type == pa.list_(pa.int8(), 6)
    assert tensor_type.shape == [2, 3]
    assert tensor_type.dim_names is None
    assert tensor_type.permutation is None

    tensor_type = pa.fixed_shape_tensor(pa.float64(), [2, 2, 3],
                                        permutation=[0, 2, 1])
    assert tensor_type.extension_name == "arrow.fixed_shape_tensor"
    assert tensor_type.storage_type == pa.list_(pa.float64(), 12)
    assert tensor_type.shape == [2, 2, 3]
    assert tensor_type.dim_names is None
    assert tensor_type.permutation == [0, 2, 1]

    tensor_type = pa.fixed_shape_tensor(pa.bool_(), [2, 2, 3],
                                        dim_names=['C', 'H', 'W'])
    assert tensor_type.extension_name == "arrow.fixed_shape_tensor"
    assert tensor_type.storage_type == pa.list_(pa.bool_(), 12)
    assert tensor_type.shape == [2, 2, 3]
    assert tensor_type.dim_names == ['C', 'H', 'W']
    assert tensor_type.permutation is None


@pytest.mark.numpy
@pytest.mark.parametrize("np_type_str", ("int8", "int64", "float32"))
def test_tensor_class_methods(np_type_str):
    from numpy.lib.stride_tricks import as_strided
    arrow_type = pa.from_numpy_dtype(np.dtype(np_type_str))

    tensor_type = pa.fixed_shape_tensor(arrow_type, [2, 3])
    storage = pa.array([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]],
                       pa.list_(arrow_type, 6))
    arr = pa.ExtensionArray.from_storage(tensor_type, storage)
    expected = np.array(
        [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]],
        dtype=np.dtype(np_type_str)
    )
    np.testing.assert_array_equal(arr.to_tensor(), expected)
    np.testing.assert_array_equal(arr.to_numpy_ndarray(), expected)

    expected = np.array([[[7, 8, 9], [10, 11, 12]]], dtype=np.dtype(np_type_str))
    result = arr[1:].to_numpy_ndarray()
    np.testing.assert_array_equal(result, expected)

    values = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
    flat_arr = np.array(values[0], dtype=np.dtype(np_type_str))
    bw = np.dtype(np_type_str).itemsize
    storage = pa.array(values, pa.list_(arrow_type, 12))

    tensor_type = pa.fixed_shape_tensor(arrow_type, [2, 2, 3], permutation=[0, 1, 2])
    result = pa.ExtensionArray.from_storage(tensor_type, storage)
    expected = np.array(
        [[[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]]],
        dtype=np.dtype(np_type_str)
    )
    np.testing.assert_array_equal(result.to_numpy_ndarray(), expected)

    result = flat_arr.reshape(1, 2, 3, 2)
    expected = np.array(
        [[[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]],
        dtype=np.dtype(np_type_str)
    )
    np.testing.assert_array_equal(result, expected)

    tensor_type = pa.fixed_shape_tensor(arrow_type, [2, 2, 3], permutation=[0, 2, 1])
    result = pa.ExtensionArray.from_storage(tensor_type, storage)
    expected = as_strided(flat_arr, shape=(1, 2, 3, 2),
                          strides=(bw * 12, bw * 6, bw, bw * 3))
    np.testing.assert_array_equal(result.to_numpy_ndarray(), expected)

    tensor_type = pa.fixed_shape_tensor(arrow_type, [2, 2, 3], permutation=[2, 0, 1])
    result = pa.ExtensionArray.from_storage(tensor_type, storage)
    expected = as_strided(flat_arr, shape=(1, 3, 2, 2),
                          strides=(bw * 12, bw, bw * 6, bw * 2))
    np.testing.assert_array_equal(result.to_numpy_ndarray(), expected)

    assert result.type.permutation == [2, 0, 1]
    assert result.type.shape == [2, 2, 3]
    assert result.to_tensor().shape == (1, 3, 2, 2)
    assert result.to_tensor().strides == (12 * bw, 1 * bw, 6 * bw, 2 * bw)


@pytest.mark.numpy
@pytest.mark.parametrize("np_type_str", ("int8", "int64", "float32"))
def test_tensor_array_from_numpy(np_type_str):
    from numpy.lib.stride_tricks import as_strided
    arrow_type = pa.from_numpy_dtype(np.dtype(np_type_str))

    arr = np.array([[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]],
                   dtype=np.dtype(np_type_str), order="C")
    tensor_array_from_numpy = pa.FixedShapeTensorArray.from_numpy_ndarray(arr)
    assert isinstance(tensor_array_from_numpy.type, pa.FixedShapeTensorType)
    assert tensor_array_from_numpy.type.value_type == arrow_type
    assert tensor_array_from_numpy.type.shape == [2, 3]
    assert tensor_array_from_numpy.type.dim_names is None

    arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]],
                   dtype=np.dtype(np_type_str), order="F")
    with pytest.raises(ValueError, match="First stride needs to be largest"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr)

    flat_arr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                        dtype=np.dtype(np_type_str))
    bw = np.dtype(np_type_str).itemsize

    arr = flat_arr.reshape(1, 3, 4)
    tensor_array_from_numpy = pa.FixedShapeTensorArray.from_numpy_ndarray(arr)
    assert tensor_array_from_numpy.type.shape == [3, 4]
    assert tensor_array_from_numpy.type.permutation == [0, 1]
    assert tensor_array_from_numpy.type.dim_names is None
    assert tensor_array_from_numpy.to_tensor() == pa.Tensor.from_numpy(arr)

    arr = as_strided(flat_arr, shape=(1, 2, 3, 2),
                     strides=(bw * 12, bw * 6, bw, bw * 3))
    tensor_array_from_numpy = pa.FixedShapeTensorArray.from_numpy_ndarray(arr)
    assert tensor_array_from_numpy.type.shape == [2, 2, 3]
    assert tensor_array_from_numpy.type.permutation == [0, 2, 1]
    assert tensor_array_from_numpy.type.dim_names is None
    assert tensor_array_from_numpy.to_tensor() == pa.Tensor.from_numpy(arr)

    arr = flat_arr.reshape(1, 2, 3, 2)
    result = pa.FixedShapeTensorArray.from_numpy_ndarray(arr)
    expected = np.array(
        [[[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]],
        dtype=np.dtype(np_type_str)
    )
    np.testing.assert_array_equal(result.to_numpy_ndarray(), expected)

    arr = np.array([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]],
                   dtype=np.dtype(np_type_str))
    expected = arr[1:]
    result = pa.FixedShapeTensorArray.from_numpy_ndarray(arr)[1:].to_numpy_ndarray()
    np.testing.assert_array_equal(result, expected)

    arr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=np.dtype(np_type_str))
    with pytest.raises(ValueError, match="Cannot convert 1D array or scalar to fixed"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr)

    arr = np.array(1, dtype=np.dtype(np_type_str))
    with pytest.raises(ValueError, match="Cannot convert 1D array or scalar to fixed"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr)

    arr = np.array([], dtype=np.dtype(np_type_str))

    with pytest.raises(ValueError, match="Cannot convert 1D array or scalar to fixed"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr.reshape((0)))

    with pytest.raises(ValueError, match="Expected a non-empty ndarray"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr.reshape((0, 3, 2)))

    with pytest.raises(ValueError, match="Expected a non-empty ndarray"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr.reshape((3, 0, 2)))

    arr = np.array([[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]],
                   dtype=np.dtype(np_type_str), order="C")
    dim_names = ["a", "b"]
    tensor_array_from_numpy = pa.FixedShapeTensorArray.from_numpy_ndarray(
        arr, dim_names=dim_names)
    assert tensor_array_from_numpy.type.value_type == arrow_type
    assert tensor_array_from_numpy.type.shape == [2, 3]
    assert tensor_array_from_numpy.type.dim_names == dim_names

    with pytest.raises(ValueError, match="The length of dim_names"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr, dim_names=['only_one'])

    with pytest.raises(TypeError, match="dim_names must be a tuple or list"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr, dim_names=123)

    with pytest.raises(TypeError, match="dim_names must be a tuple or list"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(
            arr, dim_names=(x for x in range(2)))

    with pytest.raises(TypeError, match="Each element of dim_names must be a string"):
        pa.FixedShapeTensorArray.from_numpy_ndarray(arr, dim_names=[0, 1])


@pytest.mark.numpy
@pytest.mark.parametrize("tensor_type", (
    pa.fixed_shape_tensor(pa.int8(), [2, 2, 3]),
    pa.fixed_shape_tensor(pa.int8(), [2, 2, 3], permutation=[0, 2, 1]),
    pa.fixed_shape_tensor(pa.int8(), [2, 2, 3], dim_names=['C', 'H', 'W'])
))
def test_tensor_type_ipc(tensor_type):
    storage = pa.array([[1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6]], pa.list_(pa.int8(), 12))
    arr = pa.ExtensionArray.from_storage(tensor_type, storage)
    batch = pa.RecordBatch.from_arrays([arr], ["ext"])

    # check the built array has exactly the expected clss
    tensor_class = tensor_type.__arrow_ext_class__()
    assert isinstance(arr, tensor_class)

    buf = ipc_write_batch(batch)
    del batch
    batch = ipc_read_batch(buf)

    result = batch.column(0)
    # check the deserialized array class is the expected one
    assert isinstance(result, tensor_class)
    assert result.type.extension_name == "arrow.fixed_shape_tensor"
    assert arr.storage.to_pylist() == [[1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6]]

    # we get back an actual TensorType
    assert isinstance(result.type, pa.FixedShapeTensorType)
    assert result.type.value_type == pa.int8()
    assert result.type.shape == [2, 2, 3]


def test_tensor_type_equality():
    tensor_type = pa.fixed_shape_tensor(pa.int8(), [2, 2, 3])
    assert tensor_type.extension_name == "arrow.fixed_shape_tensor"

    tensor_type2 = pa.fixed_shape_tensor(pa.int8(), [2, 2, 3])
    tensor_type3 = pa.fixed_shape_tensor(pa.uint8(), [2, 2, 3])
    assert tensor_type == tensor_type2
    assert not tensor_type == tensor_type3


def test_tensor_type_cast():
    tensor_type = pa.fixed_shape_tensor(pa.int8(), [2, 3])
    inner = pa.array(range(18), pa.int8())
    storage = pa.FixedSizeListArray.from_arrays(inner, 6)

    # cast storage -> extension type
    result = storage.cast(tensor_type)
    expected = pa.ExtensionArray.from_storage(tensor_type, storage)
    assert result.equals(expected)

    # cast extension type -> storage type
    storage_result = result.cast(storage.type)
    assert storage_result.equals(storage)


@pytest.mark.pandas
def test_extension_to_pandas_storage_type(registered_period_type):
    period_type, _ = registered_period_type
    np_arr = np.array([1, 2, 3, 4], dtype='i8')
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(period_type, storage)

    if isinstance(period_type, PeriodTypeWithToPandasDtype):
        pandas_dtype = period_type.to_pandas_dtype()
    else:
        pandas_dtype = np_arr.dtype

    # Test arrays
    result = arr.to_pandas()
    assert result.dtype == pandas_dtype

    # Test chunked arrays
    chunked_arr = pa.chunked_array([arr])
    result = chunked_arr.to_numpy()
    assert result.dtype == np_arr.dtype

    result = chunked_arr.to_pandas()
    assert result.dtype == pandas_dtype

    # Test Table.to_pandas
    data = [
        pa.array([1, 2, 3, 4]),
        pa.array(['foo', 'bar', None, None]),
        pa.array([True, None, True, False]),
        arr
    ]
    my_schema = pa.schema([('f0', pa.int8()),
                           ('f1', pa.string()),
                           ('f2', pa.bool_()),
                           ('ext', period_type)])
    table = pa.Table.from_arrays(data, schema=my_schema)
    result = table.to_pandas()
    assert result["ext"].dtype == pandas_dtype

    import pandas as pd
    # Skip tests for 2.0.x, See: GH-35821
    if (
        Version(pd.__version__) >= Version("2.1.0")
    ):
        # Check the usage of types_mapper
        result = table.to_pandas(types_mapper=pd.ArrowDtype)
        assert isinstance(result["ext"].dtype, pd.ArrowDtype)


def test_tensor_type_is_picklable(pickle_module):
    # GH-35599

    expected_type = pa.fixed_shape_tensor(pa.int32(), (2, 2))
    result = pickle_module.loads(pickle_module.dumps(expected_type))

    assert result == expected_type

    arr = [[1, 2, 3, 4], [10, 20, 30, 40], [100, 200, 300, 400]]
    storage = pa.array(arr, pa.list_(pa.int32(), 4))
    expected_arr = pa.ExtensionArray.from_storage(expected_type, storage)
    result = pickle_module.loads(pickle_module.dumps(expected_arr))

    assert result == expected_arr


@pytest.mark.parametrize(("tensor_type", "text"), [
    (
        pa.fixed_shape_tensor(pa.int8(), [2, 2, 3]),
        'fixed_shape_tensor[value_type=int8, shape=[2,2,3]]'
    ),
    (
        pa.fixed_shape_tensor(pa.int32(), [2, 2, 3], permutation=[0, 2, 1]),
        'fixed_shape_tensor[value_type=int32, shape=[2,2,3], permutation=[0,2,1]]'
    ),
    (
        pa.fixed_shape_tensor(pa.int64(), [2, 2, 3], dim_names=['C', 'H', 'W']),
        'fixed_shape_tensor[value_type=int64, shape=[2,2,3], dim_names=[C,H,W]]'
    )
])
def test_tensor_type_str(tensor_type, text):
    tensor_type_str = tensor_type.__str__()
    assert text in tensor_type_str


@pytest.mark.parametrize("storage_type,storage", [
    (pa.null(), [None] * 4),
    (pa.int64(), [1, 2, None, 4]),
    (pa.binary(), [None, b"foobar"]),
    (pa.list_(pa.int64()), [[], [1, 2], None, [3, None]]),
])
def test_opaque_type(pickle_module, storage_type, storage):
    opaque_type = pa.opaque(storage_type, "type", "vendor")
    assert opaque_type.extension_name == "arrow.opaque"
    assert opaque_type.storage_type == storage_type
    assert opaque_type.type_name == "type"
    assert opaque_type.vendor_name == "vendor"
    assert "arrow.opaque" in str(opaque_type)

    assert opaque_type == opaque_type
    assert opaque_type != storage_type
    assert opaque_type != pa.opaque(storage_type, "type2", "vendor")
    assert opaque_type != pa.opaque(storage_type, "type", "vendor2")
    assert opaque_type != pa.opaque(pa.decimal128(12, 3), "type", "vendor")

    # Pickle roundtrip
    result = pickle_module.loads(pickle_module.dumps(opaque_type))
    assert result == opaque_type

    # IPC roundtrip
    opaque_arr_class = opaque_type.__arrow_ext_class__()
    storage = pa.array(storage, storage_type)
    arr = pa.ExtensionArray.from_storage(opaque_type, storage)
    assert isinstance(arr, opaque_arr_class)

    buf = ipc_write_batch(pa.RecordBatch.from_arrays([arr], ["ext"]))
    batch = ipc_read_batch(buf)

    assert batch.column(0).type.extension_name == "arrow.opaque"
    assert isinstance(batch.column(0), opaque_arr_class)

    # cast storage -> extension type
    result = storage.cast(opaque_type)
    assert result == arr

    # cast extension type -> storage type
    inner = arr.cast(storage_type)
    assert inner == storage


def test_bool8_type(pickle_module):
    bool8_type = pa.bool8()
    storage_type = pa.int8()
    assert bool8_type.extension_name == "arrow.bool8"
    assert bool8_type.storage_type == storage_type
    assert str(bool8_type) == "extension<arrow.bool8>"

    assert bool8_type == bool8_type
    assert bool8_type == pa.bool8()
    assert bool8_type != storage_type

    # Pickle roundtrip
    result = pickle_module.loads(pickle_module.dumps(bool8_type))
    assert result == bool8_type

    # IPC roundtrip
    storage = pa.array([-1, 0, 1, 2, None], storage_type)
    arr = pa.ExtensionArray.from_storage(bool8_type, storage)
    assert isinstance(arr, pa.Bool8Array)

    # extension is registered by default
    buf = ipc_write_batch(pa.RecordBatch.from_arrays([arr], ["ext"]))
    batch = ipc_read_batch(buf)

    assert batch.column(0).type.extension_name == "arrow.bool8"
    assert isinstance(batch.column(0), pa.Bool8Array)

    # cast storage -> extension type
    result = storage.cast(bool8_type)
    assert result == arr

    # cast extension type -> storage type
    inner = arr.cast(storage_type)
    assert inner == storage


def test_bool8_to_bool_conversion():
    bool_arr = pa.array([True, False, True, True, None], pa.bool_())
    bool8_arr = pa.ExtensionArray.from_storage(
        pa.bool8(),
        pa.array([-1, 0, 1, 2, None], pa.int8()),
    )

    # cast extension type -> arrow boolean type
    assert bool8_arr.cast(pa.bool_()) == bool_arr

    # cast arrow boolean type -> extension type, expecting canonical values
    canonical_storage = pa.array([1, 0, 1, 1, None], pa.int8())
    canonical_bool8_arr = pa.ExtensionArray.from_storage(pa.bool8(), canonical_storage)
    assert bool_arr.cast(pa.bool8()) == canonical_bool8_arr


@pytest.mark.numpy
def test_bool8_to_numpy_conversion():
    arr = pa.ExtensionArray.from_storage(
        pa.bool8(),
        pa.array([-1, 0, 1, 2, None], pa.int8()),
    )

    # cannot zero-copy with nulls
    with pytest.raises(
        pa.ArrowInvalid,
        match="Needed to copy 1 chunks with 1 nulls, but zero_copy_only was True",
    ):
        arr.to_numpy()

    # nullable conversion possible with a copy, but dest dtype is object
    assert np.array_equal(
        arr.to_numpy(zero_copy_only=False),
        np.array([True, False, True, True, None], dtype=np.object_),
    )

    # zero-copy possible with non-null array
    np_arr_no_nulls = np.array([True, False, True, True], dtype=np.bool_)
    arr_no_nulls = pa.ExtensionArray.from_storage(
        pa.bool8(),
        pa.array([-1, 0, 1, 2], pa.int8()),
    )

    arr_to_np = arr_no_nulls.to_numpy()
    assert np.array_equal(arr_to_np, np_arr_no_nulls)

    # same underlying buffer
    assert arr_to_np.ctypes.data == arr_no_nulls.buffers()[1].address

    # if the user requests a writable array, a copy should be performed
    arr_to_np_writable = arr_no_nulls.to_numpy(zero_copy_only=False, writable=True)
    assert np.array_equal(arr_to_np_writable, np_arr_no_nulls)

    # different underlying buffer
    assert arr_to_np_writable.ctypes.data != arr_no_nulls.buffers()[1].address


@pytest.mark.numpy
def test_bool8_from_numpy_conversion():
    np_arr_no_nulls = np.array([True, False, True, True], dtype=np.bool_)
    canonical_bool8_arr_no_nulls = pa.ExtensionArray.from_storage(
        pa.bool8(),
        pa.array([1, 0, 1, 1], pa.int8()),
    )

    arr_from_np = pa.Bool8Array.from_numpy(np_arr_no_nulls)
    assert arr_from_np == canonical_bool8_arr_no_nulls

    # same underlying buffer
    assert arr_from_np.buffers()[1].address == np_arr_no_nulls.ctypes.data

    # conversion only valid for 1-D arrays
    with pytest.raises(
        ValueError,
        match="Cannot convert 2-D array to bool8 array",
    ):
        pa.Bool8Array.from_numpy(
            np.array([[True, False], [False, True]], dtype=np.bool_),
        )

    with pytest.raises(
        ValueError,
        match="Cannot convert 0-D array to bool8 array",
    ):
        pa.Bool8Array.from_numpy(np.bool_())

    # must use compatible storage type
    with pytest.raises(
        TypeError,
        match="Array dtype float64 incompatible with bool8 storage",
    ):
        pa.Bool8Array.from_numpy(np.array([1, 2, 3], dtype=np.float64))


def test_bool8_scalar():
    assert pa.ExtensionScalar.from_storage(pa.bool8(), -1).as_py() is True
    assert pa.ExtensionScalar.from_storage(pa.bool8(), 0).as_py() is False
    assert pa.ExtensionScalar.from_storage(pa.bool8(), 1).as_py() is True
    assert pa.ExtensionScalar.from_storage(pa.bool8(), 2).as_py() is True
    assert pa.ExtensionScalar.from_storage(pa.bool8(), None).as_py() is None

    arr = pa.ExtensionArray.from_storage(
        pa.bool8(),
        pa.array([-1, 0, 1, 2, None], pa.int8()),
    )
    assert arr[0].as_py() is True
    assert arr[1].as_py() is False
    assert arr[2].as_py() is True
    assert arr[3].as_py() is True
    assert arr[4].as_py() is None

    assert pa.scalar(-1, type=pa.bool8()).as_py() is True
    assert pa.scalar(0, type=pa.bool8()).as_py() is False
    assert pa.scalar(1, type=pa.bool8()).as_py() is True
    assert pa.scalar(2, type=pa.bool8()).as_py() is True
    assert pa.scalar(None, type=pa.bool8()).as_py() is None


@pytest.mark.parametrize("storage_type", (
    pa.string(), pa.large_string(), pa.string_view()))
def test_json(storage_type, pickle_module):
    data = ['{"a": 1}', '{"b": 2}', None]
    json_type = pa.json_(storage_type)
    storage = pa.array(data, type=storage_type)
    array = pa.array(data, type=json_type)
    json_arr_class = json_type.__arrow_ext_class__()

    assert pa.json_() == pa.json_(pa.utf8())
    assert json_type.extension_name == "arrow.json"
    assert json_type.storage_type == storage_type
    assert json_type.__class__ is pa.JsonType

    assert json_type == pa.json_(storage_type)
    assert json_type != storage_type

    assert isinstance(array, pa.JsonArray)

    assert array.to_pylist() == data
    assert array[0].as_py() == data[0]
    assert array[2].as_py() is None

    # Pickle roundtrip
    result = pickle_module.loads(pickle_module.dumps(json_type))
    assert result == json_type

    # IPC roundtrip
    buf = ipc_write_batch(pa.RecordBatch.from_arrays([array], ["ext"]))
    batch = ipc_read_batch(buf)
    reconstructed_array = batch.column(0)
    assert reconstructed_array.type == json_type
    assert reconstructed_array == array
    assert isinstance(array, json_arr_class)

    assert json_type.__arrow_ext_scalar_class__() == pa.JsonScalar
    assert isinstance(array[0], pa.JsonScalar)

    # cast storage -> extension type
    result = storage.cast(json_type)
    assert result == array

    # cast extension type -> storage type
    inner = array.cast(storage_type)
    assert inner == storage

    for storage_type in (pa.int32(), pa.large_binary(), pa.float32()):
        with pytest.raises(
                pa.ArrowInvalid,
                match=f"Invalid storage type for JsonExtensionType: {storage_type}"):
            pa.json_(storage_type)
