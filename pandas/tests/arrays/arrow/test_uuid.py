import uuid
import pytest
import pandas as pd



pa = pytest.importorskip("pyarrow")

#TEST1: first test will be basic, to see if construciton works
def test_uuid_keeps_dtype_verify():
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000") #we have our uuid in u
    arr = pa.array([u.bytes, None], type=pa.uuid()) #we put it in a function, we include None in the array to simultaneously test the case where there is a "missing value" since they cause bugs a lot so knock two birds out with one stone

    s = pd.Series(arr) #lets put it into a series and hope it does not convert to object


    assert s.dtype != object
    assert s.array._pa_array.type == pa.uuid() #we peak inside before the comparison to avoid any weird wrapper stuff inside of series. We just want to see if it maintains uuid type

#now we verified series, lets do data frame
    d = pd.DataFrame({"id": arr})
    assert d["id"].array._pa_array.type == pa.uuid()


#TEST2: lets see if it can read missing UUID values in an array
def test_uuid_isna():
    #standard, create uuid, put in array, put in series
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    arr = pa.array([u.bytes, None], type=pa.uuid())
    s = pd.Series(arr)

    result = s.isna().tolist()
    assert result == [False, True] #should return flase and true for missing value

#TEST3: for next test, maintainers talked about comparisons, and its a basic functinoality so verifying it with a test
def test_uuid_comparison_eq():
    #standard create uuid, array, then series
    u1 = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    u2 = uuid.UUID("550e8400-e29b-41d4-a716-446655440001")
    arr1 = pa.array([u1.bytes, u2.bytes, None], type=pa.uuid())
    arr2 = pa.array([u1.bytes, u1.bytes, None], type=pa.uuid())
    s1 = pd.Series(arr1)
    s2 = pd.Series(arr2)

    out = (s1 == s2).tolist()#compare
    assert out == [True, False, True]#verify results


#TEST4: we already test thge entrie array being UUID, now individual elements lets make sure bytes and does iloc[i] work
def test_uuid_getitem_scalar():
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    arr = pa.array([u.bytes, None], type=pa.uuid())
    s = pd.Series(arr)

    #extract individual elements
    v0 = s.iloc[0]
    v1 = s.iloc[1]

    assert isinstance(v0, (bytes, bytearray)) #maybe later me make uuid.UUID?
    assert v1 is None #still works?

#TEST5: Maintainer explicitly talked about this, we dont want UUID support to break pandas "semantics" so  "x in s" should check index, not value, and "x in s.array" should check balues
def test_uuid_contains_behavior():
    u = uuid.uuid4()
    arr = pa.array([u.bytes], type=pa.uuid())
    s = pd.Series(arr)

    #checks index and not value
    assert (s.iloc[0] in s) is False

    #checks values
    assert (s.iloc[0] in s.array) is True
    assert (None in s.array) is False

#TEST6: pass as chunked array to make sure nothing breaks with this since using chunked array within pandas is common in piplines
def test_series_from_pyarrow_uuid_chunkedarray():
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000")
    chunk1 = pa.array([u.bytes], type=pa.uuid())
    chunk2 = pa.array([None], type=pa.uuid())
    carr = pa.chunked_array([chunk1, chunk2])
    s = pd.Series(carr)

    #everythig should work
    assert s.dtype != object
    assert s.array._pa_array.type == pa.uuid()
    assert s.isna().tolist() == [False, True]

    
    