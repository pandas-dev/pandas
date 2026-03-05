import uuid
import pytest
import pandas as pd



pa = pytest.importorskip("pyarrow")

def test_uuid_keeps_dtype_verify():
    u = uuid.UUID("550e8400-e29b-41d4-a716-446655440000") #we have our uuid in u
    arr = pa.array([u.bytes, None], type=pa.uuid()) #we put it in a function, we include None in the array to simultaneously test the case where there is a "missing value" since they cause bugs a lot so knock two birds out with one stone

    s = pd.Series(arr) #lets put it into a series and hope it does not convert to object


    assert s.dtype != object
    assert s.array._pa_array.type == pa.uuid() #we convert it back to pyarrow  before the comparison to avoid any weird wrapper stuff inside of series. We just want to see if it maintains uuid type

#now we verified series, lets do data frame
    d = pd.DataFrame({"id": arr})
    assert d["id"].array._pa_array.type == pa.uuid()

    
    