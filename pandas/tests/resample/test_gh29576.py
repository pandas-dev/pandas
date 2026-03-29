import numpy as np
import pandas as pd

def test_resample_2qs_mar_offset():
    # GH#29576
    d = pd.Series(data=np.zeros(365), 
                  index=pd.date_range('1950-01-01', '1950-12-31', freq='D'))
    result = d.resample('2QS-MAR').mean()

    expected_index = pd.to_datetime(['1949-09-01', '1950-03-01', '1950-09-01'])
    assert result.index[0] == pd.Timestamp('1949-09-01')
