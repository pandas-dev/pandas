# Your code here


import sys
sys.path.append("C:\\Users\\CheanHui\\Desktop\\pandas\\pandas")
import pandas

from scipy.sparse import coo_matrix, eye  
import pandas as pd
from pandas.core.arrays.sparse import SparseArray, SparseDtype
import pandas._testing as tm
import numpy as np


df = pd.DataFrame.sparse.from_spmatrix(eye(5))

#print(df.loc[range(2)])                                                      


#print(df.loc[range(2)].sparse.density)


#print(df.loc[range(2)].loc[range(1)])                                          

#print(df.loc[range(2)].loc[range(1)].sparse.density)

#print(df.loc[range(2)].loc[range(1)].dtypes)


#
#expected = pd.DataFrame(mat.todense()).astype(dtype)
def test_loc():
    # Tests .loc on sparse DataFrame #34687
    df = pd.DataFrame.sparse.from_spmatrix(eye(5))
    res1 = df.loc[range(2)]
    exp1 = pd.DataFrame([[1.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0, 0.0, 0.0]], dtype=SparseDtype("float64", 0.0))
    tm.assert_frame_equal(res1, exp1)
    
    res2 = df.loc[range(2)].loc[range(1)]
    exp2 = pd.DataFrame([[1.0, 0.0, 0.0, 0.0, 0.0]], dtype=SparseDtype("float64", 0.0))
    tm.assert_frame_equal(res2, exp2)
    #m.assert_extension_array_equal(exp2, exp2)
    #print(res2.mask)
    #print(exp2.mask)

test_loc()

"""
import pandas as pd


def test_sort_index_mode_true():

    d = {'col1': [1, 2, 3], 'col2': [3, 4, 5]}
    df = pd.DataFrame(data=d)
    now = pd.Timestamp.now()
    df.index = [now, now - pd.Timedelta('1m'), now + pd.Timedelta('2m')]
    
    try:
        df1 = df.sort_index()
    except Exception as ex:
        print('test 1 failed: ', ex)

    pd.set_option('use_inf_as_na', True)

    try:
        df2 = df.sort_index()
    except Exception as ex:
        print('test 2 failed: ', ex)

test_sort_index_mode_true()
"""


import pandas as pd

df = pd.DataFrame()
df['major'] = pd.Series(list('AAAAB'))
df['minor'] = pd.Series(
    list('YXYXX'), dtype=pd.CategoricalDtype(['Y', 'X']))
df['value'] = pd.Series([1, 2, 3, 4, 5])
df.set_index(['major', 'minor'], inplace=True)

df1 = df.iloc[:2]
df2 = df.iloc[2:]

df3 = df1.join(df2, lsuffix='_left', rsuffix='_right')
#print(df)
#print("-----")
#print(df1)
#print("-----")
#print(df2)

