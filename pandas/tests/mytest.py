import pandas as pd
import numpy as np

def test_me():
    df = pd.DataFrame({'key1' : ['a', 'a', 'b', 'b', 'a'],
                    'key2' : ['one', 'two', 'one', 'two', 'one'],
                    'key3' : ['three', 'three', 'three', 'six', 'six'],
                    'data1' : np.random.randn(5),
                    'data2' : np.random.randn(5)})
    df.groupby('key1').min()
