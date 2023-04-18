import pandas as pd
import numpy as np
import pytest

def test_normgroup():
    # test normalization by l2 norm
    data = {'group': ['A', 'A', 'B', 'B', 'B', 'C', 'C'],
            'value': [1, 2, 3, 4, 5, 6, 7]}
    df = pd.DataFrame(data)
    expected = pd.DataFrame({'value': [0.447214, 0.894427, 0.424264, 0.565685, 0.707107, 0.640513, 0.767188]})
    assert df.normgroup('group', norm='l2').equals(expected)

    # test normalization by max norm
    data = {'group': ['A', 'A', 'A', 'B', 'B', 'C'],
            'value': [1, 3, 5, 2, 4, 6]}
    df = pd.DataFrame(data)
    expected = pd.DataFrame({'value': [1, 1, 1, 0.5, 1, 1]})
    assert df.normgroup('group', norm=np.inf).equals(expected)

    # test empty DataFrame
    df = pd.DataFrame(columns=['group', 'value'])
    expected = pd.DataFrame(columns=['value'])
    assert df.normgroup('group', norm='l2').equals(expected)

    # test invalid column name
    data = {'group': ['A', 'A', 'B', 'B'],
            'value': [1, 2, 3, 4]}
    df = pd.DataFrame(data)
    with pytest.raises(KeyError):
        df.normgroup('invalid_col', norm='l2')
