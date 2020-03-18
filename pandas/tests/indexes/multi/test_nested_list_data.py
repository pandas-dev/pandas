import pandas as pd
import pandas._testing as tm
import numpy as np
import pytest


'''
This code is meant to test the fix implemented for issue #14467.
https://github.com/pandas-dev/pandas/issues/14467
'''

class TestNestedListDataMultiIndex:

  def test_nested_list(self):
    # Case from issue, creating data from np.array works, and should match this case
    result = pd.DataFrame([[1, 2, 3], [4, 5, 6]],
                      index=[['gibberish']*2, [0, 1]],
                      columns=[['baldersash']*3, [10, 20, 30]])

    expected = pd.DataFrame(np.array([[1, 2, 3], [4, 5, 6]]),
                      index=[['gibberish']*2, [0, 1]],
                      columns=[['baldersash']*3, [10, 20, 30]])
    tm.assert_index_equal(result.columns, expected.columns)

  def test_nest_list_with_multiIndex(self):
    # Creating from a multiIndex should also still work
    m = pd.MultiIndex.from_arrays([['baldersash']*3, [10, 20, 30]])
    result = pd.DataFrame([[1, 2, 3], [4, 5, 6]],
                      index=[['gibberish']*2, [0, 1]],
                      columns=m)

    expected = pd.DataFrame(np.array([[1, 2, 3], [4, 5, 6]]),
                      index=[['gibberish']*2, [0, 1]],
                      columns=[['baldersash']*3, [10, 20, 30]])
    tm.assert_index_equal(result.columns, expected.columns)

  def test_wrong_length_raises_error(self):
    # Make sure the code raises an error if the nested lists have the wrong length
    with pytest.raises(ValueError):
      result = pd.DataFrame([[1, 2, 3], [4, 5, 6]],
                        index=[['gibberish']*2, [0, 1]],
                        columns=[['baldersash']*2, [10, 20]])