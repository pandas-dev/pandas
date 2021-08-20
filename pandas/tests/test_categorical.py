import pandas as pd
import unittest

def series_fill(cat, ele):
  if cat.dtype == 'object':
      raise TypeError("Given DataFrame is not of category datatype. Please change it to category datatype.")
  if ele not in cat.categories:
    raise TypeError("Element not present in categories. Cannot be filled in series.")
  ser = pd.Series(cat).fillna(ele)
  filled = cat.fillna(ser)
  return filled

class Test_series_fill(unittest.TestCase):
  
  def test_series_fill(self):
    self.assertCountEqual(series_fill(pd.Categorical(["A", "B", None]), "B"), ["A", "B", "B"])
    self.assertCountEqual(series_fill(pd.Categorical(["1", "2", "3", None]), "2"), ["1", "2", "3", "2"])

  def test_input_value(self):
    self.assertRaises(TypeError, series_fill, True)

unittest.main(argv=[''],verbosity=2, exit=False)