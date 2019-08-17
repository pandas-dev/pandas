import pandas as pd
import pytest
import unittest

class TestStringMethods(unittest.TestCase):
	
	def test1(self):
		x = pd.Series([1,2,3,4])
		xt = type(x)
		assert pd.isnull(xt)==False,"Passed"

	def test2(self):
		y = pd.DataFrame({"col": [1,2,3,4]})
		yt = type(y)
		assert pd.isnull(yt)==False,"Passed"


if __name__ == '__main__':
    unittest.main()