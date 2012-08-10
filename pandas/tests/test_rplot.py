import unittest
import pandas.tools.rplot as rplot
from pandas import read_csv
import os

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

class TestUtilityFunctions(unittest.TestCase):
	pass

class TestScaleGradient(unittest.TestCase):
	def setUp(self):
		path = os.path.join(curpath(), 'data/iris.csv')
		self.data = read_csv(path, sep=',')
		self.gradient = rplot.ScaleGradient("SepalLength", colour1=(0.2, 0.3, 0.4), colour2=(0.8, 0.7, 0.6))

	def test_gradient(self):
		for index in range(len(self.data)):
			row = self.data.irow(index)
			r, g, b = self.gradient(self.data, index)
			r1, g1, b1 = self.gradient.colour1
			r2, g2, b2 = self.gradient.colour2
			self.assertGreaterEqual(r, r1)
			self.assertGreaterEqual(g, g1)
			self.assertGreaterEqual(b, b1)
			self.assertLessEqual(r, r2)
			self.assertLessEqual(g, g2)
			self.assertLessEqual(b, b2)

class TestScaleGradient2(unittest.TestCase):
	pass

if __name__ == '__main__':
	unittest.main()