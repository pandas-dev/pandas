import unittest
import pandas.tools.rplot as rplot
from pandas import read_csv
import os

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

def between(a, b, x):
	"""Check if x is in the somewhere between a and b.

	Parameters:
	-----------
	a: float, interval start
	b: float, interval end
	x: float, value to test for

	Returns:
	--------
	True if x is between a and b, False otherwise
	"""
	if a < b:
		return x >= a and x <= b
	else:
		return x <= a and x >= b

class TestUtilityFunctions(unittest.TestCase):
	"""
	Tests for RPlot utility functions.
	"""
	def test_make_aes1(self):
		aes = rplot.make_aes()
		self.assertTrue(aes['x'] is None)
		self.assertTrue(aes['y'] is None)
		self.assertTrue(aes['size'] is None)
		self.assertTrue(aes['colour'] is None)
		self.assertTrue(aes['shape'] is None)
		self.assertTrue(aes['alpha'] is None)
		self.assertTrue(type(aes) is dict)

	def test_make_aes2(self):
		with self.assertRaises(ValueError):
			rplot.make_aes(size=rplot.ScaleShape('test'))
		with self.assertRaises(ValueError):
			rplot.make_aes(colour=rplot.ScaleShape('test'))
		with self.assertRaises(ValueError):
			rplot.make_aes(shape=rplot.ScaleSize('test'))
		with self.assertRaises(ValueError):
			rplot.make_aes(alpha=rplot.ScaleShape('test'))

	def test_dictionary_union(self):
		dict1 = {1 : 1, 2 : 2, 3 : 3}
		dict2 = {1 : 1, 2 : 2, 4 : 4}
		union = rplot.dictionary_union(dict1, dict2)
		self.assertEqual(len(union), 4)
		keys = union.keys()
		self.assertTrue(1 in keys)
		self.assertTrue(2 in keys)
		self.assertTrue(3 in keys)
		self.assertTrue(4 in keys)
		self.assertTrue(rplot.dictionary_union(dict1, {}) == dict1)
		self.assertTrue(rplot.dictionary_union({}, dict1) == dict1)
		self.assertTrue(rplot.dictionary_union({}, {}) == {})

	def test_merge_aes(self):
		layer1 = rplot.Layer(size=rplot.ScaleSize('test'))
		layer2 = rplot.Layer(shape=rplot.ScaleShape('test'))
		rplot.merge_aes(layer1, layer2)
		self.assertTrue(isinstance(layer2.aes['size'], rplot.ScaleSize))
		self.assertTrue(isinstance(layer2.aes['shape'], rplot.ScaleShape))
		self.assertTrue(layer2.aes['size'] == layer1.aes['size'])
		for key in layer2.aes.keys():
			if key != 'size' and key != 'shape':
				self.assertTrue(layer2.aes[key] is None)

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
			self.assertTrue(between(r1, r2, r))
			self.assertTrue(between(g1, g2, g))
			self.assertTrue(between(b1, b2, b))

class TestScaleGradient2(unittest.TestCase):
	def setUp(self):
		path = os.path.join(curpath(), 'data/iris.csv')
		self.data = read_csv(path, sep=',')
		self.gradient = rplot.ScaleGradient2("SepalLength", colour1=(0.2, 0.3, 0.4), colour2=(0.8, 0.7, 0.6), colour3=(0.5, 0.5, 0.5))

	def test_gradient2(self):
		for index in range(len(self.data)):
			row = self.data.irow(index)
			r, g, b = self.gradient(self.data, index)
			r1, g1, b1 = self.gradient.colour1
			r2, g2, b2 = self.gradient.colour2
			r3, g3, b3 = self.gradient.colour3
			value = row[self.gradient.column]
			a_ = min(self.data[self.gradient.column])
			b_ = max(self.data[self.gradient.column])
			scaled = (value - a_) / (b_ - a_)
			if scaled < 0.5:
				self.assertTrue(between(r1, r2, r))
				self.assertTrue(between(g1, g2, g))
				self.assertTrue(between(b1, b2, b))
			else:
				self.assertTrue(between(r2, r3, r))
				self.assertTrue(between(g2, g3, g))
				self.assertTrue(between(b2, b3, b))

if __name__ == '__main__':
	unittest.main()