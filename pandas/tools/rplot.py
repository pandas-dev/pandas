import matplotlib.pyplot as plt

class RPlot:
	def __init__(self, data):
		self.data = data
		self.ax = plt.gca()

	def __add__(self, other):
		other.plot(self.ax, self.data)

class GeomPoint:
	def __init__(self, x, y, marker='o', colour='grey', size=20, alpha=1.0):
		self.x = x
		self.y = y
		self.marker = marker
		self.colour = colour
		self.size = size
		self.alpha = alpha

	def plot(self, ax, data):
		ax.scatter(data[self.x], data[self.y], c=self.colour, marker=self.marker, s=self.size, alpha=self.alpha)