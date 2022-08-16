"""
algos/ directory is intended for individual functions from core.algorithms

In many cases these algorithms are reachable in multiple ways:
   algos.foo(x, y)
   Series(x).foo(y)
   Index(x).foo(y)
   pd.array(x).foo(y)

In most cases we profile the Series variant directly, trusting the performance
of the others to be highly correlated.
"""
