.. _indexobj.rst


*****
Index
*****

.. currentmodule:: pandas

The pandas data model is relatively simple: link data to a set of
unique labels. This labelling is implemented through the
:class:`Index` class:

.. class:: Index

   a 1-D NumPy :class:`~numpy.ndarray` subclass which represents an
   *ordered set*.

   :Parameters:
       **labels** : {array_like}

           Any data that is valid for constructing a 1-D
           :class:`~numpy.ndarray` can be used here. The order of the
           labels is preserved and they need not be sorted.

.. note::

    An Index instance will always be of **object** dtype. The reason
    for this is to avoid occasional type-casting snafus when holding
    strings in NumPy arrays.

Usage and behavior
------------------

In general, users will seldom need to create Index instances
themselves. However, understanding their internal structure will
be important for developers. Creating an Index is simple:

::

    >>> index = Index(['a', 'b', 'c', 'd', 'e'])
    >>> index
    Index([a, b, c, d, e], dtype=object)

The Index stores the labels in two ways: one as a vector, and the
other in a dict:

::

    >>> index.indexMap
    {'a': 0, 'b': 1, 'c': 2, 'd': 3, 'e': 4}

This dict allows the pandas data structures to perform fast lookups
and determine membership:

::

    >>> 'a' in index
    True
    >>> 'f' in index
    False

Slicing produces a new index with the **indexMap** field adjusted
appropriately.

::

    >>> index_slice = index[2:]
    >>> index_slice.indexMap
    {'c': 0, 'd': 1, 'e': 2}

To prevent undesired behavior, Index instances are immutable:

::

    In [11]: index[2] = 'd'
	---------------------------------------------------------------------------
	Exception                                 Traceback (most recent call last)

	pandas/core/index.pyc in __setitem__(self, key, value)
	     98     def __setitem__(self, key, value):
	     99         """Disable the setting of values."""
	--> 100         raise Exception(str(self.__class__) + ' object is immutable' )
	    101
	    102     def __getitem__(self, key):

	Exception: <class 'pandas.core.index.Index'> object is immutable

Convenience methods
-------------------

A number of methods are provided for performing common set-like
operations and comparisons:

::

    >>> index1 = Index(['a', 'b', 'c'])
    >>> index2 = Index(['c', 'd', 'e'])

    >>> index1.intersection(index2)
    Index([c], dtype=object)

    >>> index1.union(index2)
    Index([a, b, c, d, e], dtype=object)

    >>> index1.equals(index2)
    False

.. autosummary::
   :toctree: generated/

   Index.equals
   Index.intersection
   Index.union
