.. _removals:

================================
 Code to remove and other ideas
================================

Dropping Python 2 support
=========================

With Python 2.7 reaching its supported end-of-life in 2020, like some other
Python projects (e.g. IPython / Jupyter) we should seriously contemplate making
pandas 2.0 only support Python 3.5 and higher. In addition to lowering the
development burden at both the C API and pure Python level, we can also finally
look to take advantage of features (things like ``asyncio``, maybe?) only
available in Python 3.

Deprecated code to remove
=========================

* ``.ix`` indexing entirely
* ``Panel`` and ``PanelND`` classes
* Plotting?

Other ideas
===========

Here's a collection of other miscellaneous ideas that don't necessarily fit
elsewhere in these documents.

Column statistics
~~~~~~~~~~~~~~~~~

In quite a few pandas algorithms, there are characteristics of the data that
are very useful to know, such as:

* **Monotonicity**: for comparable data (e.g. numbers), is the data sorted /
  strictly increasing? In time series, this permits sorting steps to be
  skipped.

* **Null count**: for data not containing any nulls, the null handling path in
  some algorithms can be skipped entirely



Strided arrays: more trouble than they are worth?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per the general discussion around changing DataFrame's internals to contain a
list / ``std::vector`` of arrays internally, for me this begs the question of
the benefits of continuing to accommodate strided one-dimensional data.

Some pros for eliminating strided data completely:

* Guaranteeing contiguous memory internally will yield more consistent and
  predictable performance.

* Not needing to consider a stride different from 1 means simpler low-level
  array indexing code (e.g. you can work with plain C arrays). The stride is a
  complexity / overhead that leaks to every algorithm that iterates over an
  array.

* You avoid strange situations where a strided view holds onto a base ndarray
  reference to a much larger array

* **Example:** `<https://github.com/wesm/feather/issues/97>`_. Here, the
  internal orientation (column-major vs. row-major) is not clear to the user.

Some cons:

* It would not be possible to perform zero-copy computations on a strided NumPy
  array

* Relatedly, initializing a Series or DataFrame from strided memory would
  require allocating an equivalent amount of contiguous memory for each of the
  columns.

For me, at least, I don't find the cons compelling enough to warrant the code
complexity tradeoff.
