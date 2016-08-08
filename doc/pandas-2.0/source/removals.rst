.. _removals:

================================
 Code to remove and other ideas
================================

Dropping Python 2 support
=========================

Controversional ideas
=====================

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
