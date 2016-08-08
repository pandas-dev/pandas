.. _goals:

=======================
 Goals and Motivations
=======================

.. note::

  These documents are largely written by Wes McKinney, and at this point
  reflect his opinions for the time being

The pandas codebase is now over 8 years old, having grown to over 200,000 lines
of code from its original ~10,000 LOC in the original 0.1 open source release
in January 2010.

At a high level, the "pandas 2.0" effort is based on a couple of observations:

* The pandas 0.x series of releases have been primarily iterative improvements
  to the library, with new features, bug fixes, and improved
  documentation. There have also been a series of deprecations, API changes,
  and other evolutions of pandas's API to account for suboptimal design choices
  (for example: the ``.ix`` operator) made in the early days of the project
  (2010 to 2012).
* pandas's ability to support an increasingly broad set of use cases has been
  significantly constrained (as will be examined in detail in these documents)
  by its tight coupling to NumPy and therefore subject to design limitations in
  NumPy.
* Making significant functional additions (particularly filling gaps in NumPy)
  to pandas, particularly new data types, has grown increasingly complex with
  very obvious accumulations of technical debt.
* pandas is being used increasingly for very large datasets on machines with
  many cores and large amounts of RAM (100s of gigabytes to terabytes)
* Rough edges in pandas's implementation (e.g. its handling of missing data
  across data types) are being exposed to users.

These documents are largely concerned with pandas's internal design, which is
mostly invisible to average users. Advanced users of pandas are generally
familiar with some of these internal details, particular around performance and
memory use, and so the degree to which users are impacted will vary quite a
lot.

Key areas of work
=================

Possible changes or improvements to pandas's internals fall into a number of
different buckets to be explored in great detail:

* **Decoupling from NumPy while preserving interoperability**: by eliminating
  the presumption that pandas objects internally must contain data stored in
  NumPy ``ndarray`` objects, we will be able to bring more consistency to
  pandas's semantics and enable the core developers to extend pandas more
  cleanly with new data types, data structures, and computational semantics.
* **Exposing a pandas Cython and/or C/C++ API to other Python library
  developers**: the internals of Series and DataFrame are only weakly
  accessible in other developers' native code. At minimum, we wish to better
  enable developers to construct the precise data structures / memory
  representation that fill the insides of Series and DataFrame.
* **Improving user control and visibility of memory use**: pandas's memory use,
  as a result of its internal implementation, can frequently be opaque to the
  user or outright unpredictable.
* **Improving performance and system utilization**: We aim to improve both the
  micro (operations that take < 1 ms) and macro (all other operations)
  performance of pandas across the board. As part of this, we aim to make it
  easier for pandas's core developers to leverage multicore systems to
  accelerate computations (without running into any of Python's well-known
  concurrency limitations)
* **Removal of deprecated / underutilized functionality**: As the Python data
  ecosystem has grown, a number of areas of pandas (e.g. plotting and datasets
  with more than 2 dimensions) may be better served by other open source
  projects. Also, functionality that has been explicitly deprecated or
  discouraged from use (like the ``.ix`` indexing operator) would ideally be
  removed.

Non-goals / FAQ
===============

As this will be a quite nuanced discussion, especially for those not intimately
familiar with pandas's implementation details, I wanted to speak to a couple of
commonly-asked questions in brief:

1. **Will this work make it harder to use pandas with NumPy, scikit-learn,
   statsmodels, SciPy, or other libraries that depend on NumPy
   interoperability?**
  * We are not planning on it. Data that is representable without memory
    copying or conversion in NumPy arrays will continue to be 100%
    interoperable.
  * Data containing missing (NA) values may require explicit conversion where
    it is not currently required. For example: integer or boolean type arrays
    with missing data. I trust this will be seen as a positive development.
  * If anything, more performant and more precise data semantics in pandas will
    generally make production code using a downstream library like scikit-learn
    more dependable and future-proof.

2. **By decoupling from NumPy, it sounds like you are reimplementing NumPy or
   adding a new data type system**

   * Simply put: no. But it's more complicated than that because of the
     numerous interpretations of "type system".

   * pandas already contains a large amount (10s of KLOCs) of custom
     computational code (see, for example,
     `https://github.com/pydata/pandas/tree/master/pandas/src`) that implements
     functionality not present in NumPy.

   * pandas already features its own (what can be called) "logical type
     system", including things like custom data types (such as that of
     ``pandas.Categorical``), pandas-specific missing data representation, and
     implicit type casting (e.g. integer to float on introduction of missing
     data). Unfortunately, these logical data types are somewhat weakly
     expressed, and the mix of NumPy dtype objects and custom pandas types is
     problematic for many internal (implementation) and external (user API)
     reasons. I will examine in detail the difference between **physical
     types** (i.e. NumPy's dtypes) and **logical types** (i.e. what pandas
     currently has, implicitly).

Summary
=======

Overall, the goal of the pandas 2.0 project is to yield a faster, more cleanly
architected, and more future-proof library that is a drop-in replacement for
90-95% of pandas user code. There will be API / code breakages, but the intent
of any code breakage will almost always be to fix something that has been
"wrong" or inconsistent. Many advanced users will have worked around some of
these rough edges, and so their workarounds may either need to be removed or
changed to accommodate the new (and hopefully it can be agreed in each case:
better) semantics.
