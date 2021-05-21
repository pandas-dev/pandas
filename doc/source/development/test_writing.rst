.. _test_organization:

Test organization
=================
Ideally, there should be one, and only one, obvious place for a test to reside.
Until we reach that ideal, these are some rules of thumb for where a test should
be located.

1. Does your test depend only on code in ``pd._libs.tslibs``?
   This test likely belongs in one of:

   - tests.tslibs

     .. note::

          No file in ``tests.tslibs`` should import from any pandas modules
          outside of ``pd._libs.tslibs``

   - tests.scalar
   - tests.tseries.offsets

2. Does your test depend only on code in pd._libs?
   This test likely belongs in one of:

   - tests.libs
   - tests.groupby.test_libgroupby

3. Is your test for an arithmetic or comparison method?
   This test likely belongs in one of:

   - tests.arithmetic

     .. note::

         These are intended for tests that can be shared to test the behavior
         of DataFrame/Series/Index/ExtensionArray using the ``box_with_array``
         fixture.

   - tests.frame.test_arithmetic
   - tests.series.test_arithmetic

4. Is your test for a reduction method (min, max, sum, prod, ...)?
   This test likely belongs in one of:

   - tests.reductions

     .. note::

         These are intended for tests that can be shared to test the behavior
         of DataFrame/Series/Index/ExtensionArray.

   - tests.frame.test_reductions
   - tests.series.test_reductions
   - tests.test_nanops

5. Is your test for an indexing method?
   This is the most difficult case for deciding where a test belongs, because
   there are many of these tests, and many of them test more than one method
   (e.g. both ``Series.__getitem__`` and ``Series.loc.__getitem__``)

   A) Is the test specifically testing an Index method (e.g. ``Index.get_loc``,
      ``Index.get_indexer``)?
      This test likely belongs in one of:

      - tests.indexes.test_indexing
      - tests.indexes.fooindex.test_indexing

      Within that files there should be a method-specific test class e.g.
      ``TestGetLoc``.

      In most cases, neither ``Series`` nor ``DataFrame`` objects should be
      needed in these tests.

   B) Is the test for a Series or DataFrame indexing method *other* than
      ``__getitem__`` or ``__setitem__``, e.g. ``xs``, ``where``, ``take``,
      ``mask``, ``lookup``, or ``insert``?
      This test likely belongs in one of:

      - tests.frame.indexing.test_methodname
      - tests.series.indexing.test_methodname

   C) Is the test for any of ``loc``, ``iloc``, ``at``, or ``iat``?
      This test likely belongs in one of:

      - tests.indexing.test_loc
      - tests.indexing.test_iloc
      - tests.indexing.test_at
      - tests.indexing.test_iat

      Within the appropriate file, test classes correspond to either types of
      indexers (e.g. ``TestLocBooleanMask``) or major use cases
      (e.g. ``TestLocSetitemWithExpansion``).

      See the note in section D) about tests that test multiple indexing methods.

   D) Is the test for ``Series.__getitem__``, ``Series.__setitem__``,
      ``DataFrame.__getitem__``, or ``DataFrame.__setitem__``?
      This test likely belongs in one of:

      - tests.series.test_getitem
      - tests.series.test_setitem
      - tests.frame.test_getitem
      - tests.frame.test_setitem

      If many cases such a test may test multiple similar methods, e.g.

      .. code-block:: python

           import pandas as pd
           import pandas._testing as tm

           def test_getitem_listlike_of_ints():
               ser = pd.Series(range(5))

               result = ser[[3, 4]]
               expected = pd.Series([2, 3])
               tm.assert_series_equal(result, expected)

               result = ser.loc[[3, 4]]
               tm.assert_series_equal(result, expected)

    In cases like this, the test location should be based on the *underlying*
    method being tested.  Or in the case of a test for a bugfix, the location
    of the actual bug.  So in this example, we know that ``Series.__getitem__``
    calls ``Series.loc.__getitem__``, so this is *really* a test for
    ``loc.__getitem__``.  So this test belongs in ``tests.indexing.test_loc``.

6. Is your test for a DataFrame or Series method?

   A) Is the method a plotting method?
      This test likely belongs in one of:

      - tests.plotting

   B) Is the method an IO method?
      This test likely belongs in one of:

      - tests.io

   C) Otherwise
      This test likely belongs in one of:

      - tests.series.methods.test_mymethod
      - tests.frame.methods.test_mymethod

        .. note::

            If a test can be shared between DataFrame/Series using the
            ``frame_or_series`` fixture, by convention it goes in the
            ``tests.frame`` file.

7. Is your test for an Index method, not depending on Series/DataFrame?
   This test likely belongs in one of:

   - tests.indexes

8) Is your test for one of the pandas-provided ExtensionArrays (``Categorical``,
   ``DatetimeArray``, ``TimedeltaArray``, ``PeriodArray``, ``IntervalArray``,
   ``PandasArray``, ``FloatArray``, ``BoolArray``, ``StringArray``)?
   This test likely belongs in one of:

   - tests.arrays

9) Is your test for *all* ExtensionArray subclasses (the "EA Interface")?
   This test likely belongs in one of:

   - tests.extension
