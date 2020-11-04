.. _test_organization:

Test Organization
=================
Ideally, there should be one, and only one, obvious place for a test to reside.
Until we reach that ideal, these are some rules of thumb for where a test should
be located.

1) Does your test depend only on code in ``pd._libs.tslibs``?
    This test likely belongs in one of:
    - tests.tslibs
        Policy: No file in ``tests.tslibs`` should import from any pandas modules outside of ``pd._libs.tslibs``
    - tests.scalar
    - tests.tseries.offsets

2) Does your test depend only on code in pd._libs?
    This test likely belongs in one of:
    - tests.libs
    - tests.groupby.test_libgroupby

3) Is your test for an arithmetic or comparison method?
    This test likely belongs in one of:
    - tests.arithmetic
        These are intended for tests that can be shared to test the behavior of DataFrame/Series/Index/ExtensionArray using the ``box_with_array`` fixture.
    - tests.frame.test_arithmetic
    - tests.series.test_arithmetic

4) Is your test for a reduction method (min, max, sum, prod, ...)?
    This test likely belongs in one of:
    - tests.reductions
        These are intended for tests that can be shared to test the behavior of DataFrame/Series/Index/ExtensionArray.
    - tests.frame.test_reductions
    - tests.series.test_reductions
    - tests.test_nanops

5) Is your test for a DataFrame or Series method?
    A) Is the method a plotting method?
        This test likely belongs in one of:
        - tests.plotting
    B) Is the method an IO method?
        This test likely belongs in one of:
        - tests.io

    C) Is the method an indexing method?
        i) Is the method ``loc``, ``iloc``, ``at``, or ``iat``?
            This test likely belongs in one of:
            - tests.indexing.test_methodname
        ii) Otherwise
            This test likely belongs in one of:
            - tests.frame.indexing.test_methodname
            - tests.series.indexing.test_methodname

    D) Otherwise
        This test likely belongs in one of:
        - tests.series.methods.test_mymethod
        - tests.frame.methods.test_mymethod
            If a test can be shared between DataFrame/Series using the ``frame_or_series`` fixture, by convention it goes in tests.frame file.
        - tests.generic.methods.test_mymethod
            The generic/methods/ directory is only for methods with tests
            that are fully parametrized over Series+DataFrame

6) Is your test for an Index method, not depending on Series/DataFrame?
    This test likely belongs in one of:
    - tests.indexes

7) Is your test for one of the pandas-provided ExtensionArrays (Categorical, DatetimeArray, TimedeltaArray, PeriodArray, IntervalArray, PandasArray, FloatArray, BoolArray, IntervalArray, StringArray)?
    This test likely belongs in one of:
    - tests.arrays

8) Is your test for _all_ ExtensionArray subclasses (the "EA Interface")?
    This test likely belongs in one of:
    - tests.extension
