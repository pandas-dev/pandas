"""
Tests are extremely important for ndindex. All operations should produce
correct results. We test this by checking against numpy arange (the array
values do not matter, so long as they are distinct).

There are three primary types of tests that we employ to verify this:

- Exhaustive tests. These test every possible value in some range. See for
  example test_slice. This is the best type of test, but unfortunately, it is
  often impossible to do due to combinatorial explosion.

- Hypothesis tests. Hypothesis is a library that can intelligently check a
  combinatorial search space. This requires writing hypothesis strategies that
  can generate all the relevant types of indices (see helpers.py). For more
  information on hypothesis, see
  https://hypothesis.readthedocs.io/en/latest/index.html.

- Explicit tests. These are hand crafted tests that test that the output of a
  function is some exact value. These are used when a property-based test is
  difficult to write, and it is simpler to test the exact output value. These
  tests still include a correctness check (check_same()) to make sure the
  hard-coded output is actually correct.

The basic idea in both cases is the same. Take the pure index and the
ndindex(index).raw, or in the case of a transformation, the before and after
raw index, and index an arange with them. If they do not give the same output
array, or do not both produce the same error, the code is not correct.

Why bother with hypothesis if the same thing is already tested exhaustively?
The main reason is that hypothesis is much better at producing human-readable
failure examples. When an exhaustive test fails, the failure will always be
from the first set of inputs in the loop that produces a failure. Hypothesis
on the other hand attempts to "shrink" the failure input to smallest input
that still fails. For example, a failing exhaustive slice test might give
Slice(-10, -9, -10) as a the failing example, but hypothesis would shrink it
to Slice(-2, -1, -1). Another reason for the duplication is that hypothesis
can sometimes test a slightly expanded test space without any additional
consequences. For example, test_slice_reduce_hypothesis() tests all types of
array shapes, whereas test_slice_reduce_exhaustive() tests only 1-dimensional
shapes. This doesn't affect things because hypothesis will always shrink large
shapes to a 1-dimensional shape in the case of a failure. Consequently every
exhaustive test should have a corresponding hypothesis test.

For things that can only be tested with hypothesis, you can use @example, to
force certain combinations to be tested. This is useful because we require
100% test coverage, and hypothesis's randomness can cause this to be flaky
otherwise.

"""

# Variable naming conventions in the tests:

# a: numpy arange. May be reshaped to be multidimensional
# shape: a tuple of integers
# i: integer used as an integer index
# idx: generic index (Python type)
# index: generic index (ndindex type)
# s: slice (Python type)
# S: Slice (ndindex type)
# size: integer passed to arange
