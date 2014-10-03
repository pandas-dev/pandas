from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import (
    run_module_suite, assert_equal, assert_array_equal,
    assert_raises
    )
from numpy.lib.stride_tricks import as_strided, broadcast_arrays


def assert_shapes_correct(input_shapes, expected_shape):
    # Broadcast a list of arrays with the given input shapes and check the
    # common output shape.

    inarrays = [np.zeros(s) for s in input_shapes]
    outarrays = broadcast_arrays(*inarrays)
    outshapes = [a.shape for a in outarrays]
    expected = [expected_shape] * len(inarrays)
    assert_equal(outshapes, expected)


def assert_incompatible_shapes_raise(input_shapes):
    # Broadcast a list of arrays with the given (incompatible) input shapes
    # and check that they raise a ValueError.

    inarrays = [np.zeros(s) for s in input_shapes]
    assert_raises(ValueError, broadcast_arrays, *inarrays)


def assert_same_as_ufunc(shape0, shape1, transposed=False, flipped=False):
    # Broadcast two shapes against each other and check that the data layout
    # is the same as if a ufunc did the broadcasting.

    x0 = np.zeros(shape0, dtype=int)
    # Note that multiply.reduce's identity element is 1.0, so when shape1==(),
    # this gives the desired n==1.
    n = int(np.multiply.reduce(shape1))
    x1 = np.arange(n).reshape(shape1)
    if transposed:
        x0 = x0.T
        x1 = x1.T
    if flipped:
        x0 = x0[::-1]
        x1 = x1[::-1]
    # Use the add ufunc to do the broadcasting. Since we're adding 0s to x1, the
    # result should be exactly the same as the broadcasted view of x1.
    y = x0 + x1
    b0, b1 = broadcast_arrays(x0, x1)
    assert_array_equal(y, b1)


def test_same():
    x = np.arange(10)
    y = np.arange(10)
    bx, by = broadcast_arrays(x, y)
    assert_array_equal(x, bx)
    assert_array_equal(y, by)


def test_one_off():
    x = np.array([[1, 2, 3]])
    y = np.array([[1], [2], [3]])
    bx, by = broadcast_arrays(x, y)
    bx0 = np.array([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
    by0 = bx0.T
    assert_array_equal(bx0, bx)
    assert_array_equal(by0, by)


def test_same_input_shapes():
    # Check that the final shape is just the input shape.

    data = [
        (),
        (1,),
        (3,),
        (0, 1),
        (0, 3),
        (1, 0),
        (3, 0),
        (1, 3),
        (3, 1),
        (3, 3),
    ]
    for shape in data:
        input_shapes = [shape]
        # Single input.
        assert_shapes_correct(input_shapes, shape)
        # Double input.
        input_shapes2 = [shape, shape]
        assert_shapes_correct(input_shapes2, shape)
        # Triple input.
        input_shapes3 = [shape, shape, shape]
        assert_shapes_correct(input_shapes3, shape)


def test_two_compatible_by_ones_input_shapes():
    # Check that two different input shapes of the same length, but some have
    # ones, broadcast to the correct shape.

    data = [
        [[(1,), (3,)], (3,)],
        [[(1, 3), (3, 3)], (3, 3)],
        [[(3, 1), (3, 3)], (3, 3)],
        [[(1, 3), (3, 1)], (3, 3)],
        [[(1, 1), (3, 3)], (3, 3)],
        [[(1, 1), (1, 3)], (1, 3)],
        [[(1, 1), (3, 1)], (3, 1)],
        [[(1, 0), (0, 0)], (0, 0)],
        [[(0, 1), (0, 0)], (0, 0)],
        [[(1, 0), (0, 1)], (0, 0)],
        [[(1, 1), (0, 0)], (0, 0)],
        [[(1, 1), (1, 0)], (1, 0)],
        [[(1, 1), (0, 1)], (0, 1)],
    ]
    for input_shapes, expected_shape in data:
        assert_shapes_correct(input_shapes, expected_shape)
        # Reverse the input shapes since broadcasting should be symmetric.
        assert_shapes_correct(input_shapes[::-1], expected_shape)


def test_two_compatible_by_prepending_ones_input_shapes():
    # Check that two different input shapes (of different lengths) broadcast
    # to the correct shape.

    data = [
        [[(), (3,)], (3,)],
        [[(3,), (3, 3)], (3, 3)],
        [[(3,), (3, 1)], (3, 3)],
        [[(1,), (3, 3)], (3, 3)],
        [[(), (3, 3)], (3, 3)],
        [[(1, 1), (3,)], (1, 3)],
        [[(1,), (3, 1)], (3, 1)],
        [[(1,), (1, 3)], (1, 3)],
        [[(), (1, 3)], (1, 3)],
        [[(), (3, 1)], (3, 1)],
        [[(), (0,)], (0,)],
        [[(0,), (0, 0)], (0, 0)],
        [[(0,), (0, 1)], (0, 0)],
        [[(1,), (0, 0)], (0, 0)],
        [[(), (0, 0)], (0, 0)],
        [[(1, 1), (0,)], (1, 0)],
        [[(1,), (0, 1)], (0, 1)],
        [[(1,), (1, 0)], (1, 0)],
        [[(), (1, 0)], (1, 0)],
        [[(), (0, 1)], (0, 1)],
    ]
    for input_shapes, expected_shape in data:
        assert_shapes_correct(input_shapes, expected_shape)
        # Reverse the input shapes since broadcasting should be symmetric.
        assert_shapes_correct(input_shapes[::-1], expected_shape)


def test_incompatible_shapes_raise_valueerror():
    # Check that a ValueError is raised for incompatible shapes.

    data = [
        [(3,), (4,)],
        [(2, 3), (2,)],
        [(3,), (3,), (4,)],
        [(1, 3, 4), (2, 3, 3)],
    ]
    for input_shapes in data:
        assert_incompatible_shapes_raise(input_shapes)
        # Reverse the input shapes since broadcasting should be symmetric.
        assert_incompatible_shapes_raise(input_shapes[::-1])


def test_same_as_ufunc():
    # Check that the data layout is the same as if a ufunc did the operation.

    data = [
        [[(1,), (3,)], (3,)],
        [[(1, 3), (3, 3)], (3, 3)],
        [[(3, 1), (3, 3)], (3, 3)],
        [[(1, 3), (3, 1)], (3, 3)],
        [[(1, 1), (3, 3)], (3, 3)],
        [[(1, 1), (1, 3)], (1, 3)],
        [[(1, 1), (3, 1)], (3, 1)],
        [[(1, 0), (0, 0)], (0, 0)],
        [[(0, 1), (0, 0)], (0, 0)],
        [[(1, 0), (0, 1)], (0, 0)],
        [[(1, 1), (0, 0)], (0, 0)],
        [[(1, 1), (1, 0)], (1, 0)],
        [[(1, 1), (0, 1)], (0, 1)],
        [[(), (3,)], (3,)],
        [[(3,), (3, 3)], (3, 3)],
        [[(3,), (3, 1)], (3, 3)],
        [[(1,), (3, 3)], (3, 3)],
        [[(), (3, 3)], (3, 3)],
        [[(1, 1), (3,)], (1, 3)],
        [[(1,), (3, 1)], (3, 1)],
        [[(1,), (1, 3)], (1, 3)],
        [[(), (1, 3)], (1, 3)],
        [[(), (3, 1)], (3, 1)],
        [[(), (0,)], (0,)],
        [[(0,), (0, 0)], (0, 0)],
        [[(0,), (0, 1)], (0, 0)],
        [[(1,), (0, 0)], (0, 0)],
        [[(), (0, 0)], (0, 0)],
        [[(1, 1), (0,)], (1, 0)],
        [[(1,), (0, 1)], (0, 1)],
        [[(1,), (1, 0)], (1, 0)],
        [[(), (1, 0)], (1, 0)],
        [[(), (0, 1)], (0, 1)],
    ]
    for input_shapes, expected_shape in data:
        assert_same_as_ufunc(input_shapes[0], input_shapes[1],
                             "Shapes: %s %s" % (input_shapes[0], input_shapes[1]))
        # Reverse the input shapes since broadcasting should be symmetric.
        assert_same_as_ufunc(input_shapes[1], input_shapes[0])
        # Try them transposed, too.
        assert_same_as_ufunc(input_shapes[0], input_shapes[1], True)
        # ... and flipped for non-rank-0 inputs in order to test negative
        # strides.
        if () not in input_shapes:
            assert_same_as_ufunc(input_shapes[0], input_shapes[1], False, True)
            assert_same_as_ufunc(input_shapes[0], input_shapes[1], True, True)

def test_as_strided():
    a = np.array([None])
    a_view = as_strided(a)
    expected = np.array([None])
    assert_array_equal(a_view, np.array([None]))

    a = np.array([1, 2, 3, 4])
    a_view = as_strided(a, shape=(2,), strides=(2 * a.itemsize,))
    expected = np.array([1, 3])
    assert_array_equal(a_view, expected)

    a = np.array([1, 2, 3, 4])
    a_view = as_strided(a, shape=(3, 4), strides=(0, 1 * a.itemsize))
    expected = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
    assert_array_equal(a_view, expected)


if __name__ == "__main__":
    run_module_suite()
