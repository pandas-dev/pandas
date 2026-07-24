"""
Test that ndindex does not depend on NumPy, as long as no functionality
involving array indices is used.

Also test that it does not depend on SymPy (which used to be a hard
dependency, see ndindex/_crt.py).

Because the tests themselves require NumPy, we must test this in a subprocess
where we know that NumPy is not yet imported.

"""

import sys
import subprocess

import pytest

from ndindex import (ndindex, iter_indices, Integer, Slice, Tuple, ellipsis,
                     Newaxis, ChunkSize, BroadcastError, AxisError)

def _test_dependency_ndindex(mod):
    assert mod not in sys.modules

    # Run various functions that shouldn't depend on NumPy and check that they
    # don't import it.

    ndindex(1)
    assert mod not in sys.modules

    ndindex[0:1]
    assert mod not in sys.modules

    ndindex(ndindex(1))
    assert mod not in sys.modules

    ndindex(None)
    assert mod not in sys.modules

    ndindex[1, 0:1]
    assert mod not in sys.modules

    ndindex(...)
    assert mod not in sys.modules

    try:
        ndindex('')
    except IndexError:
        pass
    assert mod not in sys.modules

def _test_dependency_iter_indices(mod):
    assert mod not in sys.modules

    list(iter_indices((1, 2), (2, 1)))
    assert mod not in sys.modules

    list(iter_indices((1, 2), (2, 1), skip_axes=(1,)))
    assert mod not in sys.modules

    try:
        list(iter_indices((2, 3), (3, 2)))
    except BroadcastError:
        pass
    else:
        assert False
    assert mod not in sys.modules

    try:
        list(iter_indices((2, 3), skip_axes=(2,)))
    except AxisError:
        pass
    else:
        assert False
    assert mod not in sys.modules

def _test_dependency_indices(mod):
    for index_type, args in [
            (Integer, (1,)),
            (Slice, (2, 10, 3)),
            (ellipsis, ()),
            (Newaxis, ()),
            (Tuple, (slice(0, 3), 1)),
    ]:
        assert mod not in sys.modules, index_type

        idx = index_type(*args)
        assert mod not in sys.modules, index_type

        idx.raw
        assert mod not in sys.modules, index_type

        if hasattr(idx, '__len__'):
            len(idx)
        assert mod not in sys.modules, index_type

        idx.reduce()
        assert mod not in sys.modules, index_type

        idx.reduce((2, 3))
        assert mod not in sys.modules, index_type

        idx.newshape((2, 3))
        assert mod not in sys.modules, index_type

        try:
            idx.as_subindex(slice(1, 20, 5))
        except NotImplementedError:
            pass
        assert mod not in sys.modules, index_type

        idx.isempty((2, 3))
        assert mod not in sys.modules, index_type

        idx == idx
        assert mod not in sys.modules, index_type

        idx == 1
        assert mod not in sys.modules, index_type

        hash(idx)
        assert mod not in sys.modules, index_type

        idx.broadcast_arrays()
        assert mod not in sys.modules, index_type

        idx.expand((2, 3))
        assert mod not in sys.modules, index_type

        str(idx)
        assert mod not in sys.modules, index_type

        repr(idx)
        assert mod not in sys.modules, index_type

def _test_dependency_ChunkSize(mod):
    assert mod not in sys.modules

    c = ChunkSize((10, 10))
    assert mod not in sys.modules

    hash(c)
    assert mod not in sys.modules

    c[0]
    assert mod not in sys.modules

    len(c)
    assert mod not in sys.modules

    c.num_chunks((100, 100))
    assert mod not in sys.modules

    list(c.indices((100, 100)))
    assert mod not in sys.modules

    list(c.as_subchunks(slice(0, 50), (100, 100)))
    assert mod not in sys.modules

    c.num_subchunks(slice(0, 50), (100, 100))
    assert mod not in sys.modules

    c.containing_block(slice(0, 50), (100, 100))
    assert mod not in sys.modules

    list(c.as_subchunks(None, (100, 100)))
    assert mod not in sys.modules

    c.num_subchunks(None, (100, 100))
    assert mod not in sys.modules

    c.containing_block(None, (100, 100))
    assert mod not in sys.modules

    list(c.as_subchunks(1, (100, 100)))
    assert mod not in sys.modules

    c.num_subchunks(1, (100, 100))
    assert mod not in sys.modules

    c.containing_block(1, (100, 100))
    assert mod not in sys.modules

    list(c.as_subchunks(..., (100, 100)))
    assert mod not in sys.modules

    c.num_subchunks(..., (100, 100))
    assert mod not in sys.modules

    c.containing_block(..., (100, 100))
    assert mod not in sys.modules

    list(c.as_subchunks((0, slice(0, 50)), (100, 100)))
    assert mod not in sys.modules

    c.num_subchunks((0, slice(0, 50)), (100, 100))
    assert mod not in sys.modules

    c.containing_block((0, slice(0, 50)), (100, 100))
    assert mod not in sys.modules

@pytest.mark.no_cover
@pytest.mark.parametrize('func', [func for func in globals() if func.startswith('_test')])
def test_numpy_dependency(func):
    # This unfortunately won't go through any of the pytest machinery. We
    # reraise the exception as an AssertionError so that pytest will show it
    # in a semi-reasonable way

    try:
        subprocess.run([sys.executable, '-c', f'''\
from ndindex.tests.test_no_dependencies import {func}

{func}("numpy")'''], check=True, capture_output=True, encoding='utf-8')
    except subprocess.CalledProcessError as e:
        print(e.stdout, end='')
        raise AssertionError(e.stderr)

@pytest.mark.no_cover
@pytest.mark.parametrize('func', [func for func in globals() if func.startswith('_test')])
def test_sympy_dependency(func):
    # This unfortunately won't go through any of the pytest machinery. We
    # reraise the exception as an AssertionError so that pytest will show it
    # in a semi-reasonable way

    try:
        subprocess.run([sys.executable, '-c', f'''\
from ndindex.tests.test_no_dependencies import {func}

{func}("sympy")'''], check=True, capture_output=True, encoding='utf-8')
    except subprocess.CalledProcessError as e:
        print(e.stdout, end='')
        raise AssertionError(e.stderr)
