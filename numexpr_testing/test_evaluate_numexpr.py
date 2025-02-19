import numpy as np
import numexpr as ne
import warnings
import pytest
import sys
from operator import add, sub, mul, truediv
import math

# Branch coverage tracking
branch_coverage = {}

def record_branch(branch_id):
    """Track execution count of a specific branch."""
    if branch_id in branch_coverage:
        branch_coverage[branch_id] += 1
    else:
        branch_coverage[branch_id] = 1

# Supporting functions
def _can_use_numexpr(op, op_str, left_op, right_op, eval_type):
    """Simplified check for numexpr compatibility."""
    return True

def _bool_arith_fallback(op_str, left_op, right_op):
    """Simplified boolean arithmetic fallback."""
    return True

def _evaluate_standard(op, op_str, left_op, right_op):
    """Standard evaluation fallback."""
    return op(left_op, right_op)

_TEST_MODE = False

def _store_test_result(result):
    """Store test results (simplified)."""
    pass

def _evaluate_numexpr(op, op_str, left_op, right_op):
    record_branch(1)  # Entry point
    result = None

    if _can_use_numexpr(op, op_str, left_op, right_op, "evaluate"):
        record_branch(2)
        is_reversed = hasattr(op, '__name__') and op.__name__.startswith('r')
        
        if is_reversed:
            record_branch(3)
            left_op, right_op = right_op, left_op

        left_value = left_op
        right_value = right_op

        try:
            record_branch(4)
            result = ne.evaluate(
                f"left_value {op_str} right_value",
                local_dict={"left_value": left_value, "right_value": right_value},
                casting="safe",
            )
        except TypeError:
            record_branch(5)
            pass
        except NotImplementedError:
            record_branch(6)
            if _bool_arith_fallback(op_str, left_op, right_op):
                record_branch(7)
                pass
            else:
                record_branch(8)
                raise

        if is_reversed:
            record_branch(9)
            left_op, right_op = right_op, left_op

    if _TEST_MODE:
        record_branch(10)
        _store_test_result(result is not None)

    if result is None:
        record_branch(11)
        result = _evaluate_standard(op, op_str, left_op, right_op)

    record_branch(12)
    return result

def analyze_coverage():
    """Print branch coverage analysis"""
    total_branches = 12  # Total number of branches in the code
    covered_branches = sum(1 for count in branch_coverage.values() if count > 0)
    
    print("\nBranch Coverage Analysis:")
    print("-" * 50)
    print(f"Total branches: {total_branches}")
    print(f"Covered branches: {covered_branches}")
    coverage_percentage = (covered_branches / total_branches * 100) if total_branches > 0 else 0
    print(f"Coverage percentage: {coverage_percentage:.2f}%")
    print("\nDetailed branch execution counts:")
    for branch_id, count in sorted(branch_coverage.items()):
        print(f"Branch {branch_id}: {count} times")
    print("-" * 50)
    return covered_branches, total_branches

# Create a reversed add operation
def radd(x, y):
    radd.__name__ = 'radd'
    return add(y, x)

@pytest.fixture
def setup_arrays():
    """Fixture for common test arrays"""
    return {
        'basic': (np.array([1, 2, 3]), np.array([4, 5, 6])),
        'float': (np.array([1.5, 2.5, 3.5]), np.array([1.0, 2.0, 3.0])),
        'bool': (np.array([True, False, True]), np.array([True, True, False])),
        'small': (np.arange(10), np.arange(10)),
        'large': (np.arange(100000), np.arange(100000))
    }

# Category 1: Basic Arithmetic Operations
def test_basic_arithmetic(setup_arrays):
    left, right = setup_arrays['basic']
    
    # Addition
    result = _evaluate_numexpr(add, '+', left, right)
    assert np.array_equal(result, np.array([5, 7, 9]))
    
    # Subtraction
    result = _evaluate_numexpr(sub, '-', left, right)
    assert np.array_equal(result, np.array([-3, -3, -3]))
    
    # Multiplication
    result = _evaluate_numexpr(mul, '*', left, right)
    assert np.array_equal(result, np.array([4, 10, 18]))
    
    # Division
    result = _evaluate_numexpr(truediv, '/', left, right)
    np.testing.assert_array_almost_equal(result, np.array([0.25, 0.4, 0.5]))

# Category 2: Reversed Operations
def test_reversed_operations(setup_arrays):
    left, right = setup_arrays['basic']
    
    def reversed_sub(x, y):
        reversed_sub.__name__ = 'rsub'
        return y - x
    
    result = _evaluate_numexpr(reversed_sub, '-', left, right)
    assert np.array_equal(result, np.array([3, 3, 3]))

# Category 3: Different Data Types
def test_different_dtypes(setup_arrays):
    int_arr, float_arr = np.array([1, 2, 3]), np.array([1.5, 2.5, 3.5])
    bool_arr = np.array([True, False, True])
    
    # Integer + Float
    result = _evaluate_numexpr(add, '+', int_arr, float_arr)
    np.testing.assert_array_almost_equal(result, np.array([2.5, 4.5, 6.5]))
    
    # Boolean operations
    result = _evaluate_numexpr(add, '+', bool_arr, bool_arr)
    assert np.array_equal(result, np.array([True, False, True]))  # Changed this line
    
    # Test boolean with integers
    result = _evaluate_numexpr(add, '+', bool_arr, np.array([1, 1, 1]))
    assert np.array_equal(result, np.array([2, 1, 2]))

# Category 4: Size Thresholds
def test_size_thresholds(setup_arrays):
    small_left, small_right = setup_arrays['small']
    large_left, large_right = setup_arrays['large']
    
    # Small arrays
    result_small = _evaluate_numexpr(add, '+', small_left, small_right)
    assert np.array_equal(result_small, small_left + small_right)
    
    # Large arrays
    result_large = _evaluate_numexpr(add, '+', large_left, large_right)
    assert np.array_equal(result_large, large_left + large_right)

# Category 5: Edge Cases
def test_edge_cases():
    # Empty arrays
    result = _evaluate_numexpr(add, '+', np.array([]), np.array([]))
    assert len(result) == 0
    
    # Single element
    result = _evaluate_numexpr(add, '+', np.array([1]), np.array([2]))
    assert np.array_equal(result, np.array([3]))
    
    # NaN values
    nan_arr = np.array([np.nan, 1, 2])
    result = _evaluate_numexpr(add, '+', nan_arr, np.array([1, 2, 3]))
    assert np.isnan(result[0])
    assert np.array_equal(result[1:], np.array([3, 5]))
    
    # Infinity
    inf_arr = np.array([np.inf, -np.inf, 1])
    result = _evaluate_numexpr(add, '+', inf_arr, np.array([1, 1, 1]))
    assert np.isinf(result[0]) and result[0] > 0
    assert np.isinf(result[1]) and result[1] < 0
    assert result[2] == 2

@pytest.fixture(scope="session", autouse=True)
def final_coverage():
    yield
    covered_branches, total_branches = analyze_coverage()
    coverage_percentage = (covered_branches / total_branches * 100)
    print(f"\nFinal coverage: {coverage_percentage:.2f}%")
    assert coverage_percentage >= 90, f"Branch coverage ({coverage_percentage:.2f}%) is below 90%"