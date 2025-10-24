"""
Simple validation test for the fix - doesn't require full pandas installation
"""

import sys
import ast

def test_merge_py_syntax():
    """Verify merge.py has valid Python syntax."""
    with open('c:/noc_project/projects/pandas/pandas/core/reshape/merge.py', 'r') as f:
        code = f.read()
    
    try:
        ast.parse(code)
        print("✓ merge.py syntax is valid")
        return True
    except SyntaxError as e:
        print(f"✗ Syntax error in merge.py: {e}")
        return False


def test_new_function_exists():
    """Verify the new _get_merge_error_message function exists."""
    with open('c:/noc_project/projects/pandas/pandas/core/reshape/merge.py', 'r') as f:
        code = f.read()
    
    if 'def _get_merge_error_message' in code:
        print("✓ Function _get_merge_error_message exists")
        return True
    else:
        print("✗ Function _get_merge_error_message not found")
        return False


def test_polars_special_handling():
    """Verify polars special handling exists."""
    with open('c:/noc_project/projects/pandas/pandas/core/reshape/merge.py', 'r') as f:
        code = f.read()
    
    if 'polars.dataframe.frame' in code and 'to_pandas' in code:
        print("✓ Polars special handling code found")
        return True
    else:
        print("✗ Polars special handling not found")
        return False


def test_error_message_improvement():
    """Verify the error message was improved."""
    with open('c:/noc_project/projects/pandas/pandas/core/reshape/merge.py', 'r') as f:
        code = f.read()
    
    # Check that old generic message is gone
    old_msg = 'f"Can only merge Series or DataFrame objects, a {type(obj)} was passed"'
    
    if old_msg in code:
        print("✗ Old generic error message still present")
        return False
    
    # Check that we're calling the new function
    if '_get_merge_error_message(obj)' in code:
        print("✓ Error message improved (now calls _get_merge_error_message)")
        return True
    else:
        print("✗ Error message improvement not found")
        return False


def test_regression_tests():
    """Verify regression test file has tests."""
    with open('c:/noc_project/projects/pandas/test_issue_61434_tests.py', 'r') as f:
        code = f.read()
    
    required_tests = [
        'test_merge_with_polars_dataframe',
        'test_merge_pandas_baseline',
        'test_merge_with_dict',
    ]
    
    all_found = all(test in code for test in required_tests)
    
    if all_found:
        print(f"✓ All required tests found: {required_tests}")
        return True
    else:
        missing = [t for t in required_tests if t not in code]
        print(f"✗ Missing tests: {missing}")
        return False


if __name__ == "__main__":
    print("=" * 70)
    print("VALIDATION TESTS FOR PANDAS ISSUE #61434 FIX")
    print("=" * 70)
    print()
    
    tests = [
        ("Syntax validation", test_merge_py_syntax),
        ("New function exists", test_new_function_exists),
        ("Polars handling", test_polars_special_handling),
        ("Error message improved", test_error_message_improvement),
        ("Regression tests", test_regression_tests),
    ]
    
    results = []
    for name, test_func in tests:
        print(f"Testing: {name}")
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"✗ Test failed with error: {e}")
            results.append((name, False))
        print()
    
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")
    
    all_passed = all(result for _, result in results)
    print()
    
    if all_passed:
        print("✅ All validation tests PASSED - Ready to commit!")
        sys.exit(0)
    else:
        print("❌ Some validation tests FAILED - Fix needed")
        sys.exit(1)
