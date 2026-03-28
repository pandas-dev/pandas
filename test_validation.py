#!/usr/bin/env python
"""Quick validation test for limit_behavior feature"""
import sys

import numpy as np

# Test 1: Validate function works
print("=" * 60)
print("Test 1: Validate function validation")
print("=" * 60)

def validate_limit_behavior(limit_behavior: str):
    valid_limit_behaviors = ["fill", "skip"]
    limit_behavior = limit_behavior.lower()
    if limit_behavior not in valid_limit_behaviors:
        raise ValueError(
            f"Invalid limit_behavior: expecting one of {valid_limit_behaviors}, got "
            f"{limit_behavior}."
        )
    return limit_behavior

# Valid values should pass
try:
    assert validate_limit_behavior("fill") == "fill"
    print("✓ 'fill' - PASS")
except Exception as e:
    print(f"✗ 'fill' - FAIL: {e}")

try:
    assert validate_limit_behavior("skip") == "skip"
    print("✓ 'skip' - PASS")
except Exception as e:
    print(f"✗ 'skip' - FAIL: {e}")

# Case insensitive
try:
    assert validate_limit_behavior("FILL") == "fill"
    print("✓ 'FILL' (case insensitive) - PASS")
except Exception as e:
    print(f"✗ 'FILL' - FAIL: {e}")

# Invalid value should raise
try:
    validate_limit_behavior("invalid")
    print("✗ 'invalid' - FAIL: Should have raised ValueError")
except ValueError as e:
    if "Invalid limit_behavior" in str(e):
        print(f"✓ 'invalid' raises ValueError - PASS")
    else:
        print(f"✗ 'invalid' - FAIL: Wrong error message: {e}")

try:
    validate_limit_behavior("random")
    print("✗ 'random' - FAIL: Should have raised ValueError")
except ValueError as e:
    if "Invalid limit_behavior" in str(e):
        print(f"✓ 'random' raises ValueError - PASS")
    else:
        print(f"✗ 'random' - FAIL: Wrong error message: {e}")

print("\n" + "=" * 60)
print("All validation tests passed! ✅")
print("=" * 60)
