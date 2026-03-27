#!/usr/bin/env python
"""Quick test to verify limit_behavior logic"""
import numpy as np

# Simulate the gap detection logic from our implementation
def check_limit_behavior_skip_logic():
    """Test the gap detection logic"""
    
    # Test case 1: Single large gap
    invalid = np.array([False, True, True, True, False])
    nan_indices = np.where(invalid)[0]
    gap_ends = np.where(np.diff(nan_indices) > 1)[0]
    gap_starts = np.concatenate(([0], gap_ends + 1))
    gap_ends = np.concatenate((gap_ends, [len(nan_indices) - 1]))
    
    print("Test 1: Single large gap [1,2,3]")
    print(f"  invalid: {invalid}")
    print(f"  nan_indices: {nan_indices}")
    print(f"  gap_starts: {gap_starts}, gap_ends: {gap_ends}")
    
    revert_indices = []
    limit = 1
    for start, end in zip(gap_starts, gap_ends):
        gap_size = end - start + 1
        print(f"  Gap: indices {nan_indices[start:end+1]}, size={gap_size}, exceeds_limit={gap_size > limit}")
        if gap_size > limit:
            revert_indices.extend(nan_indices[start : end + 1])
    
    print(f"  Revert indices: {revert_indices}")
    assert revert_indices == [1, 2, 3], f"Expected [1, 2, 3], got {revert_indices}"
    print("  ✓ PASS\n")
    
    # Test case 2: Multiple gaps
    invalid = np.array([False, True, False, True, True, True, False])
    nan_indices = np.where(invalid)[0]
    gap_ends = np.where(np.diff(nan_indices) > 1)[0]
    gap_starts = np.concatenate(([0], gap_ends + 1))
    gap_ends = np.concatenate((gap_ends, [len(nan_indices) - 1]))
    
    print("Test 2: Multiple gaps [1] and [3,4,5]")
    print(f"  invalid: {invalid}")
    print(f"  nan_indices: {nan_indices}")
    print(f"  gap_starts: {gap_starts}, gap_ends: {gap_ends}")
    
    revert_indices = []
    limit = 2
    for start, end in zip(gap_starts, gap_ends):
        gap_size = end - start + 1
        print(f"  Gap: indices {nan_indices[start:end+1]}, size={gap_size}, exceeds_limit={gap_size > limit}")
        if gap_size > limit:
            revert_indices.extend(nan_indices[start : end + 1])
    
    print(f"  Revert indices: {revert_indices}")
    assert revert_indices == [3, 4, 5], f"Expected [3, 4, 5], got {revert_indices}"
    print("  ✓ PASS\n")
    
    # Test case 3: Gap exactly equal to limit
    invalid = np.array([False, True, True, False])
    nan_indices = np.where(invalid)[0]
    gap_ends = np.where(np.diff(nan_indices) > 1)[0]
    gap_starts = np.concatenate(([0], gap_ends + 1))
    gap_ends = np.concatenate((gap_ends, [len(nan_indices) - 1]))
    
    print("Test 3: Gap exactly equal to limit [1,2]")
    print(f"  invalid: {invalid}")
    print(f"  nan_indices: {nan_indices}")
    
    revert_indices = []
    limit = 2
    for start, end in zip(gap_starts, gap_ends):
        gap_size = end - start + 1
        print(f"  Gap: indices {nan_indices[start:end+1]}, size={gap_size}, exceeds_limit={gap_size > limit}")
        if gap_size > limit:
            revert_indices.extend(nan_indices[start : end + 1])
    
    print(f"  Revert indices: {revert_indices}")
    assert revert_indices == [], f"Expected [], got {revert_indices}"
    print("  ✓ PASS\n")
    
    print("✅ All logic tests passed!")

if __name__ == "__main__":
    check_limit_behavior_skip_logic()
