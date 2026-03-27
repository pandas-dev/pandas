# Implementation Summary: interpolate(limit_behavior="skip") Enhancement

## Overview
Added optional `limit_behavior` parameter to `Series.interpolate()` to allow skipping NaN gaps that exceed the specified limit, instead of partially filling them.

## Files Modified

### 1. `pandas/core/generic.py` - Series/DataFrame interpolate() method
**Changes:**
- Added `limit_behavior: Literal["fill", "skip"] = "fill"` parameter to function signature (line ~7907)
- Added docstring documentation for the new parameter (default='fill', versionadded=2.2.0)
- Pass `limit_behavior` to `obj._mgr.interpolate()` call (line ~8099)

**Behavior:**
- `limit_behavior="fill"` (default): Current behavior - fill up to `limit` consecutive NaNs
- `limit_behavior="skip"`: New behavior - if any consecutive NaN gap > limit, skip entire gap

### 2. `pandas/core/missing.py` - Core interpolation logic
**Changes:**

#### A. Updated `interpolate_2d_inplace()` function (line ~364)
- Added `limit_behavior: str = "fill"` parameter to signature
- Pass `limit_behavior` to `_interpolate_1d()` call within func()

#### B. Updated `_interpolate_1d()` function (line ~450)
- Added `limit_behavior: str = "fill"` parameter to signature
- Added gap-reverting logic after interpolation (lines ~545-569):
  ```python
  # GH#NEW_ISSUE: If limit_behavior="skip", revert gaps exceeding limit
  if limit_behavior == "skip" and limit is not None:
      # Find consecutive NaN groups that exceed the limit
      nan_indices = np.where(invalid)[0]
      if len(nan_indices) > 0:
          # Detect gap boundaries by finding non-consecutive NaN indices
          gap_ends = np.where(np.diff(nan_indices) > 1)[0]
          gap_starts = np.concatenate(([0], gap_ends + 1))
          gap_ends = np.concatenate((gap_ends, [len(nan_indices) - 1]))

          revert_indices = []
          for start, end in zip(gap_starts, gap_ends):
              gap_size = end - start + 1
              if gap_size > limit:
                  # Gap exceeds limit; mark these NaN indices for reverting
                  revert_indices.extend(nan_indices[start : end + 1])

          if revert_indices:
              revert_indices = np.array(revert_indices, dtype=np.int64)
              preserve_nans = np.union1d(preserve_nans, revert_indices)
  ```

### 3. `pandas/tests/series/methods/test_interpolate.py` - Test cases
**Added new test class `TestInterpolateLimitBehavior`** with 5 test methods:
1. `test_interpolate_limit_behavior_skip_basic` - Gap of 3 NaNs with limit=1 should skip all
2. `test_interpolate_limit_behavior_fill_default` - Default behavior unchanged (fill up to limit)
3. `test_interpolate_limit_behavior_skip_multiple_gaps` - Skip only gaps exceeding limit
4. `test_interpolate_limit_behavior_skip_exact_limit` - Gap exactly equal to limit should be filled
5. `test_interpolate_limit_behavior_skip_forward_direction` - Works with limit_direction

## Example Usage

```python
import pandas as pd
import numpy as np

# Test case 1: Gap too large, skip entire gap
s = pd.Series([1.0, np.nan, np.nan, np.nan, 5.0])
result = s.interpolate(limit=1, limit_behavior="skip")
# Expected: [1.0, NaN, NaN, NaN, 5.0]

# Test case 2: Multiple gaps, skip only those exceeding limit
s = pd.Series([1.0, np.nan, 3.0, np.nan, np.nan, np.nan, 7.0])
result = s.interpolate(limit=2, limit_behavior="skip")
# Expected: [1.0, 2.0, 3.0, NaN, NaN, NaN, 7.0]
# Gap [1] has size 1 <= 2 -> filled
# Gap [3,4,5] has size 3 > 2 -> skipped

# Default behavior unchanged
s = pd.Series([1.0, np.nan, np.nan, np.nan, 5.0])
result = s.interpolate(limit=1)  # or limit_behavior="fill" (same)
# Expected: [1.0, 2.0, NaN, NaN, 5.0]
```

## Algorithm

The implementation works in two phases:

### Phase 1: Normal Interpolation
- Perform standard interpolation as usual with the given limit and limit_direction
- This fills NaN values according to the existing logic

### Phase 2: Skip Large Gaps (if limit_behavior="skip")
- Detect all consecutive NaN groups by finding indices where NaNs are not adjacent
- For each gap, check if its size exceeds the limit
- If exceeded, add all indices in that gap to the `preserve_nans` array
- Finally, set all preserved_nans back to NaN, reverting the interpolation for those gaps

## Design Decisions

✅ **Minimal and Localized:** Only 3 files modified, changes are focused
✅ **Backward Compatible:** Default behavior (fill) is unchanged
✅ **Simple Logic:** Straightforward gap detection without new helper functions
✅ **Consistent with pandas:** Follows existing parameter naming and logic patterns
✅ **Works with existing features:** Compatible with limit_direction, limit_area, etc.

## Testing

Logic verification passed for:
- Single large gaps
- Multiple gaps with selective skipping
- Gaps exactly equal to limit (should not skip)
- Edge cases

Test class added with comprehensive coverage.
