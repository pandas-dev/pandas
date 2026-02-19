# Next Contribution Recommendations - Similar to Issue #64208

## 🎯 **PRIMARY RECOMMENDATION: Issue #64196**

**Issue Title:** BUG: For `infer_dtype` `skipna` is ignored for Period / Interval  
**Issue URL:** https://github.com/pandas-dev/pandas/issues/64196  
**Labels:** Bug, Needs Triage  
**Status:** Open, Unassigned  
**Created:** 2026-02-17  
**Comments:** 1

### Why This is an Excellent Choice:

1. ✅ **Recent Issue** - Created Feb 17, 2026 (very fresh)
2. ✅ **Unassigned** - Available for claiming
3. ✅ **Clear Bug Fix** - Not documentation, pure code fix
4. ✅ **Well-Documented** - Issue includes clear reproduction steps and expected behavior
5. ✅ **Straightforward Fix** - Issue description mentions "The wiring itself is straightforward"
6. ✅ **Similar Complexity** - Comparable to #64208 (simple parameter wiring)
7. ✅ **High Merge Probability** - Clear bug with obvious fix path

### Issue Description Summary:

The `pandas.api.types.infer_dtype` function accepts a `skipna` argument, but for `Period` and `Interval` inputs it is currently not applied:
- Interval inference does not pass `skipna` into its validation helper
- Period inference contains a `FIXME` noting that `skipna` is not actually used

### Expected Impact:

- **Functionality:** Fixes missing parameter behavior for Period/Interval dtype inference
- **API Consistency:** Makes `skipna` parameter work consistently across all dtype inference
- **User Impact:** Users can now properly control missing value handling in dtype inference

### Implementation Approach:

1. **Files to Modify:**
   - `pandas/core/dtypes/inference.py` (likely location for `infer_dtype`)
   - Possibly `pandas/core/arrays/period.py` and `pandas/core/arrays/interval.py`
   - `pandas/tests/dtypes/test_inference.py` (add test cases)

2. **Changes Required:**
   - Wire `skipna` parameter through to Period inference logic
   - Wire `skipna` parameter through to Interval inference logic
   - Remove/update FIXME comment
   - Add tests to verify `skipna=True` vs `skipna=False` behavior

### Files/Modules Requiring Modification:

- `pandas/core/dtypes/inference.py` (main fix)
- `pandas/core/arrays/period.py` (possibly)
- `pandas/core/arrays/interval.py` (possibly)
- `pandas/tests/dtypes/test_inference.py` (add tests)

### Estimated Complexity: **LOW-MEDIUM**

**Similarity to #64208:**
- ✅ Simple parameter wiring (similar to fixing test assertion)
- ✅ Clear fix path
- ✅ Well-documented issue
- ✅ Low risk of breaking changes

---

## 🥈 **SECONDARY RECOMMENDATION: Issue #64184**

**Issue Title:** BUG: Inconsistent number conversion with leading zeros in to_numeric()  
**Issue URL:** https://github.com/pandas-dev/pandas/issues/64184  
**Labels:** Bug  
**Status:** Open, Unassigned  
**Created:** 2026-02-17  
**Comments:** 0

### Why This is a Good Choice:

1. ✅ **Recent Issue** - Created Feb 17, 2026
2. ✅ **Unassigned** - Available for work
3. ✅ **Clear Bug** - Inconsistent behavior in `to_numeric()`
4. ✅ **Well-Documented** - Includes reproducible examples
5. ⚠️ **Medium Complexity** - May require understanding dtype inference logic

### Issue Description Summary:

`pd.to_numeric()` has inconsistent behavior when converting numbers with leading zeros. The dtype chosen (int64 vs float64) depends on other values in the Series, causing precision issues with large numbers containing leading zeros.

### Expected Impact:

- **Functionality:** Fixes inconsistent dtype inference in `to_numeric()`
- **User Impact:** Prevents precision loss when converting numbers with leading zeros

### Implementation Approach:

1. **Files to Modify:**
   - `pandas/core/tools/numeric.py` (likely location)
   - `pandas/tests/tools/test_to_numeric.py` (add test cases)

### Estimated Complexity: **MEDIUM**

---

## 🥉 **TERTIARY RECOMMENDATION: Issue #31142**

**Issue Title:** Series.combine with fill_value gives unexpected results  
**Issue URL:** https://github.com/pandas-dev/pandas/issues/31142  
**Labels:** Bug, good first issue, Needs Tests, combine/combine_first/update  
**Status:** Open, Unassigned  
**Created:** 2020-01-19  
**Comments:** 9

### Why This is a Good Choice:

1. ✅ **Good First Issue Label** - Maintainer approved for new contributors
2. ✅ **Unassigned** - Available for work
3. ✅ **Clear Bug** - Unexpected behavior in `Series.combine()`
4. ✅ **Well-Documented** - Includes code samples and expected output
5. ⚠️ **Older Issue** - Created in 2020, may need verification it still exists

### Issue Description Summary:

`Series.combine()` with `fill_value` gives unexpected results when index contains both string and integer labels. The issue occurs when index `0` (integer) conflicts with string index labels.

### Expected Impact:

- **Functionality:** Fixes unexpected behavior in `Series.combine()`
- **User Impact:** Correct handling of mixed-type indices

### Estimated Complexity: **MEDIUM**

---

## Other Potential Issues

### Issue #64150
- **Title:** BUG: DataFrame.unstack(0, sort=False) creates fake column names
- **Complexity:** Medium
- **Note:** Has 10 comments, may have discussion about approach

### Issue #54627
- **Title:** BUG: groupby.var() does not return arrow types with arrow backed series as input
- **Complexity:** Medium
- **Labels:** Bug, good first issue, Arrow
- **Note:** Related to Arrow dtype support

### Issue #55136
- **Title:** BUG: Conversion from datetime64[ns] to datetime does not effect .info() and probably not the DataFrame itself
- **Complexity:** Medium
- **Labels:** Bug, good first issue, Datetime
- **Note:** Older issue (2023), may need investigation

---

## Recommendation Summary

**Best Choice: Issue #64196**
- Most similar to #64208 in complexity
- Recent and unassigned
- Clear fix path
- Issue author mentions it's "straightforward"
- High likelihood of merge

**Next Steps:**
1. Review issue #64196: https://github.com/pandas-dev/pandas/issues/64196
2. Claim the issue by commenting `take`
3. Investigate the code in `pandas/core/dtypes/inference.py`
4. Implement the fix by wiring `skipna` parameter through
5. Add tests to verify the fix
6. Submit PR

---

*Analysis generated: 2026-02-18*


