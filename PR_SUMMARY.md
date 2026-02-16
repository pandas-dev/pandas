# Pandas PR #54627: groupby.var() Arrow dtype Preservation

## Summary

Working on fixing Issue #54627: `groupby.var()` does not return Arrow types when input has Arrow dtype.

## Branch
`issue-54627-arrow-var` (local only, not pushed due to fork access)

## Changes Made

### 1. Test File
`pandas/tests/groupby/test_reductions_issue_54627.py`
- Tests for Arrow dtype preservation in groupby.var()
- Tests for float64[pyarrow] and decimal[pyarrow] dtypes

### 2. Documentation
`ISSUE_54627_FIX.md`
- Root cause analysis
- Proposed fix approaches
- Alternative implementation options

## To Create PR

### Option 1: Fork and Push
```bash
# Create fork on GitHub first: https://github.com/pandas-dev/pandas/fork

# Add fork as remote
git remote add fork https://github.com/ssiweifnag/pandas.git

# Push to fork
git push fork issue-54627-arrow-var

# Create PR at: https://github.com/pandas-dev/pandas/compare/main...ssiweifnag:issue-54627-arrow-var
```

### Option 2: Work with gh CLI
```bash
# Create PR directly from branch
gh pr create --base main --head issue-54627-arrow-var \
  --title "BUG: groupby.var() does not return arrow types" \
  --body "Fixes #54627"
```

## Test Cases

### Test Case 1: Decimal Arrow Type
```python
df = DataFrame({
    "A": Series([True, True], dtype="bool[pyarrow]"),
    "B": Series([decimal.Decimal(123), decimal.Decimal(12)], 
               dtype=pd.ArrowDtype(pa.decimal128(6, 3)))
})
result = df.groupby("A").var()
# Expected: B dtype should be float64[pyarrow]
```

### Test Case 2: Float Arrow Type
```python
df = DataFrame({
    "key": Series([1, 1, 2], dtype="int64[pyarrow]"),
    "value": Series([1.0, 2.0, 3.0], dtype="float64[pyarrow]")
})
result = df.groupby("key").var()
# Expected: value dtype should be float64[pyarrow]
```

## Next Steps

1. Create GitHub fork of pandas-dev/pandas
2. Push branch to fork
3. Create PR with description
4. Address review feedback
5. Implement the actual fix (see ISSUE_54627_FIX.md)

## References
- Issue: https://github.com/pandas-dev/pandas/issues/54627
- Related: #53831, #54023
