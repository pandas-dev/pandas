# Issue #54627 Fix Proposal

## Problem

`groupby.var()` does not return Arrow types when input has Arrow dtype.

### Expected Behavior
When the input Series/DataFrame has Arrow dtype, the result should also have Arrow dtype.

### Actual Behavior
The result has numpy float64 dtype instead of Arrow dtype.

## Root Cause

The `_cython_agg_general` function uses cython operations which convert Arrow arrays to numpy arrays. The result is not converted back to Arrow dtype.

## Proposed Fix

### Location
`pandas/core/groupby/groupby.py` - `_cython_agg_general` method

### Approach

After the cython aggregation, check if the input had Arrow dtype and convert the result back to Arrow dtype.

### Code Changes

```python
# In _cython_agg_general, after getting the result:
def _cython_agg_general(
    self,
    how: str,
    alt: Callable | None = None,
    numeric_only: bool = False,
    min_count: int = -1,
    **kwargs,
):
    # ... existing code ...
    
    def array_func(values: ArrayLike) -> ArrayLike:
        # ... existing code ...
        return result
    
    new_mgr = data.grouped_reduce(array_func)
    res = self._wrap_agged_manager(new_mgr)
    
    # NEW: Check if input had Arrow dtype and convert result
    if how == "var":
        from pandas.core.arrays import ArrowExtensionArray
        from pandas.core.dtypes.dtypes import ArrowDtype
        
        # Check if the input data has Arrow dtype
        obj = self.obj
        if isinstance(obj, DataFrame):
            # For DataFrame, check all numeric columns
            has_arrow = any(
                isinstance(col.dtype, ArrowDtype) 
                for col in obj.dtypes 
                if numeric_only or col.kind in "iufc"
            )
        elif isinstance(obj, Series):
            has_arrow = isinstance(obj.dtype, ArrowDtype)
        else:
            has_arrow = False
        
        if has_arrow:
            # Convert result to Arrow dtype
            # This is a simplified version - the actual fix
            # would need more careful implementation
            pass  # TODO: Implement Arrow dtype conversion
    
    if how in ["idxmin", "idxmax"]:
        res = self._wrap_idxmax_idxmin(res, how=how, skipna=kwargs["skipna"])
    out = self._wrap_aggregated_output(res)
    return out
```

## Alternative Approach

A cleaner approach would be to follow the pattern used in `size()`:

```python
# In _wrap_aggregated_output or a new method
def _preserve_arrow_dtype(self, result):
    """Convert result to Arrow dtype if input had Arrow dtype."""
    from pandas.core.arrays import ArrowExtensionArray
    from pandas.core.dtypes.dtypes import ArrowDtype
    
    obj = self.obj
    
    # Check if any input column has Arrow dtype
    has_arrow = False
    if isinstance(obj, DataFrame):
        has_arrow = any(
            isinstance(dtype, ArrowDtype) 
            for dtype in obj.dtypes
        )
    elif isinstance(obj, Series):
        has_arrow = isinstance(obj.dtype, ArrowDtype)
    
    if not has_arrow:
        return result
    
    # Convert result columns to Arrow dtype
    if isinstance(result, DataFrame):
        for col in result.columns:
            if isinstance(result[col].dtype, np.floating):
                # Convert float64 to float64[pyarrow]
                result[col] = result[col].astype("float64[pyarrow]")
    
    return result
```

## Testing

Add tests in `pandas/tests/groupby/test_reductions.py` or create a new test file.

## Related Issues

- #53831: Related Arrow dtype preservation issue
- #54023: Arrow dtype in groupby operations
