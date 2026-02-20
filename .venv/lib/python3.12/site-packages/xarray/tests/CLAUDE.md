# Testing Guidelines for xarray

## Handling Optional Dependencies

xarray has many optional dependencies that may not be available in all testing environments. Always use the standard decorators and patterns when writing tests that require specific dependencies.

### Standard Decorators

**ALWAYS use decorators** like `@requires_dask`, `@requires_cftime`, etc. instead of conditional `if` statements.

All available decorators are defined in `xarray/tests/__init__.py` (look for `requires_*` decorators).

### DO NOT use conditional imports or skipif

❌ **WRONG - Do not do this:**

```python
def test_mean_with_cftime():
    if has_dask:  # WRONG!
        ds = ds.chunk({})
        result = ds.mean()
```

❌ **ALSO WRONG - Avoid pytest.mark.skipif in parametrize:**

```python
@pytest.mark.parametrize(
    "chunk",
    [
        pytest.param(
            True, marks=pytest.mark.skipif(not has_dask, reason="requires dask")
        ),
        False,
    ],
)
def test_something(chunk): ...
```

✅ **CORRECT - Do this instead:**

```python
def test_mean_with_cftime():
    # Test without dask
    result = ds.mean()


@requires_dask
def test_mean_with_cftime_dask():
    # Separate test for dask functionality
    ds = ds.chunk({})
    result = ds.mean()
```

✅ **OR for parametrized tests, split them:**

```python
def test_something_without_dask():
    # Test the False case
    ...


@requires_dask
def test_something_with_dask():
    # Test the True case with dask
    ...
```

### Multiple dependencies

When a test requires multiple optional dependencies:

```python
@requires_dask
@requires_scipy
def test_interpolation_with_dask(): ...
```

### Importing optional dependencies in tests

For imports within test functions, use `pytest.importorskip`:

```python
def test_cftime_functionality():
    cftime = pytest.importorskip("cftime")
    # Now use cftime
```

### Common patterns

1. **Split tests by dependency** - Don't mix optional dependency code with base functionality:

   ```python
   def test_base_functionality():
       # Core test without optional deps
       result = ds.mean()
       assert result is not None


   @requires_dask
   def test_dask_functionality():
       # Dask-specific test
       ds_chunked = ds.chunk({})
       result = ds_chunked.mean()
       assert result is not None
   ```

2. **Use fixtures for dependency-specific setup**:

   ```python
   @pytest.fixture
   def dask_array():
       pytest.importorskip("dask.array")
       import dask.array as da

       return da.from_array([1, 2, 3], chunks=2)
   ```

3. **Check available implementations**:

   ```python
   from xarray.core.duck_array_ops import available_implementations


   @pytest.mark.parametrize("implementation", available_implementations())
   def test_with_available_backends(implementation): ...
   ```

### Key Points

- CI environments intentionally exclude certain dependencies (e.g., `all-but-dask`, `bare-minimum`)
- A test failing in "all-but-dask" because it uses dask is a test bug, not a CI issue
- Look at similar existing tests for patterns to follow
