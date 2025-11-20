# Static Analysis Report

## Tool Used

**Tool:** ruff v0.14.3  
**Command:** `ruff check pandas/util/_validators.py --select ALL --ignore E501,D203,D213`  
**Target File:** `pandas/util/_validators.py`

## Key Findings

Ruff identified **88 code smells** in the target file across multiple categories:

- **EM102 (25 instances):** F-string literals in exceptions - should assign to variable first
- **TRY003 (15 instances):** Long exception messages outside exception class
- **ANN001 (23 instances):** Missing type annotations for function arguments
- **D205/D401 (12 instances):** Docstring formatting issues
- **PLC0415 (1 instance):** Import not at top-level of file
- **SIM102 (1 instance):** Nested if statements can be combined
- **FBT001/FBT002 (8 instances):** Boolean-typed positional arguments

Total: **88 code smells detected**

## Fixes Summary

### Fix #1: F-String in Exception (EM102)
**Assigned to:** Sandeep Ramavath  
**Location:** Line 292 in `pandas/util/_validators.py`  
**Issue:** Exception uses f-string literal directly instead of assigning to variable first

**Before:**
```python
def validate_bool_kwarg_nullable(value, arg_name: str) -> None:
    if (
        value is lib.no_default
        or isinstance(value, bool)
        or value is None
        or value is NA
        or (lib.is_float(value) and np.isnan(value))
    ):
        return
    raise ValueError(f"{arg_name} must be None, pd.NA, np.nan, True, or False; got {value}")
```

**After:**
```python
def validate_bool_kwarg_nullable(value, arg_name: str) -> None:
    if (
        value is lib.no_default
        or isinstance(value, bool)
        or value is None
        or value is NA
        or (lib.is_float(value) and np.isnan(value))
    ):
        return
    msg = f"{arg_name} must be None, pd.NA, np.nan, True, or False; got {value}"
    raise ValueError(msg)
```

**Rationale:** Assigning error messages to variables before raising exceptions improves code maintainability, makes testing easier, and follows Python best practices for exception handling.

---

### Fix #2: Import Not at Top-Level (PLC0415)
**Assigned to:** Nithikesh Bobbili  
**Location:** Line 314 in `pandas/util/_validators.py`  
**Issue:** Import statement inside function body instead of at module level

**Before:**
```python
def validate_fillna_kwargs(value, method, validate_scalar_dict_value: bool = True):
    """Validate the keyword arguments to 'fillna'.
    ...
    """
    from pandas.core.missing import clean_fill_method  # Import inside function

    if value is None and method is None:
        raise ValueError("Must specify a fill 'value' or 'method'.")
    if value is None and method is not None:
        method = clean_fill_method(method)
    # ...
```

**After:**
```python
# At module level (top of file, after line 20):
from pandas.core.missing import clean_fill_method

# ...

def validate_fillna_kwargs(value, method, validate_scalar_dict_value: bool = True):
    """Validate the keyword arguments to 'fillna'.
    ...
    """
    if value is None and method is None:
        raise ValueError("Must specify a fill 'value' or 'method'.")
    if value is None and method is not None:
        method = clean_fill_method(method)
    # ...
```

**Rationale:** Module-level imports are executed once at module load time rather than every function call, improving performance. It also makes dependencies more visible and follows PEP 8 conventions.

---


### Fix #3: Nested If Statements (SIM102)
**Assigned to:** Mallikarjuna  
**Location:** Lines 471-472 in `pandas/util/_validators.py`  
**Issue:** Unnecessary nested if statements can be combined with `and`

**Before:**
```python
def check_dtype_backend(dtype_backend) -> None:
    if dtype_backend is not lib.no_default:
        if dtype_backend not in ["numpy_nullable", "pyarrow"]:
            raise ValueError(
                f"dtype_backend {dtype_backend} is invalid, only 'numpy_nullable' and "
                f"'pyarrow' are allowed.",
            )
```

**After:**
```python
def check_dtype_backend(dtype_backend) -> None:
    if dtype_backend is not lib.no_default and dtype_backend not in ["numpy_nullable", "pyarrow"]:
        raise ValueError(
            f"dtype_backend {dtype_backend} is invalid, only 'numpy_nullable' and "
            f"'pyarrow' are allowed.",
        )
```

**Rationale:** Combining related conditions into a single if statement reduces nesting depth, improves readability, and makes the code more concise without losing clarity.

---

## Group Contributions

**Sandeep Ramavath:**
- Identified and fixed EM102 code smell (f-string in exception)
- Refactored exception handling to use variable assignment
- Impact: Improved exception handling best practices

**Nithikesh Bobbili:**
- Identified and fixed PLC0415 code smell (import location)
- Moved import statement to module level
- Impact: Better performance and code organization

**Mallikarjuna:**
- Identified and fixed SIM102 code smell (nested if statements)
- Simplified conditional logic by combining conditions
- Impact: Reduced code complexity and improved readability