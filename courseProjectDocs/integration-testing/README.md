# Integration Testing - Instructions to Run Tests

## Test File Location

Integration tests are located in: **`pandas/tests/util/test_integration.py`**

## Prerequisites

```bash
# Navigate to project directory
cd /Volumes/T7Shield/SWEN777/SWEN_777_Pandas

# Activate virtual environment
source venv/bin/activate
```

## How to Run Tests to Reproduce Results

### Run All Integration Tests

```bash
python -m pytest pandas/tests/util/test_integration.py -v
```

**Expected Output:**
```
collected 6 items

pandas/tests/util/test_integration.py::TestSandeepIntegration::test_series_to_dataframe_dtype_preservation PASSED
pandas/tests/util/test_integration.py::TestSandeepIntegration::test_dataframe_from_dict_mixed_series_dtypes PASSED
pandas/tests/util/test_integration.
py::TestNithikeshIntegration::test_validate_fillna_with_clean_method PASSED
pandas/tests/util/test_integration.py::TestNithikeshIntegration::test_series_fillna_integration PASSED
pandas/tests/util/test_integration.