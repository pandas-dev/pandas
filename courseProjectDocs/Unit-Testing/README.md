# Unit Testing Extension - README

## Overview

This document explains how to run the 15 additional unit test cases added to the pandas codebase as part of our SWEN 777 course project. These tests target critical edge cases and uncovered logic paths across three key pandas modules.

**Project Details:**
- **Course:** SWEN 777 - Software Architecture
- **Group Size:** 3 members  
- **Deliverable:** 15 meaningful unit test cases (5 per student)
- **Target:** Increase test coverage for uncovered or edge-case logic

## Test Files Overview

Our group added tests to three separate files to avoid interfering with the baseline test suite:

1. **`pandas/tests/test_nanops_additional.py`** (NEW FILE)
   - 5 tests for numerical operations edge cases
   - Tests empty arrays, mask scenarios, boundary conditions

2. **`pandas/tests/test_series_constructors_additional.py`** (NEW FILE)  
   - 5 tests for Series object creation edge cases
   - Tests invalid inputs, empty data, dtype inference

3. **`pandas/tests/tseries/offsets/test_offsets.py`** (MODIFIED FILE)
   - 5 tests for datetime offset edge cases (added to existing file)
   - Tests boundary timestamps, business logic, leap years

## How to Reproduce Test Results

### Prerequisites

Before running the tests, ensure you have:
- Python 3.13+ installed
- Virtual environment activated
- pandas development version (3.0.0.dev0+)
- pytest 8.4.2+ with pytest-cov 7.0.0+

### Step-by-Step Instructions

#### 1. Environment Setup
```bash
# Navigate to project directory
cd /Volumes/T7Shield/SWEN777/SWEN_777_Pandas

# Activate virtual environment
source venv/bin/activate
```

#### 2. Run All 15 Added Tests
```bash
python -m pytest \
    pandas/tests/test_nanops_additional.py \
    pandas/tests/test_series_constructors_additional.py \
    pandas/tests/tseries/offsets/test_offsets.py::test_dateoffset_boundary_values \
    pandas/tests/tseries/offsets/test_offsets.py::test_business_day_weekend_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_custom_business_hour_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_quarter_offset_leap_year \
    pandas/tests/tseries/offsets/test_offsets.py::test_offset_frequency_string_edge_cases \
    -v
```

#### 3. Generate Coverage Report
```bash
python -m pytest \
    pandas/tests/test_nanops_additional.py \
    pandas/tests/test_series_constructors_additional.py \
    pandas/tests/tseries/offsets/test_offsets.py::test_dateoffset_boundary_values \
    pandas/tests/tseries/offsets/test_offsets.py::test_business_day_weekend_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_custom_business_hour_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_quarter_offset_leap_year \
    pandas/tests/tseries/offsets/test_offsets.py::test_offset_frequency_string_edge_cases \
    --cov=pandas --cov-report=html:courseProjectDocs/Setup/htmlcov --cov-report=term
```

#### 4. Run Tests by Category (Optional)

**Nanops Tests Only:**
```bash
python -m pytest pandas/tests/test_nanops_additional.py -v
```

**Series Constructor Tests Only:**
```bash
python -m pytest pandas/tests/test_series_constructors_additional.py -v
```

**DateTime Offset Tests Only:**
```bash
python -m pytest \
    pandas/tests/tseries/offsets/test_offsets.py::test_dateoffset_boundary_values \
    pandas/tests/tseries/offsets/test_offsets.py::test_business_day_weekend_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_custom_business_hour_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_quarter_offset_leap_year \
    pandas/tests/tseries/offsets/test_offsets.py::test_offset_frequency_string_edge_cases \
    -v
```

## Expected Test Results

When you run the tests, you should see:

- **Total Tests:** 15
- **Tests Passed:** 15 
- **Tests Failed:** 0
- **Success Rate:** 100%
- **Execution Time:** ~1.04 seconds

**Sample Output:**
```
============================= test session starts ==============================
platform darwin -- Python 3.13.5, pytest-8.4.2, pluggy-1.6.0
collected 15 items

pandas/tests/test_nanops_additional.py::test_nansum_empty_array_edge_cases PASSED
pandas/tests/test_nanops_additional.py::test_nanmean_mask_edge_cases PASSED
pandas/tests/test_nanops_additional.py::test_nanvar_ddof_boundary_conditions PASSED
pandas/tests/test_nanops_additional.py::test_nanargmax_nanargmin_error_conditions PASSED
pandas/tests/test_nanops_additional.py::test_nanskew_nankurt_insufficient_samples PASSED
pandas/tests/test_series_constructors_additional.py::test_series_constructor_invalid_key_types PASSED    [  7%]
pandas/tests/test_series_constructors_additional.py::test_series_constructor_empty_edge_cases PASSED    [  8%]
pandas/tests/test_series_constructors_additional.py::test_series_constructor_mixed_dtype_edge_cases PASSED    [  9%]
pandas/tests/test_series_constructors_additional.py::test_series_constructor_memory_intensive PASSED    [ 10%]
pandas/tests/test_series_constructors_additional.py::test_series_constructor_invalid_index_length PASSED    [ 11%]
pandas/tests/tseries/offsets/test_offsets.py::test_dateoffset_boundary_values PASSED    [ 12%]
pandas/tests/tseries/offsets/test_offsets.py::test_business_day_weekend_edge_cases PASSED    [ 13%]
pandas/tests/tseries/offsets/test_offsets.py::test_custom_business_hour_edge_cases PASSED    [ 14%]
pandas/tests/tseries/offsets/test_offsets.py::test_quarter_offset_leap_year PASSED    [ 15%]
pandas/tests/tseries/offsets/test_offsets.py::test_offset_frequency_string_edge_cases PASSED    [ 16%]

============================== 15 passed in 1.04s ==============================
```

## Coverage Analysis

### Comprehensive Coverage Command
To run both baseline and additional tests for complete coverage analysis:

```bash
python -m pytest \
    pandas/tests/series/test_constructors.py \
    pandas/tests/frame/test_constructors.py \
    pandas/tests/test_nanops.py \
    pandas/tests/series/methods/test_dropna.py \
    pandas/tests/frame/methods/test_dropna.py \
    pandas/tests/test_nanops_additional.py \
    pandas/tests/test_series_constructors_additional.py \
    pandas/tests/tseries/offsets/test_offsets.py::test_dateoffset_boundary_values \
    pandas/tests/tseries/offsets/test_offsets.py::test_business_day_weekend_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_custom_business_hour_edge_cases \
    pandas/tests/tseries/offsets/test_offsets.py::test_quarter_offset_leap_year \
    pandas/tests/tseries/offsets/test_offsets.py::test_offset_frequency_string_edge_cases \
    --cov=pandas --cov-report=html:courseProjectDocs/Setup/htmlcov --cov-report=term
```

### Coverage Report Location
- **HTML Report:** `courseProjectDocs/Setup/htmlcov/index.html`
- **Terminal Output:** Displayed during test execution
- **Expected Coverage:** 11% overall (improvement from ~10% baseline)

## Troubleshooting

### Common Issues and Solutions

1. **Environment Setup**
   - Ensure virtual environment is activated
   - Verify Python 3.13+ installation
   - Check pandas development build installation

2. **Test Execution Problems**
   - Clear pytest cache: `python -m pytest --cache-clear`
   - Run tests individually if batch execution fails
   - Check for import conflicts

3. **Coverage Report Issues**
   - Ensure output directory exists: `mkdir -p courseProjectDocs/Setup/htmlcov`
   - Run with verbose output: `--cov-report=term-missing`



## Project Team Information

**Course:** SWEN 777 - Software Testing and Quality Assurance  
**Project:** Pandas Unit Testing Extension  
**Team Members:**
- Nithikesh Reddy
- Sandeep
- Malikarjuna


## Results

- **Test Execution:** All 15 tests should pass
- **Coverage Improvement:** From ~10% to 11% overall coverage
- **New Code Coverage:** 100% coverage for added test functions
- pytest-cov 7.0.0+ (for coverage analysis)
- numpy
- Virtual environment recommended

## Setting Up Test Environment

1. Create and activate a virtual environment
2. Install development dependencies from requirements-dev.txt
3. Build pandas in development mode

## Coverage Analysis

To analyze coverage improvements from these tests, use pytest with coverage flags targeting the specific modules (pandas.core.nanops, pandas.core.series, pandas.tseries.offsets) and generate both HTML and terminal reports.

## Test Design Principles
All added tests follow these principles:
1. **Edge Case Focus:** Target boundary conditions and unusual inputs
2. **Error Handling:** Test exception conditions and error paths
3. **Uncovered Logic:** Address gaps identified in coverage analysis
4. **Maintainability:** Clear naming and comprehensive documentation
5. **Integration:** Seamlessly integrate with existing test structure

## Files Modified

1. `pandas/tests/test_nanops.py` - Added 5 test functions (lines ~1280-1340)
2. `pandas/tests/series/test_constructors.py` - Added 5 test functions (lines ~890-970)
3. `pandas/tests/tseries/offsets/test_offsets.py` - Added 5 test functions (lines ~1235-1310)

## Group Members

- Member 1: Nanops module test cases (5 tests)
- Member 2: Series constructor test cases (5 tests)  
- Member 3: DateTime offset test cases (5 tests)

## Success Criteria
All 15 test cases pass successfully  
Tests cover edge cases and boundary conditions  
Tests integrate with existing pandas test suite  
Comprehensive documentation provided  
Test cases target previously uncovered code paths  

## Next Steps

For further test development:
1. Monitor coverage reports to identify additional gaps
2. Consider adding integration tests for cross-module functionality
3. Expand boundary condition testing for other pandas modules
4. Add performance benchmarks for edge case scenarios