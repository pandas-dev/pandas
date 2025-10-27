# Unit Testing Extension - README

## Overview

This document explains how to run the 15 additional unit test cases added to the pandas codebase as part of our SWEN 777 course project. These tests target critical edge cases and uncovered logic paths across three key pandas modules.

**Project Details:**
- **Course:** SWEN 777 - Software Quality Assurance
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
- **Total Tests**: 15
- **Tests Passed**: 15
- **Tests Failed**: 0
- **Success Rate**: 100%
- **Execution Time**: ~1.04 seconds

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
pandas/tests/test_series_constructors_additional.py::test_series_constructor_invalid_key_types PASSED
pandas/tests/test_series_constructors_additional.py::test_series_constructor_empty_edge_cases PASSED
pandas/tests/test_series_constructors_additional.py::test_series_constructor_mixed_dtype_edge_cases PASSED
pandas/tests/test_series_constructors_additional.py::test_series_constructor_memory_intensive PASSED
pandas/tests/test_series_constructors_additional.py::test_series_constructor_invalid_index_length PASSED
pandas/tests/tseries/offsets/test_offsets.py::test_dateoffset_boundary_values PASSED
pandas/tests/tseries/offsets/test_offsets.py::test_business_day_weekend_edge_cases PASSED
pandas/tests/tseries/offsets/test_offsets.py::test_custom_business_hour_edge_cases PASSED
pandas/tests/tseries/offsets/test_offsets.py::test_quarter_offset_leap_year PASSED
pandas/tests/tseries/offsets/test_offsets.py::test_offset_frequency_string_edge_cases PASSED

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
- **HTML Report**: `courseProjectDocs/Setup/htmlcov/index.html`
- **Terminal Output**: Displayed during test execution
- **Expected Coverage**: 11% overall (improvement from ~10% baseline)

## Unit Testing II - Mocking & Stubbing Tests

### Additional Mocking Tests (NEW)

In addition to the original 15 unit tests above, we have added 15 mocking-based tests for Unit Testing II assignment:

**New Test Files:**
1. **`pandas/tests/mocking/test_database_io.py`** - 5 tests for database I/O operations
2. **`pandas/tests/mocking/test_filesystem_io.py`** - 5 tests for file system I/O operations
3. **`pandas/tests/mocking/test_datetime.py`** - 5 tests for datetime/time-series operations

### Prerequisites for Mocking Tests
Before running the mocking tests, ensure you have:
- All prerequisites from above
- **pytest-mock 3.15.1+** (NEW REQUIREMENT)

```bash
# Install pytest-mock if not already installed
pip install pytest-mock
```

### Running Mocking Tests

#### Run All 15 Mocking Tests
```bash
# Run all mocking tests with verbose output
pytest pandas/tests/mocking/ -v
```

**Expected Output:**
```
============================= test session starts ==============================
platform darwin -- Python 3.13.5, pytest-8.4.2, pluggy-1.6.0
collected 15 items

pandas/tests/mocking/test_database_io.py::TestDatabaseIOMocking::test_read_sql_basic PASSED
pandas/tests/mocking/test_database_io.py::TestDatabaseIOMocking::test_read_sql_empty_result PASSED
pandas/tests/mocking/test_database_io.py::TestDatabaseIOMocking::test_read_sql_with_parameters PASSED
pandas/tests/mocking/test_database_io.py::TestDatabaseIOMocking::test_read_sql_dtype_handling PASSED
pandas/tests/mocking/test_database_io.py::TestDatabaseIOMocking::test_read_sql_connection_error_handling PASSED
pandas/tests/mocking/test_datetime.py::TestDateTimeOperationsMocking::test_timestamp_now_mocked PASSED
pandas/tests/mocking/test_datetime.py::TestDateTimeOperationsMocking::test_date_range_generation PASSED
pandas/tests/mocking/test_datetime.py::TestDateTimeOperationsMocking::test_time_series_resampling PASSED
pandas/tests/mocking/test_datetime.py::TestDateTimeOperationsMocking::test_rolling_window_operations PASSED
pandas/tests/mocking/test_datetime.py::TestDateTimeOperationsMocking::test_datetime_parsing_with_format PASSED
pandas/tests/mocking/test_filesystem_io.py::TestFileSystemIOMocking::test_read_csv_basic PASSED
pandas/tests/mocking/test_filesystem_io.py::TestFileSystemIOMocking::test_read_csv_with_delimiter PASSED
pandas/tests/mocking/test_filesystem_io.py::TestFileSystemIOMocking::test_read_excel_basic PASSED
pandas/tests/mocking/test_filesystem_io.py::TestFileSystemIOMocking::test_read_hdf_basic PASSED
pandas/tests/mocking/test_filesystem_io.py::TestFileSystemIOMocking::test_csv_file_not_found_handling PASSED

============================== 15 passed in 0.83s ==============================
```

#### Run Mocking Tests by Category
```bash
# Database I/O tests
pytest pandas/tests/mocking/test_database_io.py -v

# File System I/O tests
pytest pandas/tests/mocking/test_filesystem_io.py -v

# DateTime operations tests
pytest pandas/tests/mocking/test_datetime.py -v
```

#### Generate Mocking Test Coverage Report
```bash
# Generate coverage report for mocking test code
pytest pandas/tests/mocking/ --cov=pandas/tests/mocking --cov-report=term

# Expected output:
# Name                                         Stmts   Miss Branch BrPart  Cover
# --------------------------------------------------------------------------------
# pandas/tests/mocking/__init__.py                 0      0      0      0   100%
# pandas/tests/mocking/test_database_io.py        44      0      0      0   100%
# pandas/tests/mocking/test_datetime.py           70     10      8      3    81%
# pandas/tests/mocking/test_filesystem_io.py      45      1      2      1    96%
# --------------------------------------------------------------------------------
# TOTAL                                          159     11     10      4    90%
```

### Mocking Test Results Summary
- **Total Mocking Tests**: 15
- **Passed**: 15 (100%)
- **Failed**: 0
- **Execution Time**: 0.83 seconds
- **Test Code Coverage**: 90%

For detailed mocking strategy and design decisions, see: `courseProjectDocs/Unit-Testing/mocking.md`


## Project Team Information

**Course:** SWEN 777 - Software Testing and Quality Assurance  
**Team Members:**
- Nithikesh Reddy
- Sandeep
- Malikarjuna

