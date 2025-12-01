# Integration Testing Report

## Test Design Summary

This integration testing effort focuses on verifying interactions between multiple pandas modules. The tests are organized into three areas, each covering at least 2 module interactions:

### Test 1: Series-DataFrame-Dtype Integration (Sandeep Ramavath)
**Modules Integrated:**
- `pandas.core.series` (Series class)
- `pandas.core.frame` (DataFrame class)
- `pandas.core.internals` (internal data managers)
- `pandas.core.dtypes` (data type handling)

**Interactions Tested:**
1. **Series.to_frame()**: Tests dtype preservation when converting Series to DataFrame through internal manager conversion
2. **DataFrame construction from dict**: Tests how DataFrame handles multiple Series with different dtypes (int64, float32, object) during construction


### Test 2: Validation-Missing Data Integration (Nithikesh Bobbili)
**Modules Integrated:**
- `pandas.util._validators` (validation utilities)
- `pandas.core.missing` (missing data handling)
- `pandas.core.series` (Series operations)
- `pandas.core.internals` (internal data modification)

**Interactions Tested:**
1. **validate_fillna_kwargs with clean_fill_method**: Tests delegation from validator to missing data module for method normalization
2. **Series.fillna/ffill operations**: Tests complete pipeline from user API through validation to missing data handling

### Test 3: Dtype Backend-Libs Integration (Mallikarjuna)
**Modules Integrated:**
- `pandas.util._validators` (validation functions)
- `pandas._libs.lib` (C extension library with sentinel values)
- `numpy` (array handling and validation)

**Interactions Tested:**
1. **check_dtype_backend with lib.no_default**: Tests validator interaction with C library sentinel values
2. **validate_percentile with numpy arrays**: Tests pandas validation with numpy array conversion and bounds checking

## Test Data Preparation

### Input Data Generation

**Test 1 - Series/DataFrame Integration:**
- **Input**: Created Series with explicit dtype (`int32`) and sample data `[1, 2, 3]`
- **Input**: Created multiple Series with different dtypes: int64, float32, object
- **Rationale**: Different dtypes exercise type preservation logic across module boundaries

**Test 2 - Validation/Missing Data:**
- **Input**: Series with `np.nan` values: `[1.0, np.nan, 3.0, np.nan, 5.0]`
- **Input**: Method names `"pad"`, `"ffill"` and `None` values
- **Rationale**: Missing values and various method names test validation and fill method delegation

**Test 3 - Backend/Libs Validation:**
- **Input**: `lib.no_default` sentinel, valid backends (`"numpy_nullable"`, `"pyarrow"`), invalid backend string
- **Input**: Valid percentiles (`0.5`, `[0.25, 0.5, 0.75]`) and invalid (`1.5`, `[0.25, 1.5, 0.75]`)
- **Rationale**: Mix of valid/invalid inputs tests error handling across module boundaries

### Expected Output Data

All tests include explicit expected outputs:
- Series/DataFrame tests verify dtype preservation and data integrity
- Validation tests verify normalized method names and appropriate ValueError exceptions
- Backend tests verify acceptance of valid values and rejection with specific error messages

## Execution and Results

**Test File**: `pandas/tests/util/test_integration.py`

**Execution Command:**
```bash
python -m pytest pandas/tests/util/test_integration.py -v
```

**Test Results:**
```
collected 6 items

test_series_to_dataframe_dtype_preservation PASSED
test_dataframe_from_dict_mixed_series_dtypes PASSED
test_validate_fillna_with_clean_method PASSED
test_series_fillna_integration PASSED
test_check_dtype_backend_with_lib_sentinel PASSED
test_percentile_validation_with_numpy_arrays PASSED

=================================== 6 passed in 0.94s
```

**Summary:**
- **Total Tests**: 6 integration tests
- **Passed**: 6 (100%)
- **Failed**: 0
- **Execution Time**: 0.94 seconds

### Defects Discovered

**No defects were discovered during integration testing.** All module interactions functioned as expected:

- Series-to-DataFrame conversion preserves dtypes correctly
- DataFrame construction handles mixed-dtype Series properly
- Validation module correctly delegates to missing data module
- Series fillna operations integrate validation and missing data modules
- Backend validation properly handles C library sentinel values
- Percentile validation correctly integrates with NumPy array handling

All error cases (ValueError for invalid inputs) behaved as designed, raising appropriate exceptions with descriptive messages.

## Bug Reports

**No bugs identified.** All integration points between modules are functioning correctly. The following expected behaviors were verified:

1. **Type preservation across module boundaries**: Dtypes maintained through Series→DataFrame→Internals conversions
2. **Validation delegation**: Validators correctly call specialized modules (e.g., `clean_fill_method`)
3. **Error propagation**: Invalid inputs raise appropriate exceptions with clear messages
4. **Sentinel value handling**: C library sentinels (`lib.no_default`) recognized by validators

## Group Contributions

| Student | Test Cases | Modules Integrated | Coverage |
|---------|------------|-------------------|----------|
| **Sandeep Ramavath** | 2 tests | Series, DataFrame, Internals, Dtypes | Series-DataFrame conversion and construction |
| **Nithikesh Bobbili** | 2 tests | Validators, Missing Data, Series, Internals | Fillna validation and operation pipeline |
| **Mallikarjuna** | 2 tests | Validators, C Libs, NumPy | Backend validation and percentile checking |

**Total**: 6 integration tests covering 8+ distinct pandas modules with both normal and edge case scenarios.

