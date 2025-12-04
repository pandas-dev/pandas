# System Testing Report

## Executive Summary

This report documents the system-level black-box testing effort for the pandas library. Our team created 3 system test cases that validate complete end-to-end workflows through pandas' public API, treating the system as a black box without referencing internal implementation details.

---

## Test Scope and Coverage

### Testing Approach

Our system tests validate pandas' core functionality through **black-box testing**, meaning:
- Tests interact only through public APIs (DataFrame, Series, read_csv, to_csv, etc.)
- No reference to internal implementation or private methods
- Tests simulate real user workflows from start to finish
- Validation based on observable behavior and outputs

### Workflows Validated

The system tests cover three critical end-to-end workflows that represent typical pandas usage patterns:

#### 1. Data Loading and Export Workflow (Sandeep Ramavath)
**Scope:** Complete data I/O lifecycle
- **Features Tested:**
  - CSV file import (`pd.read_csv()`)
  - CSV file export (`DataFrame.to_csv()`)
  - Mixed data type handling (integers, floats, strings, dates, booleans)
  - Data persistence and round-trip integrity
  - Datetime parsing during import
  
- **User Story:** "As a data analyst, I want to load data from CSV files, work with it in pandas, and export results back to CSV so that I can share my analysis with others."

#### 2. Data Cleaning and Transformation Workflow (Nithikesh Bobbili)
**Scope:** Missing data handling and data quality
- **Features Tested:**
  - Missing value detection (`isnull()`, `sum()`)
  - Forward fill strategy (`ffill()`)
  - Backward fill strategy (`bfill()`)
  - Constant value fill (`fillna()`)
  - Data integrity preservation during cleaning
  
- **User Story:** "As a data scientist, I want to identify and handle missing values in my dataset using various filling strategies so that I can prepare clean data for analysis."

#### 3. Aggregation and Analysis Workflow (Mallikarjuna)
**Scope:** Group-by operations and statistical analysis
- **Features Tested:**
  - Categorical grouping (`groupby()`)
  - Statistical aggregations (mean, sum, count)
  - Multiple simultaneous aggregations (`agg()`)
  - Grouped data integrity
  - Result correctness verification
  
- **User Story:** "As a business analyst, I want to group data by categories and compute statistics for each group so that I can understand patterns and trends in my data."

### Coverage Metrics

| Workflow Category | Public APIs Used | Test Cases | Assertions |
|------------------|------------------|------------|------------|
| Data I/O | 2 APIs | 1 | 8 |
| Data Cleaning | 4 APIs | 1 | 11 |
| Data Aggregation | 4 APIs | 1 | 13 |
| **Total** | **10 unique APIs** | **3** | **32** |

### Out of Scope

The following are explicitly **not** tested in this system testing phase:
- Internal implementation details (block managers, internals, etc.)
- Performance benchmarks or optimization
- Edge cases requiring white-box knowledge
- Deprecated or experimental APIs
- Platform-specific behaviors

---

## Test Case Summary

### Test Case 1: CSV Data Import-Export Workflow

**Test ID:** SYS-001  
**Owner:** Sandeep Ramavath  
**Category:** Data I/O Workflow  
**Test File:** `pandas/tests/system/test_system_workflows.py::TestDataIOWorkflow::test_csv_roundtrip_workflow`

| Attribute | Details |
|-----------|---------|
| **Title** | CSV Data Import-Export Workflow |
| **Pre-conditions** | • Temporary directory available for file operations<br>• pandas library installed and functional<br>• Write permissions in test directory |
| **Test Steps** | **Step 1:** Create DataFrame with mixed data types using public API<br>&nbsp;&nbsp;- Create DataFrame with 5 columns: id (int), name (string), score (float), date (datetime), active (boolean)<br>&nbsp;&nbsp;- Use pandas constructor: `pd.DataFrame()`<br><br>**Step 2:** Export DataFrame to CSV file<br>&nbsp;&nbsp;- Call `to_csv()` method with file path<br>&nbsp;&nbsp;- Use `index=False` parameter<br>&nbsp;&nbsp;- Verify file creation on disk<br><br>**Step 3:** Import CSV file back into new DataFrame<br>&nbsp;&nbsp;- Call `pd.read_csv()` with file path<br>&nbsp;&nbsp;- Use `parse_dates` parameter for date column<br><br>**Step 4:** Verify data integrity and type preservation<br>&nbsp;&nbsp;- Check row count matches original (5 rows)<br>&nbsp;&nbsp;- Verify column names preserved<br>&nbsp;&nbsp;- Compare all values with original data<br>&nbsp;&nbsp;- Verify datetime type correctly parsed |
| **Expected Results** | • CSV file created successfully at specified path<br>• File contains 5 data rows plus header<br>• Data round-trips without any loss<br>• Integer values: [1, 2, 3, 4, 5] preserved<br>• String values: ['Alice', 'Bob', 'Charlie', 'David', 'Eve'] preserved<br>• Float values: [95.5, 87.3, 92.1, 88.7, 91.4] preserved<br>• Boolean values: [True, False, True, True, False] preserved<br>• Date column recognized as datetime64 type<br>• All assertions pass without errors |
| **Actual Results** | **PASSED** - All expected results achieved |

---

### Test Case 2: Missing Data Cleaning Workflow

**Test ID:** SYS-002  
**Owner:** Nithikesh Bobbili  
**Category:** Data Cleaning Workflow  
**Test File:** `pandas/tests/system/test_system_workflows.py::TestDataCleaningWorkflow::test_missing_data_handling_workflow`

| Attribute | Details |
|-----------|---------|
| **Title** | Missing Data Cleaning Workflow |
| **Pre-conditions** | • pandas library available<br>• numpy library available for NaN values<br>• No external dependencies or files required |
| **Test Steps** | **Step 1:** Create DataFrame with missing values using public API<br>&nbsp;&nbsp;- Create 3-column DataFrame with `np.nan` values<br>&nbsp;&nbsp;- Column A: 2 missing values at positions 1 and 3<br>&nbsp;&nbsp;- Column B: 2 missing values at positions 0 and 2<br>&nbsp;&nbsp;- Column C: 1 missing value at position 4<br><br>**Step 2:** Detect missing values using public methods<br>&nbsp;&nbsp;- Call `isnull()` to create boolean mask<br>&nbsp;&nbsp;- Call `sum()` to count missing values per column<br>&nbsp;&nbsp;- Verify counts: A=2, B=2, C=1<br><br>**Step 3:** Fill missing values using multiple strategies<br>&nbsp;&nbsp;- **Strategy 3a:** Forward fill using `ffill()`<br>&nbsp;&nbsp;&nbsp;&nbsp;- Verify propagation of last valid value<br>&nbsp;&nbsp;&nbsp;&nbsp;- Check remaining NaN count<br>&nbsp;&nbsp;- **Strategy 3b:** Backward fill using `bfill()`<br>&nbsp;&nbsp;&nbsp;&nbsp;- Verify propagation of next valid value<br>&nbsp;&nbsp;&nbsp;&nbsp;- Check remaining NaN count<br>&nbsp;&nbsp;- **Strategy 3c:** Constant fill using `fillna(0)`<br>&nbsp;&nbsp;&nbsp;&nbsp;- Verify all NaN replaced with 0<br><br>**Step 4:** Verify all missing values handled correctly<br>&nbsp;&nbsp;- Confirm no NaN values remain after constant fill<br>&nbsp;&nbsp;- Verify DataFrame shape preserved<br>&nbsp;&nbsp;- Check specific filled values match expectations |
| **Expected Results** | • Missing values correctly identified: A=2, B=2, C=1 (total 5)<br>• Forward fill leaves 1 NaN (at first position of column B)<br>• Forward fill propagates value correctly (row 1, col A = 1.0)<br>• Backward fill leaves 1 NaN (at last position of column C)<br>• Backward fill propagates value correctly (row 0, col B = 2.0)<br>• Constant fill (value=0) removes all NaN values<br>• Constant fill replaces NaN with exact value (row 1, col A = 0.0)<br>• DataFrame shape (5, 3) preserved after all operations<br>• All assertions pass without errors |
| **Actual Results** | **PASSED** - All expected results achieved |

---

### Test Case 3: Group-by Aggregation Analysis Workflow

**Test ID:** SYS-003  
**Owner:** Mallikarjuna  
**Category:** Aggregation and Analysis Workflow  
**Test File:** `pandas/tests/system/test_system_workflows.py::TestAggregationWorkflow::test_groupby_aggregation_workflow`

| Attribute | Details |
|-----------|---------|
| **Title** | Group-by Aggregation Analysis Workflow |
| **Pre-conditions** | • pandas library functional<br>• Sufficient memory for group operations<br>• No external data sources required |
| **Test Steps** | **Step 1:** Create DataFrame with categorical and numeric data<br>&nbsp;&nbsp;- Create DataFrame with 3 columns:<br>&nbsp;&nbsp;&nbsp;&nbsp;- category: ['A', 'B', 'A', 'B', 'A', 'B', 'A', 'B']<br>&nbsp;&nbsp;&nbsp;&nbsp;- value: [10, 20, 15, 25, 20, 30, 25, 35]<br>&nbsp;&nbsp;&nbsp;&nbsp;- quantity: [1, 2, 3, 4, 5, 6, 7, 8]<br>&nbsp;&nbsp;- Total 8 rows, evenly split between categories A and B<br><br>**Step 2:** Group data by category using public API<br>&nbsp;&nbsp;- Call `groupby('category')` on DataFrame<br>&nbsp;&nbsp;- Store grouped object for multiple operations<br><br>**Step 3:** Apply multiple aggregation functions<br>&nbsp;&nbsp;- **Step 3a:** Apply mean aggregation on 'value' column<br>&nbsp;&nbsp;&nbsp;&nbsp;- Calculate average for each category<br>&nbsp;&nbsp;&nbsp;&nbsp;- Verify: Category A mean = 17.5, Category B mean = 27.5<br>&nbsp;&nbsp;- **Step 3b:** Apply sum aggregation on 'value' column<br>&nbsp;&nbsp;&nbsp;&nbsp;- Calculate total for each category<br>&nbsp;&nbsp;&nbsp;&nbsp;- Verify: Category A sum = 70, Category B sum = 110<br>&nbsp;&nbsp;- **Step 3c:** Apply count aggregation<br>&nbsp;&nbsp;&nbsp;&nbsp;- Count items in each category using `size()`<br>&nbsp;&nbsp;&nbsp;&nbsp;- Verify: Category A count = 4, Category B count = 4<br><br>**Step 4:** Apply multiple aggregations simultaneously<br>&nbsp;&nbsp;- Use `agg(['mean', 'sum', 'count'])` on grouped data<br>&nbsp;&nbsp;- Create multi-column result DataFrame<br><br>**Step 5:** Verify aggregated results comprehensively<br>&nbsp;&nbsp;- Check all 6 values (2 categories × 3 aggregations)<br>&nbsp;&nbsp;- Verify result DataFrame shape is (2, 3)<br>&nbsp;&nbsp;- Confirm index contains category labels |
| **Expected Results** | • Data groups correctly into 2 categories (A and B)<br>• Category A mean aggregation: (10+15+20+25)/4 = 17.5 ✓<br>• Category B mean aggregation: (20+25+30+35)/4 = 27.5 ✓<br>• Category A sum aggregation: 10+15+20+25 = 70 ✓<br>• Category B sum aggregation: 20+25+30+35 = 110 ✓<br>• Category A count: 4 items ✓<br>• Category B count: 4 items ✓<br>• Multi-aggregation creates DataFrame with shape (2, 3)<br>• Multi-aggregation preserves all individual results<br>• Result index contains 'A' and 'B' as category labels<br>• Result columns contain 'mean', 'sum', 'count'<br>• All assertions pass without errors |
| **Actual Results** | **PASSED** - All expected results achieved |

---

## Execution and Results

### Test Environment

**Test File:** `pandas/tests/system/test_system_workflows.py`  
**Testing Framework:** pytest 8.4.2  
**Python Version:** 3.13.5  
**Pandas Version:** 3.0.0.dev0+  
**NumPy Version:** 1.26+  
**Operating System:** macOS

### Execution Command

```bash
python -m pytest pandas/tests/system/test_system_workflows.py -v
```

### Test Results Summary

```
collected 3 items

pandas/tests/system/test_system_workflows.py::TestDataIOWorkflow::test_csv_roundtrip_workflow PASSED [33%]
pandas/tests/system/test_system_workflows.py::TestDataCleaningWorkflow::test_missing_data_handling_workflow PASSED [66%]
pandas/tests/system/test_system_workflows.py::TestAggregationWorkflow::test_groupby_aggregation_workflow PASSED [100%]

=================================== 3 passed in 0.52s ===================================
```

### Detailed Test Results

| Test Case | Status | Duration | Assertions | Outcome |
|-----------|--------|----------|------------|---------|
| CSV Data Import-Export Workflow | PASSED | ~0.18s | 8 | All data round-tripped correctly |
| Missing Data Cleaning Workflow | PASSED | ~0.16s | 11 | All fill strategies worked as expected |
| Group-by Aggregation Workflow | PASSED | ~0.18s | 13 | All aggregations computed correctly |

**Summary Statistics:**
- **Total Test Cases:** 3
- **Passed:** 3 (100%)
- **Failed:** 0 (0%)
- **Skipped:** 0
- **Total Execution Time:** 0.52 seconds
- **Average Test Duration:** 0.17 seconds
- **Total Assertions:** 32
- **Assertions Passed:** 32 (100%)

### Behavioral Analysis

#### Test Case 1: CSV Roundtrip Workflow

**Expected Behavior:**
- CSV export creates valid file with proper formatting
- CSV import reconstructs DataFrame with same data
- Mixed data types (int, float, string, datetime, bool) preserved
- Datetime parsing works correctly with `parse_dates` parameter

**Actual Behavior:**
**Matches Expected** - All behaviors confirmed:
- CSV file created with proper structure (header + 5 data rows)
- All numeric values preserved exactly (no rounding errors)
- String values preserved with proper encoding
- Boolean values correctly written and read as True/False
- Datetime column parsed correctly to datetime64 dtype
- No data loss or corruption during round-trip

**Deviations:** None

#### Test Case 2: Missing Data Cleaning Workflow

**Expected Behavior:**
- `isnull()` correctly identifies NaN values
- `ffill()` propagates last valid observation forward
- `bfill()` propagates next valid observation backward
- `fillna()` replaces all NaN with specified constant
- Original DataFrame remains unchanged (immutable operations)

**Actual Behavior:**
**Matches Expected** - All behaviors confirmed:
- Missing value detection accurate (5 total NaN values identified)
- Forward fill correctly propagated values, leaving only leading NaN
- Backward fill correctly propagated values, leaving only trailing NaN
- Constant fill successfully eliminated all NaN values
- Shape and non-NaN values preserved across all operations
- Original DataFrame immutable (each operation returns new DataFrame)

**Deviations:** None

#### Test Case 3: Group-by Aggregation Workflow

**Expected Behavior:**
- `groupby()` splits data by category labels
- Aggregation functions compute correct statistics per group
- Multiple aggregations can be applied simultaneously
- Result maintains category labels as index

**Actual Behavior:**
**Matches Expected** - All behaviors confirmed:
- Data correctly split into 2 groups (A and B)
- Mean calculations accurate: A=17.5, B=27.5
- Sum calculations accurate: A=70, B=110
- Count calculations accurate: A=4, B=4
- Multi-aggregation created proper DataFrame structure
- Category labels preserved in result index
- All numeric computations precise (no floating-point errors)

**Deviations:** None

### Failures and Deviations

**Result: No failures or behavioral deviations were discovered.**

All system tests passed successfully, indicating that:
- End-to-end workflows function as designed
- Public APIs behave according to documentation
- Data integrity maintained across operations
- No unexpected errors or exceptions
- All user workflows complete successfully

### Test Coverage Analysis

The system tests successfully validated:

| Workflow Component | Validated | Evidence |
|-------------------|-----------|----------|
| File I/O operations | Yes | CSV roundtrip successful |
| Data type handling | Yes | 5 different types preserved |
| Missing value detection | Yes | All NaN values identified |
| Fill strategies | Yes | 3 strategies all worked |
| Grouping operations | Yes | Categories split correctly |
| Aggregation functions | Yes | 3 aggregations accurate |
| Multi-aggregation | Yes | Combined aggregations worked |
| Data immutability | Yes | Original data preserved |

---

## Group Contributions

### Individual Contributions

| Student | Test Cases | Workflow Validated | Assertions | LOC |
|---------|------------|-------------------|------------|-----|
| **Sandeep Ramavath** | 1 test case | Data I/O (CSV import/export) | 8 | ~50 |
| **Nithikesh Bobbili** | 1 test case | Data Cleaning (missing data handling) | 11 | ~60 |
| **Mallikarjuna** | 1 test case | Aggregation (groupby operations) | 13 | ~65 |

