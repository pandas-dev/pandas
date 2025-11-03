# Mutation Testing Setup

This document provides instructions for setting up and running mutation testing for the pandas validation functions.

---

## Tool Used

**Mutatest 3.1.0**

Mutatest is a Python mutation testing tool that generates small code changes (mutations) and runs your test suite to verify if the tests can detect these changes.

**Installation:**
```bash
pip install mutatest==3.1.0
```

**Key Features:**
- Supports substitution mutations (our primary mode)
- Random sampling of mutations for large codebases
- Detailed reporting of detected, survived, and unknown mutations

**Documentation:** https://mutatest.readthedocs.io/

---

## How to Run Mutation Tests

### Prerequisites

1. **Navigate to the repository root:**
   ```bash
   cd /Volumes/T7Shield/SWEN777/SWEN_777_Pandas
   ```

2. **Activate the virtual environment:**
   ```bash
   source venv/bin/activate
   ```

3. **Verify installations:**
   ```bash
   python --version              # Should show Python 3.13.5
   pytest --version              # Should show pytest 8.4.2
   venv/bin/mutatest --version   # Should show mutatest 3.1.0
   ```

### Step 1: Run Tests First

Before running mutation testing, verify all tests pass:

```bash
# Run all validation tests together
python -m pytest pandas/tests/util/test_validate_endpoints.py \
                 pandas/tests/util/test_validate_percentile.py \
                 pandas/tests/util/test_validate_bool_kwarg.py -v
```

**Expected Output:** All 35 tests should pass (9 + 14 + 12)

### Step 2: Run Mutation Testing

#### Student 1 (Sandeep Ramavath) - validate_endpoints

```bash
# Final run with n=40 samples
venv/bin/mutatest -s pandas/util/_validators.py \
                  -t "python -m pytest pandas/tests/util/test_validate_endpoints.py -x" \
                  -m s -n 40 --nocov
```

#### Student 2 (Nithikesh Bobbili) - validate_percentile

```bash
# Final run with n=40 samples
venv/bin/mutatest -s pandas/util/_validators.py \
                  -t "python -m pytest pandas/tests/util/test_validate_percentile.py -x" \
                  -m s -n 40 --nocov
```

#### Student 3 (Malikarjuna ) - validate_bool_kwarg

```bash
# Final run with n=40 samples
venv/bin/mutatest -s pandas/util/_validators.py \
                  -t "python -m pytest pandas/tests/util/test_validate_bool_kwarg.py -x" \
                  -m s -n 40 --nocov
```

### Command Parameters Explained

- `-s`: Source file to mutate (pandas/util/_validators.py)
- `-t`: Test command to run (pytest with specific test file)
- `-m s`: Mutation mode - substitution (changes operators, constants, etc.)
- `-n 40`: Number of mutations to sample
- `--nocov`: Disable coverage collection for faster execution
- `-x`: pytest flag to stop on first test failure

---

## Target Files

### Source File Under Test

**File:** `pandas/util/_validators.py`
- **Total Lines:** 483
- **Total Mutation Targets:** 138 identified by mutatest

### Target Functions

| Function | Line Range | Lines | Purpose | Student |
|----------|-----------|-------|---------|---------|
| validate_endpoints(closed) | 391-420 | 30 | Validates "closed" parameter for interval boundaries | Sandeep Ramavath |
| validate_percentile(q) | 339-368 | 30 | Validates percentile values in range [0, 1] | Nithikesh Bobbili |
| validate_bool_kwarg(value, arg_name) | 228-270 | 43 | Validates boolean keyword arguments | Mallikarjuna |

**Total Lines Tested:** 103 lines across 3 functions

---

## Test Files

All test files are located in `pandas/tests/util/`

### Student 1: Sandeep Ramavath - test_validate_endpoints.py

**Function Tested:** `validate_endpoints(closed)` (lines 391-420)

**Total Tests:** 9
- Initial tests: 7
- Improvement tests: 2

**Test Coverage:**
- Valid inputs: None, "left", "right"
- Invalid inputs: empty string, uppercase, integers, invalid strings
- Return type validation (tuple)
- Mutual exclusivity of left/right flags

### Student 2: Nithikesh Bobbili - test_validate_percentile.py

**Function Tested:** `validate_percentile(q)` (lines 339-368)

**Total Tests:** 14
- Initial tests: 11
- Improvement tests: 3

**Test Coverage:**
- Valid single values: 0.0, 0.5, 1.0
- Valid collections: lists, tuples, numpy arrays
- Boundary values: 0.0 and 1.0
- Invalid values: below 0, above 1
- Mixed valid/invalid in collections
- Return type validation (ndarray)
- Precise edge cases near boundaries

### Student 3: Mallikarjuna  - test_validate_bool_kwarg.py

**Function Tested:** `validate_bool_kwarg(value, arg_name)` (lines 228-270)

**Total Tests:** 12
- Initial tests: 9
- Improvement tests: 3

**Test Coverage:**
- Valid boolean values: True, False
- None handling: allowed by default, disallowed when specified
- Integer handling: disallowed by default, allowed when specified
- Invalid types: strings, lists, floats
- Parameter combinations
- Edge case: zero as integer

**Total Tests Across All Students:** 35 tests

---

## Notes

### Important Limitations

#### Random Sampling Challenge

Mutatest samples mutations randomly from the **entire source file** (pandas/util/_validators.py, 483 lines with 138 mutation targets), not just the target functions.

**Impact:**
- Target functions cover only 103 lines (~21% of file)
- With sample size n=40, expect only ~8 mutations in target functions
- Most sampled mutations fall outside target function ranges
- This causes low overall detection percentages (5-31%)

**Key Insight:** When mutations occur within target function ranges, detection rates are ~100% for all students, demonstrating excellent test quality.

#### Interpreting Mutation Scores

**Overall Scores (appear low):**
- Student 1: 16/51 detected (31%)
- Student 2: 2/41 detected (5%)
- Student 3: 4/42 detected (10%)

**Within-Function Detection (actual quality):**
- Student 1: 16/16 detected (100%) - All sampled mutations in lines 391-420 caught
- Student 2: 2/2 detected (100%) - All sampled mutations in lines 339-368 caught
- Student 3: 4/4 detected (100%) - All sampled mutations in lines 228-270 caught

**Conclusion:** Low overall percentages reflect tool limitation (random sampling), not poor test quality.

### Mutation Types Detected

The test suites successfully detect:
1. **Boolean constant mutations:** True ↔ False ↔ None
2. **Comparison operator mutations:** == ↔ != ↔ < ↔ > ↔ <= ↔ >=
3. **If statement mutations:** If_Statement ↔ If_True ↔ If_False

---

## Group Contributions

### Student 1: Sandeep Ramavath
**Function:** `validate_endpoints(closed)` (lines 391-420)

**Contributions:**
- Created initial test suite with 7 comprehensive tests
- Covered all valid values (None, "left", "right") and invalid input scenarios
- Added 2 improvement tests targeting return type validation and mutual exclusivity
- Ran initial mutation testing (n=20) and final testing (n=40)
- Analyzed mutation results and identified patterns

**Results:**
- Initial: 4 detected, 16 survived, 1 unknown, 1 timeout (22 total runs)
- Final: 16 detected, 34 survived, 1 unknown (51 total runs)
- **Achievement:** 100% detection rate for mutations within function lines 391-420

### Student 2: Nithikesh Bobbili
**Function:** `validate_percentile(q)` (lines 339-368)

**Contributions:**
- Created comprehensive initial test suite with 11 tests
- Covered single values, collections (lists, tuples, arrays), and boundary cases
- Added 3 improvement tests targeting return type and precise boundary edges
- Ran initial mutation testing (n=20) and final testing (n=40)
- Documented edge case testing strategies

**Results:**
- Initial: 3 detected, 18 survived, 1 unknown (22 total runs)
- Final: 2 detected, 38 survived, 1 timeout (41 total runs)
- **Achievement:** 100% detection rate for mutations within function lines 339-368

### Student 3: Mallikarjuna
**Function:** `validate_bool_kwarg(value, arg_name)` (lines 228-270)

**Contributions:**
- Created initial test suite with 9 tests covering boolean validation
- Tested None handling, integer handling, and invalid type scenarios
- Added 3 improvement tests targeting parameter combinations and edge cases
- Ran initial mutation testing (n=20) and final testing (n=40)
- Analyzed sampling variance effects

**Results:**
- Initial: 0 detected, 20 survived (20 total runs) - no mutations sampled in target range
- Final: 4 detected, 38 survived (42 total runs)
- **Achievement:** 100% detection rate for mutations within function lines 228-270

### Collaborative Efforts

All team members collaborated on:
- Consistent test structure using pytest class-based organization
- Following pandas testing conventions and style guidelines
- Comprehensive documentation of findings in report.md
- Analysis of mutation testing limitations and interpretation
- Understanding the impact of random sampling on results

---
