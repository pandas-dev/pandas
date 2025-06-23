
## 1. Basic Information

- **Project Name:** pandas  
- **GitHub URL:** [https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)  
- **Primary Language(s):**
  - Python (core language)
  - C / Cython (for performance-critical components)

- **What is the project used for?**  
  pandas is a powerful, open-source library used for:
  - Data manipulation and analysis
  - Working with structured data (like CSV, Excel, SQL, JSON)
  - Offering key data structures: `DataFrame` and `Series`
  - Enabling fast, flexible operations for data cleaning, filtering, grouping, merging, and more  
  It's widely used in data science, machine learning, finance, and research.

---

## 2. Contribution Guidelines

- **Are there clear steps in a CONTRIBUTING.md file?**  
  ❌ No. The project uses a `contributing.rst` file instead of `CONTRIBUTING.md`. This file provides comprehensive guidelines for contributing to pandas, including:
  1. Accepted contribution types such as bug fixes, documentation updates, feature enhancements, and suggestions.
  2. Steps to identify suitable tasks by selecting issues labeled as "good first issue" or "Docs".
  3. A version control workflow that involves: Forking the repository → Cloning it locally → Creating a new branch → Making changes → Submitting a Pull Request (PR).
  4. Instructions for setting up the development environment using conda and regularly syncing with the upstream main branch.
  5. Best practices for writing meaningful commit messages, referencing related issues, and ensuring all tests pass before submission.

- **Is there a Code of Conduct?**  
  ✅ Yes, the project follows a [Code of Conduct](https://github.com/pandas-dev/pandas/blob/main/.github/CODE_OF_CONDUCT.md) based on the Contributor Covenant to ensure a welcoming and respectful community.

- **Is a CLA (Contributor License Agreement) needed?**  
  ❌ No Contributor License Agreement is required for contributing to pandas.

- **Are first-time contributors welcomed?**  
  ✅ Yes, very much! The project:
  - Labels beginner-friendly issues (`good first issue`)
  - Offers clear contribution steps
  - Encourages community interaction on GitHub discussions and issues

## 3. Environment Setup

### Steps to Set Up Locally:

1. Fork the repository on GitHub to your account.
2. Clone the repository locally:
   ```bash
   git clone https://github.com/<your-username>/pandas.git
   cd pandas
   ```
3. Create and activate a development environment using conda:
   ```bash
   conda create -n devenv python=3.10
   conda activate devenv
   ```
4. Install development dependencies:
   ```bash
   pip install -r requirements-dev.txt
   ```
5. Build the C extensions required by pandas:
   ```bash
   python setup.py build_ext --inplace
   ```
6. (Optional but recommended) Run the test suite to validate your environment:
   ```bash
   pytest pandas/tests/
   ```

## 4. Making a Contribution

- **Open Issue Chosen:**  
  [BUG: Groupby aggregate coercion of outputs inconsistency for pyarrow dtypes #61636](https://github.com/pandas-dev/pandas/issues/61636)

- **Issue Summary:**  
  When using `groupby(...).agg()` on PyArrow-backed DataFrames, the output types are sometimes inconsistently coerced to pandas-native dtypes like `float64`, rather than preserving the original PyArrow dtypes. This leads to unexpected results and breaks downstream workflows that rely on dtype stability.

### Steps to Resolve the Issue:

1. Reproduce the issue locally by creating a DataFrame backed by PyArrow dtypes and performing `groupby(...).agg()` with functions like `'sum'`, `'first'`, etc.
2. Implement a fix to ensure aggregation on PyArrow-backed DataFrames:
   - Maintain the original PyArrow dtypes wherever applicable.
   - Avoid coercion unless required by the aggregation operation.
3. Modify the logic in the relevant pandas core modules (likely `core/groupby/aggregation.py` or `core/groupby/groupby.py`).
4. Add targeted unit tests under `pandas/tests/groupby/` to cover aggregation behavior on PyArrow-backed data.
5. Run all tests to ensure that the issue is resolved and no regressions are introduced.

## 5. Create a Pull Request Plan

### Pull Request Workflow:

1. Create a new feature branch:
   ```bash
   git checkout -b fix-groupby-coercion-pyarrow
   ```

2. Make the required code changes in the appropriate files.

3. Add and commit the changes:
   ```bash
   git add .
   git commit -m "BUG: Fix aggregate dtype coercion on pyarrow-backed GroupBy (#61636)"
   ```

4. Push the changes to your fork:
   ```bash
   git push origin fix-groupby-coercion-pyarrow
   ```

5. Open a Pull Request in GitHub from your branch to `pandas-dev/pandas:main`.

### Example PR Title:
```
BUG: Fix aggregate dtype coercion on pyarrow-backed GroupBy (#61636)
```

### PR Description:
```
This PR addresses [#61636](https://github.com/pandas-dev/pandas/issues/61636), which reports inconsistent dtype coercion during groupby aggregation on PyArrow-backed DataFrames. Specifically, aggregations like 'sum' or 'first' on columns with Arrow dtypes (e.g., int32, uint64) may return outputs with unexpected pandas-native dtypes like float64.

The fix ensures that aggregation operations on Arrow-backed columns preserve the original Arrow dtypes wherever possible, improving consistency and reliability for downstream workflows.

New unit tests have been added to validate aggregation outputs for PyArrow-backed DataFrames and confirm dtype stability.

Closes #61636.
```

### Testing the Fix:

- Run the test suite using:
   ```bash
   pytest pandas/tests/groupby/
   ```
- Confirm that all tests pass and that the new tests adequately cover the issue scenario involving Arrow dtypes.
