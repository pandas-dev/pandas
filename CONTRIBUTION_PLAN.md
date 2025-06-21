
# CONTRIBUTION_PLAN.md

## 1. Basic Information

- **Project Name:** pandas  
- **GitHub URL:** [https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)  
- **Primary Language(s):** Python, with some Cython  
- **Project Purpose:**  
  Pandas is a fast, powerful, flexible, and easy-to-use open-source data analysis and manipulation library built on top of Python. It provides data structures like `Series` (1D) and `DataFrame` (2D) that simplify working with structured data. Pandas is widely used in data science, analytics, machine learning pipelines, and scientific computing.

---

## 2. Contribution Guidelines

- **Is there a CONTRIBUTING.md file?**  
  No. The project uses a `contributing.rst` file instead of `CONTRIBUTING.md`. This file outlines detailed instructions for contributing to pandas, including:
  - Types of contributions accepted: bug fixes, documentation, enhancements, and suggestions.
  - Instructions to pick issues labeled "good first issue" or "Docs".
  - Version control workflow: Fork → Clone → Create Branch → Make Changes → Pull Request (PR).
  - Environment setup with conda and regular syncing with the main branch.
  - Guidelines for writing good commit messages, referencing issues, and ensuring tests pass.

- **Is there a Code of Conduct?**  
  Yes. The pandas project follows a [Code of Conduct](https://github.com/pandas-dev/pandas/blob/main/CODE_OF_CONDUCT.md) to ensure a harassment-free, inclusive, and welcoming environment for everyone.

- **Is a CLA (Contributor License Agreement) needed?**  
  ❌ No. Contributors are expected to comply with the open-source license (BSD 3-Clause) and follow the project's contribution and conduct guidelines.

- **Are first-time contributors welcomed?**  
  ✅ Yes. The project actively encourages contributions from first-timers and provides clear instructions for onboarding.

---

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

---

## 4. Making a Contribution

- **Open Issue Chosen:**  
  [BUG: Implicit conversion to float64 with isin() #61676](https://github.com/pandas-dev/pandas/issues/61676)

- **Issue Summary:**  
  Using `isin()` on a DataFrame column and passing a value of type `np.uint64` causes unexpected implicit conversion to `float64`, resulting in incorrect behavior.

### Steps to Resolve the Issue:
1. Reproduce the issue locally by writing a minimal test case.
2. Implement a fix by ensuring consistent data types between the DataFrame column and the values passed to `isin()`:
   - **Approach 1:** Convert the `isin()` argument to match the column's dtype.
   - **Approach 2:** Convert the DataFrame column to `uint64` if needed.
3. Add appropriate unit tests under `pandas/tests/`.
4. Ensure all existing and new tests pass using `pytest`.

---

## 5. Create a Pull Request Plan

### Pull Request Workflow:
1. Create a new feature branch:
   ```bash
   git checkout -b fix-isin-float64-conversion
   ```
2. Make the required code changes in the appropriate files.
3. Add and commit the changes:
   ```bash
   git add .
   git commit -m "BUG: Fix implicit float64 conversion in isin (#61676)"
   ```
4. Push the changes to your fork:
   ```bash
   git push origin fix-isin-float64-conversion
   ```
5. Open a Pull Request in GitHub from your branch to `pandas-dev/pandas:main`.

### Example PR Title:
```
BUG: Fix implicit float64 conversion in isin (#61676)
```

### PR Description:
```
This PR addresses issue #61676 where the use of np.uint64 with isin() leads to implicit float64 conversion, causing incorrect behavior. The fix ensures type consistency by either casting the value to match the column type or converting the column when appropriate. Additional unit tests are included to verify the fix.
```

### Testing the Fix:
- Run the test suite using `pytest pandas/tests/`.
- Confirm that all tests pass and that the new tests cover the bug scenario.
