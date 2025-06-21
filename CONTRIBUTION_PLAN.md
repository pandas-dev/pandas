# CONTRIBUTION_PLAN.md

## 1. Basic Information

- **Project Name**: pandas  
- **GitHub URL**: [https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)

- **Primary Language(s)**:  
  - Python  
  - C (some extensions)  
  - Cython (for performance-critical code)

- **What is the project used for?**  
  `pandas` is a fast, powerful, and flexible open-source data analysis and manipulation library for Python. It provides high-level data structures such as `DataFrame` and `Series`, which are essential for working with structured data, time series, and large datasets in a convenient and efficient way.

---

## 2. Contribution Guidelines

- **Are there clear steps in a `CONTRIBUTING.md` file?**  
  ✅ No

- **Is there a Code of Conduct?**  
  ✅ Yes. The project includes a [Code of Conduct](https://github.com/pandas-dev/pandas/blob/main/CODE_OF_CONDUCT.md) to ensure respectful and inclusive participation.

- **Is a CLA (Contributor License Agreement) needed?**  
  ❌ No CLA is currently required for contributions to pandas.

- **Are first-time contributors welcomed?**  
  ✅ Yes. The project welcomes first-time contributors and provides good first issues with the `good first issue` label to help newcomers get started.

---

## 3. Environment Setup

- **How do you set up the project locally?**

```bash
# Step 1: Clone your fork
git clone https://github.com/<your-username>/pandas.git
cd pandas

# Step 2: Set up a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Step 3: Upgrade pip
pip install --upgrade pip

# Step 4: Install development dependencies
pip install -r requirements-dev.txt

# Step 5: Build the C extensions (optional but recommended)
python setup.py build_ext --inplace

# Step 6: Run tests to verify setup
pytest pandas
```

- **Any dependencies or setup steps?**  
  - Python 3.8 or higher  
  - `pip`, `setuptools`, and `wheel`  
  - `pytest` for testing  
  - `numpy`, `cython`, `hypothesis`, `mypy`, `black`, and others (in `requirements-dev.txt`)  
  - Compiler tools (like GCC or MSVC) for building extensions

---

> 📘 For full and latest instructions, visit the [official contributing guide](https://pandas.pydata.org/docs/dev/development/contributing.html).
