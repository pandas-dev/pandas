
# Pandas Code Guidelines

This document consolidates all coding style and contribution guidelines for Pandas developers into one place. Please follow these rules to ensure clean, consistent, and maintainable code.

---

## ✍️ Code Formatting

### 🔹 Use `black`
- Pandas uses [Black](https://black.readthedocs.io/en/stable/) to automatically format code.
- Run `black` before committing:
  ```bash
  black pandas/
  ```

### 🔹 Line Length
- Limit lines to **88 characters** (as per Black).
- For docstrings or comments, use **80 characters** when possible.

---

## 🚫 Linting & Static Checks

### 🔹 Use `flake8`
- Run `flake8` to check for style violations:
  ```bash
  flake8 pandas/
  ```

### 🔹 Use `isort` for import sorting
- Keeps imports grouped and ordered.
- Run:
  ```bash
  isort pandas/
  ```

---

## 📦 Imports

- Standard library imports first
- Third-party packages next (e.g., `numpy`, `pytest`)
- Local pandas imports last

**Example:**
```python
import os
import numpy as np
from pandas.core.frame import DataFrame
```

---

## 🔤 Strings

- Use **f-strings** for formatting:
  ```python
  name = "Abu"
  print(f"Hello, {name}")
  ```

- Prefer `repr()` over `str()` for debugging.
- Avoid unnecessary whitespace inside `{}`:
  ```python
  f"{x}"  ✅
  f"{ x }"  ❌
  ```

---

## 🧪 Testing Guidelines

- Use `pytest` framework.
- Test files should be placed in `pandas/tests`.
- Use fixtures instead of setup/teardown methods.
- Mark known failures:
  ```python
  @pytest.mark.xfail(reason="Issue #12345")
  def test_example():
      ...
  ```

- Run all tests with:
  ```bash
  pytest pandas/
  ```

---

## 👨‍🔧 General Python Style

- Follow [PEP8](https://pep8.org/)
- Avoid `type(x) == y` — prefer `isinstance(x, y)`
- Use `is` for `None` comparison:
  ```python
  if x is None:
      ...
  ```

- Avoid importing inside functions unless necessary.
- Avoid `import *`

---

## 🧠 Function/Variable Naming

| Type         | Style        | Example         |
|--------------|--------------|-----------------|
| Function     | `snake_case` | `compute_mean()`|
| Variable     | `snake_case` | `total_count`   |
| Class        | `PascalCase` | `DataFrame`     |
| Constant     | `UPPER_CASE` | `MAX_ROWS`      |

---

## 📝 Docstrings

- Use [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) format.
- Include:
  - Short summary
  - Parameters
  - Returns
  - Examples

**Example:**
```python
def compute_mean(data):
    '''
    Compute the mean of a list.

    Parameters
    ----------
    data : list of float
        The input data.

    Returns
    -------
    float
        Mean of the input.
    '''
```

---

## 🔁 Pre-commit Hook

- Pandas uses `pre-commit` to automate linting and formatting:
  ```bash
  pip install pre-commit
  pre-commit install
  ```

---

## 💡 Git & Pull Request Guidelines

- Use meaningful commit messages.
- Format: `DOC: Add new section to code guidelines (#issue-number)`
- PRs should be small and focused.
- Always reference the issue being addressed.

---

## 📚 Reference

- [Pandas Contribution Guide](https://pandas.pydata.org/docs/development/contributing.html)
- [Pandas Code Style](https://pandas.pydata.org/docs/development/code_style.html)
- [PEP8](https://pep8.org/)
- [Black](https://black.readthedocs.io/en/stable/)
- [pytest](https://docs.pytest.org/en/stable/)

---

*This file was created to unify the code formatting and contribution expectations for new and experienced developers working on Pandas. Following these guidelines will help ensure high-quality and maintainable code.*
