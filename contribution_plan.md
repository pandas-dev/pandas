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
  No

- **Is there a Code of Conduct?**  
  ✅ Yes, the project follows a [Code of Conduct](https://github.com/pandas-dev/pandas/blob/main/.github/CODE_OF_CONDUCT.md) based on the Contributor Covenant to ensure a welcoming and respectful community.

- **Is a CLA (Contributor License Agreement) needed?**  
  ❌ No Contributor License Agreement is required for contributing to pandas.

- **Are first-time contributors welcomed?**  
  ✅ Yes, very much! The project:
  - Labels beginner-friendly issues (`good first issue`)
  - Offers clear contribution steps
  - Encourages community interaction on GitHub discussions and issues

---

## 3. Environment Setup

- **How do you set up the project locally?**

  1. **Install Anaconda**
     - Download from [https://www.anaconda.com](https://www.anaconda.com)

  2. **Create a conda environment**
     conda create -n pandas-dev python=3.10 -y
     conda activate pandas-dev
     

  3. **Clone the GitHub repository**
     git clone https://github.com/pandas-dev/pandas.git
     cd pandas

  4. **Install dependencies**
     pip install -r requirements-dev.txt

  5. **(Optional) Build pandas from source**
     python setup.py build_ext --inplace

  6. **(Optional) Run tests**
     pytest pandas

- **Any dependencies or setup steps?**  
  Yes — dependencies are managed through `requirements-dev.txt` and include:
  - `numpy`, `cython`
  - `pytest`, `mypy`, `black`, `flake8`
  - `isort`, `versioneer`, and others required for linting, testing, and building
