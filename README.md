<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://pandas.pydata.org/static/img/pandas_white.svg">
    <img alt="Pandas Logo" src="https://pandas.pydata.org/static/img/pandas.svg" width="300">
  </picture>
</p>

---

# 🐼 pandas: Python Data Analysis Made Easy

| Category    | Badges                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Testing** | [![CI - Test](https://github.com/pandas-dev/pandas/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/pandas-dev/pandas/actions/workflows/unit-tests.yml) [![Coverage](https://codecov.io/github/pandas-dev/pandas/coverage.svg?branch=main)](https://codecov.io/gh/pandas-dev/pandas)                                                                                                                                                                                                              |
| **Package** | [![PyPI Version](https://img.shields.io/pypi/v/pandas.svg)](https://pypi.org/project/pandas/) [![PyPI Downloads](https://img.shields.io/pypi/dm/pandas.svg?label=PyPI%20downloads)](https://pypi.org/project/pandas/) [![Conda Version](https://anaconda.org/conda-forge/pandas/badges/version.svg)](https://anaconda.org/conda-forge/pandas) [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/pandas.svg?label=Conda%20downloads)](https://anaconda.org/conda-forge/pandas)                      |
| **Meta**    | [![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg)](https://numfocus.org) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3509134.svg)](https://doi.org/10.5281/zenodo.3509134) [![License](https://img.shields.io/pypi/l/pandas.svg)](https://github.com/pandas-dev/pandas/blob/main/LICENSE) [![Slack](https://img.shields.io/badge/join_Slack-information-brightgreen.svg?logo=slack)](https://pandas.pydata.org/docs/dev/development/community.html#community-slack) |

---

## 🧐 What is pandas?

**pandas** is a popular Python library for working with **structured data** (like tables or spreadsheets). It's built to make **data cleaning**, **analysis**, and **manipulation** easy and powerful.

If you’ve used Excel or SQL, pandas gives you similar tools — right inside Python.

It helps you:

* Load data from CSV, Excel, or databases
* Clean up messy data
* Analyze it with stats and math
* Visualize trends and patterns (when used with other libraries like Matplotlib)

> pandas is widely used in data science, machine learning, finance, and research.

---

## 📚 Table of Contents

* [Main Features](#main-features)
* [Installation](#installation)
* [Dependencies](#dependencies)
* [Installing from Source](#installing-from-source)
* [Documentation](#documentation)
* [Getting Help](#getting-help)
* [Contributing](#contributing)

---

## 🚀 Main Features

Here’s what you can do with pandas:

* Handle **missing data** easily
* Add or remove columns in a table (called a *DataFrame*)
* Automatically align data by labels
* Group and summarize data with `groupby`
* Convert lists or dictionaries into a table format
* Select data by row, column, condition, or labels
* Merge or join datasets together
* Reshape data (pivot, stack, unstack)
* Work with **multi-level indexes** (called *MultiIndex*)
* Load data from CSV, Excel, SQL, HDF5, and more
* Work with **date/time data** like a pro

[📖 Learn more about each feature](https://pandas.pydata.org/docs/user_guide/index.html)

---

## 💻 Installation

### 📦 Using conda (recommended for beginners):

```bash
conda install -c conda-forge pandas
```

### 📦 Or using pip:

```bash
pip install pandas
```

[Release Notes & Updates](https://pandas.pydata.org/pandas-docs/stable/whatsnew/index.html)

---

## 🔗 Where to Find pandas

* **Source Code**: [GitHub - pandas-dev/pandas](https://github.com/pandas-dev/pandas)
* **PyPI Package**: [pandas on PyPI](https://pypi.org/project/pandas/)
* **Conda Package**: [pandas on conda-forge](https://anaconda.org/conda-forge/pandas)

---

## ⚙️ Dependencies

To use pandas, you’ll need:

* [NumPy](https://www.numpy.org): For fast numerical operations
* [python-dateutil](https://dateutil.readthedocs.io): For working with dates
* [tzdata](https://tzdata.readthedocs.io): Timezone support

📄 [See full dependency info here](https://pandas.pydata.org/pandas-docs/stable/install.html#dependencies)

---

## 🛠️ Installing from Source

If you want to install the latest development version:

### 1. Clone the repo:

```bash
git clone https://github.com/pandas-dev/pandas.git
cd pandas
```

### 2. Install Cython (required):

```bash
pip install cython
```

### 3. Install pandas from source:

```bash
pip install .
```

### Or, for development mode:

```bash
python -m pip install -ve . --no-build-isolation --config-settings editable-verbose=true
```

📄 [Full developer setup instructions](https://pandas.pydata.org/docs/dev/development/contributing_environment.html)

---

## 📖 Documentation

* Official docs: [pandas.pydata.org/docs](https://pandas.pydata.org/pandas-docs/stable/)
* Tutorials: [User Guide](https://pandas.pydata.org/docs/user_guide/index.html)

---

## ❓ Getting Help

If you’re stuck, here are good places to ask questions:

* [Stack Overflow](https://stackoverflow.com/questions/tagged/pandas)
* [pandas Slack](https://pandas.pydata.org/docs/dev/development/community.html#community-slack)
* [PyData mailing list](https://groups.google.com/forum/#!forum/pydata)

---

## 🤝 Contributing to pandas

We welcome **everyone** to help make pandas better!

### Ways you can contribute:

* Fix bugs
* Improve documentation
* Add new features
* Review pull requests
* Help answer questions

🧭 Start by checking out:

* [Good First Issues](https://github.com/pandas-dev/pandas/issues?q=is%3Aissue+label%3A%22good+first+issue%22)
* [Contribution Guide](https://pandas.pydata.org/docs/dev/development/contributing.html)

You can also [subscribe on CodeTriage](https://www.codetriage.com/pandas-dev/pandas) to help triage issues.

📜 [Code of Conduct](https://github.com/pandas-dev/.github/blob/master/CODE_OF_CONDUCT.md)

---

## 🏁 Background

pandas was created in 2008 by Wes McKinney while working at [AQR Capital](https://www.aqr.com/) and has grown into one of the most widely-used tools in data science and analytics.

---

## 🧭 Want to Learn More?

* [pandas Docs](https://pandas.pydata.org/pandas-docs/stable/)
* [pandas YouTube Tutorials](https://www.youtube.com/results?search_query=pandas+python+tutorial)

---

🔝 [Back to Top](#-pandas-python-data-analysis-made-easy)
