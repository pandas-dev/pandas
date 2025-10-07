# Pandas Baseline Build & Test Setup

This document provides instructions for reproducing the pandas baseline build and test results.

## Environment Setup

### Prerequisites
- Python 3.13.5
- Virtual environment support


### Step-by-Step Setup

1. **Clone the Repository**
   ```bash
   git clone https://github.com/saisandeepramavath/SWEN_777_Pandas.git
   cd SWEN_777_Pandas
   ```

2. **Create and Activate Virtual Environment**
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. **Upgrade pip**
   ```bash
   pip install --upgrade pip
   ```

4. **Install Dependencies**
   ```bash
   pip install -r requirements-dev.txt
   ```

## Running Tests

### Comprehensive Test Suite
To reproduce the test results, run the following command:

```bash
python -m pytest pandas/tests/series/test_constructors.py pandas/tests/frame/test_constructors.py pandas/tests/test_nanops.py pandas/tests/series/methods/test_dropna.py pandas/tests/frame/methods/test_dropna.py -v --cov=pandas --cov-report=html:courseProjectDocs/Setup/htmlcov --cov-report=term
```

### Individual Test Modules
You can also run individual test modules:

```bash
# Series constructors
python -m pytest pandas/tests/series/test_constructors.py -v

# DataFrame constructors  
python -m pytest pandas/tests/frame/test_constructors.py -v

# Numerical operations
python -m pytest pandas/tests/test_nanops.py -v

# Missing data handling
python -m pytest pandas/tests/series/methods/test_dropna.py pandas/tests/frame/methods/test_dropna.py -v
```

## Test Results Overview

The test suite executed includes:
- **Series Constructor Tests**: Core pandas Series creation and initialization
- **DataFrame Constructor Tests**: Core pandas DataFrame creation and initialization  
- **Numerical Operations Tests**: Mathematical operations and statistical functions
- **Missing Data Tests**: NA/NaN value handling and dropna functionality

## Coverage Report

The HTML coverage report is generated in `courseProjectDocs/Setup/htmlcov/index.html`. 
Open this file in a web browser to view detailed coverage information.



## Additional Information

- **Test Framework**: pytest with coverage reporting
- **Build System**: Meson + Ninja (pandas development build)
- **Python Version**: 3.13.5
- **Test Categories**: Unit tests focusing on core functionality