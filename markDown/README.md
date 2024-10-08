# DataFrame and Series to Markdown Converter

## Overview

This package provides utility functions to convert pandas DataFrames and Series into markdown-formatted tables. These functions are integrated into the pandas library as methods for seamless usage.

## Installation

Ensure you have pandas installed:
```bash
pip install pandas
```

1. Clone the Repository

2. Import the Module: In you Python script or Jupyter Notebook, you can import the markDown.py module like so:

```python
from markDown import markdownTable, markdownSeries
```
3. Use the Functions: Once imported, you can use the `markdownTable` and `markdownSeries` functions directly, or you can utilize the extended functionality added to Pandas DataFrame and Series.

For example, to convert a DataFrame to a markdown table:

```python
import pandas as pd

# Assuming 'df' is the DataFrame you want to convert
markdown_table = df.markdownTable()
```
Or, to convert a Series to a markdown table:

```python
import pandas as pd

# Assuming 'series' is the Series you want to convert
markdown_table = series.markdownSeries()
```

Users can also specify specific columns for DataFrame conversion:

```python
markdown_table = df.markdownTable('Column1', 'Column2', 'Column3')
```

Or customize the column names for Series conversion:

```python
markdown_table = series.markdownTable(col1='Index', col2='Value')
```
## Conclusion
These functions provide an easy way to convert pandas DataFrames and Series into markdown-formatted tables, enhancing the readability and presentation of your data in markdown-supported environments. This is simular to using `to_markdown` in pandas without having to install another library. 

pandas, the `to_markdown` method is provided by the `tabulate` library, which needs to be installed separately. The to_markdown method is available in pandas version 1.3.0 and later.

