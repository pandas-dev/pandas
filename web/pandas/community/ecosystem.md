# Ecosystem

Increasingly, packages are being built on top of pandas to address
specific needs in data preparation, analysis and visualization. This is
encouraging because it means pandas is not only helping users to handle
their data tasks but also that it provides a better starting point for
developers to build powerful and more focused data tools. The creation
of libraries that complement pandas' functionality also allows pandas
development to remain focused around its original requirements.

This is an community-maintained list of projects that build on pandas in order
to provide tools in the PyData space. The pandas core development team does not necessarily endorse any particular project on this list or have any knowledge of the maintenance status of any particular library.

For a more complete list of projects that depend on pandas, see the [libraries.io usage page for
pandas](https://libraries.io/pypi/pandas/usage) or [search pypi for
pandas](https://pypi.org/search/?q=pandas).

We'd like to make it easier for users to find these projects, if you
know of other substantial projects that you feel should be on this list,
please let us know.

## Statistics and machine learning

### [Statsmodels](https://www.statsmodels.org/)

Statsmodels is the prominent Python "statistics and econometrics
library" and it has a long-standing special relationship with pandas.
Statsmodels provides powerful statistics, econometrics, analysis and
modeling functionality that is out of pandas' scope. Statsmodels
leverages pandas objects as the underlying data container for
computation.

### [Featuretools](https://github.com/alteryx/featuretools/)

Featuretools is a Python library for automated feature engineering built
on top of pandas. It excels at transforming temporal and relational
datasets into feature matrices for machine learning using reusable
feature engineering "primitives". Users can contribute their own
primitives in Python and share them with the rest of the community.

### [Compose](https://github.com/alteryx/compose)

Compose is a machine learning tool for labeling data and prediction engineering.
It allows you to structure the labeling process by parameterizing
prediction problems and transforming time-driven relational data into
target values with cutoff times that can be used for supervised learning.

### [STUMPY](https://github.com/TDAmeritrade/stumpy)

STUMPY is a powerful and scalable Python library for modern time series analysis.
At its core, STUMPY efficiently computes something called a
[matrix profile](https://stumpy.readthedocs.io/en/latest/Tutorial_The_Matrix_Profile.html),
which can be used for a wide variety of time series data mining tasks.

## Visualization

### [Altair](https://altair-viz.github.io/)

Altair is a declarative statistical visualization library for Python.
With Altair, you can spend more time understanding your data and its
meaning. Altair's API is simple, friendly and consistent and built on
top of the powerful Vega-Lite JSON specification. This elegant
simplicity produces beautiful and effective visualizations with a
minimal amount of code. Altair works with Pandas DataFrames.

### [Bokeh](https://docs.bokeh.org)

Bokeh is a Python interactive visualization library for large datasets
that natively uses the latest web technologies. Its goal is to provide
elegant, concise construction of novel graphics in the style of
Protovis/D3, while delivering high-performance interactivity over large
data to thin clients.

[Pandas-Bokeh](https://github.com/PatrikHlobil/Pandas-Bokeh) provides a
high level API for Bokeh that can be loaded as a native Pandas plotting
backend via

```
pd.set_option("plotting.backend", "pandas_bokeh")
```

It is very similar to the matplotlib plotting backend, but provides
interactive web-based charts and maps.

### [pygwalker](https://github.com/Kanaries/pygwalker)

PyGWalker is an interactive data visualization and
exploratory data analysis tool built upon Graphic Walker
with support for visualization, cleaning, and annotation workflows.

pygwalker can save interactively created charts
to Graphic-Walker and Vega-Lite JSON.

```
import pygwalker as pyg
pyg.walk(df)
```

### [seaborn](https://seaborn.pydata.org)

Seaborn is a Python visualization library based on
[matplotlib](https://matplotlib.org). It provides a high-level,
dataset-oriented interface for creating attractive statistical graphics.
The plotting functions in seaborn understand pandas objects and leverage
pandas grouping operations internally to support concise specification
of complex visualizations. Seaborn also goes beyond matplotlib and
pandas with the option to perform statistical estimation while plotting,
aggregating across observations and visualizing the fit of statistical
models to emphasize patterns in a dataset.

```
import seaborn as sns
sns.set_theme()
```

### [plotnine](https://github.com/has2k1/plotnine/)

Hadley Wickham's [ggplot2](https://ggplot2.tidyverse.org/) is a
foundational exploratory visualization package for the R language. Based
on ["The Grammar of
Graphics"](https://www.cs.uic.edu/~wilkinson/TheGrammarOfGraphics/GOG.html)
it provides a powerful, declarative and extremely general way to
generate bespoke plots of any kind of data.
Various implementations to other languages are available.
A good implementation for Python users is [has2k1/plotnine](https://github.com/has2k1/plotnine/).

### [IPython Vega](https://github.com/vega/ipyvega)

[IPython Vega](https://github.com/vega/ipyvega) leverages [Vega](https://github.com/vega/vega) to create plots within Jupyter Notebook.

### [Plotly](https://plot.ly/python)

[Plotly's](https://plot.ly/) [Python API](https://plot.ly/python/)
enables interactive figures and web shareability. Maps, 2D, 3D, and
live-streaming graphs are rendered with WebGL and
[D3.js](https://d3js.org/). The library supports plotting directly from
a pandas DataFrame and cloud-based collaboration. Users of [matplotlib,
ggplot for Python, and
Seaborn](https://plot.ly/python/matplotlib-to-plotly-tutorial/) can
convert figures into interactive web-based plots. Plots can be drawn in
[IPython Notebooks](https://plot.ly/ipython-notebooks/) , edited with R
or MATLAB, modified in a GUI, or embedded in apps and dashboards. Plotly
is free for unlimited sharing, and has
[cloud](https://plot.ly/product/plans/),
[offline](https://plot.ly/python/offline/), or
[on-premise](https://plot.ly/product/enterprise/) accounts for private
use.

### [Lux](https://github.com/lux-org/lux)

Lux is a Python library that facilitates fast and easy experimentation with data by automating the visual data exploration process. To use Lux, simply add an extra import alongside pandas:

```python
import lux
import pandas as pd

df = pd.read_csv("data.csv")
df  # discover interesting insights!
```

By printing out a dataframe, Lux automatically [recommends a set of visualizations](https://raw.githubusercontent.com/lux-org/lux-resources/master/readme_img/demohighlight.gif) that highlights interesting trends and patterns in the dataframe. Users can leverage any existing pandas commands without modifying their code, while being able to visualize their pandas data structures (e.g., DataFrame, Series, Index) at the same time. Lux also offers a [powerful, intuitive language](https://lux-api.readthedocs.io/en/latest/source/guide/vis.html>) that allow users to create  Altair, matplotlib, or Vega-Lite visualizations without having to think at the level of code.

### [D-Tale](https://github.com/man-group/dtale)

D-Tale is a lightweight web client for visualizing pandas data structures. It
provides a rich spreadsheet-style grid which acts as a wrapper for a lot of
pandas functionality (query, sort, describe, corr...) so users can quickly
manipulate their data. There is also an interactive chart-builder using Plotly
Dash allowing users to build nice portable visualizations. D-Tale can be
invoked with the following command

```python
import dtale

dtale.show(df)
```

D-Tale integrates seamlessly with Jupyter notebooks, Python terminals, Kaggle
& Google Colab. Here are some demos of the [grid](http://alphatechadmin.pythonanywhere.com/dtale/main/1).

### [hvplot](https://hvplot.holoviz.org/index.html)

hvPlot is a high-level plotting API for the PyData ecosystem built on [HoloViews](https://holoviews.org/).
It can be loaded as a native pandas plotting backend via

```python
pd.set_option("plotting.backend", "hvplot")
```

## IDE

### [IPython](https://ipython.org/documentation.html)

IPython is an interactive command shell and distributed computing
environment. IPython tab completion works with Pandas methods and also
attributes like DataFrame columns.

### [Jupyter Notebook / Jupyter Lab](https://jupyter.org)

Jupyter Notebook is a web application for creating Jupyter notebooks. A
Jupyter notebook is a JSON document containing an ordered list of
input/output cells which can contain code, text, mathematics, plots and
rich media. Jupyter notebooks can be converted to a number of open
standard output formats (HTML, HTML presentation slides, LaTeX, PDF,
ReStructuredText, Markdown, Python) through 'Download As' in the web
interface and `jupyter convert` in a shell.

Pandas DataFrames implement `_repr_html_`and `_repr_latex` methods which
are utilized by Jupyter Notebook for displaying (abbreviated) HTML or
LaTeX tables. LaTeX output is properly escaped. (Note: HTML tables may
or may not be compatible with non-HTML Jupyter output formats.)

See [Options and Settings](https://pandas.pydata.org/docs/user_guide/options.html)
for pandas `display.` settings.

### [Spyder](https://www.spyder-ide.org/)

Spyder is a cross-platform PyQt-based IDE combining the editing,
analysis, debugging and profiling functionality of a software
development tool with the data exploration, interactive execution, deep
inspection and rich visualization capabilities of a scientific
environment like MATLAB or Rstudio.

Its [Variable
Explorer](https://docs.spyder-ide.org/current/panes/variableexplorer.html) allows
users to view, manipulate and edit pandas `Index`, `Series`, and
`DataFrame` objects like a "spreadsheet", including copying and
modifying values, sorting, displaying a "heatmap", converting data
types and more. Pandas objects can also be renamed, duplicated, new
columns added, copied/pasted to/from the clipboard (as TSV), and
saved/loaded to/from a file. Spyder can also import data from a variety
of plain text and binary files or the clipboard into a new pandas
DataFrame via a sophisticated import wizard.

Most pandas classes, methods and data attributes can be autocompleted in
Spyder's [Editor](https://docs.spyder-ide.org/current/panes/editor.html) and [IPython
Console](https://docs.spyder-ide.org/current/panes/ipythonconsole.html), and Spyder's
[Help pane](https://docs.spyder-ide.org/current/panes/help.html) can retrieve and
render Numpydoc documentation on pandas objects in rich text with Sphinx
both automatically and on-demand.

## API

### [pandas-datareader](https://github.com/pydata/pandas-datareader)

`pandas-datareader` is a remote data access library for pandas
(PyPI:`pandas-datareader`). It is based on functionality that was
located in `pandas.io.data` and `pandas.io.wb` but was split off in
v0.19. See more in the [pandas-datareader
docs](https://pandas-datareader.readthedocs.io/en/latest/):

The following data feeds are available:

- Google Finance
- Tiingo
- Morningstar
- IEX
- Robinhood
- Enigma
- Quandl
- FRED
- Fama/French
- World Bank
- OECD
- Eurostat
- TSP Fund Data
- Nasdaq Trader Symbol Definitions
- Stooq Index Data
- MOEX Data

### [pandaSDMX](https://pandasdmx.readthedocs.io)

pandaSDMX is a library to retrieve and acquire statistical data and
metadata disseminated in [SDMX](https://sdmx.org) 2.1, an
ISO-standard widely used by institutions such as statistics offices,
central banks, and international organisations. pandaSDMX can expose
datasets and related structural metadata including data flows,
code-lists, and data structure definitions as pandas Series or
MultiIndexed DataFrames.

### [fredapi](https://github.com/mortada/fredapi)

fredapi is a Python interface to the [Federal Reserve Economic Data
(FRED)](https://fred.stlouisfed.org/) provided by the Federal Reserve
Bank of St. Louis. It works with both the FRED database and ALFRED
database that contains point-in-time data (i.e. historic data
revisions). fredapi provides a wrapper in Python to the FRED HTTP API,
and also provides several convenient methods for parsing and analyzing
point-in-time data from ALFRED. fredapi makes use of pandas and returns
data in a Series or DataFrame. This module requires a FRED API key that
you can obtain for free on the FRED website.

## Domain specific

### [Geopandas](https://github.com/geopandas/geopandas)

Geopandas extends pandas data objects to include geographic information
which support geometric operations. If your work entails maps and
geographical coordinates, and you love pandas, you should take a close
look at Geopandas.

### [gurobipy-pandas](https://github.com/Gurobi/gurobipy-pandas)

gurobipy-pandas provides a convenient accessor API to connect pandas with
gurobipy. It enables users to more easily and efficiently build mathematical
optimization models from data stored in DataFrames and Series, and to read
solutions back directly as pandas objects.

### [staircase](https://github.com/staircase-dev/staircase)

staircase is a data analysis package, built upon pandas and numpy, for modelling and
manipulation of mathematical step functions. It provides a rich variety of arithmetic
operations, relational operations, logical operations, statistical operations and
aggregations for step functions defined over real numbers, datetime and timedelta domains.

### [xarray](https://github.com/pydata/xarray)

xarray brings the labeled data power of pandas to the physical sciences
by providing N-dimensional variants of the core pandas data structures.
It aims to provide a pandas-like and pandas-compatible toolkit for
analytics on multi-dimensional arrays, rather than the tabular data for
which pandas excels.

## IO

### [NTV-pandas](https://github.com/loco-philippe/ntv-pandas)

NTV-pandas provides a JSON converter with more data types than the ones supported by pandas directly.

It supports the following data types:

- pandas data types
- data types defined in the [NTV format](https://loco-philippe.github.io/ES/JSON%20semantic%20format%20(JSON-NTV).htm)
- data types defined in [Table Schema specification](http://dataprotocols.org/json-table-schema/#field-types-and-formats)

The interface is always reversible (conversion round trip) with two formats (JSON-NTV and JSON-TableSchema).

Example:

```python
import ntv_pandas as npd

jsn = df.npd.to_json(table=False)  # save df as a JSON-value (format Table Schema if table is True else format NTV )
df  = npd.read_json(jsn)  # load a JSON-value as a `DataFrame`

df.equals(npd.read_json(df.npd.to_json(df)))  # `True` in any case, whether `table=True` or not
```

### [BCPandas](https://github.com/yehoshuadimarsky/bcpandas)

BCPandas provides high performance writes from pandas to Microsoft SQL Server,
far exceeding the performance of the native ``df.to_sql`` method. Internally, it uses
Microsoft's BCP utility, but the complexity is fully abstracted away from the end user.
Rigorously tested, it is a complete replacement for ``df.to_sql``.

### [Deltalake](https://pypi.org/project/deltalake)

Deltalake python package lets you access tables stored in
[Delta Lake](https://delta.io/) natively in Python without the need to use Spark or
JVM. It provides the ``delta_table.to_pyarrow_table().to_pandas()`` method to convert
any Delta table into Pandas dataframe.

### [pandas-gbq](https://github.com/googleapis/python-bigquery-pandas)

pandas-gbq provides high performance reads and writes to and from
[Google BigQuery](https://cloud.google.com/bigquery/). Previously (before version 2.2.0),
these methods were exposed as `pandas.read_gbq` and `DataFrame.to_gbq`.
Use `pandas_gbq.read_gbq` and `pandas_gbq.to_gbq`, instead.

## Out-of-core

### [Bodo](https://bodo.ai/)

Bodo is a high-performance Python computing engine that automatically parallelizes and
optimizes your code through compilation using HPC (high-performance computing) techniques.
Designed to operate with native pandas dataframes, Bodo compiles your pandas code to execute
across multiple cores on a single machine or distributed clusters of multiple compute nodes efficiently.
Bodo also makes distributed pandas dataframes queryable with SQL.

The community edition of Bodo is free to use on up to 8 cores. Beyond that, Bodo offers a paid
enterprise edition. Free licenses of Bodo (for more than 8 cores) are available
[upon request](https://www.bodo.ai/contact) for academic and non-profit use.

### [Cylon](https://cylondata.org/)

Cylon is a fast, scalable, distributed memory parallel runtime with a pandas
like Python DataFrame API. ”Core Cylon” is implemented with C++ using Apache
Arrow format to represent the data in-memory. Cylon DataFrame API implements
most of the core operators of pandas such as merge, filter, join, concat,
group-by, drop_duplicates, etc. These operators are designed to work across
thousands of cores to scale applications. It can interoperate with pandas
DataFrame by reading data from pandas or converting data to pandas so users
can selectively scale parts of their pandas DataFrame applications.

```python
from pycylon import read_csv, DataFrame, CylonEnv
from pycylon.net import MPIConfig

# Initialize Cylon distributed environment
config: MPIConfig = MPIConfig()
env: CylonEnv = CylonEnv(config=config, distributed=True)

df1: DataFrame = read_csv('/tmp/csv1.csv')
df2: DataFrame = read_csv('/tmp/csv2.csv')

# Using 1000s of cores across the cluster to compute the join
df3: Table = df1.join(other=df2, on=[0], algorithm="hash", env=env)

print(df3)
```

### [Dask](https://docs.dask.org)

Dask is a flexible parallel computing library for analytics. Dask
provides a familiar `DataFrame` interface for out-of-core, parallel and
distributed computing.

### [Dask-ML](https://ml.dask.org)

Dask-ML enables parallel and distributed machine learning using Dask
alongside existing machine learning libraries like Scikit-Learn,
XGBoost, and TensorFlow.

### [Ibis](https://ibis-project.org/docs/)

Ibis offers a standard way to write analytics code, that can be run in multiple engines. It helps in bridging the gap between local Python environments (like pandas) and remote storage and execution systems like Hadoop components (like HDFS, Impala, Hive, Spark) and SQL databases (Postgres, etc.).


### [Koalas](https://koalas.readthedocs.io/en/latest/)

Koalas provides a familiar pandas DataFrame interface on top of Apache
Spark. It enables users to leverage multi-cores on one machine or a
cluster of machines to speed up or scale their DataFrame code.

### [Modin](https://github.com/modin-project/modin)

The ``modin.pandas`` DataFrame is a parallel and distributed drop-in replacement
for pandas. This means that you can use Modin with existing pandas code or write
new code with the existing pandas API. Modin can leverage your entire machine or
cluster to speed up and scale your pandas workloads, including traditionally
time-consuming tasks like ingesting data (``read_csv``, ``read_excel``,
``read_parquet``, etc.).

```python
# import pandas as pd
import modin.pandas as pd

df = pd.read_csv("big.csv")  # use all your cores!
```

### [Pandarallel](https://github.com/nalepae/pandarallel)

Pandarallel provides a simple way to parallelize your pandas operations on all your CPUs by changing only one line of code.
If also displays progress bars.

```python
from pandarallel import pandarallel

pandarallel.initialize(progress_bar=True)

# df.apply(func)
df.parallel_apply(func)
```

### [Vaex](https://vaex.io/docs/)

Increasingly, packages are being built on top of pandas to address
specific needs in data preparation, analysis and visualization. Vaex is
a python library for Out-of-Core DataFrames (similar to Pandas), to
visualize and explore big tabular datasets. It can calculate statistics
such as mean, sum, count, standard deviation etc, on an N-dimensional
grid up to a billion (10^9) objects/rows per second. Visualization is
done using histograms, density plots and 3d volume rendering, allowing
interactive exploration of big data. Vaex uses memory mapping, zero
memory copy policy and lazy computations for best performance (no memory
wasted).

- ``vaex.from_pandas``
- ``vaex.to_pandas_df``

### [Hail Query](https://hail.is/)

An out-of-core, preemptible-safe, distributed, dataframe library serving
the genetics community. Hail Query ships with on-disk data formats,
in-memory data formats, an expression compiler, a query planner, and a
distributed sort algorithm all designed to accelerate queries on large
matrices of genome sequencing data.

It is often easiest to use pandas to manipulate the summary statistics or
other small aggregates produced by Hail. For this reason, Hail provides
native import to and export from pandas DataFrames:

- [`Table.from_pandas`](https://hail.is/docs/latest/hail.Table.html#hail.Table.from_pandas)
- [`Table.to_pandas`](https://hail.is/docs/latest/hail.Table.html#hail.Table.to_pandas)

## Data cleaning and validation

### [pyjanitor](https://github.com/pyjanitor-devs/pyjanitor)

Pyjanitor provides a clean API for cleaning data, using method chaining.

### [Pandera](https://pandera.readthedocs.io/en/stable/)

Pandera provides a flexible and expressive API for performing data validation on dataframes
to make data processing pipelines more readable and robust.
Dataframes contain information that pandera explicitly validates at runtime. This is useful in
production-critical data pipelines or reproducible research settings.

## Extension data types

Pandas provides an interface for defining
[extension types](https://pandas.pydata.org/docs/development/extending.html#extension-types) to extend NumPy's type system.
The following libraries implement that interface to provide types not found in NumPy or pandas,
which work well with pandas' data containers.

### [awkward-pandas](https://awkward-pandas.readthedocs.io/)

Awkward-pandas provides an extension type for storing [Awkward
Arrays](https://awkward-array.org/) inside pandas' Series and
DataFrame. It also provides an accessor for using awkward functions
on Series that are of awkward type.

### [db-dtypes](https://github.com/googleapis/python-db-dtypes-pandas)

db-dtypes provides an extension types for working with types like
DATE, TIME, and JSON from database systems. This package is used
by pandas-gbq to provide natural dtypes for BigQuery data types without
a natural numpy type.

### [Pandas-Genomics](https://pandas-genomics.readthedocs.io/en/latest/)

Pandas-Genomics provides an extension type and extension array for working
 with genomics data.  It also includes `genomics` accessors for many useful properties
 and methods related to QC and analysis of genomics data.

### [Physipandas](https://github.com/mocquin/physipandas)

Physipandas provides an extension for manipulating physical quantities
 (like scalar and numpy.ndarray) in association with a physical unit
 (like meter or joule) and additional features for integration of
 `physipy` accessors with pandas Series and Dataframe.

### [Pint-Pandas](https://github.com/hgrecco/pint-pandas)

Pint-Pandas provides an extension type for storing numeric arrays with units.
These arrays can be stored inside pandas' Series and DataFrame. Operations
between Series and DataFrame columns which use pint's extension array are then
units aware.

### [Text Extensions](https://ibm.biz/text-extensions-for-pandas)

Text Extensions for Pandas provides extension types to cover common data structures for representing natural language data, plus library integrations that convert the outputs of popular natural language processing libraries into Pandas DataFrames.

## Accessors

A directory of projects providing
[extension accessors](https://pandas.pydata.org/docs/development/extending.html#registering-custom-accessors).
This is for users to discover new accessors and for library
authors to coordinate on the namespace.

  | Library                                                              | Accessor   | Classes               |
  | -------------------------------------------------------------------- | ---------- | --------------------- |
  | [awkward-pandas](https://awkward-pandas.readthedocs.io/en/latest/)   | `ak`       | `Series`              |
  | [pdvega](https://altair-viz.github.io/pdvega/)                       | `vgplot`   | `Series`, `DataFrame` |
  | [pandas-genomics](https://pandas-genomics.readthedocs.io/en/latest/) | `genomics` | `Series`, `DataFrame` |
  | [pint-pandas](https://github.com/hgrecco/pint-pandas)                | `pint`     | `Series`, `DataFrame` |
  | [physipandas](https://github.com/mocquin/physipandas)                | `physipy`  | `Series`, `DataFrame` |
  | [composeml](https://github.com/alteryx/compose)                      | `slice`    | `DataFrame`           |
  | [gurobipy-pandas](https://github.com/Gurobi/gurobipy-pandas)         | `gppd`     | `Series`, `DataFrame` |
  | [staircase](https://www.staircase.dev/)                              | `sc`       | `Series`, `DataFrame` |
  | [woodwork](https://github.com/alteryx/woodwork)                      | `slice`    | `Series`, `DataFrame` |

## Development tools

### [pandas-stubs](https://github.com/VirtusLab/pandas-stubs)

While pandas repository is partially typed, the package itself doesn't expose this information for external use.
Install pandas-stubs to enable basic type coverage of pandas API.

Learn more by reading through these issues [14468](https://github.com/pandas-dev/pandas/issues/14468),
[26766](https://github.com/pandas-dev/pandas/issues/26766), [28142](https://github.com/pandas-dev/pandas/issues/28142).

See installation and usage instructions on the [GitHub page](https://github.com/VirtusLab/pandas-stubs).

### [Hamilton](https://github.com/dagworks-inc/hamilton)

Hamilton is a declarative dataflow framework that came out of Stitch Fix. It was designed to help one manage a
Pandas code base, specifically with respect to feature engineering for machine learning models.

It prescibes an opinionated paradigm, that ensures all code is:

* unit testable
* integration testing friendly
* documentation friendly
* transformation logic is reusable, as it is decoupled from the context of where it is used.
* integratable with runtime data quality checks.

This helps one to scale your pandas code base, at the same time, keeping maintenance costs low.

For more information, see [documentation](https://hamilton.readthedocs.io/).
