.. _ecosystem:

****************
pandas Ecosystem
****************

Increasingly, packages are being built on top of pandas to address specific needs
in data preparation, analysis and visualization.
This is encouraging because it means pandas is not only helping users to handle
their data tasks but also that it provides a better starting point for developers to
build powerful and more focused data tools.
The creation of libraries that complement pandas' functionality also allows pandas
development to remain focused around it's original requirements.

This is an in-exhaustive list of projects that build on pandas in order to provide
tools in the PyData space.

We'd like to make it easier for users to find these project, if you know of other
substantial projects that you feel should be on this list, please let us know.


.. _ecosystem.stats:

Statistics and Machine Learning
-------------------------------

`Statsmodels <http://www.statsmodels.org/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Statsmodels is the prominent python "statistics and econometrics library" and it has
a long-standing special relationship with pandas. Statsmodels provides powerful statistics,
econometrics, analysis and modeling functionality that is out of pandas' scope.
Statsmodels leverages pandas objects as the underlying data container for computation.

`sklearn-pandas <https://github.com/paulgb/sklearn-pandas>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use pandas DataFrames in your `scikit-learn <http://scikit-learn.org/>`__
ML pipeline.



.. _ecosystem.visualization:

Visualization
-------------

`Bokeh <http://bokeh.pydata.org>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bokeh is a Python interactive visualization library for large datasets that natively uses
the latest web technologies. Its goal is to provide elegant, concise construction of novel
graphics in the style of Protovis/D3, while delivering high-performance interactivity over
large data to thin clients.

`seaborn <https://seaborn.pydata.org>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Seaborn is a Python visualization library based on `matplotlib
<http://matplotlib.org>`__.  It provides a high-level, dataset-oriented
interface for creating attractive statistical graphics. The plotting functions
in seaborn understand pandas objects and leverage pandas grouping operations
internally to support concise specification of complex visualizations. Seaborn
also goes beyond matplotlib and pandas with the option to perform statistical
estimation while plotting, aggregating across observations and visualizing the
fit of statistical models to emphasize patterns in a dataset.

`yhat/ggplot <https://github.com/yhat/ggplot>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hadley Wickham's `ggplot2 <http://ggplot2.org/>`__ is a foundational exploratory visualization package for the R language.
Based on `"The Grammar of Graphics" <http://www.cs.uic.edu/~wilkinson/TheGrammarOfGraphics/GOG.html>`__ it
provides a powerful, declarative and extremely general way to generate bespoke plots of any kind of data.
It's really quite incredible. Various implementations to other languages are available,
but a faithful implementation for python users has long been missing. Although still young
(as of Jan-2014), the `yhat/ggplot <https://github.com/yhat/ggplot>`__ project has been
progressing quickly in that direction.

`Vincent <https://github.com/wrobstory/vincent>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Vincent <https://github.com/wrobstory/vincent>`__ project leverages `Vega <https://github.com/trifacta/vega>`__
(that in turn, leverages `d3 <http://d3js.org/>`__) to create
plots. Although functional, as of Summer 2016 the Vincent project has not been updated
in over two years and is `unlikely to receive further updates <https://github.com/wrobstory/vincent#2015-08-12-update>`__.

`IPython Vega <https://github.com/vega/ipyvega>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like Vincent, the `IPython Vega <https://github.com/vega/ipyvega>`__ project leverages `Vega
<https://github.com/trifacta/vega>`__ to create plots, but primarily
targets the IPython Notebook environment.

`Plotly <https://plot.ly/python>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Plotly’s <https://plot.ly/>`__ `Python API <https://plot.ly/python/>`__ enables interactive figures and web shareability. Maps, 2D, 3D, and live-streaming graphs are rendered with WebGL and `D3.js <http://d3js.org/>`__. The library supports plotting directly from a pandas DataFrame and cloud-based collaboration. Users of `matplotlib, ggplot for Python, and Seaborn <https://plot.ly/python/matplotlib-to-plotly-tutorial/>`__ can convert figures into interactive web-based plots. Plots can be drawn in `IPython Notebooks <https://plot.ly/ipython-notebooks/>`__ , edited with R or MATLAB, modified in a GUI, or embedded in apps and dashboards. Plotly is free for unlimited sharing, and has `cloud <https://plot.ly/product/plans/>`__, `offline <https://plot.ly/python/offline/>`__, or `on-premise <https://plot.ly/product/enterprise/>`__ accounts for private use.

`QtPandas <https://github.com/draperjames/qtpandas>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spun off from the main pandas library, the `qtpandas <https://github.com/draperjames/qtpandas>`__
library enables DataFrame visualization and manipulation in PyQt4 and PySide applications.


.. _ecosystem.ide:

IDE
------

`IPython <http://ipython.org/documentation.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IPython is an interactive command shell and distributed computing
environment.
IPython Notebook is a web application for creating IPython notebooks.
An IPython notebook is a JSON document containing an ordered list
of input/output cells which can contain code, text, mathematics, plots
and rich media.
IPython notebooks can be converted to a number of open standard output formats
(HTML, HTML presentation slides, LaTeX, PDF, ReStructuredText, Markdown,
Python) through 'Download As' in the web interface and ``ipython nbconvert``
in a shell.

Pandas DataFrames implement ``_repr_html_`` methods
which are utilized by IPython Notebook for displaying
(abbreviated) HTML tables.  (Note: HTML tables may or may not be
compatible with non-HTML IPython output formats.)

`quantopian/qgrid <https://github.com/quantopian/qgrid>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qgrid is "an interactive grid for sorting and filtering
DataFrames in IPython Notebook" built with SlickGrid.

`Spyder <https://github.com/spyder-ide/spyder/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spyder is a cross-platform Qt-based open-source Python IDE with
editing, testing, debugging, and introspection features.
Spyder can now introspect and display Pandas DataFrames and show
both "column wise min/max and global min/max coloring."


.. _ecosystem.api:

API
-----

`pandas-datareader <https://github.com/pydata/pandas-datareader>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``pandas-datareader`` is a remote data access library for pandas (PyPI:``pandas-datareader``).
It is based on functionality that was located in ``pandas.io.data`` and ``pandas.io.wb`` but was
split off in v0.19.
See more in the  `pandas-datareader docs <https://pandas-datareader.readthedocs.io/en/latest/>`_:

The following data feeds are available:

  * Yahoo! Finance
  * Google Finance
  * FRED
  * Fama/French
  * World Bank
  * OECD
  * Eurostat
  * EDGAR Index

`quandl/Python <https://github.com/quandl/Python>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Quandl API for Python wraps the Quandl REST API to return
Pandas DataFrames with timeseries indexes.

`pydatastream <https://github.com/vfilimonov/pydatastream>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PyDatastream is a Python interface to the
`Thomson Dataworks Enterprise (DWE/Datastream) <http://dataworks.thomson.com/Dataworks/Enterprise/1.0/>`__
SOAP API to return indexed Pandas DataFrames or Panels with financial data.
This package requires valid credentials for this API (non free).

`pandaSDMX <https://pandasdmx.readthedocs.io>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pandaSDMX is a library to retrieve and acquire statistical data
and metadata disseminated in
`SDMX <http://www.sdmx.org>`_ 2.1, an ISO-standard
widely used by institutions such as statistics offices, central banks,   
and international organisations. pandaSDMX can expose datasets and related 
structural metadata including dataflows, code-lists, 
and datastructure definitions as pandas Series 
or multi-indexed DataFrames.  
   
`fredapi <https://github.com/mortada/fredapi>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fredapi is a Python interface to the `Federal Reserve Economic Data (FRED) <http://research.stlouisfed.org/fred2/>`__
provided by the Federal Reserve Bank of St. Louis. It works with both the FRED database and ALFRED database that
contains point-in-time data (i.e. historic data revisions). fredapi provides a wrapper in python to the FRED
HTTP API, and also provides several convenient methods for parsing and analyzing point-in-time data from ALFRED.
fredapi makes use of pandas and returns data in a Series or DataFrame. This module requires a FRED API key that
you can obtain for free on the FRED website.


.. _ecosystem.domain:

Domain Specific
---------------

`Geopandas <https://github.com/kjordahl/geopandas>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geopandas extends pandas data objects to include geographic information which support
geometric operations. If your work entails maps and geographical coordinates, and
you love pandas, you should take a close look at Geopandas.

`xarray <https://github.com/pydata/xarray>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xarray brings the labeled data power of pandas to the physical sciences by
providing N-dimensional variants of the core pandas data structures. It aims to
provide a pandas-like and pandas-compatible toolkit for analytics on multi-
dimensional arrays, rather than the tabular data for which pandas excels.


.. _ecosystem.out-of-core:

Out-of-core
-------------

`Dask <https://dask.readthedocs.io/en/latest/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask is a flexible parallel computing library for analytics. Dask
provides a familiar ``DataFrame`` interface for out-of-core, parallel and distributed computing.

`Dask-ML <https://dask-ml.readthedocs.io/en/latest/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask-ML enables parallel and distributed machine learning using Dask alongside existing machine learning libraries like Scikit-Learn, XGBoost, and TensorFlow.


`Blaze <http://blaze.pydata.org/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Blaze provides a standard API for doing computations with various
in-memory and on-disk backends: NumPy, Pandas, SQLAlchemy, MongoDB, PyTables,
PySpark.

`Odo <http://odo.pydata.org>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Odo provides a uniform API for moving data between different formats. It uses
pandas own ``read_csv`` for CSV IO and leverages many existing packages such as
PyTables, h5py, and pymongo to move data between non pandas formats. Its graph
based approach is also extensible by end users for custom formats that may be
too specific for the core of odo.

.. _ecosystem.data_validation:

Data validation
---------------

`Engarde <http://engarde.readthedocs.io/en/latest/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Engarde is a lightweight library used to explicitly state your assumptions abour your datasets
and check that they're *actually* true.
