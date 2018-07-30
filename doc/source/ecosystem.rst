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

This is an inexhaustive list of projects that build on pandas in order to provide
tools in the PyData space. For a list of projects that depend on pandas,
see the 
`libraries.io usage page for pandas <https://libraries.io/pypi/pandas/usage>`_
or `search pypi for pandas <https://pypi.org/search/?q=pandas>`_.

We'd like to make it easier for users to find these projects, if you know of other
substantial projects that you feel should be on this list, please let us know.


.. _ecosystem.stats:

Statistics and Machine Learning
-------------------------------

`Statsmodels <http://www.statsmodels.org/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Statsmodels is the prominent Python "statistics and econometrics library" and it has
a long-standing special relationship with pandas. Statsmodels provides powerful statistics,
econometrics, analysis and modeling functionality that is out of pandas' scope.
Statsmodels leverages pandas objects as the underlying data container for computation.

`sklearn-pandas <https://github.com/paulgb/sklearn-pandas>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use pandas DataFrames in your `scikit-learn <http://scikit-learn.org/>`__
ML pipeline.

`Featuretools <https://github.com/featuretools/featuretools/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Featuretools is a Python library for automated feature engineering built on top of pandas. It excels at transforming temporal and relational datasets into feature matrices for machine learning using reusable feature engineering "primitives". Users can contribute their own primitives in Python and share them with the rest of the community. 

.. _ecosystem.visualization:

Visualization
-------------

`Altair <https://altair-viz.github.io/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Altair is a declarative statistical visualization library for Python.
With Altair, you can spend more time understanding your data and its
meaning. Altair's API is simple, friendly and consistent and built on
top of the powerful Vega-Lite JSON specification. This elegant
simplicity produces beautiful and effective visualizations with a
minimal amount of code. Altair works with Pandas DataFrames.


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

`yhat/ggpy <https://github.com/yhat/ggpy>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hadley Wickham's `ggplot2 <http://ggplot2.org/>`__ is a foundational exploratory visualization package for the R language.
Based on `"The Grammar of Graphics" <http://www.cs.uic.edu/~wilkinson/TheGrammarOfGraphics/GOG.html>`__ it
provides a powerful, declarative and extremely general way to generate bespoke plots of any kind of data.
It's really quite incredible. Various implementations to other languages are available,
but a faithful implementation for Python users has long been missing. Although still young
(as of Jan-2014), the `yhat/ggpy <https://github.com/yhat/ggpy>`__ project has been
progressing quickly in that direction.

`IPython Vega <https://github.com/vega/ipyvega>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`IPython Vega <https://github.com/vega/ipyvega>`__ leverages `Vega
<https://github.com/trifacta/vega>`__ to create plots within Jupyter Notebook.

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
environment. IPython tab completion works with Pandas methods and also
attributes like DataFrame columns.

`Jupyter Notebook / Jupyter Lab <https://jupyter.org>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Jupyter Notebook is a web application for creating Jupyter notebooks.
A Jupyter notebook is a JSON document containing an ordered list
of input/output cells which can contain code, text, mathematics, plots
and rich media.
Jupyter notebooks can be converted to a number of open standard output formats
(HTML, HTML presentation slides, LaTeX, PDF, ReStructuredText, Markdown,
Python) through 'Download As' in the web interface and ``jupyter convert``
in a shell.

Pandas DataFrames implement ``_repr_html_``and ``_repr_latex`` methods
which are utilized by Jupyter Notebook for displaying
(abbreviated) HTML or LaTeX tables. LaTeX output is properly escaped.
(Note: HTML tables may or may not be
compatible with non-HTML Jupyter output formats.)

See :ref:`Options and Settings <options>` and :ref:`<options.available>`
for pandas ``display.`` settings.

`quantopian/qgrid <https://github.com/quantopian/qgrid>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qgrid is "an interactive grid for sorting and filtering
DataFrames in IPython Notebook" built with SlickGrid.

`Spyder <https://www.spyder-ide.org/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spyder is a cross-platform PyQt-based IDE combining the editing, analysis,
debugging and profiling functionality of a software development tool with the
data exploration, interactive execution, deep inspection and rich visualization
capabilities of a scientific environment like MATLAB or Rstudio.

Its `Variable Explorer <https://docs.spyder-ide.org/variableexplorer.html>`__
allows users to view, manipulate and edit pandas Index, DateTimeIndex, Series,
TimeSeries and DataFrame objects with a full GUI, including the ability to
sort by column, copy and modify individual values, display a "heatmap"
per-column or globally, convert column data types and more, all with user- or
automatically-adjustable column-widths and display formatting.
The Variable Explorer also allows renaming and duplicating pandas objects,
adding new columns, copying them to the clipboard as plain text (TSV), pasting
them into  different Spyder/IPython sessions, and saving and loading
one or multiple in a session to and from a file in lossless native format.
Spyder can import data from a variety of plain text and binary files
or the clipboard into a new pandas DataFrame, via a sophisticated import wizard
offering options to handle column and row separators, skipping rows, and
comments, and the ability to transpose, type-convert and preview the results.

Most pandas classes, methods and data attributes can be autocompleted in
Spyder's `Editor <https://docs.spyder-ide.org/editor.html>`__ and
`IPython Console <https://docs.spyder-ide.org/ipythonconsole.html>__,
and Spyder's `Help pane<https://docs.spyder-ide.org/help.html>__ can retrieve
and render Numpydoc documentation on pandas objects in rich text with Sphinx
both automatically and on-demand.


.. _ecosystem.api:

API
---

`pandas-datareader <https://github.com/pydata/pandas-datareader>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``pandas-datareader`` is a remote data access library for pandas (PyPI:``pandas-datareader``).
It is based on functionality that was located in ``pandas.io.data`` and ``pandas.io.wb`` but was
split off in v0.19.
See more in the  `pandas-datareader docs <https://pandas-datareader.readthedocs.io/en/latest/>`_:

The following data feeds are available:

 * Google Finance
 * Tiingo
 * Morningstar
 * IEX
 * Robinhood
 * Enigma
 * Quandl
 * FRED
 * Fama/French
 * World Bank
 * OECD
 * Eurostat
 * TSP Fund Data
 * Nasdaq Trader Symbol Definitions
 * Stooq Index Data
 * MOEX Data

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
structural metadata including data flows, code-lists,
and data structure definitions as pandas Series
or MultiIndexed DataFrames.
   
`fredapi <https://github.com/mortada/fredapi>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fredapi is a Python interface to the `Federal Reserve Economic Data (FRED) <http://research.stlouisfed.org/fred2/>`__
provided by the Federal Reserve Bank of St. Louis. It works with both the FRED database and ALFRED database that
contains point-in-time data (i.e. historic data revisions). fredapi provides a wrapper in Python to the FRED
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

`Blaze <http://blaze.pydata.org/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Blaze provides a standard API for doing computations with various
in-memory and on-disk backends: NumPy, Pandas, SQLAlchemy, MongoDB, PyTables,
PySpark.

`Dask <https://dask.readthedocs.io/en/latest/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask is a flexible parallel computing library for analytics. Dask
provides a familiar ``DataFrame`` interface for out-of-core, parallel and distributed computing.

`Dask-ML <https://dask-ml.readthedocs.io/en/latest/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask-ML enables parallel and distributed machine learning using Dask alongside existing machine learning libraries like Scikit-Learn, XGBoost, and TensorFlow.

`Odo <http://odo.pydata.org>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Odo provides a uniform API for moving data between different formats. It uses
pandas own ``read_csv`` for CSV IO and leverages many existing packages such as
PyTables, h5py, and pymongo to move data between non pandas formats. Its graph
based approach is also extensible by end users for custom formats that may be
too specific for the core of odo.

`Ray <https://ray.readthedocs.io/en/latest/pandas_on_ray.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pandas on Ray is an early stage DataFrame library that wraps Pandas and transparently distributes the data and computation. The user does not need to know how many cores their system has, nor do they need to specify how to distribute the data. In fact, users can continue using their previous Pandas notebooks while experiencing a considerable speedup from Pandas on Ray, even on a single machine. Only a modification of the import statement is needed, as we demonstrate below. Once you’ve changed your import statement, you’re ready to use Pandas on Ray just like you would Pandas.

.. code:: python

    # import pandas as pd
    import ray.dataframe as pd


`Vaex <https://docs.vaex.io/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Increasingly, packages are being built on top of pandas to address specific needs in data preparation, analysis and visualization. Vaex is a python library for Out-of-Core DataFrames (similar to Pandas), to visualize and explore big tabular datasets. It can calculate statistics such as mean, sum, count, standard deviation etc, on an N-dimensional grid up to a billion (10\ :sup:`9`) objects/rows per second. Visualization is done using histograms, density plots and 3d volume rendering, allowing interactive exploration of big data. Vaex uses memory mapping, zero memory copy policy and lazy computations for best performance (no memory wasted).

 * vaex.from_pandas
 * vaex.to_pandas_df


.. _ecosystem.data_validation:

Data validation
---------------

`Engarde <http://engarde.readthedocs.io/en/latest/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Engarde is a lightweight library used to explicitly state your assumptions about your datasets
and check that they're *actually* true.

.. _ecosystem.extensions:

Extension Data Types
--------------------

Pandas provides an interface for defining
:ref:`extension types <extending.extension-types>` to extend NumPy's type
system. The following libraries implement that interface to provide types not
found in NumPy or pandas, which work well with pandas' data containers.

`cyberpandas`_
~~~~~~~~~~~~~~

Cyberpandas provides an extension type for storing arrays of IP Addresses. These
arrays can be stored inside pandas' Series and DataFrame.

.. _ecosystem.accessors:

Accessors
---------

A directory of projects providing
:ref:`extension accessors <extending.register-accessors>`. This is for users to
discover new accessors and for library authors to coordinate on the namespace.

============== ========== =========================
Library        Accessor   Classes
============== ========== =========================
`cyberpandas`_ ``ip``     ``Series``
`pdvega`_      ``vgplot`` ``Series``, ``DataFrame``
============== ========== =========================

.. _cyberpandas: https://cyberpandas.readthedocs.io/en/latest
.. _pdvega: https://jakevdp.github.io/pdvega/
