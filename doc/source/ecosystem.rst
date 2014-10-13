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

`Statsmodels <http://statsmodels.sourceforge.net>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

`yhat/ggplot <https://github.com/yhat/ggplot>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hadley Wickham's `ggplot2 <http://ggplot2.org/>`__ is a foundational exploratory visualization package for the R language.
Based on `"The Grammer of Graphics" <http://www.cs.uic.edu/~wilkinson/TheGrammarOfGraphics/GOG.html>`__ it
provides a powerful, declarative and extremely general way to generate bespoke plots of any kind of data.
It's really quite incredible. Various implementations to other languages are available,
but a faithful implementation for python users has long been missing. Although still young
(as of Jan-2014), the `yhat/ggplot <https://github.com/yhat/ggplot>`__ project has been
progressing quickly in that direction.

`Seaborn <https://github.com/mwaskom/seaborn>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although pandas has quite a bit of "just plot it" functionality built-in, visualization and
in particular statistical graphics is a vast field with a long tradition and lots of ground
to cover. The `Seaborn <https://github.com/mwaskom/seaborn>`__ project builds on top of pandas
and `matplotlib <http://matplotlib.org>`__ to provide easy plotting of data which extends to
more advanced types of plots then those offered by pandas.

`Vincent <https://github.com/wrobstory/vincent>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Vincent <https://github.com/wrobstory/vincent>`__ project leverages `Vega <https://github.com/trifacta/vega>`__
(that in turn, leverages `d3 <http://d3js.org/>`__) to create plots . It has great support
for pandas data objects.


.. _ecosystem.ide:

IDE
------

`IPython <http://ipython.org/documentation.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qgrid is "an interactive grid for sorting and filtering
DataFrames in IPython Notebook" built with SlickGrid.

`Spyder <https://code.google.com/p/spyderlib/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spyder is a cross-platform Qt-based open-source Python IDE with
editing, testing, debugging, and introspection features.
Spyder can now introspect and display Pandas DataFrames and show
both "column wise min/max and global min/max coloring."


.. _ecosystem.api:

API
-----

`quandl/Python <https://github.com/quandl/Python>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Quandl API for Python wraps the Quandl REST API to return
Pandas DataFrames with timeseries indexes.


.. _ecosystem.domain:

Domain Specific
---------------

`Geopandas <https://github.com/kjordahl/geopandas>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geopandas extends pandas data objects to include geographic information which support
geometric operations. If your work entails maps and geographical coordinates, and
you love pandas, you should take a close look at Geopandas.

`xray <https://github.com/xray/xray>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xray brings the labeled data power of pandas to the physical sciences by
providing N-dimensional variants of the core pandas data structures. It aims to
provide a pandas-like and pandas-compatible toolkit for analytics on multi-
dimensional arrays, rather than the tabular data for which pandas excels.


.. _ecosystem.out-of-core:

Out-of-core
-------------

`Blaze <http://blaze.pydata.org/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Blaze provides a standard API for doing computations with various
in-memory and on-disk backends: NumPy, Pandas, SQLAlchemy, MongoDB, PyTables,
PySpark.
