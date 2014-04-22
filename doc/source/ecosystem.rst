.. _ecosystem:

****************
Pandas Ecosystem
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

`Statsmodels <http://statsmodels.sourceforge.net>`__
----------------------------------------------------

Statsmodels is the prominent python "statistics and econometrics library" and it has
a long-standing special relationship with pandas. Statsmodels provides powerful statistics,
econometrics, analysis and modeling functionality that is out of pandas' scope.
Statsmodels leverages pandas objects as the underlying data container for computation.

`Vincent <https://github.com/wrobstory/vincent>`__
--------------------------------------------------

The `Vincent <https://github.com/wrobstory/vincent>`__ project leverages `Vega <https://github.com/trifacta/vega>`__
(that in turn, leverages `d3 <http://d3js.org/>`__) to create plots . It has great support
for pandas data objects.

`yhat/ggplot <https://github.com/yhat/ggplot>`__
------------------------------------------------

Hadley Wickham's `ggplot2 <http://ggplot2.org/>`__ is a foundational exploratory visualization package for the R language.
Based on `"The Grammer of Graphics" <http://www.cs.uic.edu/~wilkinson/TheGrammarOfGraphics/GOG.html>`__ it
provides a powerful, declarative and extremely general way to generate bespoke plots of any kind of data.
It's really quite incredible. Various implementations to other languages are available,
but a faithful implementation for python users has long been missing. Although still young
(as of Jan-2014), the `yhat/ggplot <https://github.com/yhat/ggplot>`__ project has been
progressing quickly in that direction.


`Seaborn <https://github.com/mwaskom/seaborn>`__
------------------------------------------------

Although pandas has quite a bit of "just plot it" functionality built-in, visualization and
in particular statistical graphics is a vast field with a long tradition and lots of ground
to cover. The `Seaborn <https://github.com/mwaskom/seaborn>`__ project builds on top of pandas
and `matplotlib <http://matplotlib.org>`__ to provide easy plotting of data which extends to
more advanced types of plots then those offered by pandas.


`Geopandas <https://github.com/kjordahl/geopandas>`__
-----------------------------------------------------

Geopandas extends pandas data objects to include geographic information which support
geometric operations. If your work entails maps and geographical coordinates, and
you love pandas, you should take a close look at Geopandas.

`sklearn-pandas <https://github.com/paulgb/sklearn-pandas>`__
-------------------------------------------------------------

Use pandas DataFrames in your scikit-learn ML pipeline.

tstoolbox, wdmtoolbox, hspfbintoolbox, swmmtoolbox
--------------------------------------------------
These programs offer command line and Python libraries to manage and work with
time-series.  In the background they all use PANDAS DataFrames.  They use a
common print format at the command line and data can be piped from one command
to another for processing.  When used as Python libraries the exact same
function is used, but instead of printing returns PANDAS DataFrames.

`tstoolbox <http://pythonhosted.org//tstoolbox/>`__ is a general purpose tool
for the manipulation of time-series.

`wdmtoolbox <http://pythonhosted.org/wdmtoolbox/>`__ is a tool to create,
read, update, and delete data in Watershed Data Management (WDM) files.  WDM
files are used primarily with the Hydrological Simulation Program - FORTRAN
(HSPF).

`hspfbintoolbox <http://pythonhosted.org/hspfbintoolbox/>`__ is a tool to read
data in the HSPF binary output file.

`swmmtoolbox <http://pythonhosted.org/swmmtoolbox/>`__ is a tool to read
data in the Storm Water Management Model (SWMM) binary output file.

`QuantSoftware ToolKit <http://wiki.quantsoftware.org/>`__
----------------------------------------------------------
QSToolKit (QSTK) is a Python-based open source software framework designed to
support portfolio construction and management. We are building the QSToolKit
primarily for finance students, computing students, and quantitative analysts
with programming experience. 

`PyHIS <http://pythonhosted.org/pyhis/>`__
------------------------------------------
PyHIS (Python Hydrologic Information System) is a python module that enables
retrieval of time series water data from WaterOneFlow/WaterML web services
that are part of the national CUAHSI-HIS system.

`Ramp - Rapid Machine Learning Prototyping <http://ramp.readthedocs.org/>`__
----------------------------------------------------------------------------
Ramp is a python package for rapid machine learning prototyping. It provides a
simple, declarative syntax for exploring features, algorithms and
transformations quickly and efficiently. At its core it’s a unified
pandas-based framework for working with existing python machine learning and
statistics libraries (scikit-learn, rpy2, etc.).

`Zipline <https://github.com/quantopian/zipline>`__
---------------------------------------------------
Zipline is a Pythonic algorithmic trading library. The system is fundamentally
event-driven and a close approximation of how live-trading systems operate. 

Zipline is currently used in production as the backtesting engine powering
Quantopian (https://www.quantopian.com) -- a free, community-centered platform
that allows development and real-time backtesting of trading algorithms in the
web browser.
