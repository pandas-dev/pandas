.. Pandas documentation master file, created by


pandas: a python data analysis library
======================================

:mod:`pandas` is a python package providing convenient data structures
for time series, cross-sectional, or any other form of "labeled" data,
with tools for building statistical and econometric models.

This library was created with the following design principles:

  - Working with time series and cross-sectional data should be easy
  - The user should not have to worry (much) about handling missing data
  - Data alignment should be automatic and transparent
  - Speed matters
  - Perhaps most importantly: *things should work just like you want them to*

Many of these principles are here to address the shortcomings
frequently experienced using other languages / scientific research
environments. In MATLAB, for example, you spend a lot of time coercing
data into matrices, cleaning and aligning it, and keeping everything
homogeneous. You have to use lots of functions like **nanmean, nanstd,
repmat** (for broadcasting), and other functions which help you to
maintain reliable data. Using `NumPy <http://www.numpy.org>`__ and a
Pythonic approach, pandas helps hide the dirty details of working with
unclean data, allowing you to focus on the problem you're trying to
solve rather than the implementation.

pandas is implemented primarily using NumPy and is intended to be able
to integrate very easily with other NumPy-based scientific libraries,
such as :mod:`scikits.statsmodels`.

.. note::

   This documentation assumes general familiarity with NumPy. If you
   haven't used NumPy much or at all, please check out the `NumPy
   documentation <http://docs.scipy.org>`__ first.

See the package overview for more detail about what's in the library.

User manual
-----------

.. module:: pandas

**Date**: |today|

**Version**: |version|

**License:** BSD

**Requirements:** python 2.4 to 2.6, NumPy, and dateutil

**Code Repository:** http://pandas.googlecode.com

.. toctree::
    :maxdepth: 2

    overview
    core
    groupby
    datetools
    stats
    r_interface
    io
    examples
    missing_data
    related

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Contact
-------

Please feel free to send comments or questions directly to
wesmckinn@gmail.com or the pystatsmodels mailing list.
