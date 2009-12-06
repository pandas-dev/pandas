"""
Pandas - a library for panel, time series, or cross-sectional data analysis
===========================================================================

Main data structures (see docstrings for detailed documentation)
--------------------
Index
    Represent row or column labels in Series / DataFrame structures

Series / TimeSeries
    Represents standard 1-dimensional cross-section (resp. time series)
    As an numpy.ndarray subclass, compatible with ufuncs and other NumPy
    functions

DataFrame / DataMatrix
    Represent collections of Series objects, enable easy management
    of multiple time series / cross-sections

DateRange
    Index subclass for generating arrays of fixed frequency dates

Subpackages
-----------
core
    Implementations of core data structures, basic building blocks. Most of
    the user-relevant code is accessible through the top-level namespace
io
    Persistence, parsing, and data loading tools
lib
    C, Cython, and Fortran extensions for other components
stats
    Statistical and econometric functions
"""

