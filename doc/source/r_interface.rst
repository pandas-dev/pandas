.. currentmodule:: pandas.rpy

.. _rpy:

******************
rpy2 / R interface
******************

.. note::

   This is all highly experimental. I would like to get more people involved
   with building a nice RPy2 interface for pandas


If your computer has R and rpy2 (> 2.2) installed (which will be left to the
reader), you will be able to leverage the below functionality. On Windows,
doing this is quite an ordeal at the moment, but users on Unix-like systems
should find it quite easy. rpy2 evolves in time, and is currently reaching 
its release 2.3, while the current interface is
designed for the 2.2.x series. We recommend to use 2.2.x over other series 
unless you are prepared to fix parts of the code, yet the rpy2-2.3.0
introduces improvements such as a better R-Python bridge memory management
layer so I might be a good idea to bite the bullet and submit patches for
the few minor differences that need to be fixed.


::

    # if installing for the first time
    hg clone http://bitbucket.org/lgautier/rpy2

    cd rpy2
    hg pull
    hg update version_2.2.x
    sudo python setup.py install

.. note::

    To use R packages with this interface, you will need to install
    them inside R yourself. At the moment it cannot install them for
    you.

Once you have done installed R and rpy2, you should be able to import
``pandas.rpy.common`` without a hitch.

Transferring R data sets into Python
------------------------------------

The **load_data** function retrieves an R data set and converts it to the
appropriate pandas object (most likely a DataFrame):


.. ipython:: python

   import pandas.rpy.common as com
   infert = com.load_data('infert')

   infert.head()


Converting DataFrames into R objects
------------------------------------

.. versionadded:: 0.8

Starting from pandas 0.8, there is **experimental** support to convert
DataFrames into the equivalent R object (that is, **data.frame**):

.. ipython:: python

   from pandas import DataFrame

   df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6], 'C':[7,8,9]},
                  index=["one", "two", "three"])
   r_dataframe = com.convert_to_r_dataframe(df)

   print type(r_dataframe)
   print r_dataframe


The DataFrame's index is stored as the ``rownames`` attribute of the
data.frame instance.

You can also use **convert_to_r_matrix** to obtain a ``Matrix`` instance, but
bear in mind that it will only work with homogeneously-typed DataFrames (as
R matrices bear no information on the data type):


.. ipython:: python

   r_matrix = com.convert_to_r_matrix(df)

   print type(r_matrix)
   print r_matrix


Calling R functions with pandas objects
---------------------------------------



High-level interface to R estimators
------------------------------------
