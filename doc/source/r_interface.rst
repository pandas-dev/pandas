.. currentmodule:: pandas.rpy

.. _rpy:

.. ipython:: python
   :suppress:

   import pandas as pd
   pd.options.display.max_rows = 15


******************
rpy2 / R interface
******************

.. warning::

   In v0.16.0, the ``pandas.rpy`` interface has been **deprecated and will be
   removed in a future version**. Similar functionality can be accessed
   through the `rpy2 <https://rpy2.readthedocs.io/>`__ project.
   See the :ref:`updating <rpy.updating>` section for a guide to port your
   code from the ``pandas.rpy`` to ``rpy2`` functions.


.. _rpy.updating:

Updating your code to use rpy2 functions
----------------------------------------

In v0.16.0, the ``pandas.rpy`` module has been **deprecated** and users are
pointed to the similar functionality in ``rpy2`` itself (rpy2 >= 2.4).

Instead of importing ``import pandas.rpy.common as com``, the following imports
should be done to activate the pandas conversion support in rpy2::

    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

Converting data frames back and forth between rpy2 and pandas should be largely
automated (no need to convert explicitly, it will be done on the fly in most
rpy2 functions).

To convert explicitly, the functions are ``pandas2ri.py2ri()`` and
``pandas2ri.ri2py()``. So these functions can be used to replace the existing
functions in pandas:

- ``com.convert_to_r_dataframe(df)`` should be replaced with ``pandas2ri.py2ri(df)``
- ``com.convert_robj(rdf)`` should be replaced with ``pandas2ri.ri2py(rdf)``

Note: these functions are for the latest version (rpy2 2.5.x) and were called
``pandas2ri.pandas2ri()`` and ``pandas2ri.ri2pandas()`` previously.

Some of the other functionality in `pandas.rpy` can be replaced easily as well.
For example to load R data as done with the ``load_data`` function, the
current method::

    df_iris = com.load_data('iris')

can be replaced with::

    from rpy2.robjects import r
    r.data('iris')
    df_iris = pandas2ri.ri2py(r[name])

The ``convert_to_r_matrix`` function can be replaced by the normal
``pandas2ri.py2ri`` to convert dataframes, with a subsequent call to R
``as.matrix`` function.

.. warning::

    Not all conversion functions in rpy2 are working exactly the same as the
    current methods in pandas. If you experience problems or limitations in
    comparison to the ones in pandas, please report this at the
    `issue tracker <https://github.com/pandas-dev/pandas/issues>`_.

See also the documentation of the `rpy2 <http://rpy2.bitbucket.org/>`__ project.


R interface with rpy2
---------------------

If your computer has R and rpy2 (> 2.2) installed (which will be left to the
reader), you will be able to leverage the below functionality. On Windows,
doing this is quite an ordeal at the moment, but users on Unix-like systems
should find it quite easy. rpy2 evolves in time, and is currently reaching
its release 2.3, while the current interface is
designed for the 2.2.x series. We recommend to use 2.2.x over other series
unless you are prepared to fix parts of the code, yet the rpy2-2.3.0
introduces improvements such as a better R-Python bridge memory management
layer so it might be a good idea to bite the bullet and submit patches for
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
   :okwarning:

   import pandas.rpy.common as com
   infert = com.load_data('infert')

   infert.head()


Converting DataFrames into R objects
------------------------------------

.. versionadded:: 0.8

Starting from pandas 0.8, there is **experimental** support to convert
DataFrames into the equivalent R object (that is, **data.frame**):

.. ipython:: python

   import pandas.rpy.common as com
   df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6], 'C':[7,8,9]},
                     index=["one", "two", "three"])
   r_dataframe = com.convert_to_r_dataframe(df)

   print(type(r_dataframe))
   print(r_dataframe)


The DataFrame's index is stored as the ``rownames`` attribute of the
data.frame instance.

You can also use **convert_to_r_matrix** to obtain a ``Matrix`` instance, but
bear in mind that it will only work with homogeneously-typed DataFrames (as
R matrices bear no information on the data type):


.. ipython:: python

   import pandas.rpy.common as com
   r_matrix = com.convert_to_r_matrix(df)

   print(type(r_matrix))
   print(r_matrix)


Calling R functions with pandas objects
---------------------------------------



High-level interface to R estimators
------------------------------------
