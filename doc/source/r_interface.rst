.. _rpy:

.. ipython:: python
   :suppress:

   import pandas as pd
   pd.options.display.max_rows = 15


******************
rpy2 / R interface
******************

.. warning::

    Up to pandas 0.19, a ``pandas.rpy`` module existed with functionality to
    convert between pandas and ``rpy2`` objects. This functionality now lives in
    the `rpy2 <https://rpy2.readthedocs.io/>`__ project itself.
    See the `updating section <http://pandas.pydata.org/pandas-docs/version/0.19.0/r_interface.html#updating-your-code-to-use-rpy2-functions>`__
    of the previous documentation for a guide to port your code from the
    removed ``pandas.rpy`` to ``rpy2`` functions.


`rpy2 <http://rpy2.bitbucket.org/>`__ is an interface to R running embedded in a Python process, and also includes functionality to deal with pandas DataFrames.
Converting data frames back and forth between rpy2 and pandas should be largely
automated (no need to convert explicitly, it will be done on the fly in most
rpy2 functions).
To convert explicitly, the functions are ``pandas2ri.py2ri()`` and
``pandas2ri.ri2py()``.


See also the documentation of the `rpy2 <http://rpy2.bitbucket.org/>`__ project: https://rpy2.readthedocs.io.

In the remainder of this page, a few examples of explicit conversion is given. The pandas conversion of rpy2 needs first to be activated:

.. ipython::
    :verbatim:

    In [1]: from rpy2.robjects import pandas2ri
       ...: pandas2ri.activate()

Transferring R data sets into Python
------------------------------------

Once the pandas conversion is activated (``pandas2ri.activate()``), many conversions
of R to pandas objects will be done automatically. For example, to obtain the 'iris' dataset as a pandas DataFrame:

.. ipython::
    :verbatim:

    In [2]: from rpy2.robjects import r

    In [3]: r.data('iris')

    In [4]: r['iris'].head()
    Out[4]:
        Sepal.Length  Sepal.Width  Petal.Length  Petal.Width Species
    0           5.1          3.5           1.4          0.2  setosa
    1           4.9          3.0           1.4          0.2  setosa
    2           4.7          3.2           1.3          0.2  setosa
    3           4.6          3.1           1.5          0.2  setosa
    4           5.0          3.6           1.4          0.2  setosa

If the pandas conversion was not activated, the above could also be accomplished
by explicitly converting it with the ``pandas2ri.ri2py`` function
(``pandas2ri.ri2py(r['iris'])``).

Converting DataFrames into R objects
------------------------------------

The ``pandas2ri.py2ri`` function support the reverse operation to convert
DataFrames into the equivalent R object (that is, **data.frame**):

.. ipython::
   :verbatim:

   In [5]: df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]},
      ...:                   index=["one", "two", "three"])

   In [6]: r_dataframe = pandas2ri.py2ri(df)

   In [7]: print(type(r_dataframe))
   Out[7]: <class 'rpy2.robjects.vectors.DataFrame'>

   In [8]: print(r_dataframe)
   Out[8]:
         A B C
   one   1 4 7
   two   2 5 8
   three 3 6 9


The DataFrame's index is stored as the ``rownames`` attribute of the
data.frame instance.


..
   Calling R functions with pandas objects
   High-level interface to R estimators
