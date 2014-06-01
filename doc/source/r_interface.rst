.. currentmodule:: pandas.rpy

.. _rpy:

.. ipython:: python
   :suppress:

   from pandas import *
   import numpy as np
   np.random.seed(123456)
   import matplotlib.pyplot as plt
   plt.close('all')
   options.display.mpl_style = 'default'
   options.display.max_rows=15


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

   import pandas.rpy.common as com
   infert = com.load_data('infert')

   infert.head()


Converting DataFrames into R objects
------------------------------------

.. versionadded:: 0.8

Starting from pandas 0.8, there is **experimental** support to convert
``DataFrame`` into the equivalent R object (that is, **data.frame**) using ``convert_to_r_dataframe`` function:


.. ipython:: python

   df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6], 'C':[7,8,9]},
                  index=["one", "two", "three"])
   r_dataframe = com.convert_to_r_dataframe(df)

   type(r_dataframe)
   print(r_dataframe)

   print(r_dataframe.rownames)
   print(r_dataframe.colnames)


The ``rpy2.robjects.vectors.DataFrame`` index is stored as the ``rownames``,  and columns are stored as the
``colnames`` attributes.

You can also use ``convert_to_r_matrix`` to obtain a ``rpy2.robjects.vectors.Matrix`` instance, but
bear in mind that it will only work with homogeneously-typed DataFrames (as
R matrices bear no information on the data type):


.. ipython:: python

   r_matrix = com.convert_to_r_matrix(df)

   type(r_matrix)
   print(r_matrix)


Calling R functions with pandas objects
---------------------------------------

It is easier to use ``rpy2.robjects`` directly to call R functions.
You can retrieve R object (including R function) from R namespace by dictionary access of ``robjects.r``.

Below example shows to retrieve R's **sum** function and pass ``rpy2.robjects.vector.DataFrame``.
Note that the returned value from R **sum** is stored as ``robjects.vectors.Vectors`` type.
Thus, specify index to get raw values.

See `RPy2 documentation <http://rpy.sourceforge.net/rpy2/doc-2.2/html/index.html>`__ for more.


.. ipython:: python

   import rpy2.robjects as robjects

   rsum = robjects.r['sum']
   rsum_result = rsum(r_dataframe)

   type(rsum_result)
   rsum_result[0]


Preparing Data for R
--------------------

Load Iris dataset and convert it to R **data.frame**.
You can pass ``rpy2.robjects.vectors.DataFrame`` to R namespace using ``rpy2.robjects.r.assign``.
In following examle, `r_iris` DataFrame can be refered as `iris` on R namespace.


.. ipython:: python

   iris = com.load_data('iris')
   iris.head()

   r_iris = com.convert_to_r_dataframe(iris)
   robjects.r.assign('iris', r_iris);


You can convert each data type using R functions if required.
Function calling ``objects.r`` will execure a passed formula on R's namespace.
For example, we can check the data type using R's **str** function,
then convert "Species" column to categorical type (Factor) using R's **factor** function.


.. ipython:: python

   print(robjects.r('str(iris)'))

   robjects.r('iris$Species <- factor(iris$Species)');
   print(robjects.r('str(iris)'))


High-level interface to R estimators
------------------------------------

Use "setosa" data in iris data set to perform Linear Regression.
It is much easier to prepare and slice data on pandas side, then convert it to R **data.frame**.


.. ipython:: python

   setosa = iris[iris['Species'] == 'setosa']
   setosa.head()

   r_setosa = com.convert_to_r_dataframe(setosa)
   robjects.r.assign('setosa', r_setosa);


Once DataFrame is passed to R namespace, you can execute R formula to perform Liner Regression.


.. ipython:: python

   robjects.r('result <- lm(Sepal.Length~Sepal.Width, data=setosa)');
   print(robjects.r('summary(result)'))


You can retrieve the result from R namespace to python namespace via ``rpy2.robjects.r``.
If a returned value is R named list, you can check the list of keys via ``names`` attribute.
To get raw values, access each element specifying index.


.. ipython:: python

   result = robjects.r['result']

   print(result.names)
   print(result.rx('coefficients'))

   intercept, coef1 = result.rx('coefficients')[0]
   intercept
   coef1


``convert_robj`` function converts retrieved data to python friendly data type.
In below example, retrieved R **data.frame** of fitted values and confidence interval will be
converted to pandas ``DataFrame``.


.. ipython:: python

   robjects.r('predicted <- predict(result, setosa, interval="prediction")');
   print(robjects.r('head(predicted)'))

   predicted = robjects.r['predicted']
   type(predicted)

   predicted = com.convert_robj(predicted)
   type(predicted)
   predicted.head()


Handling Time Series
--------------------

Currently, there is no easy way to create R's built-in **ts** object from pandas time series.
Also, ``Series`` cannot be converted using ``convert_to_r_dataframe`` function.
Thus, you must create ``rpy2.robjects.vectors.Vector`` instance manually before calling ``robjects.r.assign``.

Use corresponding ``Vector`` class depending on the intended data type.
See the rpy2 documentation `Vectors and arrays <http://rpy.sourceforge.net/rpy2/doc-2.2/html/vector.html>`__ for more.

Once the ``Vector`` is passed to R's namespace, call R's **ts** function to create **ts** object.


.. ipython:: python

   idx = date_range(start='2013-01-01', freq='M', periods=48)
   vts = Series(np.random.randn(48), index=idx).cumsum()
   vts

   r_values = robjects.FloatVector(vts.values)
   robjects.r.assign('values', r_values);

   robjects.r('vts <- ts(values, start=c(2013, 1, 1), frequency=12)');
   print(robjects.r['vts'])


Below example performs Seasonal Decomposition using R's **stl** function, and get the result as `converted` ``DataFrame``.
Because R's **ts** index cannot be retrieved by ``convert_robj``, assign ``DatetimeIndex`` manually after retrieval.


.. ipython:: python

   robjects.r('result <- stl(vts, s.window=12)');
   result = robjects.r['result']

   print(result.names)

   result_ts = result.rx('time.series')[0]
   converted = com.convert_robj(result_ts)
   converted.head()

   converted.index = idx
   converted.head()


Now you have pandas ``DataFrame``, you can perform further operation easily.


.. ipython:: python

   fig, axes = plt.subplots(4, 1)

   axes[0].set_ylabel('Original');
   ax = vts.plot(ax=axes[0]);
   axes[1].set_ylabel('Trend');
   ax = converted['trend'].plot(ax=axes[1]);

   axes[2].set_ylabel('Seasonal');
   ax = converted['seasonal'].plot(ax=axes[2]);

   axes[3].set_ylabel('Residuals');
   @savefig rpy2_timeseries.png
   converted['remainder'].plot(ax=axes[3])

