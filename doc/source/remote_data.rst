.. _remote_data:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import os
   import csv
   import pandas as pd

   import numpy as np
   np.random.seed(123456)
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

   import matplotlib.pyplot as plt
   plt.close('all')

   from pandas import *
   options.display.max_rows=15
   import pandas.util.testing as tm

******************
Remote Data Access
******************

.. _remote_data.data_reader:

Functions from :mod:`pandas.io.data` extract data from various Internet
sources into a DataFrame. Currently the following sources are supported:

    - Yahoo! Finance
    - Google Finance
    - St. Louis FED (FRED)
    - Kenneth French's data library
    - World Bank

It should be noted, that various sources support different kinds of data, so not all sources implement the same methods and the data elements returned might also differ.

.. _remote_data.yahoo:

Yahoo! Finance
--------------

.. ipython:: python

    import pandas.io.data as web
    import datetime
    start = datetime.datetime(2010, 1, 1)
    end = datetime.datetime(2013, 1, 27)
    f=web.DataReader("F", 'yahoo', start, end)
    f.ix['2010-01-04']

.. _remote_data.google:

Google Finance
--------------

.. ipython:: python

    import pandas.io.data as web
    import datetime
    start = datetime.datetime(2010, 1, 1)
    end = datetime.datetime(2013, 1, 27)
    f=web.DataReader("F", 'google', start, end)
    f.ix['2010-01-04']

.. _remote_data.fred:

FRED
----

.. ipython:: python

    import pandas.io.data as web
    import datetime
    start = datetime.datetime(2010, 1, 1)
    end = datetime.datetime(2013, 1, 27)
    gdp=web.DataReader("GDP", "fred", start, end)
    gdp.ix['2013-01-01']

    # Multiple series:
    inflation = web.DataReader(["CPIAUCSL", "CPILFESL"], "fred", start, end)
    inflation.head()
.. _remote_data.ff:

Fama/French
-----------

Dataset names are listed at `Fama/French Data Library
<http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>`__.

.. ipython:: python

    import pandas.io.data as web
    ip=web.DataReader("5_Industry_Portfolios", "famafrench")
    ip[4].ix[192607]

.. _remote_data.wb:

World Bank
----------

``pandas`` users can easily access thousands of panel data series from the
`World Bank's World Development Indicators <http://data.worldbank.org>`__
by using the ``wb`` I/O functions.

For example, if you wanted to compare the Gross Domestic Products per capita in
constant dollars in North America, you would use the ``search`` function:

.. code-block:: python

    In [1]: from pandas.io import wb

    In [2]: wb.search('gdp.*capita.*const').iloc[:,:2]
    Out[2]:
                         id                                               name
    3242            GDPPCKD             GDP per Capita, constant US$, millions
    5143     NY.GDP.PCAP.KD                 GDP per capita (constant 2005 US$)
    5145     NY.GDP.PCAP.KN                      GDP per capita (constant LCU)
    5147  NY.GDP.PCAP.PP.KD  GDP per capita, PPP (constant 2005 internation...

Then you would use the ``download`` function to acquire the data from the World
Bank's servers:

.. code-block:: python

    In [3]: dat = wb.download(indicator='NY.GDP.PCAP.KD', country=['US', 'CA', 'MX'], start=2005, end=2008)

    In [4]: print(dat)
                          NY.GDP.PCAP.KD
    country       year
    Canada        2008  36005.5004978584
                  2007  36182.9138439757
                  2006  35785.9698172849
                  2005  35087.8925933298
    Mexico        2008  8113.10219480083
                  2007  8119.21298908649
                  2006  7961.96818458178
                  2005  7666.69796097264
    United States 2008  43069.5819857208
                  2007  43635.5852068142
                  2006   43228.111147107
                  2005  42516.3934699993

The resulting dataset is a properly formatted ``DataFrame`` with a hierarchical
index, so it is easy to apply ``.groupby`` transformations to it:

.. code-block:: python

    In [6]: dat['NY.GDP.PCAP.KD'].groupby(level=0).mean()
    Out[6]:
    country
    Canada           35765.569188
    Mexico            7965.245332
    United States    43112.417952
    dtype: float64

Now imagine you want to compare GDP to the share of people with cellphone
contracts around the world.

.. code-block:: python

    In [7]: wb.search('cell.*%').iloc[:,:2]
    Out[7]:
                         id                                               name
    3990  IT.CEL.SETS.FE.ZS  Mobile cellular telephone users, female (% of ...
    3991  IT.CEL.SETS.MA.ZS  Mobile cellular telephone users, male (% of po...
    4027      IT.MOB.COV.ZS  Population coverage of mobile cellular telepho...

Notice that this second search was much faster than the first one because
``pandas`` now has a cached list of available data series.

.. code-block:: python

    In [13]: ind = ['NY.GDP.PCAP.KD', 'IT.MOB.COV.ZS']
    In [14]: dat = wb.download(indicator=ind, country='all', start=2011, end=2011).dropna()
    In [15]: dat.columns = ['gdp', 'cellphone']
    In [16]: print(dat.tail())
                            gdp  cellphone
    country   year
    Swaziland 2011  2413.952853       94.9
    Tunisia   2011  3687.340170      100.0
    Uganda    2011   405.332501      100.0
    Zambia    2011   767.911290       62.0
    Zimbabwe  2011   419.236086       72.4

Finally, we use the ``statsmodels`` package to assess the relationship between
our two variables using ordinary least squares regression. Unsurprisingly,
populations in rich countries tend to use cellphones at a higher rate:

.. code-block:: python

    In [17]: import numpy as np
    In [18]: import statsmodels.formula.api as smf
    In [19]: mod = smf.ols("cellphone ~ np.log(gdp)", dat).fit()
    In [20]: print(mod.summary())
                                OLS Regression Results
    ==============================================================================
    Dep. Variable:              cellphone   R-squared:                       0.297
    Model:                            OLS   Adj. R-squared:                  0.274
    Method:                 Least Squares   F-statistic:                     13.08
    Date:                Thu, 25 Jul 2013   Prob (F-statistic):            0.00105
    Time:                        15:24:42   Log-Likelihood:                -139.16
    No. Observations:                  33   AIC:                             282.3
    Df Residuals:                      31   BIC:                             285.3
    Df Model:                           1
    ===============================================================================
                      coef    std err          t      P>|t|      [95.0% Conf. Int.]
    -------------------------------------------------------------------------------
    Intercept      16.5110     19.071      0.866      0.393       -22.384    55.406
    np.log(gdp)     9.9333      2.747      3.616      0.001         4.331    15.535
    ==============================================================================
    Omnibus:                       36.054   Durbin-Watson:                   2.071
    Prob(Omnibus):                  0.000   Jarque-Bera (JB):              119.133
    Skew:                          -2.314   Prob(JB):                     1.35e-26
    Kurtosis:                      11.077   Cond. No.                         45.8
    ==============================================================================
