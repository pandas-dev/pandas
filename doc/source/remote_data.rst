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

Functions from :mod:`pandas.io.data` and :mod:`pandas.io.ga` extract data from various Internet sources into a DataFrame. Currently the following sources are supported:

    - :ref:`Yahoo! Finance<remote_data.yahoo>`
    - :ref:`Google Finance<remote_data.google>`
    - :ref:`St.Louis FED (FRED)<remote_data.fred>`
    - :ref:`Kenneth French's data library<remote_data.ff>`
    - :ref:`World Bank<remote_data.wb>`
    - :ref:`Google Analytics<remote_data.ga>`

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

.. _remote_data.yahoo_options:

Yahoo! Finance Options
----------------------
***Experimental***

The Options class allows the download of options data from Yahoo! Finance.

The ``get_all_data`` method downloads and caches option data for all expiry months
and provides a formatted ``DataFrame`` with a hierarchical index, so its easy to get
to the specific option you want.

.. ipython:: python

      from pandas.io.data import Options
      aapl = Options('aapl', 'yahoo')
      data = aapl.get_all_data()
      data.iloc[0:5, 0:5]

      #Show the $100 strike puts at all expiry dates:
      data.loc[(100, slice(None), 'put'),:].iloc[0:5, 0:5]

      #Show the volume traded of $100 strike puts at all expiry dates:
      data.loc[(100, slice(None), 'put'),'Vol'].head()

If you don't want to download all the data, more specific requests can be made.

.. ipython:: python

      import datetime
      expiry = datetime.date(2016, 1, 1)
      data = aapl.get_call_data(expiry=expiry)
      data.iloc[0:5:, 0:5]

Note that if you call ``get_all_data`` first, this second call will happen much faster,
as the data is cached.

If a given expiry date is not available, data for the next available expiry will be
returned (January 15, 2015 in the above example).

Available expiry dates can be accessed from the ``expiry_dates`` property.

.. ipython:: python

      aapl.expiry_dates
      data = aapl.get_call_data(expiry=aapl.expiry_dates[0])
      data.iloc[0:5:, 0:5]

A list-like object containing dates can also be passed to the expiry parameter,
returning options data for all expiry dates in the list.

.. ipython:: python

      data = aapl.get_near_stock_price(expiry=aapl.expiry_dates[0:3])
      data.iloc[0:5:, 0:5]

The ``month`` and ``year`` parameters can be used to get all options data for a given month.

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

Indicators
~~~~~~~~~~

Either from exploring the World Bank site, or using the search function included,
every world bank indicator is accessible.  

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

Country Codes
~~~~~~~~~~~~~

.. versionadded:: 0.15.1

The ``country`` argument accepts a string or list of mixed 
`two <http://en.wikipedia.org/wiki/ISO_3166-1_alpha-2>`__ or `three <http://en.wikipedia.org/wiki/ISO_3166-1_alpha-3>`__ character
ISO country codes, as well as dynamic `World Bank exceptions <http://data.worldbank.org/node/18>`__ to the ISO standards.

For a list of the the hard-coded country codes (used solely for error handling logic) see ``pandas.io.wb.country_codes``.

Problematic Country Codes & Indicators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   The World Bank's country list and indicators are dynamic. As of 0.15.1, 
   :func:`wb.download()` is more flexible.  To achieve this, the warning
   and exception logic changed.
   
The world bank converts some country codes,
in their response, which makes error checking by pandas difficult.
Retired indicators still persist in the search.

Given the new flexibility of 0.15.1, improved error handling by the user
may be necessary for fringe cases.

To help identify issues:

There are at least 4 kinds of country codes:

1. Standard (2/3 digit ISO) - returns data, will warn and error properly.
2. Non-standard (WB Exceptions) - returns data, but will falsely warn.
3. Blank - silently missing from the response.
4. Bad - causes the entire response from WB to fail, always exception inducing.

There are at least 3 kinds of indicators:

1. Current - Returns data.
2. Retired - Appears in search results, yet won't return data. 
3. Bad - Will not return data.

Use the ``errors`` argument to control warnings and exceptions.  Setting
errors to ignore or warn, won't stop failed responses.  (ie, 100% bad
indicators, or a single "bad" (#4 above) country code).  

See docstrings for more info.

.. _remote_data.ga:

Google Analytics
----------------

The :mod:`~pandas.io.ga` module provides a wrapper for
`Google Analytics API <https://developers.google.com/analytics/devguides>`__
to simplify retrieving traffic data.
Result sets are parsed into a pandas DataFrame with a shape and data types
derived from the source table.

Configuring Access to Google Analytics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first thing you need to do is to setup accesses to Google Analytics API. Follow the steps below:

#. In the `Google Developers Console <https://console.developers.google.com>`__
    #. enable the Analytics API
    #. create a new project
    #. create a new Client ID for an "Installed Application" (in the "APIs & auth / Credentials section" of the newly created project)
    #. download it (JSON file)
#. On your machine
    #. rename it to ``client_secrets.json``
    #. move it to the ``pandas/io`` module directory

The first time you use the :func:`read_ga` funtion, a browser window will open to ask you to authentify to the Google API. Do proceed.

Using the Google Analytics API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following will fetch users and pageviews (metrics) data per day of the week, for the first semester of 2014, from a particular property.

.. code-block:: python

    import pandas.io.ga as ga
    ga.read_ga(
        account_id  = "2360420",
        profile_id  = "19462946",
        property_id = "UA-2360420-5",
        metrics     = ['users', 'pageviews'],
        dimensions  = ['dayOfWeek'],
        start_date  = "2014-01-01",
        end_date    = "2014-08-01",
        index_col   = 0,
        filters     = "pagePath=~aboutus;ga:country==France",
    )

The only mandatory arguments are ``metrics,`` ``dimensions`` and ``start_date``. We can only strongly recommend you to always specify the ``account_id``, ``profile_id`` and ``property_id`` to avoid accessing the wrong data bucket in Google Analytics.

The ``index_col`` argument indicates which dimension(s) has to be taken as index.

The ``filters`` argument indicates the filtering to apply to the query. In the above example, the page has URL has to contain ``aboutus`` AND the visitors country has to be France.

Detailed informations in the followings:

* `pandas & google analytics, by yhat <http://blog.yhathq.com/posts/pandas-google-analytics.html>`__
* `Google Analytics integration in pandas, by Chang She <http://quantabee.wordpress.com/2012/12/17/google-analytics-pandas/>`__
* `Google Analytics Dimensions and Metrics Reference <https://developers.google.com/analytics/devguides/reporting/core/dimsmets>`_

