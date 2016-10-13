.. _remote_data:

.. currentmodule:: pandas

******************
Remote Data Access
******************

.. _remote_data.pandas_datareader:

DataReader
----------

The sub-package ``pandas.io.data`` is removed in favor of a separately
installable `pandas-datareader package
<https://github.com/pandas-dev/pandas-datareader>`_. This will allow the data
modules to be independently updated to your pandas installation. The API for
``pandas-datareader v0.1.1`` is the same as in ``pandas v0.16.1``.
(:issue:`8961`)

   You should replace the imports of the following:

   .. code-block:: python

      from pandas.io import data, wb

   With:

   .. code-block:: python

      from pandas_datareader import data, wb


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

The first time you use the :func:`read_ga` function, a browser window will open to ask you to authentify to the Google API. Do proceed.

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

The only mandatory arguments are ``metrics,`` ``dimensions`` and ``start_date``. We strongly recommend that you always specify the ``account_id``, ``profile_id`` and ``property_id`` to avoid accessing the wrong data bucket in Google Analytics.

The ``index_col`` argument indicates which dimension(s) has to be taken as index.

The ``filters`` argument indicates the filtering to apply to the query. In the above example, the page URL has to contain ``aboutus`` AND the visitors country has to be France.

Detailed information in the following:

* `pandas & google analytics, by yhat <http://blog.yhathq.com/posts/pandas-google-analytics.html>`__
* `Google Analytics integration in pandas, by Chang She <http://quantabee.wordpress.com/2012/12/17/google-analytics-pandas/>`__
* `Google Analytics Dimensions and Metrics Reference <https://developers.google.com/analytics/devguides/reporting/core/dimsmets>`_
